! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TPM_FFTW
!   Author.
!   -------
!     George Mozdzynski
!
!   Modifications.
!   -------------- 
!     Original      October 2014
!     R. El Khatib  01-Sep-2015 More subroutines for better modularity
!     R. El Khatib  08-Jun-2023 LALL_FFTW for better flexibility
!     W. Deconinck  17-Jun-2024 Replace legacy FFTW interface with the FFTW3 interface, add documentation, and improve clarity

USE, INTRINSIC :: ISO_C_BINDING

USE PARKIND1        ,ONLY : JPIM, JPRB, JPRD
USE MPL_MODULE      ,ONLY : MPL_MYRANK, MPL_RANK, MPL_NUMPROC, &
  & MPL_SEND, MPL_RECV
USE YOMHOOK         ,ONLY : LHOOK, DR_HOOK, JPHOOK
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE OMP_LIB         ,ONLY : OMP_GET_MAX_THREADS, OMP_GET_THREAD_NUM

IMPLICIT NONE

SAVE

#ifdef __NEC__
! From NLC (NEC Numeric Library Collection)
#include "aslfftw3.f03"
#define FFTW_NO_SIMD 0
#else
#include "fftw3.f03"
#endif

! For now we use the LEGACY FFTW INTERFACE still, to be removed once we are sure no
! FFTW-pretending libraries are only implementing the legacy interface.
#define LEGACY_FFTW_INTERFACE 1

PRIVATE
PUBLIC INIT_PLANS_FFTW, DESTROY_PLANS_FFTW, FFTW_RESOL, TW, EXEC_FFTW, EXEC_EFFTW

!> @brief Cache state and execution settings for FFTW plans at one active resolution.
TYPE FFTW_TYPE
  !! Number of cached plans currently stored for each transform length `KN`.
  INTEGER(KIND=JPIM),ALLOCATABLE :: N_PLANS(:)

  !! Head nodes of the linked lists storing cached plans for each transform length.
  TYPE(FFTW_PLAN),POINTER :: FFTW_PLANS(:) => NULL()

  !! Largest transform length for which this cache has been initialised.
  INTEGER(KIND=JPIM) :: N_MAX=0

  !! N_MAX_PLANS: Maximum number of cached plans retained for any one transform length. The practical cache key is (KN, KTYPE, KLOT).
  !! EXEC_FFTW passes KN = KRLEN and KLOT = NBATCH, with NBATCH = 1 when LALL_FFTW=.FALSE. and NBATCH = KFIELDS when LALL_FFTW=.TRUE.
  !! Since KTYPE is only 1 or -1, eviction starts when one fixed KRLEN needs a fifth distinct (KTYPE, NBATCH) combination, when N_MAX_PLANS=4.
  INTEGER(KIND=JPIM) :: N_MAX_PLANS=4

  !! If true, execute all fields in one batched FFTW call; otherwise execute one field at a time.
  LOGICAL            :: LALL_FFTW=.FALSE.
END TYPE FFTW_TYPE

INTEGER(KIND=JPIM), PARAMETER :: NTYPE_C2R=1, NTYPE_R2C=-1
INTEGER(KIND=JPIM), PARAMETER :: TPM_FFTW_PLAN_FLAGS = FFTW_ESTIMATE + FFTW_NO_SIMD
INTEGER(KIND=JPIM), PARAMETER :: NPLAN_ID_UNINITIALISED = 123456
INTEGER(KIND=JPIM), PARAMETER :: NPLAN_ID_DESTROYED     = 999999

!> @brief Metadata for one cached FFTW plan together with its linked-list successor.
TYPE FFTW_PLAN
  !! Sentinel used to detect uninitialised or recycled plan records.
  INTEGER(KIND=JPIM) :: NPLAN_ID=NPLAN_ID_UNINITIALISED

  !! FFTW plan handle returned by the modern FFTW interface.
  TYPE(C_PTR)        :: NPLAN=C_NULL_PTR

  !! Batch size used when the plan was created.
  INTEGER(KIND=JPIM) :: NLOT

  !! Transform direction selector associated with the plan (`1` or `-1`).
  INTEGER(KIND=JPIM) :: NTYPE

  !! Pointer to the next cached plan for the same transform length.
  TYPE(FFTW_PLAN),POINTER :: NEXT_PLAN => NULL()
END TYPE FFTW_PLAN

!> @brief Array of per-resolution FFTW caches indexed by the active transform handle.
TYPE(FFTW_TYPE),ALLOCATABLE,TARGET :: FFTW_RESOL(:)
!> @brief Pointer to the FFTW cache associated with the currently active resolution.
TYPE(FFTW_TYPE),POINTER     :: TW

! Redefine JPCB as some FFTW implementations already define it in their fftw3.f03
#define JPCB TPM_FFTW_JPCB
INTEGER, PARAMETER :: JPCB = MERGE(C_DOUBLE_COMPLEX, C_FLOAT_COMPLEX, JPRB == JPRD)

! Per-thread persistent FFTW scratch buffers, registered once in INIT_PLANS_FFTW,
! grown on demand in ENSURE_FFT_BUFFER, with slots freed in DESTROY_PLANS_FFTW.
! Replaces the per-call TPM_FFTW_ALLOC_COMPLEX/TPM_FFTW_FREE pair in EXEC_FFTW_IMPL.
! EXEC_FFTW_IMPL runs inside OMP parallel regions; slot index = OMP_GET_THREAD_NUM()+1.
TYPE(C_PTR),       ALLOCATABLE :: ZFFT_PERSIST_PTRS(:)
INTEGER(C_SIZE_T), ALLOCATABLE :: ZFFT_PERSIST_CAPS(:)

! FFTW planning + wisdom handling, controlled with ECTRANS_FFTW_WISDOM environment variable:
!   * CREATE_PLAN_FFTW uses FFTW_MEASURE (slow first planning, faster FFTs subsequently)
!   * INIT_PLANS_FFTW imports wisdom from <ECTRANS_FFTW_WISDOM>.{dp,sp} if the file is present
!   * DESTROY_PLANS_FFTW exports the unified wisdom to <ECTRANS_FFTW_WISDOM>.{dp,sp}
! Default (ECTRANS_FFTW_WISDOM unset) behaviour: TPM_FFTW_PLAN_FLAGS
! (FFTW_ESTIMATE+FFTW_NO_SIMD), no wisdom I/O, bit-identical to the unoptimised path.
LOGICAL :: LUSE_MEASURE     = .FALSE.
LOGICAL :: LWISDOM_IMPORTED = .FALSE.

! Bytes-per-word for the dense packing used by the wisdom transport code below.
INTEGER(KIND=JPIM), PARAMETER :: IBPW_WISDOM = 4

! libc free() for releasing the malloc'd string returned by
! fftw{,f}_export_wisdom_to_string
INTERFACE
  SUBROUTINE C_FREE(PTR) BIND(C, NAME='free')
    IMPORT :: C_PTR
    TYPE(C_PTR), VALUE :: PTR
  END SUBROUTINE C_FREE
END INTERFACE


INTERFACE TPM_FFTW_PLAN_MANY_DFT_C2R
#if LEGACY_FFTW_INTERFACE
  PROCEDURE TPM_FFTW_PLAN_MANY_DFT_C2R_LEGACY
  PROCEDURE TPM_FFTWF_PLAN_MANY_DFT_C2R_LEGACY
#else
  PROCEDURE FFTW_PLAN_MANY_DFT_C2R   ! defined in fftw3.f03
  PROCEDURE FFTWF_PLAN_MANY_DFT_C2R  ! defined in fftw3.f03
#endif
END INTERFACE TPM_FFTW_PLAN_MANY_DFT_C2R

INTERFACE TPM_FFTW_PLAN_MANY_DFT_R2C
#if LEGACY_FFTW_INTERFACE
  PROCEDURE TPM_FFTW_PLAN_MANY_DFT_R2C_LEGACY
  PROCEDURE TPM_FFTWF_PLAN_MANY_DFT_R2C_LEGACY
#else
  PROCEDURE FFTW_PLAN_MANY_DFT_R2C  ! defined in fftw3.f03
  PROCEDURE FFTWF_PLAN_MANY_DFT_R2C ! defined in fftw3.f03
#endif
END INTERFACE TPM_FFTW_PLAN_MANY_DFT_R2C

INTERFACE TPM_FFTW_EXECUTE_DFT_C2R
#if LEGACY_FFTW_INTERFACE
  PROCEDURE TPM_FFTW_EXECUTE_DFT_C2R_LEGACY
  PROCEDURE TPM_FFTWF_EXECUTE_DFT_C2R_LEGACY
#else
  PROCEDURE FFTW_EXECUTE_DFT_C2R    ! defined in fftw3.f03
  PROCEDURE FFTWF_EXECUTE_DFT_C2R   ! defined in fftw3.f03
#endif
  PROCEDURE TPM_FFTW_EXECUTE_DFT_C2R_RANK2
END INTERFACE TPM_FFTW_EXECUTE_DFT_C2R

INTERFACE TPM_FFTW_EXECUTE_DFT_R2C
#if LEGACY_FFTW_INTERFACE
  PROCEDURE TPM_FFTW_EXECUTE_DFT_R2C_LEGACY
  PROCEDURE TPM_FFTWF_EXECUTE_DFT_R2C_LEGACY
#else
  PROCEDURE FFTW_EXECUTE_DFT_R2C   ! defined in fftw3.f03
  PROCEDURE FFTWF_EXECUTE_DFT_R2C  ! defined in fftw3.f03
#endif
  PROCEDURE TPM_FFTW_EXECUTE_DFT_R2C_RANK2
END INTERFACE TPM_FFTW_EXECUTE_DFT_R2C

! ------------------------------------------------------------------
CONTAINS
! ------------------------------------------------------------------

FUNCTION TPM_FFTW_ALLOC_COMPLEX(N) RESULT(RES)
!> @brief Allocate a complex work buffer through the precision-specific FFTW backend.
!!
!! @param[in] N Number of complex elements to allocate.
IMPLICIT NONE
INTEGER(KIND=C_SIZE_T),INTENT(IN) :: N
TYPE(C_PTR) :: RES
IF (JPRB == JPRD) THEN
  RES=FFTW_ALLOC_COMPLEX(N)
ELSE
  RES=FFTWF_ALLOC_COMPLEX(N)
END IF
END FUNCTION TPM_FFTW_ALLOC_COMPLEX


SUBROUTINE TPM_FFTW_FREE(PTR)
!> @brief Release a complex work buffer allocated through the FFTW precision wrapper.
!!
!! @param[in] PTR C pointer to a work buffer allocated by `TPM_FFTW_ALLOC_COMPLEX`.
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: PTR
IF (JPRB == JPRD) THEN
  CALL FFTW_FREE(PTR)
ELSE
  CALL FFTWF_FREE(PTR)
END IF
END SUBROUTINE TPM_FFTW_FREE


SUBROUTINE WISDOM_FILENAME(CDPATH, KSTAT)
!> @brief Resolve the FFTW wisdom file path from `ECTRANS_FFTW_WISDOM`.
!!
!! @param[out] CDPATH Wisdom path with `.dp`/`.sp` suffix appended by precision.
!! @param[out] KSTAT Zero on success, non-zero when the env var is unset or empty.
IMPLICIT NONE
CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: CDPATH
INTEGER(KIND=JPIM),            INTENT(OUT) :: KSTAT
CHARACTER(LEN=4096) :: ZBUF
INTEGER :: ILEN, ISTAT

CALL GET_ENVIRONMENT_VARIABLE('ECTRANS_FFTW_WISDOM', VALUE=ZBUF, &
  & LENGTH=ILEN, STATUS=ISTAT)
IF (ISTAT /= 0 .OR. ILEN == 0 .OR. ILEN > LEN(ZBUF)) THEN
  KSTAT = 1
  RETURN
ENDIF
IF (JPRB == JPRD) THEN
  CDPATH = ZBUF(1:ILEN) // '.dp'
ELSE
  CDPATH = ZBUF(1:ILEN) // '.sp'
ENDIF
KSTAT = 0
END SUBROUTINE WISDOM_FILENAME


SUBROUTINE STR_TO_CSTR(CDIN, CDOUT)
!> @brief Convert a Fortran string to a NUL-terminated `C_CHAR` array for `const char *` dummies.
!!
!! @param[in] CDIN Fortran character string to convert.
!! @param[out] CDOUT Converted NUL-terminated C string.
IMPLICIT NONE
CHARACTER(LEN=*),                    INTENT(IN)  :: CDIN
CHARACTER(KIND=C_CHAR), ALLOCATABLE, INTENT(OUT) :: CDOUT(:)
INTEGER(KIND=JPIM) :: IL, J
IL = LEN_TRIM(CDIN)
ALLOCATE(CDOUT(IL+1))
DO J = 1, IL
  CDOUT(J) = CDIN(J:J)
END DO
CDOUT(IL+1) = C_NULL_CHAR
END SUBROUTINE STR_TO_CSTR


SUBROUTINE IMPORT_FFTW_WISDOM
!> @brief Import FFTW wisdom from `<ECTRANS_FFTW_WISDOM>.{dp,sp}`; all ranks read the file.
!!
!! @note No-op when the env var is unset, and only attempted once per run.
IMPLICIT NONE
CHARACTER(LEN=:),       ALLOCATABLE :: ZPATH
CHARACTER(KIND=C_CHAR), ALLOCATABLE :: ZCPATH(:)
INTEGER(KIND=JPIM) :: ISTAT
INTEGER(C_INT)     :: IRET
LOGICAL            :: LLPRINT

IF (LWISDOM_IMPORTED) RETURN
CALL WISDOM_FILENAME(ZPATH, ISTAT)
IF (ISTAT /= 0) RETURN
CALL STR_TO_CSTR(ZPATH, ZCPATH)

IF (JPRB == JPRD) THEN
  IRET = FFTW_IMPORT_WISDOM_FROM_FILENAME(ZCPATH)
ELSE
  IRET = FFTWF_IMPORT_WISDOM_FROM_FILENAME(ZCPATH)
ENDIF
! Mark imported even on failure: missing/malformed file should not be
! retried on subsequent resolution setups.
LWISDOM_IMPORTED = .TRUE.
LLPRINT = (MPL_NUMPROC < 1) .OR. (MPL_RANK == 1)
IF (LLPRINT) THEN
  IF (IRET == 1) THEN
    WRITE(0,'(A,A)') 'TPM_FFTW: imported FFTW wisdom from ', TRIM(ZPATH)
  ELSE
    WRITE(0,'(A,A,A)') &
      & 'TPM_FFTW: WARNING failed to import FFTW wisdom from ', &
      & TRIM(ZPATH), ' (continuing with empty wisdom)'
  ENDIF
ENDIF
END SUBROUTINE IMPORT_FFTW_WISDOM


SUBROUTINE PACK_LOCAL_WISDOM(IBUF, IWORDS, ILEN_BYTES)
!> @brief Export the local FFTW planner wisdom to a packed JPIM word buffer.
!!
!! @param[out] IBUF Packed wisdom bytes, `IBPW_WISDOM` bytes per word, zero-padded tail.
!! @param[out] IWORDS Number of words stored in `IBUF`.
!! @param[out] ILEN_BYTES Actual NUL-terminated byte length (bytes past it are padding).
!! @note The FFTW C string is malloc'd and released here via `C_FREE`.
IMPLICIT NONE
INTEGER(KIND=JPIM), ALLOCATABLE, INTENT(OUT) :: IBUF(:)
INTEGER(KIND=JPIM),              INTENT(OUT) :: IWORDS
INTEGER(KIND=JPIM),              INTENT(OUT) :: ILEN_BYTES
CHARACTER(KIND=C_CHAR), POINTER :: ZLOC(:)
TYPE(C_PTR)        :: IPTR
INTEGER(KIND=JPIM) :: JJ, IHUGE, IWORD, IBPOS

IF (JPRB == JPRD) THEN
  IPTR = FFTW_EXPORT_WISDOM_TO_STRING()
ELSE
  IPTR = FFTWF_EXPORT_WISDOM_TO_STRING()
ENDIF

ILEN_BYTES = 0
IF (C_ASSOCIATED(IPTR)) THEN
  IHUGE = 16 * 1024 * 1024 !! 16 MB
  CALL C_F_POINTER(IPTR, ZLOC, [IHUGE])
  DO JJ = 1, IHUGE
    IF (ZLOC(JJ) == C_NULL_CHAR) EXIT
    ILEN_BYTES = ILEN_BYTES + 1
  ENDDO
ENDIF

IWORDS = (ILEN_BYTES + IBPW_WISDOM - 1) / IBPW_WISDOM
ALLOCATE(IBUF(MAX(IWORDS, 1)))
IBUF(:) = 0
DO JJ = 1, ILEN_BYTES
  IWORD = (JJ - 1) / IBPW_WISDOM + 1
  IBPOS = MOD(JJ - 1, IBPW_WISDOM)
  IBUF(IWORD) = IOR(IBUF(IWORD), &
    & ISHFT(IAND(ICHAR(ZLOC(JJ)), 255), 8 * IBPOS))
ENDDO

IF (C_ASSOCIATED(IPTR)) CALL C_FREE(IPTR)
NULLIFY(ZLOC)
END SUBROUTINE PACK_LOCAL_WISDOM


SUBROUTINE UNPACK_AND_IMPORT_WISDOM(IBUF, ILEN_BYTES, LDOK)
!> @brief Unpack a packed word buffer and import it into the local FFTW planner.
!!
!! @param[in] IBUF Packed wisdom bytes (`IBPW_WISDOM` bytes per word).
!! @param[in] ILEN_BYTES Number of significant bytes in `IBUF`.
!! @param[out] LDOK True when the import succeeded (or there was nothing to import).
!! @note FFTW wisdom is additive, so this merges the incoming plans with local knowledge.
IMPLICIT NONE
INTEGER(KIND=JPIM), INTENT(IN)  :: IBUF(:)
INTEGER(KIND=JPIM), INTENT(IN)  :: ILEN_BYTES
LOGICAL,            INTENT(OUT) :: LDOK
CHARACTER(KIND=C_CHAR), ALLOCATABLE, TARGET :: ZSTR(:)
INTEGER(KIND=JPIM) :: JJ, IWORD, IBPOS, IBYTE
INTEGER(C_INT)     :: IRET

LDOK = .TRUE.
IF (ILEN_BYTES <= 0) RETURN

ALLOCATE(ZSTR(ILEN_BYTES + 1))
DO JJ = 1, ILEN_BYTES
  IWORD = (JJ - 1) / IBPW_WISDOM + 1
  IBPOS = MOD(JJ - 1, IBPW_WISDOM)
  IBYTE = IAND(ISHFT(IBUF(IWORD), -8 * IBPOS), 255)
  ZSTR(JJ) = ACHAR(IBYTE)
ENDDO
ZSTR(ILEN_BYTES + 1) = C_NULL_CHAR

IF (JPRB == JPRD) THEN
  IRET = FFTW_IMPORT_WISDOM_FROM_STRING(ZSTR)
ELSE
  IRET = FFTWF_IMPORT_WISDOM_FROM_STRING(ZSTR)
ENDIF
DEALLOCATE(ZSTR)
LDOK = (IRET == 1)
END SUBROUTINE UNPACK_AND_IMPORT_WISDOM


SUBROUTINE EXPORT_FFTW_WISDOM
!> @brief Export FFTW wisdom to `<ECTRANS_FFTW_WISDOM>.{dp,sp}`.
!!
!! In MPI runs each rank only measures plans for its own latitudes, so any single
!! rank's wisdom is insufficient. FFTW wisdom is additive (importing A then B yields
!! the union), which lets us reduce it in a butterfly tree instead of flat-gathering:
!! at step k (stride 2^k), ranks with bit k clear receive from the peer with bit k
!! set, import the peer's wisdom locally, and re-export the accumulated union upward.
!! Ranks with bit k set send once and then drop out. After ceil(log2(NPROC)) steps
!! rank 1 holds the full union and writes it to disk. Per-hop message size is bounded
!! by the number of *unique* plans in the accumulated union, not by the sum over ranks.
!! @note No-op when the env var is unset. Uses `MPL_SEND`/`MPL_RECV` only.
IMPLICIT NONE
CHARACTER(LEN=:),       ALLOCATABLE :: ZPATH
CHARACTER(KIND=C_CHAR), ALLOCATABLE :: ZCPATH(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IBUF_LOC(:), IBUF_RECV(:)
INTEGER(KIND=JPIM) :: ISTAT, ILEN_LOC, IWORDS_LOC
INTEGER(KIND=JPIM) :: IHEADER(2)
INTEGER(KIND=JPIM) :: IR0, ISTRIDE, IPEER_R0, IPEER_1B
INTEGER(C_INT)     :: IRET
LOGICAL            :: LLROOT, LLMPI, LLACTIVE, LLOK
INTEGER(KIND=JPIM), PARAMETER :: ITAG_HDR = 12341
INTEGER(KIND=JPIM), PARAMETER :: ITAG_BUF = 12342

! this does nothing if ECTRANS_FFTW_WISDOM environment variable not set
CALL WISDOM_FILENAME(ZPATH, ISTAT)
IF (ISTAT /= 0) RETURN

LLMPI  = (MPL_NUMPROC >= 1)
LLROOT = (.NOT. LLMPI) .OR. (MPL_RANK == 1)

! Snapshot the local planner state as a packed buffer.
CALL PACK_LOCAL_WISDOM(IBUF_LOC, IWORDS_LOC, ILEN_LOC)

IF (LLMPI .AND. MPL_NUMPROC > 1) THEN
  ! Butterfly reduce-to-root (rank 1 is the root). Ranks are 1-based
  ! in MPL; work with the 0-based rank IR0 for the bit arithmetic.
  IR0      = MPL_RANK - 1
  ISTRIDE  = 1
  LLACTIVE = .TRUE.
  DO WHILE (ISTRIDE < MPL_NUMPROC .AND. LLACTIVE)
    IF (MOD(IR0, 2*ISTRIDE) == 0) THEN
      ! Receiver: pull from IR0 + ISTRIDE if that rank exists.
      IPEER_R0 = IR0 + ISTRIDE
      IF (IPEER_R0 < MPL_NUMPROC) THEN
        IPEER_1B = IPEER_R0 + 1
        CALL MPL_RECV(IHEADER, KSOURCE=IPEER_1B, KTAG=ITAG_HDR, &
          & CDSTRING='TPM_FFTW:WISDOM_TREE:HDR')
        IF (IHEADER(2) > 0) THEN
          ALLOCATE(IBUF_RECV(IHEADER(2)))
          CALL MPL_RECV(IBUF_RECV, KSOURCE=IPEER_1B, KTAG=ITAG_BUF, &
            & CDSTRING='TPM_FFTW:WISDOM_TREE:BUF')
          CALL UNPACK_AND_IMPORT_WISDOM(IBUF_RECV, IHEADER(1), LLOK)
          IF (.NOT. LLOK) THEN
            WRITE(0,'(A,I0)') &
              & 'TPM_FFTW: WARNING failed to merge wisdom from rank ', IPEER_1B
          ENDIF
          DEALLOCATE(IBUF_RECV)
          ! Re-export the (now larger) union so subsequent steps
          ! forward the accumulated wisdom upward.
          DEALLOCATE(IBUF_LOC)
          CALL PACK_LOCAL_WISDOM(IBUF_LOC, IWORDS_LOC, ILEN_LOC)
        ENDIF
      ENDIF
    ELSE
      ! Sender: push to IR0 - ISTRIDE and drop out.
      IPEER_R0 = IR0 - ISTRIDE
      IPEER_1B = IPEER_R0 + 1
      IHEADER(1) = ILEN_LOC
      IHEADER(2) = IWORDS_LOC
      CALL MPL_SEND(IHEADER, KDEST=IPEER_1B, KTAG=ITAG_HDR, &
        & CDSTRING='TPM_FFTW:WISDOM_TREE:HDR')
      IF (IWORDS_LOC > 0) THEN
        CALL MPL_SEND(IBUF_LOC(1:IWORDS_LOC), KDEST=IPEER_1B, &
          & KTAG=ITAG_BUF, CDSTRING='TPM_FFTW:WISDOM_TREE:BUF')
      ENDIF
      LLACTIVE = .FALSE.
    ENDIF
    ISTRIDE = ISTRIDE * 2
  ENDDO
ENDIF

DEALLOCATE(IBUF_LOC)

! Rank 1's planner now holds the accumulated union; write it out.
IF (LLROOT) THEN
  CALL STR_TO_CSTR(ZPATH, ZCPATH)
  IF (JPRB == JPRD) THEN
    IRET = FFTW_EXPORT_WISDOM_TO_FILENAME(ZCPATH)
  ELSE
    IRET = FFTWF_EXPORT_WISDOM_TO_FILENAME(ZCPATH)
  ENDIF
  IF (IRET == 1) THEN
    WRITE(0,'(A,A)') 'TPM_FFTW: exported unified FFTW wisdom to ', TRIM(ZPATH)
  ELSE
    WRITE(0,'(A,A)') 'TPM_FFTW: WARNING failed to export FFTW wisdom to ', &
      & TRIM(ZPATH)
  ENDIF
ENDIF

END SUBROUTINE EXPORT_FFTW_WISDOM


SUBROUTINE ENSURE_FFT_BUFFER(KSIZE_COMPLEX, PTR_OUT)
!> @brief Provide the calling OMP thread's persistent FFTW scratch buffer.
!!
!! @param[in] KSIZE_COMPLEX Minimum number of complex elements the buffer must hold.
!! @param[out] PTR_OUT C pointer to the thread-owned scratch buffer (FFTW-aligned).
!! @note Grows on demand by freeing the old FFTW allocation and re-allocating via
!!   `TPM_FFTW_ALLOC_COMPLEX`, preserving the SIMD alignment required by MEASURE plans;
!!   never shrinks. Each OMP thread exclusively owns slot `OMP_GET_THREAD_NUM()+1`, so the
!!   buffer grow path needs no locking; the registry itself is grown under the
!!   `FFTW_SCRATCH_GROW` critical region if more threads appear than were registered.
IMPLICIT NONE
INTEGER(C_SIZE_T), INTENT(IN)  :: KSIZE_COMPLEX
TYPE(C_PTR),       INTENT(OUT) :: PTR_OUT
INTEGER(KIND=JPIM) :: ITID
INTEGER(C_SIZE_T)  :: ISIZE

ISIZE = MAX(KSIZE_COMPLEX, 1_C_SIZE_T)
ITID = OMP_GET_THREAD_NUM() + 1

IF (.NOT. ALLOCATED(ZFFT_PERSIST_PTRS)) THEN
  CALL ABORT_TRANS('ENSURE_FFT_BUFFER: scratch registry not allocated, call INIT_PLANS_FFTW first')
ENDIF

IF (ITID > SIZE(ZFFT_PERSIST_PTRS)) THEN
  !$OMP CRITICAL (FFTW_SCRATCH_GROW)
  ! Recheck under lock: another thread may have grown the registry already.
  IF (ITID > SIZE(ZFFT_PERSIST_PTRS)) CALL GROW_FFT_REGISTRY(ITID)
  !$OMP END CRITICAL (FFTW_SCRATCH_GROW)
ENDIF

! From here on the calling thread exclusively owns slot ITID: no lock needed.
IF (ZFFT_PERSIST_CAPS(ITID) < ISIZE) THEN
  IF (C_ASSOCIATED(ZFFT_PERSIST_PTRS(ITID))) CALL TPM_FFTW_FREE(ZFFT_PERSIST_PTRS(ITID))
  ZFFT_PERSIST_PTRS(ITID) = TPM_FFTW_ALLOC_COMPLEX(ISIZE)
  ZFFT_PERSIST_CAPS(ITID) = ISIZE
ENDIF
PTR_OUT = ZFFT_PERSIST_PTRS(ITID)
END SUBROUTINE ENSURE_FFT_BUFFER


SUBROUTINE GROW_FFT_REGISTRY(KMIN)
!> @brief Grow the persistent-scratch registry to hold at least slot `KMIN`, preserving slots.
!!
!! @param[in] KMIN Minimum number of registry slots required (1-based thread index).
!! @note Caller must hold the `FFTW_SCRATCH_GROW` critical region.
IMPLICIT NONE
INTEGER(KIND=JPIM), INTENT(IN) :: KMIN
INTEGER(KIND=JPIM) :: NOLD, NNEW, JN
TYPE(C_PTR), ALLOCATABLE :: ZNEW_PTRS(:)
INTEGER(C_SIZE_T), ALLOCATABLE :: ZNEW_CAPS(:)

NOLD = SIZE(ZFFT_PERSIST_PTRS)
NNEW = MAX(KMIN, MAX(1_JPIM, 2_JPIM*NOLD))
ALLOCATE(ZNEW_PTRS(NNEW))
ALLOCATE(ZNEW_CAPS(NNEW))
ZNEW_PTRS(1:NOLD) = ZFFT_PERSIST_PTRS(:)
ZNEW_CAPS(1:NOLD) = ZFFT_PERSIST_CAPS(:)
DO JN=NOLD+1,NNEW
  ZNEW_PTRS(JN) = C_NULL_PTR
  ZNEW_CAPS(JN) = 0_C_SIZE_T
ENDDO
DEALLOCATE(ZFFT_PERSIST_PTRS)
DEALLOCATE(ZFFT_PERSIST_CAPS)
ALLOCATE(ZFFT_PERSIST_PTRS(NNEW))
ALLOCATE(ZFFT_PERSIST_CAPS(NNEW))
ZFFT_PERSIST_PTRS(:) = ZNEW_PTRS(:)
ZFFT_PERSIST_CAPS(:) = ZNEW_CAPS(:)
DEALLOCATE(ZNEW_PTRS)
DEALLOCATE(ZNEW_CAPS)
END SUBROUTINE GROW_FFT_REGISTRY


SUBROUTINE TPM_FFTW_DESTROY_PLAN(PLAN)
!> @brief Destroy an FFTW plan using the precision-specific backend.
!!
!! @param[in] PLAN FFTW plan handle to destroy.
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: PLAN
#if LEGACY_FFTW_INTERFACE == 1
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY
PLAN_LEGACY=TPM_FFTW_LEGACY_PLAN_FROM_CPTR(PLAN)
IF (JPRB == JPRD) THEN
  CALL DFFTW_DESTROY_PLAN(PLAN_LEGACY)
ELSE
  CALL SFFTW_DESTROY_PLAN(PLAN_LEGACY)
END IF
#else
IF (JPRB == JPRD) THEN
  CALL FFTW_DESTROY_PLAN(PLAN)
ELSE
  CALL FFTWF_DESTROY_PLAN(PLAN)
END IF
#endif
END SUBROUTINE TPM_FFTW_DESTROY_PLAN


SUBROUTINE TPM_FFTW_EXECUTE_DFT_C2R_RANK2(PLAN,IN,OUT)
!> @brief Execute a complex-to-real FFTW plan on rank-2 arrays via flattened views.
!!
!! @param[in] PLAN FFTW plan handle to execute.
!! @param[inout] IN Rank-2 complex work array passed to FFTW as a contiguous rank-1 view.
!! @param[inout] OUT Rank-2 real work array passed to FFTW as a contiguous rank-1 view.
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: PLAN
COMPLEX(KIND=JPCB),INTENT(INOUT),CONTIGUOUS,TARGET :: IN(:,:)
REAL(KIND=JPRB),INTENT(INOUT),CONTIGUOUS,TARGET :: OUT(:,:)
COMPLEX(KIND=JPCB),POINTER :: IN_FLAT(:)
REAL(KIND=JPRB),POINTER :: OUT_FLAT(:)
IN_FLAT(1:SIZE(IN)) => IN
OUT_FLAT(1:SIZE(OUT)) => OUT
CALL TPM_FFTW_EXECUTE_DFT_C2R(PLAN,IN_FLAT,OUT_FLAT)
END SUBROUTINE TPM_FFTW_EXECUTE_DFT_C2R_RANK2


SUBROUTINE TPM_FFTW_EXECUTE_DFT_R2C_RANK2(PLAN,IN,OUT)
!> @brief Execute a real-to-complex FFTW plan on rank-2 arrays via flattened views.
!!
!! @param[in] PLAN FFTW plan handle to execute.
!! @param[inout] IN Rank-2 real work array passed to FFTW as a contiguous rank-1 view.
!! @param[inout] OUT Rank-2 complex work array passed to FFTW as a contiguous rank-1 view.
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: PLAN
REAL(KIND=JPRB),INTENT(INOUT),CONTIGUOUS,TARGET :: IN(:,:)
COMPLEX(KIND=JPCB),INTENT(INOUT),CONTIGUOUS,TARGET :: OUT(:,:)
REAL(KIND=JPRB),POINTER :: IN_FLAT(:)
COMPLEX(KIND=JPCB),POINTER :: OUT_FLAT(:)
IN_FLAT(1:SIZE(IN)) => IN
OUT_FLAT(1:SIZE(OUT)) => OUT
CALL TPM_FFTW_EXECUTE_DFT_R2C(PLAN,IN_FLAT,OUT_FLAT)
END SUBROUTINE TPM_FFTW_EXECUTE_DFT_R2C_RANK2


SUBROUTINE INIT_PLANS_FFTW(KDLON)
!> @brief Allocate plan bookkeeping for the active resolution.
!!
!! @param[in] KDLON Upper bound on the transform lengths indexed in the plan cache.
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: KDLON

INTEGER(KIND=JPIM) :: INTH
CHARACTER(LEN=4096) :: ZBUF
INTEGER :: IELEN, IESTAT

TW%N_MAX=KDLON
ALLOCATE(TW%FFTW_PLANS(TW%N_MAX))
ALLOCATE(TW%N_PLANS(TW%N_MAX))
TW%N_PLANS(:)=0

! Pre-size the per-thread persistent FFT scratch registry (first call only;
! the slot arrays persist across resolution destroy/re-create cycles).
IF (.NOT. ALLOCATED(ZFFT_PERSIST_PTRS)) THEN
  INTH = OMP_GET_MAX_THREADS()
  IF (INTH < 1) INTH = 1
  ALLOCATE(ZFFT_PERSIST_PTRS(INTH))
  ALLOCATE(ZFFT_PERSIST_CAPS(INTH))
  ZFFT_PERSIST_PTRS(:) = C_NULL_PTR
  ZFFT_PERSIST_CAPS(:) = 0_C_SIZE_T
ENDIF

! Decide whether to use FFTW "measure" mode with wisdom import/export.
CALL GET_ENVIRONMENT_VARIABLE('ECTRANS_FFTW_WISDOM', VALUE=ZBUF, &
  & LENGTH=IELEN, STATUS=IESTAT)
LUSE_MEASURE = (IESTAT == 0 .AND. IELEN > 0 .AND. IELEN <= LEN(ZBUF))

! Import previous FFTW wisdom if present; no-op if the env var is unset.
CALL IMPORT_FFTW_WISDOM

RETURN  
END SUBROUTINE INIT_PLANS_FFTW


SUBROUTINE CREATE_PLAN_FFTW(KPLAN,KTYPE,KN,KLOT)
!> @brief Reuse or create a cached FFTW plan for a given transform shape.
!!
!! @param[out] KPLAN FFTW plan handle matching the requested transform layout.
!! @param[in] KTYPE Transform direction selector: `1` for complex-to-real and `-1` for real-to-complex.
!! @param[in] KN Real transform length used as the one-dimensional FFT extent.
!! @param[in] KLOT Number of transforms to execute in the batched FFTW plan.
!! @note Access to the shared plan cache is serialized by the named OpenMP critical region `FFTW_CREATE`.
!!   This includes plan lookup, eviction, and creation, so heavy multi-threaded contention here can limit
!!   scalability unless the needed plans are created early and then reused.
IMPLICIT NONE
TYPE(C_PTR),INTENT(OUT) :: KPLAN
INTEGER(KIND=JPIM),INTENT(IN) :: KTYPE,KN,KLOT

TYPE(C_PTR) :: IPLAN
INTEGER(KIND=JPIM) :: IRANK, ISTRIDE
INTEGER(KIND=JPIM) :: JL
INTEGER(KIND=JPIM) :: IRDIST,ICDIST,IN(1),IEMBED(1)
INTEGER(KIND=JPIM) :: CEMBED(1)
INTEGER(KIND=JPIM) :: IFLAG
REAL(KIND=JPRB), POINTER :: ZDUM(:)
COMPLEX(KIND=JPCB), POINTER :: CDUM(:)
TYPE(C_PTR) :: ZDUMP
LOGICAL :: LLFOUND
LOGICAL, PARAMETER :: LLRESTRICT_PLANS=.TRUE.
TYPE(FFTW_PLAN),POINTER :: CURR_FFTW_PLAN,START_FFTW_PLAN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE, ZHOOK_HANDLE2
IF (LHOOK) CALL DR_HOOK('CREATE_PLAN_FFTW',0,ZHOOK_HANDLE)

IF( KN > TW%N_MAX )THEN
  CALL ABORT_TRANS('CREATE_PLAN_FFTW: KN > N_MAX THAT WAS INITIALISED IN INIT_PLANS_FFTW')
ENDIF

IRANK=1
ISTRIDE=1
IN(1)=KN
IEMBED(1)=IN(1)
ICDIST=KN/2+1
CEMBED(1)=ICDIST
IRDIST=ICDIST*2

!$OMP CRITICAL (FFTW_CREATE)
LLFOUND=.FALSE.
IF( TW%FFTW_PLANS(KN)%NPLAN_ID /= NPLAN_ID_UNINITIALISED )THEN
  WRITE(*,'("CREATE_PLAN_FFTW.1: PLAN_ID=",I10)')TW%FFTW_PLANS(KN)%NPLAN_ID
  CALL ABORT_TRANS('CREATE_PLAN_FFTW.1: NPLAN_ID /= NPLAN_ID_UNINITIALISED')
ENDIF
CURR_FFTW_PLAN=>TW%FFTW_PLANS(KN)
IF( CURR_FFTW_PLAN%NPLAN_ID /= NPLAN_ID_UNINITIALISED )THEN
  WRITE(*,'("CREATE_PLAN_FFTW.2: PLAN_ID=",I10)')CURR_FFTW_PLAN%NPLAN_ID
  CALL ABORT_TRANS('CREATE_PLAN_FFTW.2: NPLAN_ID /= NPLAN_ID_UNINITIALISED')
ENDIF
! search for plan in existing plans
DO JL=1,TW%N_PLANS(KN)
  IF( KLOT == CURR_FFTW_PLAN%NLOT .AND. KTYPE == CURR_FFTW_PLAN%NTYPE )THEN
    LLFOUND=.TRUE.
    IPLAN=CURR_FFTW_PLAN%NPLAN
    EXIT
  ELSEIF( JL /= TW%N_PLANS(KN) )THEN
    CURR_FFTW_PLAN=>CURR_FFTW_PLAN%NEXT_PLAN
    IF( CURR_FFTW_PLAN%NPLAN_ID /= NPLAN_ID_UNINITIALISED )THEN
      WRITE(*,'("CREATE_PLAN_FFTW.3: PLAN_ID=",I10)')CURR_FFTW_PLAN%NPLAN_ID
      CALL ABORT_TRANS('CREATE_PLAN_FFTW.3: NPLAN_ID /= NPLAN_ID_UNINITIALISED')
    ENDIF
  ENDIF
ENDDO
IF( .NOT.LLFOUND )THEN
  IF( LLRESTRICT_PLANS )THEN
    IF( TW%N_PLANS(KN) == TW%N_MAX_PLANS )THEN
      ! destroy the plan at the start of the list
      CALL TPM_FFTW_DESTROY_PLAN(TW%FFTW_PLANS(KN)%NPLAN)
      TW%FFTW_PLANS(KN)%NPLAN_ID=NPLAN_ID_DESTROYED
       ! mark the plan record as destroyed but keep it in the list to be recycled by the next plan creation with the same KN,
       ! which avoids costly deallocation and reallocation of the FFTW_PLAN record and preserves the linked list structure
       ! without needing to update any pointers
      START_FFTW_PLAN=>TW%FFTW_PLANS(KN)
      TW%FFTW_PLANS(KN)=TW%FFTW_PLANS(KN)%NEXT_PLAN
      ! DEALLOCATE(START_FFTW_PLAN)
      TW%N_PLANS(KN)=TW%N_PLANS(KN)-1
    ENDIF
  ENDIF
  ! Dummy buffer sized for the real plan (ICDIST*KLOT complex = IRDIST*KLOT reals,
  ! in-place). FFTW_ESTIMATE does not touch it, but FFTW_MEASURE (used when
  ! ECTRANS_FFTW_WISDOM is set) does, so the correct sizing must be in place.
  ZDUMP=TPM_FFTW_ALLOC_COMPLEX(INT(MAX(ICDIST*KLOT,1),C_SIZE_T))
  CALL C_F_POINTER(ZDUMP,ZDUM,[MAX(IRDIST*KLOT,1)])
  CALL C_F_POINTER(ZDUMP,CDUM,[MAX(ICDIST*KLOT,1)])
  ! FFTW_MEASURE (changes numerical results vs baseline) when ECTRANS_FFTW_WISDOM
  ! is set, otherwise the default ESTIMATE+NO_SIMD path (bit-identical).
  IF (LUSE_MEASURE) THEN
    IFLAG=FFTW_MEASURE
  ELSE
    IFLAG=TPM_FFTW_PLAN_FLAGS
  ENDIF
  IF( KTYPE==NTYPE_C2R )THEN
    IF (LHOOK) CALL DR_HOOK('FFTW_PLAN_MANY_DFT_C2R',0,ZHOOK_HANDLE2)
    IPLAN=TPM_FFTW_PLAN_MANY_DFT_C2R(IRANK,IN,KLOT,CDUM,CEMBED,ISTRIDE,ICDIST,&
         & ZDUM,IEMBED,ISTRIDE,IRDIST,IFLAG)
    IF (LHOOK) CALL DR_HOOK('FFTW_PLAN_MANY_DFT_C2R',1,ZHOOK_HANDLE2)
  ELSEIF( KTYPE==NTYPE_R2C )THEN
    IF (LHOOK) CALL DR_HOOK('FFTW_PLAN_MANY_DFT_R2C',0,ZHOOK_HANDLE2)
    IPLAN=TPM_FFTW_PLAN_MANY_DFT_R2C(IRANK,IN,KLOT,ZDUM,IEMBED,ISTRIDE,IRDIST,&
         & CDUM,CEMBED,ISTRIDE,ICDIST,IFLAG)
    IF (LHOOK) CALL DR_HOOK('FFTW_PLAN_MANY_DFT_R2C',1,ZHOOK_HANDLE2)
  ELSE
    CALL ABORT_TRANS('FFTW_PLAN: INVALID KTYPE')
  ENDIF
  CALL TPM_FFTW_FREE(ZDUMP)
  KPLAN=IPLAN
  TW%N_PLANS(KN)=TW%N_PLANS(KN)+1
  IF( TW%N_PLANS(KN) /= 1 )THEN
    ALLOCATE(CURR_FFTW_PLAN%NEXT_PLAN)
    CURR_FFTW_PLAN=>CURR_FFTW_PLAN%NEXT_PLAN
  ENDIF
  IF( CURR_FFTW_PLAN%NPLAN_ID /= NPLAN_ID_UNINITIALISED )THEN
    WRITE(*,'("CREATE_PLAN_FFTW.4: PLAN_ID=",I10)')CURR_FFTW_PLAN%NPLAN_ID
    CALL ABORT_TRANS('CREATE_PLAN_FFTW.4: NPLAN_ID /= NPLAN_ID_UNINITIALISED')
  ENDIF
  CURR_FFTW_PLAN%NPLAN=IPLAN
  CURR_FFTW_PLAN%NLOT=KLOT
  CURR_FFTW_PLAN%NTYPE=KTYPE
  CURR_FFTW_PLAN%NEXT_PLAN=>NULL()
ELSE
  KPLAN=IPLAN
ENDIF
!$OMP END CRITICAL (FFTW_CREATE)

IF (LHOOK) CALL DR_HOOK('CREATE_PLAN_FFTW',1,ZHOOK_HANDLE)
RETURN
END SUBROUTINE CREATE_PLAN_FFTW


SUBROUTINE DESTROY_PLAN_FFTW(KPLAN)
!> @brief Destroy one cached FFTW plan inside the module-level OpenMP critical region.
!!
!! @param[in] KPLAN FFTW plan handle to destroy.
!! @note Plan destruction is serialized by the named OpenMP critical region `FFTW_DESTROY`.
!!   This is usually a teardown path, but repeated destruction from many threads will not run concurrently.
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: KPLAN
!$OMP CRITICAL (FFTW_DESTROY)
CALL TPM_FFTW_DESTROY_PLAN(KPLAN)
!$OMP END CRITICAL (FFTW_DESTROY)
RETURN
END SUBROUTINE DESTROY_PLAN_FFTW


SUBROUTINE DESTROY_PLANS_FFTW
!> @brief Destroy all cached FFTW plans and release the resolution-local plan tables.
!!
!! @note Plans are destroyed one by one through `DESTROY_PLAN_FFTW`, so teardown repeatedly enters the
!!   `FFTW_DESTROY` critical region and is serialized across threads.
IMPLICIT NONE
INTEGER(KIND=JPIM) :: JL, JN
INTEGER(KIND=JPIM) :: JT
TYPE(FFTW_PLAN),POINTER :: CURR_FFTW_PLAN, NEXT_FFTW_PLAN
DO JN=1,TW%N_MAX
  CURR_FFTW_PLAN=>TW%FFTW_PLANS(JN)
  DO JL=1,TW%N_PLANS(JN)
    CALL DESTROY_PLAN_FFTW(CURR_FFTW_PLAN%NPLAN)
    NEXT_FFTW_PLAN=>CURR_FFTW_PLAN%NEXT_PLAN
    IF( JL /= 1 ) THEN
      DEALLOCATE( CURR_FFTW_PLAN )
    ENDIF
    CURR_FFTW_PLAN => NEXT_FFTW_PLAN
  ENDDO
ENDDO
IF( ASSOCIATED(TW) ) THEN
  IF( ASSOCIATED(TW%FFTW_PLANS) ) DEALLOCATE(TW%FFTW_PLANS)
  IF( ALLOCATED(TW%N_PLANS) )     DEALLOCATE(TW%N_PLANS)
  TW%N_MAX=0
ENDIF

! Free per-thread persistent FFT scratch slots. The slot arrays themselves
! remain allocated so a later INIT_PLANS_FFTW + EXEC_FFTW_IMPL cycle
! (e.g. another resolution) can reuse the registry.
IF (ALLOCATED(ZFFT_PERSIST_PTRS)) THEN
  DO JT=1,SIZE(ZFFT_PERSIST_PTRS)
    IF (C_ASSOCIATED(ZFFT_PERSIST_PTRS(JT))) THEN
      CALL TPM_FFTW_FREE(ZFFT_PERSIST_PTRS(JT))
      ZFFT_PERSIST_PTRS(JT) = C_NULL_PTR
    ENDIF
    ZFFT_PERSIST_CAPS(JT) = 0_C_SIZE_T
  ENDDO
ENDIF

! Export unified wisdom if ECTRANS_FFTW_WISDOM is set; no-op otherwise.
CALL EXPORT_FFTW_WISDOM

RETURN
END SUBROUTINE DESTROY_PLANS_FFTW


SUBROUTINE COPY_PREEL_TO_ZFFT(KLEN,KOFF,LD_TRANSPOSED,PREEL,ZFFT)
!> @brief Copy a batched `PREEL` slab into the rank-2 FFT work array.
!!
!! @param[in] KLEN Number of values copied per field.
!! @param[in] KOFF One-based offset of the transform segment inside `PREEL`.
!! @param[in] LD_TRANSPOSED Selects `PREEL(point,field)` when true and `PREEL(field,point)` when false.
!! @param[in] PREEL Source work array.
!! @param[inout] ZFFT Destination rank-2 FFT work buffer.
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: KLEN
INTEGER(KIND=JPIM),INTENT(IN) :: KOFF
LOGICAL,INTENT(IN) :: LD_TRANSPOSED
REAL(KIND=JPRB),INTENT(IN) :: PREEL(:,:)
REAL(KIND=JPRB),INTENT(INOUT) :: ZFFT(:,:)

INTEGER(KIND=JPIM) :: JJ,JF

IF( LD_TRANSPOSED )THEN
  DO JF=1,SIZE(ZFFT,2)
    DO JJ=1,KLEN
      ZFFT(JJ,JF)=PREEL(KOFF+JJ-1,JF)
    ENDDO
  ENDDO
ELSE
  DO JF=1,SIZE(ZFFT,2)
    DO JJ=1,KLEN
      ZFFT(JJ,JF)=PREEL(JF,KOFF+JJ-1)
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE COPY_PREEL_TO_ZFFT

SUBROUTINE COPY_ZFFT_TO_PREEL(KLEN,KOFF,LD_TRANSPOSED,ZFFT,PREEL,PSCALE)
!> @brief Copy a rank-2 FFT work array back into `PREEL`, optionally applying a scale factor.
!!
!! @param[in] KLEN Number of values copied per field.
!! @param[in] KOFF One-based offset of the transform segment inside `PREEL`.
!! @param[in] LD_TRANSPOSED Selects `PREEL(point,field)` when true and `PREEL(field,point)` when false.
!! @param[in] ZFFT Source rank-2 FFT work buffer.
!! @param[inout] PREEL Destination work array.
!! @param[in] PSCALE Optional scale factor applied to each copied value when present.
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: KLEN
INTEGER(KIND=JPIM),INTENT(IN) :: KOFF
LOGICAL,INTENT(IN) :: LD_TRANSPOSED
REAL(KIND=JPRB),INTENT(IN) :: ZFFT(:,:)
REAL(KIND=JPRB),INTENT(INOUT) :: PREEL(:,:)
REAL(KIND=JPRB), OPTIONAL,INTENT(IN) :: PSCALE

INTEGER(KIND=JPIM) :: JJ,JF

IF( LD_TRANSPOSED )THEN
  IF (PRESENT(PSCALE)) THEN
     DO JF=1,SIZE(ZFFT,2)
       DO JJ=1,KLEN
         PREEL(KOFF+JJ-1,JF)=ZFFT(JJ,JF)*PSCALE
       ENDDO
     ENDDO
  ELSE
    DO JF=1,SIZE(ZFFT,2)
      DO JJ=1,KLEN
        PREEL(KOFF+JJ-1,JF)=ZFFT(JJ,JF)
      ENDDO
    ENDDO
  ENDIF
ELSE
  IF (PRESENT(PSCALE)) THEN
     DO JF=1,SIZE(ZFFT,2)
       DO JJ=1,KLEN
         PREEL(JF,KOFF+JJ-1)=ZFFT(JJ,JF)*PSCALE
       ENDDO
     ENDDO
  ELSE
    DO JF=1,SIZE(ZFFT,2)
      DO JJ=1,KLEN
        PREEL(JF,KOFF+JJ-1)=ZFFT(JJ,JF)
      ENDDO
    ENDDO
  ENDIF
ENDIF
END SUBROUTINE COPY_ZFFT_TO_PREEL


SUBROUTINE COPY_PREEL_JF_TO_ZFFT_1(KLEN,KOFF,LD_TRANSPOSED,PREEL,KFIELD,ZFFT1)
!> @brief Copy one field from `PREEL` into the rank-1 FFT work array.
!!
!! @param[in] KLEN Number of values copied for the selected field.
!! @param[in] KOFF One-based offset of the transform segment inside `PREEL`.
!! @param[in] LD_TRANSPOSED Selects `PREEL(point,field)` when true and `PREEL(field,point)` when false.
!! @param[in] PREEL Source work array.
!! @param[in] KFIELD Field index copied into the work array.
!! @param[inout] ZFFT1 Destination rank-1 FFT work buffer.
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: KLEN
INTEGER(KIND=JPIM),INTENT(IN) :: KOFF
LOGICAL,INTENT(IN) :: LD_TRANSPOSED
REAL(KIND=JPRB),INTENT(IN) :: PREEL(:,:)
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELD
REAL(KIND=JPRB),INTENT(INOUT) :: ZFFT1(:)

INTEGER(KIND=JPIM) :: JJ

IF( LD_TRANSPOSED )THEN
  DO JJ=1,KLEN
    ZFFT1(JJ)=PREEL(KOFF+JJ-1,KFIELD)
  ENDDO
ELSE
  DO JJ=1,KLEN
    ZFFT1(JJ)=PREEL(KFIELD,KOFF+JJ-1)
  ENDDO
ENDIF
END SUBROUTINE COPY_PREEL_JF_TO_ZFFT_1


SUBROUTINE COPY_ZFFT_1_TO_PREEL_JF(KLEN,KOFF,LD_TRANSPOSED,ZFFT1,PREEL,KFIELD,PSCALE)
!> @brief Copy one FFT work vector back into `PREEL`, optionally applying a scale factor.
!!
!! @param[in] KLEN Number of values copied for the selected field.
!! @param[in] KOFF One-based offset of the transform segment inside `PREEL`.
!! @param[in] LD_TRANSPOSED Selects `PREEL(point,field)` when true and `PREEL(field,point)` when false.
!! @param[in] ZFFT1 Source rank-1 FFT work buffer.
!! @param[inout] PREEL Destination work array.
!! @param[in] KFIELD Field index updated in `PREEL`.
!! @param[in] PSCALE Optional scale factor applied to each copied value when present.
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: KLEN
INTEGER(KIND=JPIM),INTENT(IN) :: KOFF
LOGICAL,INTENT(IN) :: LD_TRANSPOSED
REAL(KIND=JPRB),INTENT(IN) :: ZFFT1(:)
REAL(KIND=JPRB),INTENT(INOUT) :: PREEL(:,:)
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELD
REAL(KIND=JPRB), OPTIONAL,INTENT(IN) :: PSCALE

INTEGER(KIND=JPIM) :: JJ

IF( LD_TRANSPOSED )THEN
  IF (PRESENT(PSCALE)) THEN
    DO JJ=1,KLEN
      PREEL(KOFF+JJ-1,KFIELD)=ZFFT1(JJ)*PSCALE
    ENDDO
  ELSE
    DO JJ=1,KLEN
      PREEL(KOFF+JJ-1,KFIELD)=ZFFT1(JJ)
    ENDDO
  ENDIF
ELSE
  IF (PRESENT(PSCALE)) THEN
    DO JJ=1,KLEN
      PREEL(KFIELD,KOFF+JJ-1)=ZFFT1(JJ)*PSCALE
    ENDDO
  ELSE
    DO JJ=1,KLEN
      PREEL(KFIELD,KOFF+JJ-1)=ZFFT1(JJ)
    ENDDO
  ENDIF
ENDIF
END SUBROUTINE COPY_ZFFT_1_TO_PREEL_JF


SUBROUTINE EXEC_FFTW_IMPL(CDNAME,KTYPE,KRLEN,KCLEN,KOFF,KFIELDS,LD_ALL,PREEL,LD_TRANSPOSED)
!> @brief Execute the shared FFTW path for both PREEL memory layouts and batching modes.
!!
!! @param[in] CDNAME Name passed to `DR_HOOK` and used in runtime error messages.
!! @param[in] KTYPE Transform direction selector: `1` for complex-to-real and `-1` for real-to-complex.
!! @param[in] KRLEN Number of real values in each transform.
!! @param[in] KCLEN Storage length of the packed FFT segment in real-valued form.
!! @param[in] KOFF One-based offset of the transform segment inside `PREEL`.
!! @param[in] KFIELDS Number of fields to transform.
!! @param[in] LD_ALL If true, execute all fields in one batched FFTW call; otherwise transform one field at a time.
!! @param[inout] PREEL Work array holding the transform input on entry and the transformed output on return.
!! @param[in] LD_TRANSPOSED Selects the `PREEL(point,field)` layout when true and `PREEL(field,point)` when false.
!! @note The internal call to `CREATE_PLAN_FFTW` accesses the shared FFTW plan cache inside the `FFTW_CREATE`
!!   critical region. Reusing warmed-up plans avoids repeated planning work, but cache access itself remains
!!   serialized across threads.
!! @note Scratch comes from the per-thread persistent buffers (`ENSURE_FFT_BUFFER`); no allocation
!!   happens here and reuse changes no values since the buffer is fully overwritten before any read.
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN) :: CDNAME
INTEGER(KIND=JPIM),INTENT(IN) :: KTYPE
INTEGER(KIND=JPIM),INTENT(IN) :: KRLEN
INTEGER(KIND=JPIM),INTENT(IN) :: KCLEN
INTEGER(KIND=JPIM),INTENT(IN) :: KOFF
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
LOGICAL,INTENT(IN) :: LD_ALL
REAL(KIND=JPRB),INTENT(INOUT) :: PREEL(:,:)
LOGICAL,INTENT(IN) :: LD_TRANSPOSED

REAL(KIND=JPRB), POINTER :: ZFFT(:,:)
COMPLEX(KIND=JPCB), POINTER :: CFFT(:,:)
TYPE(C_PTR) :: ZFFTP
TYPE(C_PTR) :: IPLAN
INTEGER(KIND=JPIM) :: JF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE, ZHOOK_HANDLE2
INTEGER(KIND=JPIM) :: NBATCH
REAL(KIND=JPRB) :: ZSCALE

IF (LHOOK) CALL DR_HOOK(CDNAME,0,ZHOOK_HANDLE)

IF ( (KTYPE /= NTYPE_R2C) .AND. (KTYPE /= NTYPE_C2R) ) THEN
  CALL ABORT_TRANS('TPM_FFTW:'//TRIM(CDNAME)//' : WRONG VALUE KTYPE')
ENDIF

NBATCH = MERGE(KFIELDS,1,LD_ALL)
CALL CREATE_PLAN_FFTW(IPLAN, KTYPE, KRLEN, NBATCH)

! Per-thread persistent work array large enough to hold the complex FFT data for the
! entire batch of fields. The size is KCLEN/2 complex numbers per field, and there
! are NBATCH fields in the batch. This size is also sufficient for the real data.
! The buffer grows on demand in ENSURE_FFT_BUFFER and its slots are freed in
! DESTROY_PLANS_FFTW, so no per-call allocation/free happens here. Each OMP thread
! owns its slot and the buffer is fully overwritten below before any read, so reuse
! changes no values.
CALL ENSURE_FFT_BUFFER(INT(KCLEN/2*NBATCH,C_SIZE_T), ZFFTP)

! It is chosen to alias the complex and real views of the work array through C_F_POINTER,
! so that FFTW performs in-place transforms.
! It should be investigated if out-of-place transforms with separate work arrays for real and complex data can improve performance
! by allowing more flexible memory access patterns in FFTW.
CALL C_F_POINTER(ZFFTP,ZFFT,[KCLEN,  NBATCH])
CALL C_F_POINTER(ZFFTP,CFFT,[KCLEN/2,NBATCH])

IF( LD_ALL ) THEN
  ! All fields are transformed together in a single FFTW call, so the work array is laid out as ZFFT(transform_point,field).
  IF (KTYPE==NTYPE_C2R) THEN
    CALL COPY_PREEL_TO_ZFFT(KCLEN,KOFF,LD_TRANSPOSED,PREEL,ZFFT)
    IF (LHOOK) CALL DR_HOOK('FFTW_EXECUTE_DFT_C2R',0,ZHOOK_HANDLE2)
    CALL TPM_FFTW_EXECUTE_DFT_C2R(IPLAN,CFFT,ZFFT)
    IF (LHOOK) CALL DR_HOOK('FFTW_EXECUTE_DFT_C2R',1,ZHOOK_HANDLE2)
    CALL COPY_ZFFT_TO_PREEL(KRLEN,KOFF,LD_TRANSPOSED,ZFFT,PREEL)
  ELSE
    CALL COPY_PREEL_TO_ZFFT(KRLEN,KOFF,LD_TRANSPOSED,PREEL,ZFFT)
    IF (LHOOK) CALL DR_HOOK('FFTW_EXECUTE_DFT_R2C',0,ZHOOK_HANDLE2)
    CALL TPM_FFTW_EXECUTE_DFT_R2C(IPLAN,ZFFT,CFFT)
    IF (LHOOK) CALL DR_HOOK('FFTW_EXECUTE_DFT_R2C',1,ZHOOK_HANDLE2)
    ! Real-to-complex transforms require scaling by 1/KRLEN, which can be applied in the copy back to PREEL to save cycles in the FFT work array.
    ZSCALE = 1.0_JPRB/REAL(KRLEN,JPRB)
    CALL COPY_ZFFT_TO_PREEL(KCLEN,KOFF,LD_TRANSPOSED,ZFFT,PREEL,PSCALE=ZSCALE)
  ENDIF
ELSE
  ! All fields are transformed separately in a loop over `JF`, so the work array is laid out as ZFFT(transform_point,1) and only one field is copied in and out at a time.
  IF (KTYPE==NTYPE_C2R) THEN
    DO JF=1,KFIELDS
      CALL COPY_PREEL_JF_TO_ZFFT_1(KCLEN,KOFF,LD_TRANSPOSED,PREEL,JF,ZFFT(:,1))
      IF (LHOOK) CALL DR_HOOK('FFTW_EXECUTE_DFT_C2R',0,ZHOOK_HANDLE2)
      CALL TPM_FFTW_EXECUTE_DFT_C2R(IPLAN,CFFT(:,1),ZFFT(:,1))
      IF (LHOOK) CALL DR_HOOK('FFTW_EXECUTE_DFT_C2R',1,ZHOOK_HANDLE2)
      CALL COPY_ZFFT_1_TO_PREEL_JF(KRLEN,KOFF,LD_TRANSPOSED,ZFFT(:,1),PREEL,JF)
    ENDDO
  ELSE
    ! Real-to-complex transforms require scaling by 1/KRLEN, which can be applied in the copy back to PREEL to save cycles in the FFT work array.
    ZSCALE = 1.0_JPRB/REAL(KRLEN,JPRB)
    DO JF=1,KFIELDS
      CALL COPY_PREEL_JF_TO_ZFFT_1(KRLEN,KOFF,LD_TRANSPOSED,PREEL,JF,ZFFT(:,1))
      IF (LHOOK) CALL DR_HOOK('FFTW_EXECUTE_DFT_R2C',0,ZHOOK_HANDLE2)
      CALL TPM_FFTW_EXECUTE_DFT_R2C(IPLAN,ZFFT(:,1),CFFT(:,1))
      IF (LHOOK) CALL DR_HOOK('FFTW_EXECUTE_DFT_R2C',1,ZHOOK_HANDLE2)
      CALL COPY_ZFFT_1_TO_PREEL_JF(KCLEN,KOFF,LD_TRANSPOSED,ZFFT(:,1),PREEL,JF,PSCALE=ZSCALE)
    ENDDO
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK(CDNAME,1,ZHOOK_HANDLE)
END SUBROUTINE EXEC_FFTW_IMPL


SUBROUTINE EXEC_FFTW(KTYPE,KRLEN,KCLEN,KOFF,KFIELDS,LD_ALL,PREEL)
!> @brief Execute FFTW transforms for arrays stored as `PREEL(field,point)`.
!!
!! @param[in] KTYPE Transform direction selector: `1` for complex-to-real and `-1` for real-to-complex.
!! @param[in] KRLEN Number of real values in each transform.
!! @param[in] KCLEN Storage length of the packed FFT segment in real-valued form.
!! @param[in] KOFF One-based offset of the transform segment inside `PREEL`.
!! @param[in] KFIELDS Number of fields to transform.
!! @param[in] LD_ALL If true, execute all fields in one batched FFTW call; otherwise transform one field at a time.
!! @param[inout] PREEL Work array laid out as `PREEL(field,point)` containing the input on entry and the output on return.
!! @note This routine shares the module-wide FFTW plan cache, so calls may serialize briefly in `CREATE_PLAN_FFTW`
!!   when accessing the `FFTW_CREATE` critical region.
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)   :: KTYPE
INTEGER(KIND=JPIM),INTENT(IN)   :: KRLEN
INTEGER(KIND=JPIM),INTENT(IN)   :: KCLEN
INTEGER(KIND=JPIM),INTENT(IN)   :: KOFF
INTEGER(KIND=JPIM),INTENT(IN)   :: KFIELDS
LOGICAL           ,INTENT(IN)   :: LD_ALL
REAL(KIND=JPRB), INTENT(INOUT)  :: PREEL(:,:)

CALL EXEC_FFTW_IMPL('EXEC_FFTW',KTYPE,KRLEN,KCLEN,KOFF,KFIELDS,LD_ALL,PREEL,LD_TRANSPOSED=.FALSE.)
END SUBROUTINE EXEC_FFTW

SUBROUTINE EXEC_EFFTW(KTYPE,KRLEN,KCLEN,KOFF,KFIELDS,LD_ALL,PREEL)
!> @brief Execute FFTW transforms for arrays stored as `PREEL(point,field)`.
!!
!! @param[in] KTYPE Transform direction selector: `1` for complex-to-real and `-1` for real-to-complex.
!! @param[in] KRLEN Number of real values in each transform.
!! @param[in] KCLEN Storage length of the packed FFT segment in real-valued form.
!! @param[in] KOFF One-based offset of the transform segment inside `PREEL`.
!! @param[in] KFIELDS Number of fields to transform.
!! @param[in] LD_ALL If true, execute all fields in one batched FFTW call; otherwise transform one field at a time.
!! @param[inout] PREEL Work array laid out as `PREEL(point,field)` containing the input on entry and the output on return.
!! @note This routine shares the module-wide FFTW plan cache, so calls may serialize briefly in `CREATE_PLAN_FFTW`
!!   when accessing the `FFTW_CREATE` critical region. Reusing warmed-up plans avoids repeated planning work,
!!   but cache access itself remains serialized across threads.
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)   :: KTYPE
INTEGER(KIND=JPIM),INTENT(IN)   :: KRLEN
INTEGER(KIND=JPIM),INTENT(IN)   :: KCLEN
INTEGER(KIND=JPIM),INTENT(IN)   :: KOFF
INTEGER(KIND=JPIM),INTENT(IN)   :: KFIELDS
LOGICAL           ,INTENT(IN)   :: LD_ALL
REAL(KIND=JPRB), INTENT(INOUT)  :: PREEL(:,:)

CALL EXEC_FFTW_IMPL('EXEC_EFFTW',KTYPE,KRLEN,KCLEN,KOFF,KFIELDS,LD_ALL,PREEL,LD_TRANSPOSED=.TRUE.)
END SUBROUTINE EXEC_EFFTW

! -----------------------------------------------------------------------------
! Following routines are for legacy FFTW interface support and can be removed
! once the legacy interface is no longer needed.
! -----------------------------------------------------------------------------

#if LEGACY_FFTW_INTERFACE
FUNCTION TPM_FFTW_CPTR_FROM_LEGACY_PLAN(PLAN_LEGACY) RESULT(PLAN)
IMPLICIT NONE
INTEGER(KIND=C_INTPTR_T),INTENT(IN) :: PLAN_LEGACY
TYPE(C_PTR) :: PLAN
PLAN=TRANSFER(PLAN_LEGACY,PLAN)
END FUNCTION TPM_FFTW_CPTR_FROM_LEGACY_PLAN


FUNCTION TPM_FFTW_LEGACY_PLAN_FROM_CPTR(PLAN) RESULT(PLAN_LEGACY)
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: PLAN
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY
PLAN_LEGACY=TRANSFER(PLAN,PLAN_LEGACY)
END FUNCTION TPM_FFTW_LEGACY_PLAN_FROM_CPTR


FUNCTION TPM_FFTW_PLAN_MANY_DFT_C2R_LEGACY(RANK,N,HOWMANY,IN,INEMBED,ISTRIDE,IDIST,OUT,ONEMBED,OSTRIDE,ODIST,FLAGS) RESULT(PLAN)
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: RANK, HOWMANY
INTEGER(KIND=JPIM),INTENT(IN) :: ISTRIDE, IDIST, OSTRIDE, ODIST, FLAGS
INTEGER(KIND=JPIM),INTENT(IN) :: N(*), INEMBED(*), ONEMBED(*)
COMPLEX(KIND=C_DOUBLE_COMPLEX),INTENT(INOUT) :: IN(*)
REAL(KIND=C_DOUBLE),INTENT(INOUT) :: OUT(*)
TYPE(C_PTR) :: PLAN
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY

CALL DFFTW_PLAN_MANY_DFT_C2R(PLAN_LEGACY,RANK,N,HOWMANY,IN,INEMBED,ISTRIDE,IDIST,&
  & OUT,ONEMBED,OSTRIDE,ODIST,FLAGS)
PLAN=TPM_FFTW_CPTR_FROM_LEGACY_PLAN(PLAN_LEGACY)
END FUNCTION TPM_FFTW_PLAN_MANY_DFT_C2R_LEGACY


FUNCTION TPM_FFTWF_PLAN_MANY_DFT_C2R_LEGACY(RANK,N,HOWMANY,IN,INEMBED,ISTRIDE,IDIST,OUT,ONEMBED,OSTRIDE,ODIST,FLAGS) RESULT(PLAN)
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: RANK, HOWMANY
INTEGER(KIND=JPIM),INTENT(IN) :: ISTRIDE, IDIST, OSTRIDE, ODIST, FLAGS
INTEGER(KIND=JPIM),INTENT(IN) :: N(*), INEMBED(*), ONEMBED(*)
COMPLEX(KIND=C_FLOAT_COMPLEX),INTENT(INOUT) :: IN(*)
REAL(KIND=C_FLOAT),INTENT(INOUT) :: OUT(*)
TYPE(C_PTR) :: PLAN
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY

CALL SFFTW_PLAN_MANY_DFT_C2R(PLAN_LEGACY,RANK,N,HOWMANY,IN,INEMBED,ISTRIDE,IDIST,&
  & OUT,ONEMBED,OSTRIDE,ODIST,FLAGS)
PLAN=TPM_FFTW_CPTR_FROM_LEGACY_PLAN(PLAN_LEGACY)
END FUNCTION TPM_FFTWF_PLAN_MANY_DFT_C2R_LEGACY


FUNCTION TPM_FFTW_PLAN_MANY_DFT_R2C_LEGACY(RANK,N,HOWMANY,IN,INEMBED,ISTRIDE,IDIST,OUT,ONEMBED,OSTRIDE,ODIST,FLAGS) RESULT(PLAN)
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: RANK, HOWMANY
INTEGER(KIND=JPIM),INTENT(IN) :: ISTRIDE, IDIST, OSTRIDE, ODIST, FLAGS
INTEGER(KIND=JPIM),INTENT(IN) :: N(*), INEMBED(*), ONEMBED(*)
REAL(KIND=C_DOUBLE),INTENT(INOUT) :: IN(*)
COMPLEX(KIND=C_DOUBLE_COMPLEX),INTENT(INOUT) :: OUT(*)
TYPE(C_PTR) :: PLAN
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY

CALL DFFTW_PLAN_MANY_DFT_R2C(PLAN_LEGACY,RANK,N,HOWMANY,IN,INEMBED,ISTRIDE,IDIST,&
  & OUT,ONEMBED,OSTRIDE,ODIST,FLAGS)
PLAN=TPM_FFTW_CPTR_FROM_LEGACY_PLAN(PLAN_LEGACY)
END FUNCTION TPM_FFTW_PLAN_MANY_DFT_R2C_LEGACY


FUNCTION TPM_FFTWF_PLAN_MANY_DFT_R2C_LEGACY(RANK,N,HOWMANY,IN,INEMBED,ISTRIDE,IDIST,OUT,ONEMBED,OSTRIDE,ODIST,FLAGS) RESULT(PLAN)
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: RANK, HOWMANY
INTEGER(KIND=JPIM),INTENT(IN) :: ISTRIDE, IDIST, OSTRIDE, ODIST, FLAGS
INTEGER(KIND=JPIM),INTENT(IN) :: N(*), INEMBED(*), ONEMBED(*)
REAL(KIND=C_FLOAT),INTENT(INOUT) :: IN(*)
COMPLEX(KIND=C_FLOAT_COMPLEX),INTENT(INOUT) :: OUT(*)
TYPE(C_PTR) :: PLAN
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY

CALL SFFTW_PLAN_MANY_DFT_R2C(PLAN_LEGACY,RANK,N,HOWMANY,IN,INEMBED,ISTRIDE,IDIST,&
  & OUT,ONEMBED,OSTRIDE,ODIST,FLAGS)
PLAN=TPM_FFTW_CPTR_FROM_LEGACY_PLAN(PLAN_LEGACY)
END FUNCTION TPM_FFTWF_PLAN_MANY_DFT_R2C_LEGACY


SUBROUTINE TPM_FFTW_EXECUTE_DFT_C2R_LEGACY(PLAN,IN,OUT)
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: PLAN
COMPLEX(KIND=C_DOUBLE_COMPLEX),INTENT(INOUT) :: IN(*)
REAL(KIND=C_DOUBLE),INTENT(INOUT) :: OUT(*)
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY

PLAN_LEGACY=TPM_FFTW_LEGACY_PLAN_FROM_CPTR(PLAN)
CALL DFFTW_EXECUTE_DFT_C2R(PLAN_LEGACY,IN,OUT)
END SUBROUTINE TPM_FFTW_EXECUTE_DFT_C2R_LEGACY


SUBROUTINE TPM_FFTWF_EXECUTE_DFT_C2R_LEGACY(PLAN,IN,OUT)
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: PLAN
COMPLEX(KIND=C_FLOAT_COMPLEX),INTENT(INOUT) :: IN(*)
REAL(KIND=C_FLOAT),INTENT(INOUT) :: OUT(*)
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY

PLAN_LEGACY=TPM_FFTW_LEGACY_PLAN_FROM_CPTR(PLAN)
CALL SFFTW_EXECUTE_DFT_C2R(PLAN_LEGACY,IN,OUT)
END SUBROUTINE TPM_FFTWF_EXECUTE_DFT_C2R_LEGACY


SUBROUTINE TPM_FFTW_EXECUTE_DFT_R2C_LEGACY(PLAN,IN,OUT)
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: PLAN
REAL(KIND=C_DOUBLE),INTENT(INOUT) :: IN(*)
COMPLEX(KIND=C_DOUBLE_COMPLEX),INTENT(INOUT) :: OUT(*)
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY

PLAN_LEGACY=TPM_FFTW_LEGACY_PLAN_FROM_CPTR(PLAN)
CALL DFFTW_EXECUTE_DFT_R2C(PLAN_LEGACY,IN,OUT)
END SUBROUTINE TPM_FFTW_EXECUTE_DFT_R2C_LEGACY


SUBROUTINE TPM_FFTWF_EXECUTE_DFT_R2C_LEGACY(PLAN,IN,OUT)
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: PLAN
REAL(KIND=C_FLOAT),INTENT(INOUT) :: IN(*)
COMPLEX(KIND=C_FLOAT_COMPLEX),INTENT(INOUT) :: OUT(*)
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY

PLAN_LEGACY=TPM_FFTW_LEGACY_PLAN_FROM_CPTR(PLAN)
CALL SFFTW_EXECUTE_DFT_R2C(PLAN_LEGACY,IN,OUT)
END SUBROUTINE TPM_FFTWF_EXECUTE_DFT_R2C_LEGACY
#endif

END MODULE TPM_FFTW
