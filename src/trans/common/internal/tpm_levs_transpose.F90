! (C) Copyright 2026- ECMWF.
! (C) Copyright 2026- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TPM_LEVS_TRANSPOSE

! Helper subroutines for the "collar==waveset" (Option A) LEVS-communicator variant of
! the TRGTOL and TRLTOG transposes.
!
! When the grid-point NS collars are aligned 1:1 with the Fourier wave-sets
! (N_REGIONS_NS==NPRTRW and every N_REGIONS(:)==NPRTRV, LEQ_REGIONS on), every
! rank has MYSETW==NS and MYSETV==EW. The 32-rank sub-communicator
! MPL_ALL_LEVS_COMM (fixed wave-set W, varying V-set) then coincides exactly
! with the set of grid-point EW-sectors of one collar. All non-zero G<->L
! transpose traffic is confined to that sub-communicator, so the alltoallv can
! be issued on MPL_ALL_LEVS_COMM instead of the global world communicator --
! removing the NPROC-wide synchronisation and shrinking the counts/displs
! arrays from NPROC to NPRTRV.
!
! LEVS_TR_ALIGNED() is the runtime predicate the call sites AND with the
! LUSE_LEVS_TR opt-in flag: enabling the flag on a non-aligned decomposition
! is therefore a safe no-op (world-comm path is used).
!
! VERIFY_LEVS_ON_NODE() is the placement diagnostic: it checks that the
! NPRTRV ranks of MPL_ALL_LEVS_COMM physically share one node -- the condition
! under which the LEVS alltoallv is genuinely intra-node.

USE EC_PARKIND, ONLY : JPIM, JPIB, JPRD

IMPLICIT NONE

PRIVATE
PUBLIC :: LEVS_TR_ALIGNED, VERIFY_LEVS_ON_NODE
PUBLIC :: ACCOUNT_TR_BW, DUMP_TR_SKEW

LOGICAL, SAVE :: LALIGN_DONE = .FALSE.
LOGICAL, SAVE :: LALIGN_CACHE = .FALSE.

! --- G<->L and M<->L transpose bandwidth accounting (LREPORT_TR_BW) ---------
! Two independent counters: index 1  = TRLTOG (L->G, g157),
!                           index 2  = TRGTOL (G->L, g158),
!                           index 11 = TRMTOL (M->L, g807 alltoallv span),
!                           index 12 = TRLTOM (L->M, g806 alltoallv span).
! Skip the first NBW_WARMUP steady calls (plan/first-touch/allocation noise),
! then accumulate element-count+time over the next NBW_SAMPLE calls and report
! once. Payload is carried as an INTEGER(JPIB) count of JPRB elements (not
! real bytes) so the cross-rank reduction can use an integer SUM: MPL_ALLREDUCE 
! aborts on a real-mode SUM unless LDREPROD is supplied, whereas integer SUM
! is exact.  Indices 1/2 measure the whole transpose (pack+comm+unpack); 3/4 the
! non-blocking alltoallv issue only; 5/6 the MPL_WAIT span (comm completion not
! hidden by the deferred self-copy); 7/8 the deferred self-copy; 9/10 a pre-issue 
! barrier (load-imbalance skew feeding the alltoallv); 11/12 the blocking M<->L 
! alltoallv span (no pack/unpack: the Fourier buffers are pre-laid-out by LTINV/LTDIR);
! 13/14 a pre-issue barrier (LTINV/LTDIR load-imbalance skew feeding the M<->L alltoallv).
! WAIT vs SKEW separates pure transfer from pack-time imbalance; WAIT+SELF
! separate whether the region is comm-bound or copy-bound.
INTEGER(KIND=JPIM), PARAMETER :: NBW_IDX    = 14
INTEGER(KIND=JPIM), PARAMETER :: NBW_WARMUP = 4
INTEGER(KIND=JPIM), PARAMETER :: NBW_SAMPLE = 12

INTEGER(KIND=JPIM), SAVE :: NBW_CALLS(NBW_IDX) = 0
INTEGER(KIND=JPIB), SAVE :: IBW_BYTES(NBW_IDX) = 0_JPIB
REAL(KIND=JPRD),    SAVE :: ZBW_TIME(NBW_IDX)  = 0.0_JPRD

! Straggler-profiling (LDUMP_TR_SKEW): separate per-index call counter and a
! saved unit for the world-gathered per-rank barrier-wait dump (tr_skew.csv).
INTEGER(KIND=JPIM), SAVE :: NSK_CALLS(NBW_IDX) = 0
INTEGER(KIND=JPIM), SAVE :: ISK_UNIT           = -1
LOGICAL,            SAVE :: LSK_OPEN           = .FALSE.

CONTAINS

! ---------------------------------------------------------------------------

LOGICAL FUNCTION LEVS_TR_ALIGNED() RESULT(LLALIGNED)

! Return .TRUE. if the current decomposition aligns grid-point collars 1:1
! with Fourier wave-sets, so that MPL_ALL_LEVS_COMM equals the set of EW
! sectors of a collar. Result is cached: the decomposition is fixed after
! setup.

USE TPM_DISTR,      ONLY : LEQ_REGIONS, NPRTRW, NPRTRV
USE EQ_REGIONS_MOD, ONLY : N_REGIONS, N_REGIONS_NS

IF (LALIGN_DONE) THEN
  LLALIGNED = LALIGN_CACHE
  RETURN
ENDIF

LLALIGNED = .FALSE.
IF (LEQ_REGIONS) THEN
  IF (N_REGIONS_NS == NPRTRW) THEN
    IF (ALL(N_REGIONS(1:N_REGIONS_NS) == NPRTRV)) THEN
      LLALIGNED = .TRUE.
    ENDIF
  ENDIF
ENDIF

LALIGN_CACHE = LLALIGNED
LALIGN_DONE  = .TRUE.

END FUNCTION LEVS_TR_ALIGNED

! ---------------------------------------------------------------------------

SUBROUTINE VERIFY_LEVS_ON_NODE()

! Verify that the NPRTRV ranks of MPL_ALL_LEVS_COMM are co-located on a single physical node.
! Each rank hashes its hostname to an integer and does a MIN/MAX allreduce over
! MPL_ALL_LEVS_COMM; if MIN==MAX the whole sub-communicator shares one hostname.
! The LEVS-comm root (MYSETV==1) prints the result for its wave-set. Runs once if requested.

USE TPM_DISTR,      ONLY : MYPROC, MYSETW, MYSETV, NPRTRV
USE TPM_GEN,        ONLY : NOUT
USE MPL_GROUPS,     ONLY : MPL_ALL_LEVS_COMM
USE MPL_MODULE,     ONLY : MPL_ALLREDUCE

LOGICAL, SAVE :: LDONE = .FALSE.

CHARACTER(LEN=256) :: CLHOST
INTEGER(KIND=JPIM) :: ISTATUS, IHASH, IMIN, IMAX, JC

IF (LDONE) RETURN
LDONE = .TRUE.

CLHOST = ' '
CALL GET_ENVIRONMENT_VARIABLE('HOSTNAME', CLHOST, STATUS=ISTATUS)

! Simple, order-sensitive rolling hash of the hostname string.
IHASH = 5381
DO JC = 1, LEN_TRIM(CLHOST)
  IHASH = MOD(IHASH * 33 + ICHAR(CLHOST(JC:JC)), 1000000007_JPIM)
ENDDO

IMIN = IHASH
IMAX = IHASH
CALL MPL_ALLREDUCE(IMIN, 'MIN', KCOMM=MPL_ALL_LEVS_COMM, CDSTRING='VERIFY_LEVS_ON_NODE:MIN')
CALL MPL_ALLREDUCE(IMAX, 'MAX', KCOMM=MPL_ALL_LEVS_COMM, CDSTRING='VERIFY_LEVS_ON_NODE:MAX')

IF (MYSETV == 1) THEN
  IF (IMIN == IMAX) THEN
    WRITE(NOUT,'(A,I5,A,I5,A,A)') 'VERIFY_LEVS_ON_NODE: wave-set ', MYSETW, &
      & ' : all ', NPRTRV, ' LEVS ranks ON ONE NODE, host=', TRIM(CLHOST)
  ELSE
    WRITE(NOUT,'(A,I5,A,I5,A)') 'VERIFY_LEVS_ON_NODE: wave-set ', MYSETW, &
      & ' : LEVS ranks SPAN MULTIPLE NODES (', NPRTRV, ' ranks) -- LEVS_TR not intra-node!'
  ENDIF
ENDIF

END SUBROUTINE VERIFY_LEVS_ON_NODE

! ---------------------------------------------------------------------------

SUBROUTINE ACCOUNT_TR_BW(KIDX, CDLABEL, KBYTES, PTIME, KCOMM, KCOMMSIZE, LPRINT)

! Accumulate the payload byte count KBYTES and wall time PTIME of one transpose
! comm region (KIDX: 1=TRLTOG L->G, 2=TRGTOL G->L, 11=TRMTOL M->L,
! 12=TRLTOM L->M) into the module counters. The first NBW_WARMUP counted calls
! are discarded; the next NBW_SAMPLE are summed. On the call that completes the
! sample window, reduce across the transpose's communicator and print achieved
! bandwidth and its fraction of ZTR_BW_ROOFLINE, once from the reduction root.
! G<->L regions reduce over MPL_ALL_LEVS_COMM (whose ranks share one node under
! Option A); the M<->L regions reduce over MPL_ALL_MS_COMM (one rank per node,
! so the aggregate is the full transpose traffic). Pass KCOMM to override the
! reduction communicator, KCOMMSIZE its rank count (printed for the report),
! and LPRINT the rank that prints (e.g. MYPROC==1 for one summary line instead
! of one per sub-comm).
!
! KBYTES is passed as an INTEGER(JPIB) (computed at the call site, where the
! payload kind JPRB is known) so the cross-rank reduction can use an integer
! SUM: exact/reproducible, avoiding FIAT's real-mode SUM reproducibility abort.
! Time is reduced with MAX (order-independent, so no reproducibility guard).
!
! Aggregate bandwidth = SUM(bytes over the comm ranks) divided by MAX(time over
! those ranks): the slowest rank bounds the collective, and the summed bytes
! are the total moved on that comm's memory/network system.

USE TPM_DISTR,        ONLY : MYSETW, MYSETV, NPRTRV
USE TPM_GEN,          ONLY : NOUT
USE TPM_ECTRANS_OPTS, ONLY : LREPORT_TR_BW, ZTR_BW_ROOFLINE
USE MPL_GROUPS,       ONLY : MPL_ALL_LEVS_COMM
USE MPL_MODULE,       ONLY : MPL_ALLREDUCE, MPL_ALL_MS_COMM

INTEGER(KIND=JPIM), INTENT(IN) :: KIDX
CHARACTER(LEN=*),   INTENT(IN) :: CDLABEL
INTEGER(KIND=JPIB), INTENT(IN) :: KBYTES
REAL(KIND=JPRD),    INTENT(IN) :: PTIME
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KCOMM
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KCOMMSIZE
LOGICAL,            OPTIONAL, INTENT(IN) :: LPRINT

INTEGER(KIND=JPIB) :: IBYTES_SUM
REAL(KIND=JPRD)    :: ZTIME_MAX, ZBYTES_SUM, ZGBPS, ZPCT, ZGB
INTEGER(KIND=JPIM) :: ISTEADY, KCOMM_RESOLVED, KCOMMSIZE_RESOLVED
LOGICAL            :: LPRINT_RESOLVED

IF (.NOT. LREPORT_TR_BW) RETURN
IF (KIDX < 1 .OR. KIDX > NBW_IDX) RETURN

KCOMM_RESOLVED  = MPL_ALL_LEVS_COMM
IF (PRESENT(KCOMM)) KCOMM_RESOLVED = KCOMM
KCOMMSIZE_RESOLVED = NPRTRV
IF (PRESENT(KCOMMSIZE)) KCOMMSIZE_RESOLVED = KCOMMSIZE
LPRINT_RESOLVED = (MYSETV == 1)
IF (PRESENT(LPRINT)) LPRINT_RESOLVED = LPRINT

NBW_CALLS(KIDX) = NBW_CALLS(KIDX) + 1

! Discard warmup, and ignore anything after the sample window (already printed).
IF (NBW_CALLS(KIDX) <= NBW_WARMUP) RETURN
IF (NBW_CALLS(KIDX) >  NBW_WARMUP + NBW_SAMPLE) RETURN

IBW_BYTES(KIDX) = IBW_BYTES(KIDX) + KBYTES
ZBW_TIME(KIDX)  = ZBW_TIME(KIDX)  + PTIME

IF (NBW_CALLS(KIDX) /= NBW_WARMUP + NBW_SAMPLE) RETURN

! Sample window complete for this index: reduce over the transpose comm.
ISTEADY    = NBW_SAMPLE
IBYTES_SUM = IBW_BYTES(KIDX)
ZTIME_MAX  = ZBW_TIME(KIDX)
CALL MPL_ALLREDUCE(IBYTES_SUM, 'SUM', KCOMM=KCOMM_RESOLVED, CDSTRING='ACCOUNT_TR_BW:BYTES')
CALL MPL_ALLREDUCE(ZTIME_MAX,  'MAX', KCOMM=KCOMM_RESOLVED, CDSTRING='ACCOUNT_TR_BW:TIME')

IF (LPRINT_RESOLVED) THEN
  ZBYTES_SUM = REAL(IBYTES_SUM, JPRD)
  ZGBPS = 0.0_JPRD
  IF (ZTIME_MAX > 0.0_JPRD) ZGBPS = ZBYTES_SUM / ZTIME_MAX / 1.0E9_JPRD
  ZPCT = 0.0_JPRD
  IF (ZTR_BW_ROOFLINE > 0.0) ZPCT = 100.0_JPRD * ZGBPS / REAL(ZTR_BW_ROOFLINE, JPRD)
  ZGB  = ZBYTES_SUM / 1.0E9_JPRD
  WRITE(NOUT,'(A)')          '=== ecTrans transpose bandwidth (LREPORT_TR_BW) ==='
  WRITE(NOUT,'(A,A)')        '  region                  = ', TRIM(CDLABEL)
  WRITE(NOUT,'(A,I5)')       '  wave-set (MYSETW)       = ', MYSETW
  WRITE(NOUT,'(A,I5,A,I5)')  '  ranks in reduction comm = ', KCOMMSIZE_RESOLVED, '   samples summed  = ', ISTEADY
  WRITE(NOUT,'(A,F12.3)')    '  node payload / window   [GB]   = ', ZGB
  WRITE(NOUT,'(A,F12.6)')    '  node comm time / window [s]    = ', ZTIME_MAX
  WRITE(NOUT,'(A,F12.2)')    '  achieved node bandwidth [GB/s] = ', ZGBPS
  WRITE(NOUT,'(A,F12.1,A,F8.2,A)') '  roofline [GB/s]         = ', REAL(ZTR_BW_ROOFLINE, JPRD), &
    & '   -> ', ZPCT, ' % of peak'
ENDIF

END SUBROUTINE ACCOUNT_TR_BW

! ---------------------------------------------------------------------------

SUBROUTINE DUMP_TR_SKEW(KIDX, PWAIT)

! Straggler profiler (LDUMP_TR_SKEW). Collective over the WORLD communicator:
! all-gather every rank's pre-issue barrier wait PWAIT and, on global rank 1,
! append one wide CSV row (idx, call, then NPROC waits in microseconds) to
! tr_skew.csv for each sampled call. Offline this reveals whether the skew is
! caused by the same ranks/nodes arriving last every step (systematic, => a
! decomposition/placement/slow-node issue) or by a random set (OS/MPI jitter).
!
! The barrier wait is the natural per-rank imbalance proxy: the last-arriving
! rank has 0 wait, the earliest has wait~=arrival spread. A systematically 
! slow rank therefore shows a persistently small wait.
!
! Same warmup/sample window as ACCOUNT_TR_BW; must be called by ALL world ranks
! in lockstep so the all-gather and the window guard stay aligned.

USE TPM_DISTR,        ONLY : MYPROC, NPROC
USE TPM_ECTRANS_OPTS, ONLY : LREPORT_TR_BW, LDUMP_TR_SKEW
USE MPL_MODULE,       ONLY : MPL_ALLGATHERV

INTEGER(KIND=JPIM), INTENT(IN) :: KIDX
REAL(KIND=JPRD),    INTENT(IN) :: PWAIT

REAL(KIND=JPRD)    :: ZSEND(1)
REAL(KIND=JPRD),    ALLOCATABLE :: ZALL(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: ICOUNTS(:)
INTEGER(KIND=JPIM) :: JR

IF (.NOT. LREPORT_TR_BW) RETURN
IF (.NOT. LDUMP_TR_SKEW) RETURN
IF (KIDX < 1 .OR. KIDX > NBW_IDX) RETURN

NSK_CALLS(KIDX) = NSK_CALLS(KIDX) + 1
IF (NSK_CALLS(KIDX) <= NBW_WARMUP) RETURN
IF (NSK_CALLS(KIDX) >  NBW_WARMUP + NBW_SAMPLE) RETURN

ALLOCATE(ZALL(NPROC))
ALLOCATE(ICOUNTS(NPROC))
ICOUNTS(:) = 1
ZSEND(1)   = PWAIT
CALL MPL_ALLGATHERV(ZSEND(1:1), ZALL, ICOUNTS, CDSTRING='DUMP_TR_SKEW')

IF (MYPROC == 1) THEN
  IF (.NOT. LSK_OPEN) THEN
    ! NEWUNIT returns a negative unit number, so guard reopen with a flag.
    OPEN(NEWUNIT=ISK_UNIT, FILE='tr_skew.csv', STATUS='REPLACE', ACTION='WRITE')
    LSK_OPEN = .TRUE.
    WRITE(ISK_UNIT,'(A,I0)') '# per-rank pre-issue barrier wait [microseconds]; nproc=', NPROC
    WRITE(ISK_UNIT,'(A)')    '# idx(9=TRLTOG,10=TRGTOL,13=TRMTOL,14=TRLTOM),call,wait_rank1..wait_rankNPROC'
  ENDIF
  WRITE(ISK_UNIT,'(I0,A,I0)',ADVANCE='NO') KIDX, ',', NSK_CALLS(KIDX) - NBW_WARMUP
  DO JR = 1, NPROC
    WRITE(ISK_UNIT,'(A,I0)',ADVANCE='NO') ',', NINT(ZALL(JR) * 1.0E6_JPRD, JPIM)
  ENDDO
  WRITE(ISK_UNIT,'(A)') ''
ENDIF

DEALLOCATE(ZALL, ICOUNTS)

END SUBROUTINE DUMP_TR_SKEW

! ---------------------------------------------------------------------------

END MODULE TPM_LEVS_TRANSPOSE
