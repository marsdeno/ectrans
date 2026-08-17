! (C) Copyright 1995- ECMWF.
! (C) Copyright 1995- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TRLTOG_MOD

IMPLICIT NONE

PUBLIC TRLTOG
PRIVATE TRLTOG_COMM

CONTAINS

SUBROUTINE TRLTOG(PGLAT, KF_FS, KF_GP, KF_SCALARS_G, KVSET, KPTRGP, PGP, PGPUV, PGP3A, PGP3B, PGP2)

!**** *TRLTOG * - head routine for transposition of grid point data from latitudinal
!                 to column structure (this takes place between inverse
!                 FFT and grid point calculations)
!                 TRLTOG is the inverse of TRGTOL

!**   Interface.
!     ----------
!        *call* *TRLTOG(...)

!        Explicit arguments :
!        --------------------
!           PGLAT    -  Latitudinal data ready for direct FFT (input)
!           PGP    -  Blocked grid point data    (output)
!           KVSET    - "v-set" for each field      (input)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        R. El Khatib *Meteo-France*

!     Modifications.
!     --------------
!        Original  : 18-Aug-2014 from trltog
!        R. El Khatib 09-Sep-2020 NSTACK_MEMORY_TR
!     ------------------------------------------------------------------

USE PARKIND1 ,ONLY : JPIM     ,JPRB
USE YOMHOOK  ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DISTR,ONLY : D
USE TRGL_MOD ,ONLY: TRGL_BUFFERS, ALLOCATE_BUFFERS_CST, TRGL_PROLOG, ALLOCATE_BUFFERS_SR

IMPLICIT NONE

REAL(KIND=JPRB),   INTENT(IN) :: PGLAT(KF_FS,D%NLENGTF)
INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS,KF_GP
INTEGER(KIND=JPIM),INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM),INTENT(IN) :: KVSET(KF_GP)
INTEGER(KIND=JPIM),OPTIONAL, INTENT(IN) :: KPTRGP(:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)     :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)     :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)     :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)     :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)     :: PGP2(:,:,:)

TYPE(TRGL_BUFFERS) :: YDBUFS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TRLTOG',0,ZHOOK_HANDLE)

YDBUFS%LLTRGTOL = .FALSE.
CALL ALLOCATE_BUFFERS_CST(YDBUFS)
CALL GSTATS(1806, 0)
CALL TRGL_PROLOG(KF_FS, KF_GP, KVSET, YDBUFS)
CALL GSTATS(1806, 1)
CALL ALLOCATE_BUFFERS_SR(YDBUFS, KF_GP)

CALL TRLTOG_COMM(PGLAT, KF_FS, KF_GP, KF_SCALARS_G, KVSET, KPTRGP, PGP, PGPUV, PGP3A, PGP3B, PGP2, &
  &              YDBUFS)

IF (LHOOK) CALL DR_HOOK('TRLTOG',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE TRLTOG

SUBROUTINE TRLTOG_COMM(PGLAT, KF_FS, KF_GP, KF_SCALARS_G, KVSET, KPTRGP, PGP, PGPUV, PGP3A, PGP3B, &
  &                    PGP2,YDBUFS)


!**** *trltog * - transposition of grid point data from latitudinal
!                 to column structure. This takes place between inverse
!                 FFT and grid point calculations.
!                 TRLTOG_COMM is the inverse of TRGTOL

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *call* *trltog(...)

!        Explicit arguments :
!        --------------------
!           PGLAT    -  Latitudinal data ready for direct FFT (input)
!           PGP    -  Blocked grid point data    (output)
!           KVSET    - "v-set" for each field      (input)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        MPP Group *ECMWF*

!     Modifications.
!     --------------
!        Original  : 95-10-01
!        D.Dent    : 97-08-04 Reorganisation to allow NPRTRV
!                             to differ from NPRGPEW
!        =99-03-29= Mats Hamrud and Deborah Salmond
!                   JUMP in FFT's changed to 1
!                   KINDEX introduced and PCOMBUF not used for same PE
!         01-11-23  Deborah Salmond and John Hague
!                   LIMP_NOOLAP Option for non-overlapping message passing
!                               and buffer packing
!         01-12-18  Peter Towers
!                   Improved vector performance of LTOG_PACK,LTOG_UNPACK
!         03-0-02   G. Radnoti: Call barrier always when nproc>1
!         08-01-01  G.Mozdzynski: cleanup
!         09-01-02  G.Mozdzynski: use non-blocking recv and send
!        R. El Khatib 09-Sep-2020 64 bits addressing for PGLAT
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPIB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE MPL_MODULE  ,ONLY : MPL_RECV, MPL_SEND, MPL_WAIT, JP_NON_BLOCKING_STANDARD, MPL_WAITANY, &
  &                     JP_BLOCKING_STANDARD, MPL_BARRIER, JP_BLOCKING_BUFFERED, MPL_ALLTOALLV
USE MPL_GROUPS  ,ONLY : MPL_ALL_LEVS_COMM

USE TPM_GEN         ,ONLY : NTRANS_SYNC_LEVEL, NSTACK_MEMORY_TR, NOUT
USE TPM_DISTR       ,ONLY : D, MTAGLG, NPRCIDS, MYPROC, NPROC, NPRTRV
USE TPM_ECTRANS_OPTS,ONLY : LUSE_OPT2, LUSE_OPT3, LUSE_LEVS_TR, LVERIFY_LEVS_NODE, LREPORT_TR_BW
USE TPM_LEVS_TRANSPOSE,ONLY : LEVS_TR_ALIGNED, VERIFY_LEVS_ON_NODE, ACCOUNT_TR_BW, DUMP_TR_SKEW

USE TRGL_MOD, ONLY: TRGL_BUFFERS, TRGL_VARS, TRGL_ALLOCATE_VARS, TRGL_ALLOCATE_HEAP_BUFFER, &
  &                 TRGL_INIT_VARS, TRGL_INIT_OFF_VARS, TGRL_COPY_ZCOMBUF, TGRL_COPY_PGLAT, &
  &                 TGRL_INIT_PACKING_VARS

IMPLICIT NONE


INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS,KF_GP
REAL(KIND=JPRB),INTENT(IN)     :: PGLAT(KF_FS,D%NLENGTF)
INTEGER(KIND=JPIM), INTENT(IN) :: KVSET(KF_GP)
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPTRGP(:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PGP2(:,:,:)

TYPE (TRGL_BUFFERS), INTENT(INOUT), TARGET :: YDBUFS
! LOCAL VARIABLES
TYPE(TRGL_VARS) :: YLVARS

! OPT2: single-alltoallv counts/displs (world-comm indexed 1:NPROC).
! LUSE_OPT2=.FALSE. executes origianl per-peer MPL_SEND/MPL_RECV path
INTEGER(KIND=JPIM) :: ILENS(NPROC), IOFFS(NPROC), ILENR(NPROC), IOFFR(NPROC)
INTEGER(KIND=JPIM) :: JP
! LEVS-comm variant (Option A): counts/displs indexed over the NPRTRV ranks
! of MPL_ALL_LEVS_COMM. LLEVS gates the intra-node alltoallv path.
INTEGER(KIND=JPIM) :: ILENS_L(NPRTRV), IOFFS_L(NPRTRV), ILENR_L(NPRTRV), IOFFR_L(NPRTRV)
INTEGER(KIND=JPIM) :: JPL
LOGICAL :: LLEVS
! LLOVERLAP gates the OPT3 nonblocking self-copy overlap.
LOGICAL :: LLOVERLAP
LOGICAL, SAVE :: LLEVS_REPORTED = .FALSE.
! OPT3: request handle for the non-blocking MPL_IALLTOALLV overlap.
INTEGER(KIND=JPIM) :: IREQ_A2A

INTEGER(KIND=JPIM) :: IREQ_SEND(NPROC)
INTEGER(KIND=JPIM) :: IREQ_RECV(NPROC)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IRECV
INTEGER(KIND=JPIM) :: ISEND, ITAG, JL, JFLD, INS, INR, JNR
INTEGER(KIND=JPIM) :: II,ILEN
INTEGER(KIND=JPIM) :: ISEND_FLD_START,ISEND_FLD_END
! hoist KNDOFF(ISEND) out of the JL-parallel inner loop so each
! thread reads a scalar rather than re-derefing YDBUFS%INDOFF per JL.
INTEGER(KIND=JPIM) :: I_KNDOFF_ISEND

!     LOCAL ARRAYS
REAL(KIND=JPRB), TARGET :: ZCOMBUFS_STACK(-1:YDBUFS%ISENDCOUNT,MERGE (YDBUFS%INSEND,0,NSTACK_MEMORY_TR/=0))
REAL(KIND=JPRB), TARGET :: ZCOMBUFR_STACK(-1:YDBUFS%IRECVCOUNT,MERGE (YDBUFS%INRECV,0,NSTACK_MEMORY_TR/=0))

REAL(KIND=JPRB), ALLOCATABLE, TARGET, SAVE :: ZCOMBUFS_HEAP(:,:)
REAL(KIND=JPRB), ALLOCATABLE, TARGET, SAVE :: ZCOMBUFR_HEAP(:,:)

REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZCOMBUFS(:,:)
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZCOMBUFR(:,:)

! OPT2: 1D rank-remap pointers for MPL_ALLTOALLV zero-copy binding.
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZSEND_1D(:), ZRECV_1D(:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE_BAR
! OPT2: DR_HOOK TRLTOG_COMM phases, split TRLTOG_COMM wall time into pack/alltoallv/self-copy/unpack.
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE_PACK, ZHOOK_HANDLE_A2A
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE_SELF, ZHOOK_HANDLE_UNPACK

! LREPORT_TR_BW: wall-clock timing + payload byte count for this transpose.
INTEGER(KIND=JPIB) :: ITR_CLK0, ITR_CLK1, ITR_CLK_RATE, ITR_BYTES
INTEGER(KIND=JPIB) :: ICM_CLK0, ICM_CLK1, IPH_CLK0, IPH_CLK1
! ZCM_DT = non-blocking issue span (idx 3); ZWT_DT = MPL_WAIT span, the comm
! completion not hidden by the self-copy (idx 5); ZSC_DT = deferred self-copy
! span, the OPT3 overlap budget (idx 7); ZBR_DT = pre-issue barrier span, the
! load-imbalance skew across the alltoallv comm (idx 9, LREPORT_TR_BW only).
REAL(KIND=JPRD)    :: ZTR_DT, ZCM_DT, ZWT_DT, ZSC_DT, ZBR_DT

!     ------------------------------------------------------------------

!*       0.    Some initializations
!              --------------------
ASSOCIATE(KNSEND=>YDBUFS%INSEND, KNRECV=>YDBUFS%INRECV, KSENDTOT=>YDBUFS%ISENDTOT, &
  &       KRECVTOT=>YDBUFS%IRECVTOT, KSEND=>YDBUFS%ISEND, KRECV=>YDBUFS%IRECV, &
  &       KINDEX=>YDBUFS%IINDEX, KNDOFF=>YDBUFS%INDOFF)

IF (NSTACK_MEMORY_TR == 0) THEN
  CALL TRGL_ALLOCATE_HEAP_BUFFER(ZCOMBUFS_HEAP, YDBUFS%ISENDCOUNT, YDBUFS%INSEND)
  CALL TRGL_ALLOCATE_HEAP_BUFFER(ZCOMBUFR_HEAP, YDBUFS%IRECVCOUNT, YDBUFS%INRECV)

! Now, force the OS to allocate this shared array right now, not when it starts to be used which is
! an OPEN-MP loop, that would cause a threads synchronization lock :
  IF (YDBUFS%INSEND > 0 .AND. YDBUFS%ISENDCOUNT >=-1) ZCOMBUFS_HEAP(-1,1)=HUGE(1._JPRB)
  ZCOMBUFS (-1:,1:) => ZCOMBUFS_HEAP
  ZCOMBUFR (-1:,1:) => ZCOMBUFR_HEAP
ELSE
  ZCOMBUFS (-1:,1:) => ZCOMBUFS_STACK
  ZCOMBUFR (-1:,1:) => ZCOMBUFR_STACK
ENDIF

! LEVS-comm intra-node transpose selection (needs collar==waveset alignment).
LLEVS     = LUSE_LEVS_TR .AND. LEVS_TR_ALIGNED()
LLOVERLAP = LUSE_OPT3

ITAG = MTAGLG

IF (LHOOK) CALL DR_HOOK('TRLTOG_BAR',0,ZHOOK_HANDLE_BAR)
CALL GSTATS_BARRIER(762)
IF (LHOOK) CALL DR_HOOK('TRLTOG_BAR',1,ZHOOK_HANDLE_BAR)

IF (LREPORT_TR_BW) CALL SYSTEM_CLOCK(COUNT=ITR_CLK0, COUNT_RATE=ITR_CLK_RATE)

CALL GSTATS(805,0)

! Pre-post IRECV: only needed for the per-task path
! (LUSE_OPT2=.FALSE.), skip in the single-alltoallv path.
IF (.NOT. LUSE_OPT2) THEN
  IF (NTRANS_SYNC_LEVEL <= 0) THEN
    !...Receive loop.........................................................
    DO INR=1,KNRECV
      IRECV=KRECV(INR)
      CALL MPL_RECV(ZCOMBUFR(-1:KRECVTOT(IRECV),INR), KSOURCE=NPRCIDS(IRECV), &
        &           KMP_TYPE=JP_NON_BLOCKING_STANDARD, KREQUEST=IREQ_RECV(INR), KTAG=ITAG, &
        &           CDSTRING='TRLTOG_COMM: NON-BLOCKING IRECV' )
    ENDDO
  ENDIF
ENDIF

CALL GSTATS(805,1)

CALL GSTATS(1806,0)
YDBUFS%LLINDER = PRESENT(KPTRGP)
YDBUFS%LLPGPONLY = PRESENT(PGP)
CALL TRGL_ALLOCATE_VARS(YLVARS, KF_GP,KF_FS)
CALL TRGL_INIT_VARS(YLVARS, KF_SCALARS_G, PGP, PGPUV, PGP3A, PGP3B, PGP2)
CALL GSTATS(1806,1)

! Copy local contribution

! Copy local contribution
!
! OPT3: when LUSE_OPT3, defer the self-copy until after the nonblocking
! alltoallv has been posted so it can overlap with the collective. The self-
! copy touches only PGLAT (read) and PGP/PGPUV/PGP3A/PGP3B/PGP2 (write); the
! alltoallv touches only ZCOMBUFS (read) and ZCOMBUFR (write)

IF(KRECVTOT(MYPROC) > 0 .AND. .NOT. LLOVERLAP) THEN
  CALL TRGL_INIT_OFF_VARS(YDBUFS,YLVARS,KVSET,KPTRGP,KF_GP)
  IF (LHOOK) CALL DR_HOOK('TRLTOG_SELF',0,ZHOOK_HANDLE_SELF)
  CALL GSTATS(1604,0)
  CALL TGRL_COPY_PGLAT(PGLAT, YDBUFS, YLVARS, PGP, PGPUV,PGP3A, PGP3B,PGP2)
  CALL GSTATS(1604,1)
  IF (LHOOK) CALL DR_HOOK('TRLTOG_SELF',1,ZHOOK_HANDLE_SELF)
ENDIF
!
! loop over the number of processors we need to communicate with.
! NOT MYPROC
!
! Now overlapping buffer packing/unpacking with sends/waits
! Time as if all communications to avoid double accounting

CALL GSTATS(805,0)

IF (LUSE_OPT2) THEN
  ! OPT2: single MPL_ALLTOALLV replaces the per-peer ISEND / IRECV / WAITANY
  ! Buffer layout unchanged, so TGRL_COPY_ZCOMBUF unpack works
  ! bit-identically

  !....Pack loop..........................................................

  IF (LHOOK) CALL DR_HOOK('TRLTOG_PACK',0,ZHOOK_HANDLE_PACK)
  ISEND_FLD_START = 1
  ISEND_FLD_END   = KF_FS
  DO INS=1,KNSEND
    ISEND=KSEND(INS)
    ILEN = KSENDTOT(ISEND)/KF_FS
    I_KNDOFF_ISEND = KNDOFF(ISEND)
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,JL,II)
    DO JL=1,ILEN
      II = KINDEX(I_KNDOFF_ISEND+JL)
      DO JFLD=ISEND_FLD_START,ISEND_FLD_END
        ZCOMBUFS((JFLD-ISEND_FLD_START)*ILEN+JL,INS) = PGLAT(JFLD,II)
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    ZCOMBUFS(-1,INS) = 1
    ZCOMBUFS(0,INS)  = KF_FS
  ENDDO
  IF (LHOOK) CALL DR_HOOK('TRLTOG_PACK',1,ZHOOK_HANDLE_PACK)

<<<<<<< HEAD
  ! Build alltoallv counts/displs. MYPROC is excluded from KSEND/KRECV lists
  ! by TRGL_PROLOG so ILENS/ILENR are 0 on MYPROC entries; 
  ! the self-transfer is handled by TGRL_COPY_PGLAT
  ILENS = 0
  IOFFS = 0
  ILENR = 0
  IOFFR = 0
  DO INS=1,KNSEND
    JP = KSEND(INS)
    ILENS(JP) = KSENDTOT(JP) + 2
    IOFFS(JP) = (INS-1) * (YDBUFS%ISENDCOUNT + 2)
  ENDDO
  DO INR=1,KNRECV
    JP = KRECV(INR)
    ILENR(JP) = KRECVTOT(JP) + 2
    IOFFR(JP) = (INR-1) * (YDBUFS%IRECVCOUNT + 2)
  ENDDO
=======
  ! Build alltoallv counts/displs (world-comm indexed 1:NPROC). MYPROC is
  ! excluded from KSEND/KRECV lists by TRGL_PROLOG so ILENS/ILENR are 0
  ! on MYPROC entries; the self-transfer is handled by TGRL_COPY_PGLAT
  ! above.
  !
  ! LEVS-comm variant (Option A alignment): see trgtol_mod.F90 for the full
  ! rationale -- under collar==waveset alignment every peer JP lives in
  ! MPL_ALL_LEVS_COMM at local rank MOD(JP-1,NPRTRV), so the alltoallv can be
  ! issued on the 32-rank intra-node communicator.
  ! (LLEVS was resolved before the pack; do not recompute here.)

  ! (#4) placement diagnostic -- runs at most once, independent of LLEVS.
  IF (LVERIFY_LEVS_NODE) CALL VERIFY_LEVS_ON_NODE()

  ! Report the branch actually taken (runtime LLEVS), once, from the LEVS
  ! root, so the intra-node alltoallv routing is unambiguous.
  IF (LVERIFY_LEVS_NODE .AND. MYPROC == 1 .AND. .NOT. LLEVS_REPORTED) THEN
    LLEVS_REPORTED = .TRUE.
    WRITE(NOUT,'(A,L1,A,L1,A)') 'TRLTOG LEVS branch: LLEVS=', LLEVS, &
      & '  (LUSE_LEVS_TR=', LUSE_LEVS_TR, ')'
  ENDIF

  IF (LLEVS) THEN
    ILENS_L = 0
    IOFFS_L = 0
    ILENR_L = 0
    IOFFR_L = 0
    DO INS=1,KNSEND
      JP  = KSEND(INS)
      JPL = MOD(JP-1, NPRTRV) + 1
      ILENS_L(JPL) = KSENDTOT(JP) + 2
      IOFFS_L(JPL) = (INS-1) * (YDBUFS%ISENDCOUNT + 2)
    ENDDO
    DO INR=1,KNRECV
      JP  = KRECV(INR)
      JPL = MOD(JP-1, NPRTRV) + 1
      ILENR_L(JPL) = KRECVTOT(JP) + 2
      IOFFR_L(JPL) = (INR-1) * (YDBUFS%IRECVCOUNT + 2)
    ENDDO
  ELSE
    ILENS = 0
    IOFFS = 0
    ILENR = 0
    IOFFR = 0
    DO INS=1,KNSEND
      JP = KSEND(INS)
      ILENS(JP) = KSENDTOT(JP) + 2
      IOFFS(JP) = (INS-1) * (YDBUFS%ISENDCOUNT + 2)
    ENDDO
    DO INR=1,KNRECV
      JP = KRECV(INR)
      ILENR(JP) = KRECVTOT(JP) + 2
      IOFFR(JP) = (INR-1) * (YDBUFS%IRECVCOUNT + 2)
    ENDDO
  ENDIF
>>>>>>> 8ea5cef (G<->L and M<->L transposes: opt-in LEVS-comm routing + bandwidth reporting)

  ZSEND_1D(1:SIZE(ZCOMBUFS)) => ZCOMBUFS
  ZRECV_1D(1:SIZE(ZCOMBUFR)) => ZCOMBUFR

  ! Skip the collective when NPROC=1
  !
  ! OPT3: with LUSE_OPT3, MPL_ALLTOALLV is issued in non-blocking mode
  ! (KREQUEST=IREQ_A2A) and MPL_WAIT is called AFTER the self-copy so the two
  ! can overlap. When Opt3 is off, this is a plain blocking alltoallv.
  IF (NPROC > 1) THEN
    IF (LHOOK) CALL DR_HOOK('TRLTOG_A2A',0,ZHOOK_HANDLE_A2A)
    IF (LREPORT_TR_BW) CALL SYSTEM_CLOCK(COUNT=ICM_CLK0)
    ZWT_DT = 0.0_JPRD
    ZSC_DT = 0.0_JPRD
    ZBR_DT = 0.0_JPRD
    ! Diagnostic (LREPORT_TR_BW only): barrier over the same comm the alltoallv
    ! synchronizes on, timed as idx 9. Ranks enter the alltoallv staggered by
    ! their (spectrally imbalanced) pack times; this barrier charges that skew
    ! separately so the following MPL_WAIT (idx 5) reflects ~pure transfer.
    IF (LREPORT_TR_BW) THEN
      CALL SYSTEM_CLOCK(COUNT=IPH_CLK0)
      IF (LLEVS) THEN
        CALL MPL_BARRIER(KCOMM=MPL_ALL_LEVS_COMM, CDSTRING='TRLTOG_COMM: SKEW BARRIER')
      ELSE
        CALL MPL_BARRIER(CDSTRING='TRLTOG_COMM: SKEW BARRIER')
      ENDIF
      CALL SYSTEM_CLOCK(COUNT=IPH_CLK1)
      ZBR_DT = REAL(IPH_CLK1 - IPH_CLK0, JPRD) / REAL(ITR_CLK_RATE, JPRD)
      CALL SYSTEM_CLOCK(COUNT=ICM_CLK0)
    ENDIF
    IF (LLEVS) THEN
      ! Intra-node collective over the 32-rank LEVS communicator.
      IF (LUSE_OPT3) THEN
        CALL MPL_ALLTOALLV(PSENDBUF=ZSEND_1D, KSENDCOUNTS=ILENS_L, KSENDDISPL=IOFFS_L, &
          &                PRECVBUF=ZRECV_1D, KRECVCOUNTS=ILENR_L, KRECVDISPL=IOFFR_L, &
          &                KCOMM=MPL_ALL_LEVS_COMM, &
          &                KMP_TYPE=JP_NON_BLOCKING_STANDARD, KREQUEST=IREQ_A2A, &
          &                CDSTRING='TRLTOG_COMM: IALLTOALLV(LEVS)')
      ELSE
        CALL MPL_ALLTOALLV(PSENDBUF=ZSEND_1D, KSENDCOUNTS=ILENS_L, KSENDDISPL=IOFFS_L, &
          &                PRECVBUF=ZRECV_1D, KRECVCOUNTS=ILENR_L, KRECVDISPL=IOFFR_L, &
          &                KCOMM=MPL_ALL_LEVS_COMM, &
          &                KMP_TYPE=JP_BLOCKING_STANDARD, &
          &                CDSTRING='TRLTOG_COMM: ALLTOALLV(LEVS)')
      ENDIF
    ELSE IF (LUSE_OPT3) THEN
      CALL MPL_ALLTOALLV(PSENDBUF=ZSEND_1D, KSENDCOUNTS=ILENS, KSENDDISPL=IOFFS, &
        &                PRECVBUF=ZRECV_1D, KRECVCOUNTS=ILENR, KRECVDISPL=IOFFR, &
        &                KMP_TYPE=JP_NON_BLOCKING_STANDARD, KREQUEST=IREQ_A2A, &
        &                CDSTRING='TRLTOG_COMM: IALLTOALLV')
    ELSE
      CALL MPL_ALLTOALLV(PSENDBUF=ZSEND_1D, KSENDCOUNTS=ILENS, KSENDDISPL=IOFFS, &
        &                PRECVBUF=ZRECV_1D, KRECVCOUNTS=ILENR, KRECVDISPL=IOFFR, &
        &                KMP_TYPE=JP_BLOCKING_STANDARD, &
        &                CDSTRING='TRLTOG_COMM: ALLTOALLV')
    ENDIF
    IF (LREPORT_TR_BW) THEN
      CALL SYSTEM_CLOCK(COUNT=ICM_CLK1)
      ZCM_DT = REAL(ICM_CLK1 - ICM_CLK0, JPRD) / REAL(ITR_CLK_RATE, JPRD)
    ENDIF
    IF (LHOOK) CALL DR_HOOK('TRLTOG_A2A',1,ZHOOK_HANDLE_A2A)
  ENDIF

<<<<<<< HEAD
  ! deferred self-copy (runs in the shadow of IALLTOALLV when LUSE_OPT3).
  IF (LUSE_OPT3 .AND. KRECVTOT(MYPROC) > 0) THEN
=======
  ! Deferred self-copy (runs in the shadow of IALLTOALLV when LUSE_OPT3).
  ! Timer 1604 accounting is unchanged; only the wall-clock position moves.
  IF (LLOVERLAP .AND. KRECVTOT(MYPROC) > 0) THEN
>>>>>>> 8ea5cef (G<->L and M<->L transposes: opt-in LEVS-comm routing + bandwidth reporting)
    CALL TRGL_INIT_OFF_VARS(YDBUFS,YLVARS,KVSET,KPTRGP,KF_GP)
    IF (LHOOK) CALL DR_HOOK('TRLTOG_SELF',0,ZHOOK_HANDLE_SELF)
    IF (LREPORT_TR_BW) CALL SYSTEM_CLOCK(COUNT=IPH_CLK0)
    CALL GSTATS(1604,0)
    CALL TGRL_COPY_PGLAT(PGLAT, YDBUFS, YLVARS, PGP, PGPUV,PGP3A, PGP3B,PGP2)
    CALL GSTATS(1604,1)
    IF (LREPORT_TR_BW) THEN
      CALL SYSTEM_CLOCK(COUNT=IPH_CLK1)
      ZSC_DT = REAL(IPH_CLK1 - IPH_CLK0, JPRD) / REAL(ITR_CLK_RATE, JPRD)
    ENDIF
    IF (LHOOK) CALL DR_HOOK('TRLTOG_SELF',1,ZHOOK_HANDLE_SELF)
  ENDIF

  ! Wait for the nonblocking collective to complete before touching ZCOMBUFR.
  IF (LLOVERLAP .AND. NPROC > 1) THEN
    IF (LHOOK) CALL DR_HOOK('TRLTOG_A2A',0,ZHOOK_HANDLE_A2A)
    IF (LREPORT_TR_BW) CALL SYSTEM_CLOCK(COUNT=IPH_CLK0)
    CALL MPL_WAIT(KREQUEST=IREQ_A2A, CDSTRING='TRLTOG_COMM: WAIT FOR IALLTOALLV')
    IF (LREPORT_TR_BW) THEN
      CALL SYSTEM_CLOCK(COUNT=IPH_CLK1)
      ZWT_DT = REAL(IPH_CLK1 - IPH_CLK0, JPRD) / REAL(ITR_CLK_RATE, JPRD)
    ENDIF
    IF (LHOOK) CALL DR_HOOK('TRLTOG_A2A',1,ZHOOK_HANDLE_A2A)
  ENDIF

  !  Unpack loop..........................................................

  IF (LHOOK) CALL DR_HOOK('TRLTOG_UNPACK',0,ZHOOK_HANDLE_UNPACK)
  CALL TGRL_INIT_PACKING_VARS(YDBUFS,YLVARS, KVSET, KF_GP)

  DO INR=1,KNRECV
    CALL TGRL_COPY_ZCOMBUF(YDBUFS, YLVARS, INR, ZCOMBUFR, KPTRGP, PGP, PGPUV, PGP3A, PGP3B, PGP2)
  ENDDO
  IF (LHOOK) CALL DR_HOOK('TRLTOG_UNPACK',1,ZHOOK_HANDLE_UNPACK)

ELSE
  ! original point-to-point path

  !  Pack+send loop.........................................................

  ISEND_FLD_START = 1
  ISEND_FLD_END   = KF_FS
  DO INS=1,KNSEND
    ISEND=KSEND(INS)
    ILEN = KSENDTOT(ISEND)/KF_FS
    I_KNDOFF_ISEND = KNDOFF(ISEND)
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,JL,II)
    DO JL=1,ILEN
      II = KINDEX(I_KNDOFF_ISEND+JL)
      DO JFLD=ISEND_FLD_START,ISEND_FLD_END
        ZCOMBUFS((JFLD-ISEND_FLD_START)*ILEN+JL,INS) = PGLAT(JFLD,II)
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    ZCOMBUFS(-1,INS) = 1
    ZCOMBUFS(0,INS)  = KF_FS
    IF (NTRANS_SYNC_LEVEL <= 1) THEN
      CALL MPL_SEND(ZCOMBUFS(-1:KSENDTOT(ISEND),INS), KDEST=NPRCIDS(ISEND), &
        &           KMP_TYPE=JP_NON_BLOCKING_STANDARD, KREQUEST=IREQ_SEND(INS), KTAG=ITAG, &
        &           CDSTRING='TRLTOG_COMM: NON-BLOCKING ISEND')
    ELSE
      CALL MPL_SEND(ZCOMBUFS(-1:KSENDTOT(ISEND),INS), KDEST=NPRCIDS(ISEND), &
        &           KMP_TYPE=JP_BLOCKING_BUFFERED, KTAG=ITAG, &
        &           CDSTRING='TRLTOG_COMM: BLOCKING BUFFERED BSEND')
    ENDIF
  ENDDO

  !  Unpack loop.........................................................

  CALL TGRL_INIT_PACKING_VARS(YDBUFS,YLVARS, KVSET, KF_GP)

  DO JNR=1,KNRECV

    IF (NTRANS_SYNC_LEVEL <= 0) THEN
      CALL MPL_WAITANY(KREQUEST=IREQ_RECV(1:KNRECV), KINDEX=INR, &
        &              CDSTRING='TRLTOG_COMM: WAIT FOR ANY RECEIVES')
    ELSE
      INR = JNR
      IRECV=KRECV(INR)
      CALL MPL_RECV(ZCOMBUFR(-1:KRECVTOT(IRECV),INR), KSOURCE=NPRCIDS(IRECV), &
            & KMP_TYPE=JP_BLOCKING_STANDARD, KTAG=ITAG, CDSTRING='TRLTOG_COMM: BLOCKING RECV')
    ENDIF

    CALL TGRL_COPY_ZCOMBUF(YDBUFS, YLVARS, INR, ZCOMBUFR, KPTRGP, PGP, PGPUV, PGP3A, PGP3B, PGP2)
  ENDDO

  IF (NTRANS_SYNC_LEVEL <= 1) THEN
    IF(KNSEND > 0) THEN
      CALL MPL_WAIT(KREQUEST=IREQ_SEND(1:KNSEND),CDSTRING='TRLTOG_COMM: WAIT FOR ISENDS')
    ENDIF
  ENDIF

  IF (NTRANS_SYNC_LEVEL >= 1) THEN
    CALL MPL_BARRIER(CDSTRING='TRLTOG_COMM: BARRIER AT END')
  ENDIF
ENDIF

CALL GSTATS(805,1)

! LREPORT_TR_BW: close the timer and account this rank's payload. KSENDTOT is
! already the element count sent to a peer (= latitudes*KF_FS, cf. ISENDTOT =
! IPOS*KF_FS in trgl_mod, and the MPL_SEND uses ZCOMBUFS(-1:KSENDTOT(ISEND))).
! So bytes = sum over the KNSEND peers of KSENDTOT(KSEND(INS)) * sizeof(JPRB).
IF (LREPORT_TR_BW) THEN
  CALL SYSTEM_CLOCK(COUNT=ITR_CLK1)
  ZTR_DT = REAL(ITR_CLK1 - ITR_CLK0, JPRD) / REAL(ITR_CLK_RATE, JPRD)
  ITR_BYTES = 0_JPIB
  DO INS=1,KNSEND
    ITR_BYTES = ITR_BYTES + INT(KSENDTOT(KSEND(INS)), JPIB)
  ENDDO
  ITR_BYTES = ITR_BYTES * INT(STORAGE_SIZE(1.0_JPRB)/8, JPIB)
  CALL ACCOUNT_TR_BW(1, 'TRLTOG (L->G, g157)', ITR_BYTES, ZTR_DT)
  ! Comm exchange only (LEVS/global alltoallv), same payload.
  IF (NPROC > 1) CALL ACCOUNT_TR_BW(3, 'TRLTOG-COMM (g157 issue)', ITR_BYTES, ZCM_DT)
  ! Opt3 decomposition: WAIT = comm completion not hidden by the self-copy;
  ! SELF = the deferred self-copy overlapped against the alltoallv. Only the
  ! time field is meaningful (payload is the whole-transpose payload).
  IF (NPROC > 1) CALL ACCOUNT_TR_BW(5, 'TRLTOG-WAIT (g157 unhidden comm)', ITR_BYTES, ZWT_DT)
  IF (NPROC > 1) CALL ACCOUNT_TR_BW(7, 'TRLTOG-SELF (g157 overlap copy)', ITR_BYTES, ZSC_DT)
  ! Pre-issue barrier = load-imbalance skew feeding the alltoallv (time only).
  IF (NPROC > 1) CALL ACCOUNT_TR_BW(9, 'TRLTOG-SKEW (g157 pre-issue barrier)', ITR_BYTES, ZBR_DT)
  ! World-gather the per-rank barrier wait for straggler profiling (opt-in).
  IF (NPROC > 1) CALL DUMP_TR_SKEW(9, ZBR_DT)
ENDIF

CALL GSTATS_BARRIER2(762)

END ASSOCIATE

END SUBROUTINE TRLTOG_COMM

END MODULE TRLTOG_MOD
