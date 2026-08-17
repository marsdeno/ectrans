! (C) Copyright 2026- ECMWF.
! (C) Copyright 2026- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TPM_ECTRANS_OPTS

! Runtime toggles for ecTrans CPU-path optimisations. Each flag defaults
! to .TRUE. (optimisation enabled). Setting the corresponding
! environment variable ECTRANS_DISABLE_OPT<n>=1 before the run switches
! to the pre-optimisation code path in the same binary. 
!
! Flags initialised by INIT_ECTRANS_OPTS, invoked from SETUP_TRANS (CPU-side)
!
! All four OPT flags are declared in this commit even if implementation comes in later commit.

IMPLICIT NONE

PRIVATE
PUBLIC :: INIT_ECTRANS_OPTS
PUBLIC :: LUSE_OPT1, LUSE_OPT2, LUSE_OPT3, LUSE_OPT4, LUSE_OPT5, LUSE_OPT6
PUBLIC :: LUSE_OPT7, LUSE_OPT8
PUBLIC :: LUSE_LEVS_TR, LVERIFY_LEVS_NODE
PUBLIC :: LREPORT_TR_BW, ZTR_BW_ROOFLINE, LDUMP_TR_SKEW
PUBLIC :: LUSE_RECTANGULAR_DECOMP
PUBLIC :: LREPORT_FLT_TIME
PUBLIC :: LREPORT_STAGE_TIME

SAVE

! FTINV: fused FOURIER_IN + FSC per-JM copy (FOURIER_IN_MOD / FTINV_CTL_MOD).
LOGICAL :: LUSE_OPT1  = .TRUE.

! TRLTOG: single MPL_ALLTOALLV + fused self-copy/pack (TRLTOG_MOD).
LOGICAL :: LUSE_OPT2  = .TRUE.

! TRLTOG: non-blocking alltoallv overlapping the deferred self-copy.
! Only meaningful when LUSE_OPT2 is also .TRUE.; forced OFF otherwise.
LOGICAL :: LUSE_OPT3 = .TRUE.

! LTINV: fused LEINV + ASRE1B (LEINV_ASRE_FUSED_MOD, called from LTINV_MOD).
LOGICAL :: LUSE_OPT4 = .TRUE.

! FFT (tpm_fftw EXEC_FFTW LD_ALL=FALSE): KLOT=N_FFT_BLK batched block phase.
! Default .FALSE. because the batched FFTW plan is not bit-id with non-batched
! Enable by setting env variable ECTRANS_ENABLE_OPT5=1
LOGICAL :: LUSE_OPT5    = .FALSE.

! FTDIR (ftdir_ctl_mod): fused FFT + FOURIER_OUT (EXEC_FFTW_R2C_TO_FOUBUF).
! Skips one PREEL write + one PREEL read + a zero-fill of dead PREEL columns
! by unpacking FFT output straight from per-thread ZFFT into FOUBUF_IN.
! Default .FALSE. because the batched block phase (shared with OPT5)
! is not bit-identical to the KLOT=1 baseline
! Enable with ECTRANS_ENABLE_OPT6=1.
LOGICAL :: LUSE_OPT6   = .FALSE.
! OPT7: TRLTOG: hoist KNDOFF(ISEND) out of the JL-parallel send-pack loop
! (both the LUSE_OPT2 single-A2A branch and the legacy per-peer branch), so
! each thread reads a scalar rather than re-dereferencing YDBUFS%INDOFF on
! every JL iteration. Bit-identical.
! Enable ECTRANS_DISABLE_OPT7=1 to revert to the un-hoisted form.
LOGICAL :: LUSE_OPT7 = .TRUE.
! OPT8: LTINV: per-kernel DR_HOOK regions inside LTINV_MOD (PRFI1B/VDTUV/
! SPNSDE/LEINV/ASRE1B) so drhook.prof profiles attribute cost to individual
! routines instead of the enclosing LTINV_MOD region. Instrumentation-only:
! no functional or performance change when DR_HOOK is disabled (LHOOK is
! compile-time constant folded false in production builds without DR_HOOK=1
! at runtime), and no change to the numerics in either case.
! Enable ECTRANS_DISABLE_OPT8=1 to omit the extra regions.
LOGICAL :: LUSE_OPT8 = .TRUE.
! TRGTOL/TRLTOG: route the G<->L transpose alltoallv through the 32-rank
! MPL_ALL_LEVS_COMM (fixed wave-set, varying V-set) instead of the global
! world communicator. VALID ONLY when the grid-point collars are aligned
! 1:1 with the Fourier wave-sets (Option A: N_REGIONS_NS==NPRTRW and every
! N_REGIONS(:)==NPRTRV). Under that alignment all non-zero G<->L traffic is
! confined to the LEVS sub-communicator, so the collective no longer
! synchronises across all NPROC ranks. The call sites AND this flag with a
! runtime alignment predicate, so enabling it on a non-aligned decomposition
! falls back to the world-comm path.
! OPT-IN: default .FALSE. Enable with ECTRANS_ENABLE_LEVS_TR=1.
LOGICAL :: LUSE_LEVS_TR = .FALSE.

! Diagnostic: on the first TRGTOL, gather hostnames within MPL_ALL_LEVS_COMM
! and report whether the 32-rank group is co-located on a single node (the
! precondition for LUSE_LEVS_TR to yield intra-node communication). Prints
! once per LEVS-comm root. Enable with ECTRANS_VERIFY_LEVS_NODE=1.
LOGICAL :: LVERIFY_LEVS_NODE = .FALSE.

! Diagnostic: instrument the G<->L transpose comm regions (TRLTOG L->G and
! TRGTOL G->L) with a wall-clock timer and a payload byte counter, then
! report achieved node bandwidth (GB/s) and its fraction of a configurable
! memory-bandwidth roofline. Accumulates over a fixed count of steady-state
! calls (skipping warmup) and prints once from the LEVS-comm root.
! Enable with ECTRANS_REPORT_TR_BW=1.
LOGICAL :: LREPORT_TR_BW = .FALSE.

! Roofline denominator (GB/s per node) for the LREPORT_TR_BW %-of-peak
! figure. Defaults to 350.0 (roughly STREAM for 2x EPYC 7742).
! Override with ECTRANS_TR_BW_ROOFLINE=<real>.
REAL :: ZTR_BW_ROOFLINE = 350.0

! Straggler profiling: when set (and LREPORT_TR_BW on), the G<->L transpose
! pre-issue barrier all-gathers every rank's barrier-wait over the world comm
! for the sample window and dumps it to tr_skew.csv (rank 1 only), so the skew
! can be attributed to systematic slow ranks/nodes vs random per-step jitter.
! Enable with ECTRANS_DUMP_TR_SKEW=1.
LOGICAL :: LDUMP_TR_SKEW = .FALSE.

! Grid decomposition selector. DEFAULT .FALSE. keeps the stock Leopardi
! equal-regions partition (variable sectors per collar, polar caps) with the
! standard SUMPLATBEQ latitude split -- the traditional, always-supported
! decomposition used throughout the rest of IFS. Setting
! ECTRANS_RECTANGULAR_DECOMP=1 opts in to a "rectangular" decomposition
! instead: SUMP_TRANS0 overrides the equal-regions output so that the NS
! collar count equals NPRTRW and every collar holds exactly NPROC/NPRTRW=
! NPRTRV ranks, and SUMPLAT cuts grid latitudes on the Fourier wave-set
! boundaries (LUSE_FOURIER_BOUNDARIES). That alignment is the precondition
! for LUSE_LEVS_TR (on-node transpose); it is opt-in because it changes grid
! load balance vs the native decomposition and is only useful when actually
! exercising LUSE_LEVS_TR.
! Known issue : for at least some truncation/NPRTRW combinations the Fourier-
! boundary-derived collars come out degenerate (an empty collar absorbing
! zero grid points, with another collar absorbing the rest), which aborts
! in SUSTAONL when LDSPLIT=.TRUE. and silently leaves a rank with zero grid
! points otherwise. Not yet fixed -- treat as experimental.
! Cached here by INIT_ECTRANS_OPTS, called from SETUP_TRANS0 (before
! SUMP_TRANS0 runs) so the flag is available in time for both SUMP_TRANS0
! and the later SUMP_TRANS call sites.
LOGICAL :: LUSE_RECTANGULAR_DECOMP = .FALSE.

! Diagnostic: instrument the FLT/butterfly (MULT_BUTM) OMP-parallel dispatch
! in LTINV_CTL/LTDIR_CTL (the per-wavenumber SCHEDULE(DYNAMIC,1) loop -- the
! active call path in this build, since HAVE_BATCHED_BLAS is off; the
! LEINV_BATCHED/LEINVAD_BATCHED equivalents are only reachable when that
! macro is defined) with wall-clock timing of the parallel region plus a
! per-OMP-thread busy-time breakdown (min/max/avg across threads), to
! directly quantify load imbalance across the coarse per-wavenumber task
! dispatch (each butterfly-eligible wavenumber's whole MULT_BUTM call runs
! on a single thread with no internal parallelism).
! Prints once per rank-1 call (first few steady-state calls only).
! Only meaningful when S%LUSEFLT is on (--flt / LDUSEFLT=.TRUE.).
! Enable with ECTRANS_REPORT_FLT_TIME=1.
LOGICAL :: LREPORT_FLT_TIME = .FALSE.

! Diagnostic: per-stage wall-clock timing breakdown (TIMEF, not GSTATS) of the
! transform stages in FTDIR_CTL/FTINV_CTL/LTDIR_CTL/LTINV_CTL and their
! LTDIR_BATCHED_CTL/LTINV_BATCHED_CTL counterparts (transpose / Legendre /
! Fourier regions), accumulated and averaged over the first 10 calls on rank
! 1 only. Also gates FTDIR_CTL's one-time rank-1 dump of the grid
! decomposition (grid_decomp.csv: latitude, NLOEN, owning wave-set), which is
! informational context for the same report.
! Enable with ECTRANS_REPORT_STAGE_TIME=1.
LOGICAL :: LREPORT_STAGE_TIME = .FALSE.

LOGICAL, PRIVATE :: LINITIALISED = .FALSE.

CONTAINS

SUBROUTINE INIT_ECTRANS_OPTS()

! Read all optimisation env vars once and cache them as module LOGICALs

USE TPM_GEN, ONLY : NOUT, NPRINTLEV

IF (LINITIALISED) RETURN

LUSE_OPT1  = .NOT. ENV_DISABLE('ECTRANS_DISABLE_OPT1')
LUSE_OPT2  = .NOT. ENV_DISABLE('ECTRANS_DISABLE_OPT2')
! OPT3 is only wired inside the LUSE_OPT2 branch of trltog_mod,
! so disabling OPT2 implicitly disables OPT3 too.
LUSE_OPT3 = LUSE_OPT2 .AND. .NOT. ENV_DISABLE('ECTRANS_DISABLE_OPT3')
LUSE_OPT4 = .NOT. ENV_DISABLE('ECTRANS_DISABLE_OPT4')
LUSE_OPT5   = ENV_ENABLE('ECTRANS_ENABLE_OPT5')
LUSE_OPT6  = ENV_ENABLE('ECTRANS_ENABLE_OPT6')
LUSE_OPT7 = .NOT. ENV_DISABLE('ECTRANS_DISABLE_OPT7')
LUSE_OPT8 = .NOT. ENV_DISABLE('ECTRANS_DISABLE_OPT8')
LUSE_LEVS_TR       = ENV_ENABLE('ECTRANS_ENABLE_LEVS_TR')
LVERIFY_LEVS_NODE  = ENV_ENABLE('ECTRANS_VERIFY_LEVS_NODE')
LREPORT_TR_BW      = ENV_ENABLE('ECTRANS_REPORT_TR_BW')
LDUMP_TR_SKEW      = ENV_ENABLE('ECTRANS_DUMP_TR_SKEW')
LUSE_RECTANGULAR_DECOMP = ENV_ENABLE('ECTRANS_RECTANGULAR_DECOMP')
CALL ENV_REAL('ECTRANS_TR_BW_ROOFLINE', ZTR_BW_ROOFLINE)
LREPORT_FLT_TIME = ENV_ENABLE('ECTRANS_REPORT_FLT_TIME')
LREPORT_STAGE_TIME = ENV_ENABLE('ECTRANS_REPORT_STAGE_TIME')

IF (NPRINTLEV > 0) THEN
  WRITE(NOUT,'(A)')       '=== ecTrans CPU optimisation flags ==='
  WRITE(NOUT,'(A,L1,A)')  '  LUSE_OPT1  (fused FOURIER_IN+FSC)        = ', LUSE_OPT1, &
   &                      '   (ECTRANS_DISABLE_OPT1=1 to disable)'
  WRITE(NOUT,'(A,L1,A)')  '  LUSE_OPT2  (TRLTOG single alltoallv)     = ', LUSE_OPT2, &
   &                      '   (ECTRANS_DISABLE_OPT2=1 to disable)'
  WRITE(NOUT,'(A,L1,A)')  '  LUSE_OPT3 (TRLTOG nonblocking overlap)  = ', LUSE_OPT3, &
   &                      '   (ECTRANS_DISABLE_OPT3=1 to disable)'
  WRITE(NOUT,'(A,L1,A)')  '  LUSE_OPT4 (fused LEINV+ASRE1B)          = ', LUSE_OPT4, &
   &                      '   (ECTRANS_DISABLE_OPT4=1 to disable)'
  WRITE(NOUT,'(A,L1,A)')  '  LUSE_OPT5   (FFT KLOT=N_FFT_BLK block)    = ', LUSE_OPT5, &
   &                      '   (ECTRANS_ENABLE_OPT5=1 to enable; NOT bit-id)'
  WRITE(NOUT,'(A,L1,A)')  '  LUSE_OPT6  (fused FTDIR->FOUBUF)          = ', LUSE_OPT6, &
   &                      '   (ECTRANS_ENABLE_OPT6=1 to enable; NOT bit-id)'
  WRITE(NOUT,'(A,L1,A)')  '  LUSE_OPT7 (TRLTOG KNDOFF hoist)           = ', LUSE_OPT7, &
   &                      '   (ECTRANS_DISABLE_OPT7=1 to disable)'
  WRITE(NOUT,'(A,L1,A)')  '  LUSE_OPT8 (LTINV per-kernel DR_HOOK)      = ', LUSE_OPT8, &
   &                      '   (ECTRANS_DISABLE_OPT8=1 to disable)'
  WRITE(NOUT,'(A,L1,A)')  '  LUSE_LEVS_TR (G<->L transpose on LEVS comm)= ', LUSE_LEVS_TR, &
   &                      '   (ECTRANS_ENABLE_LEVS_TR=1; needs collar==waveset)'
  WRITE(NOUT,'(A,L1,A)')  '  LVERIFY_LEVS_NODE (LEVS on-node diag)      = ', LVERIFY_LEVS_NODE, &
   &                      '   (ECTRANS_VERIFY_LEVS_NODE=1 to enable)'
  WRITE(NOUT,'(A,L1,A)')  '  LREPORT_TR_BW (G<->L transpose BW report)  = ', LREPORT_TR_BW, &
   &                      '   (ECTRANS_REPORT_TR_BW=1 to enable)'
  WRITE(NOUT,'(A,F8.1,A)')'  ZTR_BW_ROOFLINE (GB/s per node)            = ', ZTR_BW_ROOFLINE, &
   &                      '   (ECTRANS_TR_BW_ROOFLINE=<real> to override)'
  WRITE(NOUT,'(A,L1,A)')  '  LDUMP_TR_SKEW (per-rank skew dump)         = ', LDUMP_TR_SKEW, &
   &                      '   (ECTRANS_DUMP_TR_SKEW=1 to enable)'
  WRITE(NOUT,'(A,L1,A)')  '  LUSE_RECTANGULAR_DECOMP (opt-in, needs LEVS_TR)= ', LUSE_RECTANGULAR_DECOMP, &
   &                      '   (ECTRANS_RECTANGULAR_DECOMP=1; default=native eq-regions)'
  WRITE(NOUT,'(A,L1,A)')  '  LREPORT_FLT_TIME (MULT_BUTM load balance)  = ', LREPORT_FLT_TIME, &
   &                      '   (ECTRANS_REPORT_FLT_TIME=1 to enable)'
  WRITE(NOUT,'(A,L1,A)')  '  LREPORT_STAGE_TIME (FTDIR/FTINV/LTDIR/LTINV)= ', LREPORT_STAGE_TIME, &
   &                      '   (ECTRANS_REPORT_STAGE_TIME=1 to enable)'
ENDIF

LINITIALISED = .TRUE.

END SUBROUTINE INIT_ECTRANS_OPTS

! -----------------------------------------------------------------------

LOGICAL FUNCTION ENV_DISABLE(CDNAME) RESULT(LLDIS)

! Return .TRUE. if env variable CDNAME is set to "1", .FALSE. otherwise
! (including if env variable CDNAME is unset)

CHARACTER(LEN=*), INTENT(IN) :: CDNAME

CHARACTER(LEN=8) :: CLVAL
INTEGER          :: ISTATUS

LLDIS = .FALSE.
CALL GET_ENVIRONMENT_VARIABLE(CDNAME, CLVAL, STATUS=ISTATUS)
IF (ISTATUS == 0) THEN
  IF (TRIM(CLVAL) == '1') LLDIS = .TRUE.
ENDIF

END FUNCTION ENV_DISABLE

! -----------------------------------------------------------------------

LOGICAL FUNCTION ENV_ENABLE(CDNAME) RESULT(LLEN)

! Return .TRUE. if env variable CDNAME is set to "1", .FALSE. otherwise
! (including if env variable CDNAME is unset)

CHARACTER(LEN=*), INTENT(IN) :: CDNAME

CHARACTER(LEN=8) :: CLVAL
INTEGER          :: ISTATUS

LLEN = .FALSE.
CALL GET_ENVIRONMENT_VARIABLE(CDNAME, CLVAL, STATUS=ISTATUS)
IF (ISTATUS == 0) THEN
  IF (TRIM(CLVAL) == '1') LLEN = .TRUE.
ENDIF

END FUNCTION ENV_ENABLE

! -----------------------------------------------------------------------

SUBROUTINE ENV_REAL(CDNAME, PVAL)

! If the environment variable CDNAME is set to a parseable real value,
! overwrite PVAL with it. Otherwise leave PVAL unchanged 

CHARACTER(LEN=*), INTENT(IN)    :: CDNAME
REAL,             INTENT(INOUT) :: PVAL

CHARACTER(LEN=32) :: CLVAL
INTEGER           :: ISTATUS, IOS
REAL              :: ZTMP

CALL GET_ENVIRONMENT_VARIABLE(CDNAME, CLVAL, STATUS=ISTATUS)
IF (ISTATUS == 0 .AND. LEN_TRIM(CLVAL) > 0) THEN
  READ(CLVAL, *, IOSTAT=IOS) ZTMP
  IF (IOS == 0) PVAL = ZTMP
ENDIF

END SUBROUTINE ENV_REAL

END MODULE TPM_ECTRANS_OPTS
