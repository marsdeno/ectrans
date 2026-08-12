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
PUBLIC :: LREPORT_FLT_TIME

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
LREPORT_FLT_TIME = ENV_ENABLE('ECTRANS_REPORT_FLT_TIME')

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
  WRITE(NOUT,'(A,L1,A)')  '  LREPORT_FLT_TIME (MULT_BUTM load balance)  = ', LREPORT_FLT_TIME, &
   &                      '   (ECTRANS_REPORT_FLT_TIME=1 to enable)'
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

END MODULE TPM_ECTRANS_OPTS
