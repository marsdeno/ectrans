! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTDIR_CTL_MOD
CONTAINS
SUBROUTINE LTDIR_CTL(KF_FS,KF_UV,KF_SCALARS, &
 & PSPVOR,PSPDIV,PSPSCALAR, &
 & PSPSC3A,PSPSC3B,PSPSC2, &
 & KFLDPTRUV,KFLDPTRSC)

!**** *LTDIR_CTL* - Control routine for direct Legendre transform

!     Purpose.
!     --------
!        Direct Legendre transform

!**   Interface.
!     ----------
!     CALL LTDIR_CTL(...)

!     Explicit arguments :
!     --------------------
!     KF_FS      - number of fields in Fourier space
!     KF_UV      - local number of spectral u-v fields
!     KF_SCALARS - local number of scalar spectral fields
!     PSPVOR(:,:) - spectral vorticity (output)
!     PSPDIV(:,:) - spectral divergence (output)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (output)
!     KFLDPTRUV(:) - field pointer for vorticity and divergence (input)
!     KFLDPTRSC(:) - field pointer for scalarvalued fields (input)

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPRD

USE TPM_GEN         ,ONLY : LALLOPERM
USE TPM_TRANS       ,ONLY : FOUBUF, FOUBUF_IN
USE TPM_DISTR       ,ONLY : D, MYPROC
USE TPM_ECTRANS_OPTS, ONLY : LREPORT_FLT_TIME, LREPORT_STAGE_TIME
USE OMP_LIB         ,ONLY : OMP_GET_WTIME, OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS

USE LTDIR_MOD       ,ONLY : LTDIR
USE TRLTOM_MOD      ,ONLY : TRLTOM
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS,KF_UV,KF_SCALARS
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSC3B(:,:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSC2(:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRSC(:)

INTEGER(KIND=JPIM) :: JM,IM,IBLEN,ILED2
REAL(KIND=JPRD) :: ZT_S153, ZT_S103, ZT_S1645
REAL(KIND=JPRD), EXTERNAL :: TIMEF
REAL(KIND=JPRD), SAVE :: ZT_ACC_153 = 0.0_JPRD, ZT_ACC_103 = 0.0_JPRD
REAL(KIND=JPRD), SAVE :: ZT_ACC_1645 = 0.0_JPRD
INTEGER(KIND=JPIM), SAVE :: ZT_NCALLS = 0

! LREPORT_FLT_TIME: per-OMP-thread busy-time breakdown for the per-wavenumber
! LTDIR dispatch (SCHEDULE(DYNAMIC,1) over D%NUMP) -- see ltinv_ctl_mod.F90
! for the symmetric instrumentation on the inverse-transform side.
REAL(KIND=JPRD) :: ZFLT_WALL0, ZFLT_WALL1, ZFLT_T0
REAL(KIND=JPRD), ALLOCATABLE, SAVE :: ZFLT_THREAD_BUSY(:)
INTEGER(KIND=JPIM) :: IFLT_TID, IFLT_NTHREADS
INTEGER(KIND=JPIM), SAVE :: IFLT_CALLS_REPORTED = 0

!     ------------------------------------------------------------------

! Transposition from Fourier space distribution to spectral space distribution

IBLEN = D%NLENGT0B*2*KF_FS
IF (ALLOCATED(FOUBUF)) THEN
  IF (MAX(1,IBLEN) > SIZE(FOUBUF)) THEN
    DEALLOCATE(FOUBUF)
    ALLOCATE(FOUBUF(MAX(1,IBLEN)))
  ENDIF
ELSE
  ALLOCATE(FOUBUF(MAX(1,IBLEN)))
ENDIF

IF (LREPORT_STAGE_TIME) ZT_S153 = TIMEF()
CALL GSTATS(153,0)
CALL TRLTOM(FOUBUF_IN,FOUBUF,2*KF_FS)
CALL GSTATS(153,1)
IF (LREPORT_STAGE_TIME) ZT_ACC_153 = ZT_ACC_153 + (TIMEF() - ZT_S153)
IF (.NOT.LALLOPERM) DEALLOCATE(FOUBUF_IN)

! Direct Legendre transform

IF (LREPORT_STAGE_TIME) ZT_S103 = TIMEF()
CALL GSTATS(103,0)
ILED2 = 2*KF_FS
IF (LREPORT_STAGE_TIME) ZT_S1645 = TIMEF()
CALL GSTATS(1645,0)
IF(KF_FS>0) THEN
  IF (LREPORT_FLT_TIME) THEN
    IFLT_NTHREADS = OMP_GET_MAX_THREADS()
    IF (.NOT. ALLOCATED(ZFLT_THREAD_BUSY)) THEN
      ALLOCATE(ZFLT_THREAD_BUSY(0:IFLT_NTHREADS-1))
    ELSEIF (SIZE(ZFLT_THREAD_BUSY) /= IFLT_NTHREADS) THEN
      DEALLOCATE(ZFLT_THREAD_BUSY)
      ALLOCATE(ZFLT_THREAD_BUSY(0:IFLT_NTHREADS-1))
    ENDIF
    ZFLT_THREAD_BUSY(:) = 0.0_JPRD
    ZFLT_WALL0 = OMP_GET_WTIME()
  ENDIF

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JM,IM,ZFLT_T0,IFLT_TID)
  DO JM=1,D%NUMP
    IF (LREPORT_FLT_TIME) ZFLT_T0 = OMP_GET_WTIME()
    IM = D%MYMS(JM)
    CALL LTDIR(IM,JM,KF_FS,KF_UV,KF_SCALARS,ILED2, &
     & PSPVOR,PSPDIV,PSPSCALAR,&
     & PSPSC3A,PSPSC3B,PSPSC2 , &
     & KFLDPTRUV,KFLDPTRSC)
    IF (LREPORT_FLT_TIME) THEN
      IFLT_TID = OMP_GET_THREAD_NUM()
      ZFLT_THREAD_BUSY(IFLT_TID) = ZFLT_THREAD_BUSY(IFLT_TID) + (OMP_GET_WTIME() - ZFLT_T0)
    ENDIF
  ENDDO
!$OMP END PARALLEL DO
  IF (LREPORT_FLT_TIME) THEN
    ZFLT_WALL1 = OMP_GET_WTIME()
    IF (MYPROC == 1 .AND. IFLT_CALLS_REPORTED <= 10) THEN
      IFLT_CALLS_REPORTED = IFLT_CALLS_REPORTED + 1
      WRITE(0,'(A,I0,A,I0,A,F10.6,A,F10.6,A,F10.6,A,F10.6,A,F6.2)') &
       & 'LTDIR FLT: NUMP=', D%NUMP, ' nthreads=', IFLT_NTHREADS, &
       & ' wall(s)=', ZFLT_WALL1-ZFLT_WALL0, &
       & ' busy_min=', MINVAL(ZFLT_THREAD_BUSY(0:IFLT_NTHREADS-1)), &
       & ' busy_max=', MAXVAL(ZFLT_THREAD_BUSY(0:IFLT_NTHREADS-1)), &
       & ' busy_avg=', SUM(ZFLT_THREAD_BUSY(0:IFLT_NTHREADS-1))/REAL(IFLT_NTHREADS,JPRD), &
       & ' imbalance(max/avg)=', MAXVAL(ZFLT_THREAD_BUSY(0:IFLT_NTHREADS-1)) / &
       &    MAX(1.0E-9_JPRD, SUM(ZFLT_THREAD_BUSY(0:IFLT_NTHREADS-1))/REAL(IFLT_NTHREADS,JPRD))
    ENDIF
  ENDIF
ENDIF
CALL GSTATS(1645,1)
IF (LREPORT_STAGE_TIME) ZT_ACC_1645 = ZT_ACC_1645 + (TIMEF() - ZT_S1645)

IF (.NOT.LALLOPERM) DEALLOCATE(FOUBUF)
CALL GSTATS(103,1)
IF (LREPORT_STAGE_TIME) THEN
  ZT_ACC_103 = ZT_ACC_103 + (TIMEF() - ZT_S103)
  ZT_NCALLS = ZT_NCALLS + 1
  IF (MYPROC == 1 .AND. ZT_NCALLS <= 10) THEN
    WRITE(0,'(A,F12.3,A,F12.3,A,F12.3)') 'LTDIR wall: g153=', ZT_ACC_153/ZT_NCALLS, &
     & 'g103=', ZT_ACC_103/ZT_NCALLS, 'g1645=', ZT_ACC_1645/ZT_NCALLS
  ENDIF
ENDIF

!     -----------------------------------------------------------------

END SUBROUTINE LTDIR_CTL
END MODULE LTDIR_CTL_MOD
