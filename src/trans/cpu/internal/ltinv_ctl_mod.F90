! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTINV_CTL_MOD
CONTAINS
SUBROUTINE LTINV_CTL(KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,&
 & PSPVOR,PSPDIV,PSPSCALAR,&
 & PSPSC3A,PSPSC3B,PSPSC2,&
 & KFLDPTRUV,KFLDPTRSC,FSPGL_PROC)

!**** *LTINV_CTL* - Control routine for inverse Legandre transform.

!     Purpose.
!     --------
!        Control routine for the inverse LEGENDRE transform

!**   Interface.
!     ----------
!     CALL INV_TRANS_CTL(...)
!     KF_OUT_LT    - number of fields coming out from inverse LT
!     KF_UV        - local number of spectral u-v fields
!     KF_SCALARS   - local number of scalar spectral fields
!     KF_SCDERS    - local number of derivatives of scalar spectral fields
!     PSPVOR(:,:)  - spectral vorticity (input)
!     PSPDIV(:,:)  - spectral divergence (input)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
!     KFLDPTRUV(:) - field pointer array for vor./div.
!     KFLDPTRSC(:) - field pointer array for PSPSCALAR
!     FSPGL_PROC  - external procedure to be executed in fourier space
!                   before transposition

!     Method.
!     -------

!     Externals.
!     ----------
!

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-06-03

!     ------------------------------------------------------------------

USE PARKIND1   ,ONLY : JPIM     ,JPRB ,JPRD

USE TPM_GEN    ,ONLY : LALLOPERM
USE TPM_TRANS  ,ONLY : FOUBUF, FOUBUF_IN
USE TPM_DISTR  ,ONLY : D, MYPROC
USE TPM_FLT    ,ONLY : S
USE TPM_ECTRANS_OPTS, ONLY : LREPORT_FLT_TIME, LREPORT_STAGE_TIME
USE OMP_LIB    ,ONLY : OMP_GET_WTIME, OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS

USE LTINV_MOD  ,ONLY : LTINV
USE TRMTOL_MOD ,ONLY : TRMTOL

IMPLICIT NONE

INTEGER(KIND=JPIM) :: JM,IM,IBLEN,ILEI2,IDIM1
REAL(KIND=JPRD) :: ZT_S102, ZT_S1647, ZT_S152
REAL(KIND=JPRD), EXTERNAL :: TIMEF
REAL(KIND=JPRD), SAVE :: ZT_102_ACC = 0.0_JPRD, ZT_1647_ACC = 0.0_JPRD
REAL(KIND=JPRD), SAVE :: ZT_152_ACC = 0.0_JPRD
INTEGER(KIND=JPIM), SAVE :: ZT_CALLS = 0

! LREPORT_FLT_TIME: per-OMP-thread busy-time breakdown for the per-wavenumber
! LTINV dispatch (SCHEDULE(DYNAMIC,1) over D%NUMP), to quantify load
! imbalance -- most relevant when S%LUSEFLT is on, since FLT/butterfly cost
! varies strongly with wavenumber while dense-GEMM cost is comparatively flat.
REAL(KIND=JPRD) :: ZFLT_WALL0, ZFLT_WALL1, ZFLT_T0
REAL(KIND=JPRD), ALLOCATABLE, SAVE :: ZFLT_THREAD_BUSY(:)
INTEGER(KIND=JPIM) :: IFLT_TID, IFLT_NTHREADS
INTEGER(KIND=JPIM), SAVE :: IFLT_CALLS_REPORTED = 0

INTEGER(KIND=JPIM),INTENT(IN) :: KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS
REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPVOR(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPDIV(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPSCALAR(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPSC3A(:,:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPSC3B(:,:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPSC2(:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRSC(:)
EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC

!     ------------------------------------------------------------------

IF (LREPORT_STAGE_TIME) ZT_S102 = TIMEF()
CALL GSTATS(102,0)
ILEI2 = 8*KF_UV + 2*KF_SCALARS + 2*KF_SCDERS
IDIM1 = 2*KF_OUT_LT
IBLEN = D%NLENGT0B*2*KF_OUT_LT
IF (ALLOCATED(FOUBUF)) THEN
  IF (MAX(1,IBLEN) > SIZE(FOUBUF)) THEN
    DEALLOCATE(FOUBUF)
    ALLOCATE(FOUBUF(MAX(1,IBLEN)))
  ENDIF
ELSE
  ALLOCATE(FOUBUF(MAX(1,IBLEN)))
ENDIF
IF (ALLOCATED(FOUBUF_IN)) THEN
  IF (MAX(1,IBLEN) > SIZE(FOUBUF_IN)) THEN
    DEALLOCATE(FOUBUF_IN)
    ALLOCATE(FOUBUF_IN(MAX(1,IBLEN)))
  ENDIF
ELSE
  ALLOCATE(FOUBUF_IN(MAX(1,IBLEN)))
  FOUBUF_IN(:) = 0
ENDIF

! Following switch necessary when latlon grids are used with different increments in NS and EW direction.
! Otherwise unassigned values will appear in output. This is very likely a bug (ATLAS-149)
IF (S%LDLL) THEN
  FOUBUF_IN(:) = 0
ENDIF

IF(KF_OUT_LT > 0) THEN
  IF (LREPORT_STAGE_TIME) ZT_S1647 = TIMEF()
  CALL GSTATS(1647,0)

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

  !!!WARNING!!! Duplication of code besides the FSPGL_PROC argument.
              ! It seems that gfortran 10 does not retain the value
              ! of FSPGL_PROC within the OMP region.
  IF( PRESENT(FSPGL_PROC) ) THEN
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JM,IM,ZFLT_T0,IFLT_TID)
    DO JM=1,D%NUMP
      IF (LREPORT_FLT_TIME) ZFLT_T0 = OMP_GET_WTIME()
      IM = D%MYMS(JM)
      CALL LTINV(IM,JM,KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,ILEI2,IDIM1,&
       & PSPVOR,PSPDIV,PSPSCALAR ,&
       & PSPSC3A,PSPSC3B,PSPSC2 , &
       & KFLDPTRUV,KFLDPTRSC,FSPGL_PROC)
      IF (LREPORT_FLT_TIME) THEN
        IFLT_TID = OMP_GET_THREAD_NUM()
        ZFLT_THREAD_BUSY(IFLT_TID) = ZFLT_THREAD_BUSY(IFLT_TID) + (OMP_GET_WTIME() - ZFLT_T0)
      ENDIF
    ENDDO
    !$OMP END PARALLEL DO
  ELSE
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JM,IM,ZFLT_T0,IFLT_TID)
    DO JM=1,D%NUMP
      IF (LREPORT_FLT_TIME) ZFLT_T0 = OMP_GET_WTIME()
      IM = D%MYMS(JM)
      CALL LTINV(IM,JM,KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,ILEI2,IDIM1,&
       & PSPVOR,PSPDIV,PSPSCALAR ,&
       & PSPSC3A,PSPSC3B,PSPSC2 , &
       & KFLDPTRUV,KFLDPTRSC)
      IF (LREPORT_FLT_TIME) THEN
        IFLT_TID = OMP_GET_THREAD_NUM()
        ZFLT_THREAD_BUSY(IFLT_TID) = ZFLT_THREAD_BUSY(IFLT_TID) + (OMP_GET_WTIME() - ZFLT_T0)
      ENDIF
    ENDDO
    !$OMP END PARALLEL DO
  ENDIF
  IF (LREPORT_FLT_TIME) THEN
    ZFLT_WALL1 = OMP_GET_WTIME()
    IF (MYPROC == 1 .AND. IFLT_CALLS_REPORTED <= 10) THEN
      IFLT_CALLS_REPORTED = IFLT_CALLS_REPORTED + 1
      WRITE(0,'(A,I0,A,I0,A,F10.6,A,F10.6,A,F10.6,A,F10.6,A,F6.2)') &
       & 'LTINV FLT: NUMP=', D%NUMP, ' nthreads=', IFLT_NTHREADS, &
       & ' wall(s)=', ZFLT_WALL1-ZFLT_WALL0, &
       & ' busy_min=', MINVAL(ZFLT_THREAD_BUSY(0:IFLT_NTHREADS-1)), &
       & ' busy_max=', MAXVAL(ZFLT_THREAD_BUSY(0:IFLT_NTHREADS-1)), &
       & ' busy_avg=', SUM(ZFLT_THREAD_BUSY(0:IFLT_NTHREADS-1))/REAL(IFLT_NTHREADS,JPRD), &
       & ' imbalance(max/avg)=', MAXVAL(ZFLT_THREAD_BUSY(0:IFLT_NTHREADS-1)) / &
       &    MAX(1.0E-9_JPRD, SUM(ZFLT_THREAD_BUSY(0:IFLT_NTHREADS-1))/REAL(IFLT_NTHREADS,JPRD))
    ENDIF
  ENDIF
  CALL GSTATS(1647,1)
  IF (LREPORT_STAGE_TIME) ZT_1647_ACC = ZT_1647_ACC + (TIMEF() - ZT_S1647)
ENDIF

CALL GSTATS(102,1)
IF (LREPORT_STAGE_TIME) ZT_102_ACC = ZT_102_ACC + (TIMEF() - ZT_S102)

IF (LREPORT_STAGE_TIME) ZT_S152 = TIMEF()
CALL GSTATS(152,0)
CALL TRMTOL(FOUBUF_IN,FOUBUF,2*KF_OUT_LT)
CALL GSTATS(152,1)
IF (LREPORT_STAGE_TIME) THEN
  ZT_152_ACC = ZT_152_ACC + (TIMEF() - ZT_S152)
  ZT_CALLS = ZT_CALLS + 1

  IF (MYPROC == 1 .AND. ZT_CALLS <= 10) THEN
    WRITE(0,'(A,F12.3,A,F12.3,A,F12.3)') 'LTINV wall: g102=', ZT_102_ACC/ZT_CALLS, &
     & 'g1647=', ZT_1647_ACC/ZT_CALLS, 'g152=', ZT_152_ACC/ZT_CALLS
  ENDIF
ENDIF
IF (.NOT.LALLOPERM) DEALLOCATE(FOUBUF_IN)
!     ------------------------------------------------------------------

END SUBROUTINE LTINV_CTL
END MODULE LTINV_CTL_MOD
