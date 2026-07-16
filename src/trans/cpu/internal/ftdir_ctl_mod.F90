! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTDIR_CTL_MOD
USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

! persistent ZGTF buffer (grow-only flat 1D, pointer-remapped to 2D).
REAL(KIND=JPRB), ALLOCATABLE, TARGET, SAVE :: ZGTF_PERSIST_FLAT(:)
INTEGER(KIND=JPIM), SAVE :: ZGTF_PERSIST_CAP = 0

CONTAINS
SUBROUTINE FTDIR_CTL(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS, &
 & KVSETUV,KVSETSC,KPTRGP,&
 & KVSETSC3A,KVSETSC3B,KVSETSC2,&
 & PGP,PGPUV,PGP3A,PGP3B,PGP2)


!**** *FTDIR_CTL - Direct Fourier transform control

!     Purpose. Control routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!     CALL FTDIR_CTL(..)

!     Explicit arguments :
!     --------------------
!     KF_UV_G      - global number of spectral u-v fields
!     KF_SCALARS_G - global number of scalar spectral fields
!     KF_GP        - total number of output gridpoint fields
!     KF_FS        - total number of fields in fourier space
!     PGP     -  gridpoint array
!     KVSETUV - "B" set in spectral/fourier space for
!                u and v variables
!     KVSETSC - "B" set in spectral/fourier space for
!                scalar variables
!     KPTRGP  -  pointer array to fields in gridpoint space

!     Method.
!     -------

!     Externals.  TRGTOL      - transposition routine
!     ----------  FOURIER_OUT - copy fourier data to Fourier buffer
!                 FTDIR       - fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        R. El Khatib 09-Sep-2020 NSTACK_MEMORY_TR
!      R. El Khatib 01-Jun-2022 contiguous pointer
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPIB     ,JPRB

USE TPM_GEN   ,ONLY : NSTACK_MEMORY_TR
USE TPM_TRANS       ,ONLY : FOUBUF_IN
USE TPM_DISTR       ,ONLY : D, MYPROC, MYSETW, NPROC
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_ECTRANS_OPTS,ONLY : LUSE_OPT6
USE TPM_DIM         ,ONLY : R
USE TPM_FFTW        ,ONLY : TW, CREATE_PLAN_FFTW, EXEC_FFTW_R2C_TO_FOUBUF, N_FFT_BLK

USE TRGTOL_MOD      ,ONLY : TRGTOL
USE FOURIER_OUT_MOD ,ONLY : FOURIER_OUT
USE FTDIR_MOD       ,ONLY : FTDIR
!

IMPLICIT NONE

! Dummy arguments

INTEGER(KIND=JPIM),INTENT(IN) :: KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPTRGP(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3A(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3B(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC2(:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP2(:,:,:)

! Local variables
REAL(KIND=JPRB),TARGET  :: ZGTF_STACK(KF_FS*MIN(1,MAX(0,NSTACK_MEMORY_TR)),D%NLENGTF)
REAL(KIND=JPRB),POINTER, CONTIGUOUS :: ZGTF(:,:)

INTEGER(KIND=JPIM) :: IST,JGL,IBLEN
INTEGER(KIND=JPIM) :: IVSETUV(KF_UV_G)
INTEGER(KIND=JPIM) :: IVSETSC(KF_SCALARS_G)
INTEGER(KIND=JPIM) :: IVSET(KF_GP)
INTEGER(KIND=JPIM) :: IFGP2,IFGP3A,IFGP3B,IOFF,J3
INTEGER(KIND=JPIM) :: IGLG,IRLEN,IPLAN_KLOT
INTEGER(KIND=JPIM) :: ICLEN,IKOFF  ! OPT6: fused-FFT args
INTEGER(KIND=JPIB) :: IDUMMY_PLAN  ! OPT6: for warming plan cache
INTEGER(KIND=JPIM) :: INEED, J  ! persistent ZGTF buffer sizing
INTEGER(KIND=JPIB),ALLOCATABLE :: IPLAN_FFTW(:)

!     ------------------------------------------------------------------

! Field distribution in Spectral/Fourier space

IF(PRESENT(KVSETUV)) THEN
  IVSETUV(:) = KVSETUV(:)
ELSE
  IVSETUV(:) = -1
ENDIF
IVSETSC(:) = -1
IF(PRESENT(KVSETSC)) THEN
  IVSETSC(:) = KVSETSC(:)
ELSE
  IOFF=0
  IF(PRESENT(KVSETSC2)) THEN
    IFGP2=UBOUND(KVSETSC2,1)
    IVSETSC(1:IFGP2)=KVSETSC2(:)
    IOFF=IOFF+IFGP2
  ENDIF
  IF(PRESENT(KVSETSC3A)) THEN
    IFGP3A=UBOUND(KVSETSC3A,1)
    DO J3=1,UBOUND(PGP3A,3)
      IVSETSC(IOFF+1:IOFF+IFGP3A)=KVSETSC3A(:)
      IOFF=IOFF+IFGP3A
    ENDDO
  ENDIF
  IF(PRESENT(KVSETSC3B)) THEN
    IFGP3B=UBOUND(KVSETSC3B,1)
    DO J3=1,UBOUND(PGP3B,3)
      IVSETSC(IOFF+1:IOFF+IFGP3B)=KVSETSC3B(:)
      IOFF=IOFF+IFGP3B
    ENDDO
  ENDIF
ENDIF

IST = 1
IF(KF_UV_G > 0) THEN
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
ENDIF
IF(KF_SCALARS_G > 0) THEN
  IVSET(IST:IST+KF_SCALARS_G-1) = IVSETSC(:)
  IST = IST+KF_SCALARS_G
ENDIF

IF (NSTACK_MEMORY_TR == 1) THEN
  ZGTF => ZGTF_STACK(:,:)
ELSE
  ! persistent flat buffer + 2D pointer rank-remap, grows to maxiumum needed size
  INEED = MAX(1, KF_FS) * MAX(1, D%NLENGTF)
  IF (.NOT.ALLOCATED(ZGTF_PERSIST_FLAT) .OR. INEED > ZGTF_PERSIST_CAP) THEN
    IF (ALLOCATED(ZGTF_PERSIST_FLAT)) DEALLOCATE(ZGTF_PERSIST_FLAT)
    ALLOCATE(ZGTF_PERSIST_FLAT(INEED))
    ZGTF_PERSIST_CAP = INEED
    !$OMP PARALLEL DO SCHEDULE(STATIC)
    DO J = 1, INEED
      ZGTF_PERSIST_FLAT(J) = 0.0_JPRB
    ENDDO
    !$OMP END PARALLEL DO
  ENDIF
  ZGTF(1:MAX(1,KF_FS), 1:MAX(1,D%NLENGTF)) => ZGTF_PERSIST_FLAT(1:INEED)
ENDIF

! Transposition

CALL GSTATS(158,0)
CALL TRGTOL(ZGTF,KF_FS,KF_GP,KF_SCALARS_G,IVSET,KPTRGP,&
 &PGP,PGPUV,PGP3A,PGP3B,PGP2)
CALL GSTATS(158,1)
CALL GSTATS(106,0)

! Fourier transform

IBLEN=D%NLENGT0B*2*KF_FS
IF (ALLOCATED(FOUBUF_IN)) THEN
  IF (MAX(1,IBLEN) > SIZE(FOUBUF_IN)) THEN
    DEALLOCATE(FOUBUF_IN)
    ALLOCATE(FOUBUF_IN(MAX(1,IBLEN)))
  ENDIF
ELSE
  ALLOCATE(FOUBUF_IN(MAX(1,IBLEN)))
ENDIF

CALL GSTATS(1640, 0)
! If this rank has any Fourier fields, Fourier transform them
IF (KF_FS > 0) THEN
  ! pre-resolve FFTW plan IDs serially before the OMP loop.
  ! removes !$OMP CRITICAL(FFTW_CREATE) from the parallel FFT hot path.
  ALLOCATE(IPLAN_FFTW(D%NDGL_FS))
  IPLAN_FFTW(:) = 0_JPIB
  IPLAN_KLOT = MERGE(KF_FS, 1, TW%LALL_FFTW)
  DO JGL = 1, D%NDGL_FS
    IGLG = D%NPTRLS(MYSETW) + JGL - 1
    IRLEN = G%NLOEN(IGLG) + R%NNOEXTZL
    IF (G%NLOEN(IGLG) > 1) THEN
      CALL CREATE_PLAN_FFTW(IPLAN_FFTW(JGL), -1, IRLEN, IPLAN_KLOT)
      ! OPT6: warm the plan cache serially for the two KLOTs used inside
      ! EXEC_FFTW_R2C_TO_FOUBUF (block + tail). FFTW plan creation is not
      ! thread-safe; creating them here avoids races in the parallel loop.
      IF (LUSE_OPT6) THEN
        IF (KF_FS >= N_FFT_BLK) THEN
          CALL CREATE_PLAN_FFTW(IDUMMY_PLAN, -1, IRLEN, N_FFT_BLK)
        ENDIF
        IF (MOD(KF_FS, N_FFT_BLK) /= 0 .OR. KF_FS < N_FFT_BLK) THEN
          CALL CREATE_PLAN_FFTW(IDUMMY_PLAN, -1, IRLEN, 1_JPIM)
        ENDIF
      ENDIF
    ENDIF
  ENDDO

  ! Loop over latitudes
  ! OPT6 : fused FFT + FOURIER_OUT writing straight into FOUBUF_IN, replacing
  ! the FTDIR + FOURIER_OUT split when NLOEN>1.   NLOEN==1 latitudes fall back
  ! to the split pipeline. 
  ! Original code unchanged when off (default)
  !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JGL, IGLG, IRLEN, ICLEN, IKOFF)
  DO JGL = 1, D%NDGL_FS
    IGLG = D%NPTRLS(MYSETW) + JGL - 1
    IF (LUSE_OPT6 .AND. G%NLOEN(IGLG) > 1) THEN
      IRLEN = G%NLOEN(IGLG) + R%NNOEXTZL
      ICLEN = (IRLEN/2 + 1) * 2
      IKOFF = D%NSTAGTF(JGL) + 1
      CALL EXEC_FFTW_R2C_TO_FOUBUF(IRLEN, ICLEN, IKOFF, KF_FS, ZGTF, JGL)
    ELSE
      ! Fourier transform
      CALL FTDIR(ZGTF, KF_FS, JGL, IPLAN_FFTW(JGL))

      ! Save Fourier data in FOUBUF_IN
      CALL FOURIER_OUT(ZGTF, KF_FS, JGL)
    ENDIF
  ENDDO
  !$OMP END PARALLEL DO

  DEALLOCATE(IPLAN_FFTW)
ENDIF
CALL GSTATS(1640, 1)

CALL GSTATS(106,1)

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR_CTL
END MODULE FTDIR_CTL_MOD



