! (C) Copyright 2026- ECMWF.
! (C) Copyright 2026- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LEINV_ASRE_FUSED_MOD
CONTAINS
SUBROUTINE LEINV_ASRE_FUSED(KM,KMLOC,KFC,KIFC,KF_OUT_LT,KSL,KDGLU,PIA)

!**** *LEINV_ASRE_FUSED* - Fused inverse Legendre transform + ASRE1B recombination.
!
!     Purpose.
!     --------
!        Fuses LEINV with ASRE1B so the transposed intermediate arrays 
!        PAOA1/PSOA1 are eliminated. The two ZC output buffers of the
!        GEMMs feed directly into FOUBUF_IN using the same ISTAN/ISTAS 
!        offset scheme that ASRE1B computes.
!
!     Interface.
!     ----------
!        CALL LEINV_ASRE_FUSED(KM,KMLOC,KFC,KIFC,KF_OUT_LT,KSL,KDGLU,PIA)
!
!        Explicit arguments :
!        -------------------
!         KM        - zonal wavenumber
!         KMLOC     - local zonal wavenumber
!         KFC       - full Fourier field count (2*KF_OUT_LT)
!         KIFC      - packed field count (KFC for KM/=0, KFC/2 for KM==0)
!         KF_OUT_LT - number of output fields
!         KSL       - starting latitude index (ISL)
!         KDGLU     - number of latitudes in this hemisphere for KM
!         PIA       - spectral coefficients
!
!     ------------------------------------------------------------------

USE PARKIND1         ,ONLY : JPRD, JPRM, JPIM     ,JPRB
USE YOMHOOK          ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM          ,ONLY : R
USE TPM_FLT          ,ONLY : S
USE TPM_TRANS        ,ONLY : FOUBUF_IN
USE TPM_DISTR        ,ONLY : D
USE BUTTERFLY_ALG_MOD,ONLY : MULT_BUTM
USE ECTRANS_BLAS_MOD ,ONLY : GEMM

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KM
INTEGER(KIND=JPIM), INTENT(IN)  :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN)  :: KFC
INTEGER(KIND=JPIM), INTENT(IN)  :: KIFC
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_OUT_LT
INTEGER(KIND=JPIM), INTENT(IN)  :: KSL
INTEGER(KIND=JPIM), INTENT(IN)  :: KDGLU
REAL(KIND=JPRB),    INTENT(IN)  :: PIA(:,:)

!     LOCAL
INTEGER(KIND=JPIM) :: IA, ILA, ILS, IS, ISKIP, ISL, IEND
INTEGER(KIND=JPIM) :: IFLD, JGL, JK, J, JI
INTEGER(KIND=JPIM) :: IGLS, IPROC, IPROCS, ISTAN, ISTAS
INTEGER(KIND=JPIM) :: ITHRESHOLD
REAL(KIND=JPRB)    :: ZBA((R%NSMAX-KM+2)/2,KIFC), ZBS((R%NSMAX-KM+3)/2,KIFC)
REAL(KIND=JPRB)    :: ZC_ANTI(KDGLU,KIFC), ZC_SYMM(KDGLU,KIFC)
REAL(KIND=JPRB)    :: ZA, ZS
CHARACTER(LEN=1) :: CLX
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

CLX = 'S'
IF (JPRB == JPRD) CLX = 'D'

ISL  = KSL
IEND = KSL + KDGLU - 1

IA  = 1+MOD(R%NSMAX-KM+2,2)
IS  = 1+MOD(R%NSMAX-KM+1,2)
ILA = (R%NSMAX-KM+2)/2
ILS = (R%NSMAX-KM+3)/2

IF (KM == 0) THEN
  ISKIP = 2
ELSE
  ISKIP = 1
ENDIF

IF (KDGLU > 0) THEN

  ITHRESHOLD = S%ITHRESHOLD

  ! 1. +++++++++++++ anti-symmetric

  IF (LHOOK) CALL DR_HOOK('LEINV_FUSED_'//CLX//'PACK_1',0,ZHOOK_HANDLE)
  IFLD = 0
  DO JK=1,KFC,ISKIP
    IFLD = IFLD + 1
    DO J=1,ILA
      ZBA(J,IFLD) = PIA(IA+1+(J-1)*2, JK)
    ENDDO
  ENDDO
  IF (LHOOK) CALL DR_HOOK('LEINV_FUSED_'//CLX//'PACK_1',1,ZHOOK_HANDLE)

  IF (ILA <= ITHRESHOLD .OR. .NOT. S%LUSEFLT) THEN
    IF (LHOOK) CALL DR_HOOK('LEINV_FUSED_'//CLX//'GEMM_1',0,ZHOOK_HANDLE)
    CALL GEMM('N','N',KDGLU,KIFC,ILA,1.0_JPRB,S%FA(KMLOC)%RPNMA,KDGLU,&
     &        ZBA,ILA,0._JPRB,ZC_ANTI,KDGLU)
    IF (LHOOK) CALL DR_HOOK('LEINV_FUSED_'//CLX//'GEMM_1',1,ZHOOK_HANDLE)
  ELSE
    IF (LHOOK) CALL DR_HOOK('LEINV_FUSED_'//CLX//'BUTM_1',0,ZHOOK_HANDLE)
    CALL MULT_BUTM('N',S%FA(KMLOC)%YBUT_STRUCT_A,KIFC,ZBA,ZC_ANTI)
    IF (LHOOK) CALL DR_HOOK('LEINV_FUSED_'//CLX//'BUTM_1',1,ZHOOK_HANDLE)
  ENDIF

  ! 2. +++++++++++++ symmetric

  IF (LHOOK) CALL DR_HOOK('LEINV_FUSED_'//CLX//'PACK_2',0,ZHOOK_HANDLE)
  IFLD = 0
  DO JK=1,KFC,ISKIP
    IFLD = IFLD + 1
    DO J=1,ILS
      ZBS(J,IFLD) = PIA(IS+1+(J-1)*2, JK)
    ENDDO
  ENDDO
  IF (LHOOK) CALL DR_HOOK('LEINV_FUSED_'//CLX//'PACK_2',1,ZHOOK_HANDLE)

  IF (ILS <= ITHRESHOLD .OR. .NOT. S%LUSEFLT) THEN
    IF (LHOOK) CALL DR_HOOK('LEINV_FUSED_'//CLX//'GEMM_2',0,ZHOOK_HANDLE)
    CALL GEMM('N','N',KDGLU,KIFC,ILS,1.0_JPRB,S%FA(KMLOC)%RPNMS,KDGLU,&
     &        ZBS,ILS,0._JPRB,ZC_SYMM,KDGLU)
    IF (LHOOK) CALL DR_HOOK('LEINV_FUSED_'//CLX//'GEMM_2',1,ZHOOK_HANDLE)
  ELSE
    IF (LHOOK) CALL DR_HOOK('LEINV_FUSED_'//CLX//'BUTM_2',0,ZHOOK_HANDLE)
    CALL MULT_BUTM('N',S%FA(KMLOC)%YBUT_STRUCT_S,KIFC,ZBS,ZC_SYMM)
    IF (LHOOK) CALL DR_HOOK('LEINV_FUSED_'//CLX//'BUTM_2',1,ZHOOK_HANDLE)
  ENDIF

  ! 3. +++++++++++++ FUSED write: transpose ZC_ANTI/ZC_SYMM and recombine
  !    directly into FOUBUF_IN. Replaces the two explicit transpose loops in
  !    LEINV (writes to PAOA1/PSOA1) plus the ASRE1B reader/adder loop.
  !    Loop order: JI outer, IFLD inner -> FOUBUF_IN writes are stride-1
  !    inside a contiguous ISTAN..ISTAN+KFC block; ZC reads are stride-KDGLU
  !    but the ZC working set fits in L2 (2*KDGLU*KIFC*8B ~ few MB)

  IF (LHOOK) CALL DR_HOOK('LEINV_FUSED_'//CLX//'WRITE',0,ZHOOK_HANDLE)

  IF (KM == 0) THEN
    ! For KM==0, only odd-indexed slots receive non-zero values (matching
    ! the pre-zero of PAOA1(2::2,:)/PSOA1(2::2,:) in reference LEINV and the
    ! JFLD=2,4,... = 0 outputs of reference ASRE1B). ISKIP=2, KIFC=KFC/2.
    DO JI=1,KDGLU
      JGL   = ISL + JI - 1
      IGLS  = R%NDGL + 1 - JGL
      IPROC  = D%NPROCL(JGL)
      IPROCS = D%NPROCL(IGLS)
      ISTAN = (D%NSTAGT0B(IPROC)  + D%NPNTGTB1(KMLOC, JGL))  * 2*KF_OUT_LT
      ISTAS = (D%NSTAGT0B(IPROCS) + D%NPNTGTB1(KMLOC, IGLS)) * 2*KF_OUT_LT
      DO IFLD=1,KIFC
        JK = 2*IFLD - 1
        ZA = ZC_ANTI(JI, IFLD)
        ZS = ZC_SYMM(JI, IFLD)
        FOUBUF_IN(ISTAN + JK  )   = ZA + ZS
        FOUBUF_IN(ISTAN + JK+1)   = 0.0_JPRB
        FOUBUF_IN(ISTAS + JK  )   = ZS - ZA
        FOUBUF_IN(ISTAS + JK+1)   = 0.0_JPRB
      ENDDO
    ENDDO
  ELSE
    ! KM/=0: ISKIP=1, KIFC=KFC, no zero-fill needed.
    DO JI=1,KDGLU
      JGL   = ISL + JI - 1
      IGLS  = R%NDGL + 1 - JGL
      IPROC  = D%NPROCL(JGL)
      IPROCS = D%NPROCL(IGLS)
      ISTAN = (D%NSTAGT0B(IPROC)  + D%NPNTGTB1(KMLOC, JGL))  * 2*KF_OUT_LT
      ISTAS = (D%NSTAGT0B(IPROCS) + D%NPNTGTB1(KMLOC, IGLS)) * 2*KF_OUT_LT
      DO IFLD=1,KIFC
        ZA = ZC_ANTI(JI, IFLD)
        ZS = ZC_SYMM(JI, IFLD)
        FOUBUF_IN(ISTAN + IFLD) = ZA + ZS
        FOUBUF_IN(ISTAS + IFLD) = ZS - ZA
      ENDDO
    ENDDO
  ENDIF

  IF (LHOOK) CALL DR_HOOK('LEINV_FUSED_'//CLX//'WRITE',1,ZHOOK_HANDLE)

ENDIF
!     ------------------------------------------------------------------

END SUBROUTINE LEINV_ASRE_FUSED
END MODULE LEINV_ASRE_FUSED_MOD
