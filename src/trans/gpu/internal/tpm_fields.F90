! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TPM_FIELDS

USE PARKIND_ECTRANS  ,ONLY : JPIM, JPRB, JPRBT, JPRD, JPRL

IMPLICIT NONE

SAVE

TYPE FIELDS_TYPE
REAL(KIND=JPRD) ,ALLOCATABLE :: RPNM(:,:) ! Legendre polynomials
REAL(KIND=JPRD) ,ALLOCATABLE :: RMU(:)    ! sin(theta) for Gaussian latitudes
REAL(KIND=JPRBT) ,ALLOCATABLE :: RW(:)     ! Weights of the Gaussian quadrature
REAL(KIND=JPRBT) ,ALLOCATABLE :: R1MU2(:)  ! 1.-MU*MU, cos(theta)**2
REAL(KIND=JPRBT) ,ALLOCATABLE :: RACTHE(:) ! 1./SQRT(R1MU2), 1/(cos(theta))

REAL(KIND=JPRBT) ,ALLOCATABLE :: REPSNM(:) ! eps(n,m) used in the Legendre transforms
REAL(KIND=JPRBT) ,ALLOCATABLE :: RN(:)     ! n (to avoid integer to real conversion)
REAL(KIND=JPRBT) ,ALLOCATABLE :: RLAPIN(:) ! eigen-values of the inverse Laplace operator
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NLTN(:) ! R%NTMAX+2-JN

REAL(KIND=JPRBT) ,ALLOCATABLE :: RMU2(:)    ! sin(theta) for dual input/output latitudes
REAL(KIND=JPRBT) ,ALLOCATABLE :: RACTHE2(:) ! 1./SQRT(R1MU2), 1/(cos(theta)) dual input/output latitudes
END TYPE FIELDS_TYPE

!flat copies of the above
REAL(KIND=JPRBT) ,ALLOCATABLE :: F_RW(:)     ! Weights of the Gaussian quadrature

TYPE(FIELDS_TYPE),ALLOCATABLE,TARGET :: FIELDS_RESOL(:)
TYPE(FIELDS_TYPE),POINTER     :: F

! scratch arrays for ltinv and ltdir and associated dimension variables

REAL(KIND=JPRBT),ALLOCATABLE :: ZAA(:,:,:)  !! JPRL for 1/2
REAL(KIND=JPRBT),ALLOCATABLE :: ZAS(:,:,:)  !! JPRL for 1/2

REAL(KIND=JPRBT), POINTER :: IZBA(:,:,:)    !! JPRL for 1/2
REAL(KIND=JPRBT),ALLOCATABLE,TARGET  :: IZBS(:,:,:) !! JPRL for 1/2
REAL(KIND=JPRBT),ALLOCATABLE :: IZCA(:,:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: IZCS(:,:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: IZCAT(:,:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: IZCST(:,:,:)

REAL(KIND=JPRBT),ALLOCATABLE :: DZBA(:,:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: DZBS(:,:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: DZBAT(:,:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: DZBST(:,:,:) !! JPRL for 1/2
REAL(KIND=JPRBT),ALLOCATABLE :: DZCA(:,:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: DZCS(:,:,:)
!REAL(KIND=JPRBT),POINTER :: DZCAT(:,:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: DZCAT(:,:,:)
REAL(KIND=JPRBT),ALLOCATABLE,TARGET :: DZCST(:,:,:)

! Arrays used for rescaling to allow half-precision Legende transforms
REAL(KIND=JPRBT), ALLOCATABLE :: ZAMAX(:,:)
REAL(KIND=JPRBT), ALLOCATABLE :: ZSMAX(:,:)

INTEGER(KIND=JPIM) :: LDZAA
INTEGER(KIND=JPIM) :: LDZAS
INTEGER(KIND=JPIM) :: TDZAA
INTEGER(KIND=JPIM) :: TDZAS

INTEGER(KIND=JPIM) :: ILDZBA
INTEGER(KIND=JPIM) :: ILDZBS
INTEGER(KIND=JPIM) :: ITDZBA
INTEGER(KIND=JPIM) :: ITDZBS
INTEGER(KIND=JPIM) :: ILDZCA
INTEGER(KIND=JPIM) :: ILDZCS
INTEGER(KIND=JPIM) :: ITDZCA
INTEGER(KIND=JPIM) :: ITDZCS



INTEGER(KIND=JPIM) :: DLDZBA
INTEGER(KIND=JPIM) :: DLDZBS
INTEGER(KIND=JPIM) :: DTDZBA
INTEGER(KIND=JPIM) :: DTDZBS
INTEGER(KIND=JPIM) :: DLDZCA
INTEGER(KIND=JPIM) :: DLDZCS
INTEGER(KIND=JPIM) :: DTDZCA
INTEGER(KIND=JPIM) :: DTDZCS

REAL(KIND=JPRB),ALLOCATABLE, TARGET :: ZIA(:,:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: ZEPSNM(:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: ZSOA1(:,:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: ZAOA1(:,:,:)
INTEGER(KIND=JPIM),ALLOCATABLE :: ISTAN(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE :: ISTAS(:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: ZSIA(:,:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: ZAIA(:,:,:)
REAL(KIND=JPRBT),ALLOCATABLE, TARGET :: ZOA1(:,:,:)
REAL(KIND=JPRBT),ALLOCATABLE, TARGET :: ZOA2(:,:,:)

END MODULE TPM_FIELDS