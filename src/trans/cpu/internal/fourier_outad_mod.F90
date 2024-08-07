! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FOURIER_OUTAD_MOD
CONTAINS
SUBROUTINE FOURIER_OUTAD(PREEL, KFIELDS, KGL)

!**** *FOURIER_OUTAD* - Copy fourier data from local array to buffer - adjoint

!     Purpose.
!     --------
!        Routine for copying fourier data from local array to buffer

!**   Interface.
!     ----------
!     CALL FOURIER_OUTAD(...)

!     Explicit arguments :  PREEL - local fourier/GP array
!     --------------------  KFIELDS - number of fields
!                           KGL - local index of latitude we are currently on
!
!     Externals.  None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01

!     ------------------------------------------------------------------

USE PARKIND1,     ONLY : JPIM, JPRB
USE TPM_DISTR,    ONLY : D, MYSETW
USE TPM_TRANS,    ONLY : FOUBUF_IN
USE TPM_GEOMETRY, ONLY : G

IMPLICIT NONE

REAL(KIND=JPRB),    INTENT(OUT) :: PREEL(:,:)
INTEGER(KIND=JPIM), INTENT(IN) :: KFIELDS
INTEGER(KIND=JPIM), INTENT(IN) :: KGL

INTEGER(KIND=JPIM) :: JM, JF, IGLG, IPROC, IR, II, ISTA

!     ------------------------------------------------------------------

! Determine global latitude index corresponding to local latitude index KGL
IGLG = D%NPTRLS(MYSETW) + KGL - 1

! Loop over all zonal wavenumbers relevant for this latitude
DO JM = 0, G%NMEN(IGLG)
  ! Get the member of the W-set responsible for this zonal wavenumber in the "m" representation
  IPROC = D%NPROCM(JM)

  ! Compute offset in FFT work array PREEL corresponding to wavenumber JM and latitude KGL
  IR = 2 * JM + 1 + D%NSTAGTF(KGL)
  II = 2 * JM + 2 + D%NSTAGTF(KGL)

  ! Compute offset for extraction of the fields from the l-to-m transposition buffer, FOUBUF, IN
  ISTA = (D%NSTAGT1B(D%MSTABF(IPROC)) + D%NPNTGTB0(JM,KGL)) * 2 * KFIELDS

  ! Copy all fields from l-to-m transposition buffer to FFT work array
  DO JF = 1, KFIELDS
    PREEL(JF,IR) = FOUBUF_IN(ISTA+2*JF-1)
    PREEL(JF,II) = FOUBUF_IN(ISTA+2*JF)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE FOURIER_OUTAD
END MODULE FOURIER_OUTAD_MOD