! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FOURIER_IN_MOD
CONTAINS
SUBROUTINE FOURIER_IN(PREEL, KFIELDS, KGL)

!**** *FOURIER_IN* - Copy fourier data from buffer to local array

!     Purpose.
!     --------
!        Routine for copying fourier data from buffer to local array

!**   Interface.
!     ----------
!     CALL FOURIER_IN(...)

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
USE TPM_TRANS,    ONLY : FOUBUF
USE TPM_GEOMETRY, ONLY : G

IMPLICIT NONE

REAL(KIND=JPRB),    INTENT(OUT) :: PREEL(:,:)
INTEGER(KIND=JPIM), INTENT(IN)  :: KFIELDS
INTEGER(KIND=JPIM), INTENT(IN)  :: KGL

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

  ! Compute offset for extraction of the fields from the m-to-l transposition buffer, FOUBUF
  ISTA = (D%NSTAGT0B(D%MSTABF(IPROC)) + D%NPNTGTB0(JM,KGL)) * 2 * KFIELDS

  ! Copy all fields from m-to-l transposition buffer to FFT work array
  ! SIMD-vectorise the per-field copy
  !$OMP SIMD
  DO JF = 1, KFIELDS
    PREEL(JF,IR) = FOUBUF(ISTA+2*JF-1)
    PREEL(JF,II) = FOUBUF(ISTA+2*JF)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE FOURIER_IN

! =====================================================================
! FOURIER_IN + FSC fusion (CPU inverse-transform hot path)
! ---------------------------------------------------------------------
! Combines the FOUBUF->PREEL copy (FOURIER_IN) with FSC's per-latitude
! Fourier-space operations:
!   * scale U and V by 1/(a*cos(theta))
!   * scale N-S derivatives by 1/(a*cos(theta))
!   * compute scalar E-W derivatives on the fly
! into a single JM loop, eliminating 3 extra passes over PREEL that FSC
! would otherwise do after FOURIER_IN.
!
! Not implemented for LATLON and LUVDER cases, so fall back to split FOURIER_IN + FSC.
!
! =====================================================================
SUBROUTINE FOURIER_IN_FSC(PREEL, KF_OUT_LT, KGL, &
 & KF_UV, IST_UV, KF_SCALARS, IST_SC, KF_SCDERS, IST_NS, IST_EW)

USE PARKIND1,     ONLY : JPIM, JPRB
USE TPM_DISTR,    ONLY : D, MYSETW
USE TPM_TRANS,    ONLY : FOUBUF
USE TPM_GEOMETRY, ONLY : G
USE TPM_FIELDS,   ONLY : F

IMPLICIT NONE

REAL(KIND=JPRB),    INTENT(INOUT) :: PREEL(:,:)
INTEGER(KIND=JPIM), INTENT(IN)    :: KF_OUT_LT   ! # rows loaded from FOUBUF
INTEGER(KIND=JPIM), INTENT(IN)    :: KGL
INTEGER(KIND=JPIM), INTENT(IN)    :: KF_UV, IST_UV
INTEGER(KIND=JPIM), INTENT(IN)    :: KF_SCALARS, IST_SC
INTEGER(KIND=JPIM), INTENT(IN)    :: KF_SCDERS, IST_NS, IST_EW

INTEGER(KIND=JPIM) :: JM, JF, IGLG, IPROC, IR, II, ISTA, IMEN, ISTAGTF, IEND_COPY
REAL(KIND=JPRB)    :: ZACHTE2, ZMUL, ZRE, ZIM

! ------------------------------------------------------------------

IGLG    = D%NPTRLS(MYSETW) + KGL - 1
IMEN    = G%NMEN(IGLG)
ISTAGTF = D%NSTAGTF(KGL)
! LATLON+LDLL branch is excluded by caller, so ZACHTE == ZACHTE2 always
ZACHTE2 = REAL(F%RACTHE(IGLG), JPRB)

! Copy-only leading rows: 1..IEND_COPY (Vor and/or Div, or empty)
IF (KF_UV > 0) THEN
  IEND_COPY = IST_UV - 1
ELSE IF (KF_SCALARS > 0) THEN
  IEND_COPY = IST_SC - 1
ELSE
  IEND_COPY = KF_OUT_LT
ENDIF

DO JM = 0, IMEN
  IPROC = D%NPROCM(JM)
  IR    = 2*JM + 1 + ISTAGTF
  II    = 2*JM + 2 + ISTAGTF
  ISTA  = (D%NSTAGT0B(D%MSTABF(IPROC)) + D%NPNTGTB0(JM,KGL)) * 2 * KF_OUT_LT
  ZMUL  = ZACHTE2 * REAL(JM, JPRB)

  ! pure copy (Vor + Div, or empty)
  !$OMP SIMD
  DO JF = 1, IEND_COPY
    PREEL(JF, IR) = FOUBUF(ISTA + 2*JF - 1)
    PREEL(JF, II) = FOUBUF(ISTA + 2*JF)
  ENDDO

  ! U and V -- load + scale by 1/(a*cos(theta))
  IF (KF_UV > 0) THEN
    !$OMP SIMD
    DO JF = IST_UV, IST_UV + 2*KF_UV - 1
      PREEL(JF, IR) = FOUBUF(ISTA + 2*JF - 1) * ZACHTE2
      PREEL(JF, II) = FOUBUF(ISTA + 2*JF)     * ZACHTE2
    ENDDO
  ENDIF

  ! scalars -- copy, and if derivatives requested also compute
  ! E-W derivatives on the fly from the just-loaded scalar values.
  IF (KF_SCALARS > 0) THEN
    IF (KF_SCDERS > 0) THEN
      !$OMP SIMD PRIVATE(ZRE, ZIM)
      DO JF = 0, KF_SCALARS - 1
        ZRE = FOUBUF(ISTA + 2*(IST_SC + JF) - 1)
        ZIM = FOUBUF(ISTA + 2*(IST_SC + JF))
        PREEL(IST_SC + JF, IR) =  ZRE
        PREEL(IST_SC + JF, II) =  ZIM
        PREEL(IST_EW + JF, IR) = -ZIM * ZMUL
        PREEL(IST_EW + JF, II) =  ZRE * ZMUL
      ENDDO
    ELSE
      !$OMP SIMD
      DO JF = IST_SC, IST_SC + KF_SCALARS - 1
        PREEL(JF, IR) = FOUBUF(ISTA + 2*JF - 1)
        PREEL(JF, II) = FOUBUF(ISTA + 2*JF)
      ENDDO
    ENDIF
  ENDIF

  ! N-S derivatives -- load + scale by 1/(a*cos(theta))
  IF (KF_SCDERS > 0) THEN
    !$OMP SIMD
    DO JF = IST_NS, IST_NS + KF_SCDERS - 1
      PREEL(JF, IR) = FOUBUF(ISTA + 2*JF - 1) * ZACHTE2
      PREEL(JF, II) = FOUBUF(ISTA + 2*JF)     * ZACHTE2
    ENDDO
  ENDIF
ENDDO

! ------------------------------------------------------------------

END SUBROUTINE FOURIER_IN_FSC
END MODULE FOURIER_IN_MOD
