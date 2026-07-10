! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE LEINVAD_BATCHED_MOD
CONTAINS
SUBROUTINE LEINVAD_BATCHED(KF_OUT_LT, NBATCH, &
    KFC_VALS, KIFC_VALS, KDGLU_VALS, ILA_VALS, ILS_VALS, &
    ISKIP_VALS, ISL_VALS, KMLOC_VALS, &
    RPNMA_PTRS, RPNMS_PTRS, &
    ZWORK_BA, ZWORK_BS, ZWORK_CA, ZWORK_CS)

! Batched LEINVAD kernel: calls batched GEMM('T','N') for all wavenumbers
! (both anti- and symmetric, fused into one BLAS dispatch).
!
! Input:  ZWORK_B* packed as (KDGLU x KIFC) per wavenumber
! Output: ZWORK_C* packed as (ILA/ILS x KIFC) per wavenumber
!
! GEMM('T','N', ILA, KIFC, KDGLU):
!   C(ILA, KIFC) = A(KDGLU, ILA)^T * B(KDGLU, KIFC)
!
! FLT (Fast Legendre Transform) support: when S%LUSEFLT and ILA/ILS exceeds
! S%ITHRESHOLD for a given wavenumber, the dense matrix RPNMA/RPNMS is not
! used; instead the butterfly structure S%FA(KMLOC)%YBUT_STRUCT_A/S is used
! via MULT_BUTM('T'). Wavenumbers are partitioned into:
!   - GEMM subset: contributes to the fused batched GEMM('T','N') call
!   - Butterfly subset: per-JM MULT_BUTM('T') call

USE PARKIND1  ,ONLY : JPIM, JPRB, JPRD, JPRM
USE ECTRANS_BLAS_MOD, ONLY : GEMM
USE ECTRANS_BLAS_BATCHED_MOD, ONLY : GEMM_BATCHED
USE TPM_FLT, ONLY : S
USE BUTTERFLY_ALG_MOD, ONLY : MULT_BUTM
USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_LOC, C_F_POINTER

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KF_OUT_LT
INTEGER(KIND=JPIM), INTENT(IN) :: NBATCH
INTEGER(KIND=JPIM), INTENT(IN) :: KFC_VALS(*), KIFC_VALS(*), KDGLU_VALS(*)
INTEGER(KIND=JPIM), INTENT(IN) :: ILA_VALS(*), ILS_VALS(*)
INTEGER(KIND=JPIM), INTENT(IN) :: ISKIP_VALS(*), ISL_VALS(*), KMLOC_VALS(*)
TYPE(C_PTR), INTENT(IN) :: RPNMA_PTRS(*), RPNMS_PTRS(*)
REAL(KIND=JPRB), TARGET, INTENT(IN) :: ZWORK_BA(*), ZWORK_BS(*)
REAL(KIND=JPRB), TARGET, INTENT(INOUT) :: ZWORK_CA(*), ZWORK_CS(*)

INTEGER(KIND=JPIM) :: I, KMLOC, KIFC, KDGLU, ILA, ILS, ITHRESHOLD
INTEGER(KIND=JPIM) :: IOFF_BA, IOFF_BS, IOFF_C
INTEGER(KIND=JPIM) :: NG    ! combined anti + sym GEMM batch count
INTEGER(KIND=JPIM) :: NBUT  ! combined anti + sym butterfly count
INTEGER(KIND=JPIM) :: M(2*NBATCH), N(2*NBATCH), K(2*NBATCH)
INTEGER(KIND=JPIM) :: LDA(2*NBATCH), LDB(2*NBATCH), LDC(2*NBATCH)
TYPE(C_PTR) :: A_PTRS(2*NBATCH), B_PTRS(2*NBATCH), C_PTRS(2*NBATCH)

INTEGER(KIND=JPIM) :: BUT_KMLOC(2*NBATCH)
INTEGER(KIND=JPIM) :: BUT_OFF_B(2*NBATCH), BUT_OFF_C(2*NBATCH)
INTEGER(KIND=JPIM) :: BUT_BROWS(2*NBATCH)   ! input rows (= KDGLU)
INTEGER(KIND=JPIM) :: BUT_CROWS(2*NBATCH)   ! output rows (= ILA or ILS)
INTEGER(KIND=JPIM) :: BUT_KIFC(2*NBATCH)
LOGICAL            :: BUT_IS_SYM(2*NBATCH)
REAL(KIND=JPRB), POINTER :: ZB_BUF(:,:), ZC_BUF(:,:)

ITHRESHOLD = S%ITHRESHOLD
NG = 0
NBUT = 0

! Antisymmetric GEMM('T','N') -> ZWORK_CA, or queue butterfly call
!   GEMM:       ZWORK_CA(ILA × KIFC) = RPNMA(KDGLU × ILA)^T · ZWORK_BA(KDGLU × KIFC)
!   Butterfly:  ZWORK_CA(ILA × KIFC) = transpose-mult of butterfly · ZWORK_BA(KDGLU × KIFC)
IOFF_BA = 1
IOFF_C = 1
DO I = 1, NBATCH
  KMLOC = KMLOC_VALS(I)
  ILA = ILA_VALS(KMLOC)
  KIFC = KIFC_VALS(I)
  KDGLU = KDGLU_VALS(I)

  IF (ILA <= ITHRESHOLD .OR. .NOT. S%LUSEFLT) THEN
    NG = NG + 1
    M(NG) = ILA
    N(NG) = KIFC
    K(NG) = KDGLU
    LDA(NG) = KDGLU
    LDB(NG) = KDGLU
    LDC(NG) = ILA

    A_PTRS(NG) = RPNMA_PTRS(KMLOC)
    B_PTRS(NG) = C_LOC(ZWORK_BA(IOFF_BA))
    C_PTRS(NG) = C_LOC(ZWORK_CA(IOFF_C))
  ELSE
    ! T3.4: queue butterfly call (anti, adjoint)
    NBUT = NBUT + 1
    BUT_KMLOC(NBUT)  = KMLOC
    BUT_OFF_B(NBUT)  = IOFF_BA
    BUT_OFF_C(NBUT)  = IOFF_C
    BUT_BROWS(NBUT)  = KDGLU
    BUT_CROWS(NBUT)  = ILA
    BUT_KIFC(NBUT)   = KIFC
    BUT_IS_SYM(NBUT) = .FALSE.
  ENDIF

  IOFF_BA = IOFF_BA + KDGLU * KIFC
  IOFF_C = IOFF_C + ILA * KIFC
ENDDO

! Symmetric GEMM('T','N') -> ZWORK_CS, or queue butterfly call
! (continues accumulating into the same combined batch)
IOFF_BS = 1
IOFF_C = 1
DO I = 1, NBATCH
  KMLOC = KMLOC_VALS(I)
  ILS = ILS_VALS(KMLOC)
  KIFC = KIFC_VALS(I)
  KDGLU = KDGLU_VALS(I)

  IF (ILS <= ITHRESHOLD .OR. .NOT. S%LUSEFLT) THEN
    NG = NG + 1
    M(NG) = ILS
    N(NG) = KIFC
    K(NG) = KDGLU
    LDA(NG) = KDGLU
    LDB(NG) = KDGLU
    LDC(NG) = ILS

    A_PTRS(NG) = RPNMS_PTRS(KMLOC)
    B_PTRS(NG) = C_LOC(ZWORK_BS(IOFF_BS))
    C_PTRS(NG) = C_LOC(ZWORK_CS(IOFF_C))
  ELSE
    ! queue butterfly call (sym, adjoint)
    NBUT = NBUT + 1
    BUT_KMLOC(NBUT)  = KMLOC
    BUT_OFF_B(NBUT)  = IOFF_BS
    BUT_OFF_C(NBUT)  = IOFF_C
    BUT_BROWS(NBUT)  = KDGLU
    BUT_CROWS(NBUT)  = ILS
    BUT_KIFC(NBUT)   = KIFC
    BUT_IS_SYM(NBUT) = .TRUE.
  ENDIF

  IOFF_BS = IOFF_BS + KDGLU * KIFC
  IOFF_C = IOFF_C + ILS * KIFC
ENDDO

IF (NBUT > 0) THEN
  !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) DEFAULT(SHARED) &
  !$OMP   PRIVATE(I, KMLOC, ZB_BUF, ZC_BUF)
  DO I = 1, NBUT
    KMLOC = BUT_KMLOC(I)
    IF (BUT_IS_SYM(I)) THEN
      CALL C_F_POINTER(C_LOC(ZWORK_BS(BUT_OFF_B(I))), ZB_BUF, &
                       [BUT_BROWS(I), BUT_KIFC(I)])
      CALL C_F_POINTER(C_LOC(ZWORK_CS(BUT_OFF_C(I))), ZC_BUF, &
                       [BUT_CROWS(I), BUT_KIFC(I)])
      CALL MULT_BUTM('T', S%FA(KMLOC)%YBUT_STRUCT_S, BUT_KIFC(I), &
                     ZB_BUF, ZC_BUF)
    ELSE
      CALL C_F_POINTER(C_LOC(ZWORK_BA(BUT_OFF_B(I))), ZB_BUF, &
                       [BUT_BROWS(I), BUT_KIFC(I)])
      CALL C_F_POINTER(C_LOC(ZWORK_CA(BUT_OFF_C(I))), ZC_BUF, &
                       [BUT_CROWS(I), BUT_KIFC(I)])
      CALL MULT_BUTM('T', S%FA(KMLOC)%YBUT_STRUCT_A, BUT_KIFC(I), &
                     ZB_BUF, ZC_BUF)
    ENDIF
  ENDDO
  !$OMP END PARALLEL DO
ENDIF

IF (NG > 0) THEN
  CALL GEMM_BATCHED('T', 'N', M, N, K, 1.0_JPRB, 0.0_JPRB, &
      A_PTRS, LDA, B_PTRS, LDB, C_PTRS, LDC, NG)
ENDIF

END SUBROUTINE LEINVAD_BATCHED
END MODULE LEINVAD_BATCHED_MOD
