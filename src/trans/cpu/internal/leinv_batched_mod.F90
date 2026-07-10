! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE LEINV_BATCHED_MOD
CONTAINS
SUBROUTINE LEINV_BATCHED(KF_OUT_LT, NBATCH, &
    KFC_VALS, KIFC_VALS, KDGLU_VALS, ILA_VALS, ILS_VALS, &
    ISKIP_VALS, ISL_VALS, KMLOC_VALS, &
    RPNMA_PTRS, RPNMS_PTRS, &
    ZWORK_BA, ZWORK_BS, ZWORK_CA, ZWORK_CS)

! Batched LEINV kernel: packs B matrices from PIA and calls batched GEMM
! for all wavenumbers (both anti- and symmetric, fused), then unpacks C to
! ZAOA1/ZSOA1.
!
! Antisymmetric: GEMM('N','N', KDGLU, KIFC, ILA, 1.0, RPNMA, KDGLU, ZBA, ILA, 0.0, ZC, KDGLU)
!   then ZAOA1(JK, ISL+JI-1) = ZC(JI, IFLD)   (transpose)
!
! Symmetric:     GEMM('N','N', KDGLU, KIFC, ILS, 1.0, RPNMS, KDGLU, ZBS, ILS, 0.0, ZC, KDGLU)
!   then PSOA1(JK, ISL+JI-1) = ZC(JI, IFLD)   (transpose)
!
! anti and sym share TRANSA='N', TRANSB='N', ALPHA=1, BETA=0 and
! write to disjoint output buffers (ZWORK_CA vs ZWORK_CS), so 
! both subsets combined into one batch with single GEMM_BATCHED call.
!
! FLT support: when S%LUSEFLT and ILA/ILS exceeds
! S%ITHRESHOLD for a given wavenumber, the dense matrix RPNMA/RPNMS is not
! used; instead the butterfly structure S%FA(KMLOC)%YBUT_STRUCT_A/S is used
! via MULT_BUTM. Wavenumbers are partitioned into:
!   - GEMM subset: contributes to the fused batched GEMM call
!   - Butterfly subset: per-JM MULT_BUTM call

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

! butterfly call metadata (anti + sym fused for OMP load balance)
! Sized at 2*NBATCH max (worst case: every KMLOC is butterfly for both subsets)
INTEGER(KIND=JPIM) :: BUT_KMLOC(2*NBATCH)
INTEGER(KIND=JPIM) :: BUT_OFF_B(2*NBATCH), BUT_OFF_C(2*NBATCH)
INTEGER(KIND=JPIM) :: BUT_BROWS(2*NBATCH)   ! input rows (= ILA or ILS)
INTEGER(KIND=JPIM) :: BUT_CROWS(2*NBATCH)   ! output rows (= KDGLU)
INTEGER(KIND=JPIM) :: BUT_KIFC(2*NBATCH)
LOGICAL            :: BUT_IS_SYM(2*NBATCH)
REAL(KIND=JPRB), POINTER :: ZB_BUF(:,:), ZC_BUF(:,:)

ITHRESHOLD = S%ITHRESHOLD
NG = 0
NBUT = 0

! Antisymmetric: accumulate GEMM triplets into the combined batch,
! or queue butterfly call metadata for the OMP-parallel butterfly pass
IOFF_BA = 1
IOFF_C = 1
DO I = 1, NBATCH
  KMLOC = KMLOC_VALS(I)
  ILA = ILA_VALS(KMLOC)
  KIFC = KIFC_VALS(I)
  KDGLU = KDGLU_VALS(I)

  IF (ILA <= ITHRESHOLD .OR. .NOT. S%LUSEFLT) THEN
    NG = NG + 1
    M(NG) = KDGLU
    N(NG) = KIFC
    K(NG) = ILA
    LDA(NG) = KDGLU
    LDB(NG) = ILA
    LDC(NG) = KDGLU

    A_PTRS(NG) = RPNMA_PTRS(KMLOC)
    B_PTRS(NG) = C_LOC(ZWORK_BA(IOFF_BA))
    C_PTRS(NG) = C_LOC(ZWORK_CA(IOFF_C))
  ELSE
    ! queue butterfly call (anti)
    NBUT = NBUT + 1
    BUT_KMLOC(NBUT)  = KMLOC
    BUT_OFF_B(NBUT)  = IOFF_BA
    BUT_OFF_C(NBUT)  = IOFF_C
    BUT_BROWS(NBUT)  = ILA
    BUT_CROWS(NBUT)  = KDGLU
    BUT_KIFC(NBUT)   = KIFC
    BUT_IS_SYM(NBUT) = .FALSE.
  ENDIF

  IOFF_BA = IOFF_BA + ILA * KIFC
  IOFF_C = IOFF_C + KDGLU * KIFC
ENDDO

! Symmetric: continue accumulating into the SAME combined batch,
! or queue butterfly call metadata
IOFF_BS = 1
IOFF_C = 1
DO I = 1, NBATCH
  KMLOC = KMLOC_VALS(I)
  ILS = ILS_VALS(KMLOC)
  KIFC = KIFC_VALS(I)
  KDGLU = KDGLU_VALS(I)

  IF (ILS <= ITHRESHOLD .OR. .NOT. S%LUSEFLT) THEN
    NG = NG + 1
    M(NG) = KDGLU
    N(NG) = KIFC
    K(NG) = ILS
    LDA(NG) = KDGLU
    LDB(NG) = ILS
    LDC(NG) = KDGLU

    A_PTRS(NG) = RPNMS_PTRS(KMLOC)
    B_PTRS(NG) = C_LOC(ZWORK_BS(IOFF_BS))
    C_PTRS(NG) = C_LOC(ZWORK_CS(IOFF_C))
  ELSE
    ! T3.4: queue butterfly call (sym)
    NBUT = NBUT + 1
    BUT_KMLOC(NBUT)  = KMLOC
    BUT_OFF_B(NBUT)  = IOFF_BS
    BUT_OFF_C(NBUT)  = IOFF_C
    BUT_BROWS(NBUT)  = ILS
    BUT_CROWS(NBUT)  = KDGLU
    BUT_KIFC(NBUT)   = KIFC
    BUT_IS_SYM(NBUT) = .TRUE.
  ENDIF

  IOFF_BS = IOFF_BS + ILS * KIFC
  IOFF_C = IOFF_C + KDGLU * KIFC
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
      CALL MULT_BUTM('N', S%FA(KMLOC)%YBUT_STRUCT_S, BUT_KIFC(I), &
                     ZB_BUF, ZC_BUF)
    ELSE
      CALL C_F_POINTER(C_LOC(ZWORK_BA(BUT_OFF_B(I))), ZB_BUF, &
                       [BUT_BROWS(I), BUT_KIFC(I)])
      CALL C_F_POINTER(C_LOC(ZWORK_CA(BUT_OFF_C(I))), ZC_BUF, &
                       [BUT_CROWS(I), BUT_KIFC(I)])
      CALL MULT_BUTM('N', S%FA(KMLOC)%YBUT_STRUCT_A, BUT_KIFC(I), &
                     ZB_BUF, ZC_BUF)
    ENDIF
  ENDDO
  !$OMP END PARALLEL DO
ENDIF

! single batched GEMM for combined anti + sym contributions
IF (NG > 0) THEN
  CALL GEMM_BATCHED('N', 'N', M, N, K, 1.0_JPRB, 0.0_JPRB, &
      A_PTRS, LDA, B_PTRS, LDB, C_PTRS, LDC, NG)
ENDIF

END SUBROUTINE LEINV_BATCHED
END MODULE LEINV_BATCHED_MOD
