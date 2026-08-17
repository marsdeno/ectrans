! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SUMP_TRANS0_MOD
CONTAINS
SUBROUTINE SUMP_TRANS0

! Set up distributed environment for the transform package (part 0)

USE EC_PARKIND      ,ONLY : JPIM
USE MPL_MODULE      ,ONLY : MPL_GROUPS_CREATE, MPL_MYRANK, MPL_NPROC
USE TPM_GEN         ,ONLY : NOUT, NERR, LMPOFF, NPRINTLEV
USE TPM_DISTR       ,ONLY : LEQ_REGIONS, MTAGDISTGP, MTAGDISTSP, MTAGGL, &
     &                      MTAGLETR, MTAGLG, MTAGLM, MTAGML, MTAGPART,  &
     &                      MYSETV, MYSETW, NPRCIDS,                     &
     &                      NPRGPEW, NPRGPNS, NPRTRNS, NPRTRV, NPRTRW,   &
     &                      MYPROC, NPROC

USE EQ_REGIONS_MOD  ,ONLY : EQ_REGIONS, MY_REGION_EW, MY_REGION_NS,      &
     &                      N_REGIONS, N_REGIONS_EW, N_REGIONS_NS
USE PE2SET_MOD      ,ONLY : PE2SET
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE TPM_ECTRANS_OPTS,ONLY : LUSE_RECTANGULAR_DECOMP

IMPLICIT NONE

LOGICAL :: LLP1,LLP2
INTEGER(KIND=JPIM) :: IPROC,JJ

!     ------------------------------------------------------------------

LLP1 = NPRINTLEV>0
LLP2 = NPRINTLEV>1
IF(LLP1) WRITE(NOUT,*) '=== ENTER ROUTINE SUMP_TRANS0 ==='


NPROC = NPRGPNS*NPRGPEW
NPRTRNS = NPRTRW
IF(MOD(NPROC,NPRTRW) /= 0 .OR. NPRTRW > NPROC) THEN
  CALL ABORT_TRANS('SUMP_TRANS0: NPROC INCONSISTENT WITH NPRTRW')
ENDIF
NPRTRV = NPROC/NPRTRW
IF(LLP1) WRITE(NOUT,*)'NPROC =',NPROC,' NPRGPNS=',NPRGPNS,' NPRGPEW=',&
 & NPRGPEW,' NPRTRW=',NPRTRW,' NPRTRV=',NPRTRV

IF(NPROC > 1 ) THEN
  IPROC = MPL_NPROC()
  IF(IPROC /= NPROC) THEN
    WRITE(NOUT,*) 'SUMP_TRANS0: NPROC=',NPROC,' BUT MPL_NPROC RETURNS',&
                 &  IPROC
    CALL ABORT_TRANS('SUMP_TRANS0: NPROC INCONSISTENT WITH MPL_NPROC')
  ENDIF
  MYPROC = MPL_MYRANK()
ELSE
  MYPROC = 1
ENDIF

IF (MYPROC > NPROC) THEN
  CALL ABORT_TRANS('SUMP_TRANS0: INCONSISTENCY IN NUMBER OF PROCESSORS USED')
ENDIF

IF( LEQ_REGIONS )THEN
  ALLOCATE(N_REGIONS(NPROC+2))
  N_REGIONS(:)=0
  CALL EQ_REGIONS(NPROC)
  ! Decomposition selector is cached by INIT_ECTRANS_OPTS, called from
  ! SETUP_TRANS0 before SUMP_TRANS0 runs (see setup_trans0.F90).
  IF( LUSE_RECTANGULAR_DECOMP )THEN
    ! Override (opt-in only): make NS collar count match Fourier wave-set
    ! count (NPRTRW) so that SUMPLAT can use Fourier boundaries directly with
    ! KPROCA=NPRTRW. This forces the rectangular decomposition that aligns
    ! grid collars 1:1 with wave-sets (precondition for LUSE_LEVS_TR).
    N_REGIONS_NS = NPRTRW
    N_REGIONS_EW = NPROC / NPRTRW
    DEALLOCATE(N_REGIONS)
    ALLOCATE(N_REGIONS(N_REGIONS_NS))
    N_REGIONS(:) = N_REGIONS_EW
  ENDIF
  ! DEFAULT (ECTRANS_RECTANGULAR_DECOMP not set): keep the stock EQ_REGIONS
  ! output untouched (N_REGIONS_NS, N_REGIONS_EW, N_REGIONS set inside
  ! EQ_REGIONS) -- the Leopardi equal-regions partition.
ELSE
  N_REGIONS_NS=NPRGPNS
  ALLOCATE(N_REGIONS(N_REGIONS_NS))
  N_REGIONS(:)=NPRGPEW
  N_REGIONS_EW=NPRGPEW
ENDIF
CALL PE2SET(MYPROC,MY_REGION_NS,MY_REGION_EW,MYSETW,MYSETV)
IF(LLP1) WRITE(NOUT,*)'MYPROC=',MYPROC,'MY_REGION_NS =',MY_REGION_NS,&
 & ' MY_REGION_EW=',MY_REGION_EW,' MYSETW=',MYSETW,' MYSETV=',MYSETV

! Rank-1 diagnostic (always on): print which grid decomposition is live 
IF( LEQ_REGIONS .AND. MYPROC == 1 )THEN
  IF( .NOT. LUSE_RECTANGULAR_DECOMP )THEN
    WRITE(NOUT,'(A,I0,A,I0,A,20(1X,I0))') &
     & 'DECOMP: NATIVE_EQREGIONS  N_REGIONS_NS=',N_REGIONS_NS, &
     & ' N_REGIONS_EW(max)=',N_REGIONS_EW, &
     & '  N_REGIONS(1:min(20,ns))=', N_REGIONS(1:MIN(20,N_REGIONS_NS))
    WRITE(NERR,'(A,I0,A,I0,A,20(1X,I0))') &
     & 'DECOMP: NATIVE_EQREGIONS  N_REGIONS_NS=',N_REGIONS_NS, &
     & ' N_REGIONS_EW(max)=',N_REGIONS_EW, &
     & '  N_REGIONS(1:min(20,ns))=', N_REGIONS(1:MIN(20,N_REGIONS_NS))
  ELSE
    WRITE(NOUT,'(A,I0,A,I0,A,I0)') &
     & 'DECOMP: RECTANGULAR  N_REGIONS_NS=NPRTRW=',N_REGIONS_NS, &
     & '  N_REGIONS_EW=NPRTRV=',N_REGIONS_EW,'  (all collars uniform=',N_REGIONS(1)
    WRITE(NERR,'(A,I0,A,I0,A,I0)') &
     & 'DECOMP: RECTANGULAR  N_REGIONS_NS=NPRTRW=',N_REGIONS_NS, &
     & '  N_REGIONS_EW=NPRTRV=',N_REGIONS_EW,'  (all collars uniform=',N_REGIONS(1)
  ENDIF
ENDIF


ALLOCATE(NPRCIDS(NPROC))
IF(LLP2)WRITE(NOUT,9) 'NPRCIDS   ',SIZE(NPRCIDS   ),SHAPE(NPRCIDS   )
DO JJ=1,NPROC
  NPRCIDS(JJ) = JJ
ENDDO

! Message passing tags

MTAGLETR = 18000
MTAGML   = 19000
MTAGLG   = 20000
MTAGPART = 21000
MTAGDISTSP = 22000
MTAGGL   = 23000
MTAGLM   = 24000
MTAGDISTGP = 25000

! Create communicators for MPI groups

IF (.NOT.LMPOFF) THEN
  CALL MPL_GROUPS_CREATE(NPRTRW, NPRTRV)
ENDIF

! Setup labels for timing package (gstats)

! CF ifs/utility    GSTATS_OUTPUT_IFS

!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)
 
END SUBROUTINE SUMP_TRANS0
END MODULE SUMP_TRANS0_MOD
