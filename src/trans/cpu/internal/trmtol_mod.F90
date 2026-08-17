! (C) Copyright 1995- ECMWF.
! (C) Copyright 1995- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TRMTOL_MOD

CONTAINS
SUBROUTINE TRMTOL(PFBUF_IN,PFBUF,KFIELD)

!**** *trmtol * - transposition in Fourier space

!     Purpose.
!     --------
!              Transpose Fourier buffer data from partitioning
!              over wave numbers to partitioning over latitudes.
!              It is called between direct FFT and direct Legendre
!              transform.
!              This routine is the inverse of TRLTOM.


!**   Interface.
!     ----------
!        *call* *trmtol(...)*

!        Explicit arguments : PFBUF  - Fourier coefficient buffer. It is
!        --------------------          used for both input and output.
!                             KFIELD - Number of fields communicated

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        MPP Group *ECMWF*

!     Modifications.
!     --------------
!        Original : 95-10-01
!        Modified : 97-06-17 G. Mozdzynski - control MPI mailbox use
!                                            (NCOMBFLEN) for nphase.eq.1
!        Modified : 99-05-28  D.Salmond - Optimise copies.
!        Modified : 00-02-02  M.Hamrud  - Remove NPHASE
!        D.Salmond : 01-11-23 LIMP_NOOLAP Option for non-overlapping message
!                             passing and buffer packing
!        G.Mozdzynski: 08-01-01 Cleanup
!        Y.Seity   : 07-08-31 add barrien synchronisation under LSYNC_TRANS
!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPIB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE MPL_MODULE  ,ONLY : MPL_ALLTOALLV, MPL_BARRIER, MPL_ALL_MS_COMM, MPL_WAIT, JP_NON_BLOCKING_STANDARD

USE TPM_DISTR       ,ONLY : D, MTAGML, MYSETW, NPRTRW, NPROC, MYPROC
USE TPM_ECTRANS_OPTS,ONLY : LREPORT_TR_BW
USE TPM_LEVS_TRANSPOSE,ONLY : ACCOUNT_TR_BW, DUMP_TR_SKEW
!USE TPM_GEN         ,ONLY : LSYNC_TRANS


IMPLICIT NONE


INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELD
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFBUF(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFBUF_IN(:)

INTEGER(KIND=JPIM) :: ILENS(NPRTRW),IOFFS(NPRTRW),ILENR(NPRTRW),IOFFR(NPRTRW)

INTEGER(KIND=JPIM) :: ITAG, J, ILEN, ISTA

! LREPORT_TR_BW: wall-clock timing + payload for the M<->L transpose.
! ZCM_DT = blocking alltoallv span (idx 11); ZBR_DT = pre-issue skew-barrier
! span (idx 13, load-imbalance from LTINV feeding the alltoallv).
INTEGER(KIND=JPIB) :: ICM_CLK0, ICM_CLK1, IPH_CLK0, IPH_CLK1, ITR_CLK_RATE, ITR_BYTES
REAL(KIND=JPRD)    :: ZCM_DT, ZBR_DT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE_BAR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE_BAR2

!INTEGER(KIND=JPIM) :: IREQ


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TRMTOL',0,ZHOOK_HANDLE)


ITAG = MTAGML

DO J=1,NPRTRW
  ILENS(J) = D%NLTSFTB(J)*KFIELD
  IOFFS(J) = D%NSTAGT0B(J)*KFIELD
  ILENR(J) = D%NLTSGTB(J)*KFIELD
  IOFFR(J) = D%NSTAGT0B(D%MSTABF(J))*KFIELD
ENDDO

IF(NPROC > 1) THEN
  IF (LHOOK) CALL DR_HOOK('TRMTOL_BAR',0,ZHOOK_HANDLE_BAR)
  CALL GSTATS_BARRIER(764)
  IF (LHOOK) CALL DR_HOOK('TRMTOL_BAR',1,ZHOOK_HANDLE_BAR)

  ! Diagnostic (LREPORT_TR_BW only): barrier over the same comm the alltoallv
  ! synchronizes on, timed as idx 13. Ranks enter the alltoallv staggered by
  ! their (spectrally imbalanced) LTINV times; this barrier charges that skew
  ! separately so the alltoallv span (idx 11) reflects ~pure transfer.
  IF (LREPORT_TR_BW) THEN
    CALL SYSTEM_CLOCK(COUNT=IPH_CLK0, COUNT_RATE=ITR_CLK_RATE)
    CALL MPL_BARRIER(KCOMM=MPL_ALL_MS_COMM, CDSTRING='TRMTOL: SKEW BARRIER')
    CALL SYSTEM_CLOCK(COUNT=IPH_CLK1)
    ZBR_DT = REAL(IPH_CLK1 - IPH_CLK0, JPRD) / REAL(ITR_CLK_RATE, JPRD)
    CALL SYSTEM_CLOCK(COUNT=ICM_CLK0)
  ENDIF

  CALL GSTATS(807,0)
  CALL MPL_ALLTOALLV(PSENDBUF=PFBUF_IN,KSENDCOUNTS=ILENS,&
   & PRECVBUF=PFBUF,KRECVCOUNTS=ILENR,KSENDDISPL=IOFFS,KRECVDISPL=IOFFR,&
   & KCOMM=MPL_ALL_MS_COMM,CDSTRING='TRMTOL:')

  ! LREPORT_TR_BW: close the timer and account this rank's payload. ILENS is
  ! already the element count sent to a peer (= F-space latitudes * KFIELD),
  ! so bytes = sum over the NPRTRW peers of ILENS * sizeof(JPRB). Both the
  ! alltoallv span and the skew barrier reduce over MPL_ALL_MS_COMM (one rank
  ! per node), the same comm the transpose runs on.
  IF (LREPORT_TR_BW) THEN
    CALL SYSTEM_CLOCK(COUNT=ICM_CLK1)
    ZCM_DT = REAL(ICM_CLK1 - ICM_CLK0, JPRD) / REAL(ITR_CLK_RATE, JPRD)
    ITR_BYTES = 0_JPIB
    DO J=1,NPRTRW
      ITR_BYTES = ITR_BYTES + INT(ILENS(J), JPIB)
    ENDDO
    ITR_BYTES = ITR_BYTES * INT(STORAGE_SIZE(1.0_JPRB)/8, JPIB)
    CALL ACCOUNT_TR_BW(11, 'TRMTOL-COMM (M->L, g807 alltoallv span)', ITR_BYTES, ZCM_DT, &
      & KCOMM=MPL_ALL_MS_COMM, KCOMMSIZE=NPRTRW, LPRINT=(MYPROC==1))
    CALL ACCOUNT_TR_BW(13, 'TRMTOL-SKEW (g764 pre-issue barrier)', ITR_BYTES, ZBR_DT, &
      & KCOMM=MPL_ALL_MS_COMM, KCOMMSIZE=NPRTRW, LPRINT=(MYPROC==1))
    ! World-gather the per-rank barrier wait for straggler profiling (opt-in).
    CALL DUMP_TR_SKEW(13, ZBR_DT)
  ENDIF
!Faster on Cray - because of peculiarity of their MPICH
! CALL MPL_ALLTOALLV(PSENDBUF=PFBUF_IN,KSENDCOUNTS=ILENS,&
!  & PRECVBUF=PFBUF,KRECVCOUNTS=ILENR,KSENDDISPL=IOFFS,KRECVDISPL=IOFFR,&
!  & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ,&
!  & KCOMM=MPL_ALL_MS_COMM,CDSTRING='TRMTOL:')
! CALL MPL_WAIT(KREQUEST=IREQ,CDSTRING='TRMTOL: WAIT')

  CALL GSTATS(807,1)
  IF (LHOOK) CALL DR_HOOK('TRMTOL_BAR2',0,ZHOOK_HANDLE_BAR2)
  CALL GSTATS_BARRIER2(764)
  IF (LHOOK) CALL DR_HOOK('TRMTOL_BAR2',1,ZHOOK_HANDLE_BAR2)
ELSE
  ILEN = D%NLTSGTB(MYSETW)*KFIELD
  ISTA = D%NSTAGT0B(MYSETW)*KFIELD+1
  CALL GSTATS(1608,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(J)
  DO J=ISTA,ISTA+ILEN-1
    PFBUF(J) = PFBUF_IN(J)
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1608,1)
ENDIF


IF (LHOOK) CALL DR_HOOK('TRMTOL',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE TRMTOL
END MODULE TRMTOL_MOD
