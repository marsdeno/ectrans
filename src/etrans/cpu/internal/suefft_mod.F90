! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE SUEFFT_MOD
CONTAINS
SUBROUTINE SUEFFT

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPM_GEN         ,ONLY : NOUT, NPRINTLEV
USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_GEOMETRY    ,ONLY : G
#ifdef WITH_FFT992
USE TPMALD_FFT      ,ONLY : TALD
USE TPM_FFT992      ,ONLY : T992
#endif
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, INIT_PLANS_FFTW
#endif
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

!

IMPLICIT NONE

INTEGER(KIND=JPIM) :: JGL,IGLG, ILATS
LOGICAL :: LLP1,LLP2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUEFFT_MOD:SUEFFT',0,ZHOOK_HANDLE)

IF(.NOT.D%LGRIDONLY) THEN
        
  LLP1 = NPRINTLEV>0
  LLP2 = NPRINTLEV>1
  IF(LLP1) WRITE(NOUT,*) '=== ENTER ROUTINE SUEFFT ==='

#ifdef WITH_FFT992
  IF( TALD%LFFT992 )THEN

    ALLOCATE(T992%TRIGS(R%NDLON+R%NNOEXTZL,D%NDGL_FS))
    IF(LLP2)WRITE(NOUT,9) 'T992%TRIGS    ',SIZE(T992%TRIGS),SHAPE(T992%TRIGS)
    ALLOCATE(T992%NFAX(19,D%NDGL_FS))
    IF(LLP2)WRITE(NOUT,9) 'T992%NFAX     ',SIZE(T992%NFAX),SHAPE(T992%NFAX)
    ALLOCATE(T992%LUSEFFT992(D%NDGL_FS))
    IF(LLP2)WRITE(NOUT,9) 'T992%LUSEFFT992',SIZE(T992%LUSEFFT992),SHAPE(T992%LUSEFFT992)

    !
    ! create TRIGS and NFAX for latitude lengths supported by FFT992,
    ! that is just with factors 2, 3 or 5
    !

    T992%LFFT992=.TRUE.
    ILATS=0
    DO JGL=1,D%NDGL_FS
      IGLG = D%NPTRLS(MYSETW)+JGL-1
      IF (G%NLOEN(IGLG)>1) THEN
        CALL SET99B(T992%TRIGS(1,JGL),T992%NFAX(1,JGL),G%NLOEN(IGLG)+R%NNOEXTZL,T992%LUSEFFT992(JGL))
        IF( .NOT.T992%LUSEFFT992(JGL) )THEN
          ILATS=ILATS+1
        ENDIF
      ENDIF
    ENDDO

    ALLOCATE(TALD%TRIGSE(R%NDGL+R%NNOEXTZG))
    IF(LLP2)WRITE(NOUT,9) 'TALD%TRIGSE    ',SIZE(TALD%TRIGSE),SHAPE(TALD%TRIGSE)
    ALLOCATE(TALD%NFAXE(19))
    IF(LLP2)WRITE(NOUT,9) 'TALD%NFAXE    ',SIZE(TALD%NFAXE),SHAPE(TALD%NFAXE)
    CALL SET99(TALD%TRIGSE,TALD%NFAXE,R%NDGL+R%NNOEXTZG)


  ELSE
#ifndef WITH_FFTW
    CALL ABORT_TRANS('SUEFFT: FFTW backend requested but FFTW is not enabled in this build')
#endif
  ENDIF
#endif
#ifdef WITH_FFTW
    ! trans FTINV/FTDIR always use FFTW when available, so FFTW plans are
    ! required even when the LAM FFTs use FFT992
    CALL INIT_PLANS_FFTW(MAX(R%NDLON+R%NNOEXTZL,R%NDGL+R%NNOEXTZG))
#else
#ifndef WITH_FFT992
    CALL ABORT_TRANS('SUEFFT: no Fourier backend available in this build')
#endif
#endif

ENDIF

IF (LHOOK) CALL DR_HOOK('SUEFFT_MOD:SUEFFT',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

END SUBROUTINE SUEFFT
END MODULE SUEFFT_MOD
