MODULE EFTINV_MOD
CONTAINS
SUBROUTINE EFTINV(PREEL,KFIELDS)

!**** *FTINV - Inverse Fourier transform

!     Purpose. Routine for Fourier to Grid-point transform
!     --------

!**   Interface.
!     ----------
!        CALL FTINV(..)

!        Explicit arguments :  PREEL   - Fourier/grid-point array
!        --------------------  KFIELDS - number of fields

!     Method.
!     -------

!     Externals.  FFT992 - FFT routine
!     ----------
!

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        G. Radnoti 01-04-24 : 2D model (NLOEN=1)
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        G. Mozdzynski (Oct 2014): support for FFTW transforms
!        G. Mozdzynski (Jun 2015): Support alternative FFTs to FFTW
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRB
USE PARKIND_ECTRANS, ONLY : JPRBT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DISTR       ,ONLY : D, MYSETW,  MYPROC, NPROC
USE TPM_GEOMETRY    ,ONLY : G
use tpm_gen, only: nout
USE TPM_FFT         ,ONLY : T !, TB
!USE BLUESTEIN_MOD   ,ONLY : BLUESTEIN_FFT
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, EXEC_FFTW
#endif
#ifdef HAVE_CUFFT
USE TPM_FFTC        ,ONLY : CREATE_PLAN_FFT, destroy_plan_fft
USE CUDA_DEVICE_MOD
#endif
USE TPM_DIM         ,ONLY : R
USE TPM_FFTH        ,ONLY : CREATE_PLAN_FFT, EXECUTE_PLAN_FFT
USE ISO_C_BINDING

IMPLICIT NONE

INTEGER (KIND=JPIM), INTENT(IN)    :: KFIELDS
REAL (KIND=JPRBT),   INTENT(INOUT) :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IRLEN, ICLEN
!INTEGER(KIND=JPIM) :: IPLAN_C2R
TYPE(C_PTR) :: IPLAN_C2R
integer :: istat
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EFTINV_MOD:EFTINV',0,ZHOOK_HANDLE)

IRLEN=R%NDLON+R%NNOEXTZG
ICLEN=D%NLENGTF/D%NDGL_FS

CALL CREATE_PLAN_FFT(IPLAN_C2R,1,IRLEN,KFIELDS*D%NDGL_FS,LDNONSTRIDED=.TRUE.)



! !$acc data present(PREEL)
! !$acc update host(PREEL)
! !$acc end data
! write (*,*) __FILE__, __LINE__
! write (*,*) 'N = ',IRLEN
! write (*,*) 'lot = ',KFIELDS*D%NDGL_FS
! write (*,*) 'shape(PREEL) = ',shape(PREEL)
! write (*,*) 'FFTH INPUT:'
! write (*,*) PREEL
! call flush(6)


CALL EXECUTE_PLAN_FFT(1,IRLEN,PREEL(1,1),PREEL(1,1),IPLAN_C2R)


! !$acc data present(PREEL)
! !$acc update host(PREEL)
! !$acc end data
! write (*,*) __FILE__, __LINE__
! write (*,*) 'FFTH OUTPUT:'
! write (*,*) PREEL
! call flush(6)


#ifdef HAVE_CUFFT
CALL CREATE_PLAN_FFT (IPLAN_C2R, +1, KN=IRLEN, KLOT=KFIELDS*D%NDGL_FS, &
                    & KISTRIDE=1, KIDIST=ICLEN/2, KOSTRIDE=1, KODIST=ICLEN)
!$acc host_data use_device(PREEL)
CALL EXECUTE_PLAN_FFTC_INPLACE (IPLAN_C2R, +1, PREEL (1, 1))
!$acc end host_data
istat = cuda_Synchronize()
#endif

IF (LHOOK) CALL DR_HOOK('EFTINV_MOD:EFTINV',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE EFTINV
END MODULE EFTINV_MOD
