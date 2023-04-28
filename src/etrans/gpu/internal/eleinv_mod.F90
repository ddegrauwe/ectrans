MODULE ELEINV_MOD
CONTAINS
SUBROUTINE ELEINV(KFC,KF_OUT_LT,PFFT)

!**** *LEINV* - Inverse Legendre transform.

!     Purpose.
!     --------
!        Inverse Legendre tranform of all variables(kernel).

!**   Interface.
!     ----------
!        CALL LEINV(...)

!        Explicit arguments :  KM - zonal wavenumber (input-c)
!        --------------------  KFC - number of fields to tranform (input-c)
!                              PIA - spectral fields
!                              for zonal wavenumber KM (input)
!                              PAOA1 - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM (output)
!                              PSOA1 - symmetric part of Fourier
!                              fields for zonal wavenumber KM (output)
!                              PLEPO - Legendre polonomials for zonal
!                              wavenumber KM (input-c)

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   MXMAOP - calls SGEMVX (matrix multiply)
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From LEINV in IFS CY22R1
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        R. El Khatib 01-Sep-2015 support for FFTW transforms

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DISTR       ,ONLY : D, D_NUMP
USE TPM_DIM         ,ONLY : R
#ifdef HAVE_CUFFT
USE TPM_FFTC        ,ONLY : CREATE_PLAN_FFT
USE CUDA_DEVICE_MOD
#endif
USE TPMALD_DIM      ,ONLY : RALD
USE TPMALD_FFT      ,ONLY : TALD
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

USE TPM_FFTH        ,ONLY : CREATE_PLAN_FFT, EXECUTE_PLAN_FFT
USE ISO_C_BINDING

!

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KFC
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_OUT_LT
REAL(KIND=JPRB),    INTENT(INOUT)  :: PFFT(:,:,:)

INTEGER(KIND=JPIM) :: IRLEN, ICLEN, JLOT
!INTEGER(KIND=JPIM) :: IPLAN_C2R
TYPE(C_PTR) :: IPLAN_C2R
REAL (KIND=JPRB)   :: ZSCAL

integer :: istat

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

IF (LHOOK) CALL DR_HOOK('ELEINV_MOD:ELEINV',0,ZHOOK_HANDLE)

IRLEN=R%NDGL+R%NNOEXTZG
ICLEN=RALD%NDGLSUR+R%NNOEXTZG
JLOT=UBOUND(PFFT,2)*UBOUND (PFFT,3)

IF (JLOT==0) THEN
  write (*,*) __FILE__, __LINE__
  write (*,*) 'N = ',IRLEN
  write (*,*) 'lot = ',JLOT
  IF (LHOOK) CALL DR_HOOK('ELEINV_MOD:ELEINV',1,ZHOOK_HANDLE)
  RETURN
ENDIF

CALL CREATE_PLAN_FFT(IPLAN_C2R,1,IRLEN,JLOT,LDNONSTRIDED=.TRUE.)

! !$acc data present(PFFT)
! !$acc update host(PFFT)
! !$acc end data
! write (*,*) __FILE__, __LINE__
! write (*,*) 'N = ',IRLEN
! write (*,*) 'lot = ',JLOT
! write (*,*) 'shape(PFFT) = ',shape(PFFT)
! write (*,*) 'FFTH INPUT:'
! write (*,*) PFFT
! call flush(6)

CALL EXECUTE_PLAN_FFT(1,IRLEN,PFFT(1,1,1),PFFT(1,1,1),IPLAN_C2R)

! !$acc data present(PFFT)
! !$acc update host(PFFT)
! !$acc end data
! write (*,*) __FILE__, __LINE__
! write (*,*) 'FFTH OUTPUT:'
! write (*,*) PFFT
! call flush(6)


#ifdef HAVE_CUFFT
CALL CREATE_PLAN_FFT (IPLAN_C2R, +1, KN=IRLEN, KLOT=UBOUND (PFFT,2)*UBOUND (PFFT, 3), &
                    & KISTRIDE=1, KIDIST=ICLEN/2, KOSTRIDE=1, KODIST=ICLEN)
					
!$acc host_data use_device (PFFT) 
CALL EXECUTE_PLAN_FFTC_INPLACE(IPLAN_C2R, +1, PFFT (1, 1, 1))
!$acc end host_data
istat = cuda_Synchronize()
#endif

IF (LHOOK) CALL DR_HOOK('ELEINV_MOD:ELEINV',1,ZHOOK_HANDLE)

END SUBROUTINE ELEINV
END MODULE ELEINV_MOD
