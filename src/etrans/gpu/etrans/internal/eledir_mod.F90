MODULE ELEDIR_MOD
CONTAINS
SUBROUTINE ELEDIR(ALLOCATOR,PFFT)

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
USE TPMALD_DIM      ,ONLY : RALD
USE TPMALD_FFT      ,ONLY : TALD

USE TPM_HICFFT      ,ONLY : EXECUTE_DIR_FFT

USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

USE ISO_C_BINDING
USE BUFFERED_ALLOCATOR_MOD

!

IMPLICIT NONE

TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR
REAL(KIND=JPRB),    INTENT(INOUT)  :: PFFT(:,:,:)

INTEGER(KIND=JPIM) :: IRLEN, ICLEN, JLOT, JJ
!INTEGER(KIND=JPIM) :: IPLAN_C2R
TYPE(C_PTR) :: IPLAN_C2R
REAL (KIND=JPRB)   :: ZSCAL
REAL (KIND=JPRB), POINTER :: ZFFT_L(:)  ! 1D copy
INTEGER(KIND=JPIM) :: OFFSETS(2)   ! daand: why isn't OFFSETS(1) not enough?
INTEGER(KIND=JPIM) :: LOENS(1)
integer :: istat
character(len=32) :: cfrmt

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

IF (LHOOK) CALL DR_HOOK('ELEDIR_MOD:ELEDIR',0,ZHOOK_HANDLE)

IRLEN=R%NDGL+R%NNOEXTZG
LOENS(1)=IRLEN
JLOT=UBOUND(PFFT,2)*UBOUND (PFFT,3)

! compute offsets; TODO: avoid recomputing/putting on device every time.
DO JJ=1,SIZE(OFFSETS)
  OFFSETS(JJ)=(JJ-1)*(IRLEN+2)
ENDDO

CALL C_F_POINTER(C_LOC(PFFT), ZFFT_L, (/UBOUND(PFFT,1)*UBOUND(PFFT,2)*UBOUND(PFFT,3)/) )

IF (JLOT==0) THEN
  IF (LHOOK) CALL DR_HOOK('ELEDIR_MOD:ELEDIR',1,ZHOOK_HANDLE)
  RETURN
ENDIF

!write (6,*) __FILE__, __LINE__; call flush(6)



!$ACC DATA PRESENT(PFFT) COPYIN(LOENS,OFFSETS)

#ifdef gnarls
!$acc update host(pfft)
write (6,*) __FILE__, __LINE__; call flush(6)
write (*,*) 'performing FFT with batch size ',JLOT,' on data with shape ',shape(PFFT)
write (*,*) 'input:'
write (cfrmt,*) '(4X,',UBOUND(PFFT,1),'F10.5)'
write (*,cfrmt) PFFT
call flush(6)
#endif

CALL EXECUTE_DIR_FFT(ZFFT_L(:),ZFFT_L(:),-JLOT, &    ! -JLOT to have hicfft make distinction between zonal and meridional direction. Don't worry, abs(JLOT) is used internally ...
    & LOENS=LOENS, &
    & OFFSETS=OFFSETS,ALLOC=ALLOCATOR%PTR)


#ifdef gnarls
!$acc update host(pfft)
write (6,*) __FILE__, __LINE__; call flush(6)
write (*,*) 'output:'
write (cfrmt,*) '(4X,',UBOUND(PFFT,1),'F10.5)'
write (*,cfrmt) PFFT
call flush(6)
#endif

!$ACC END DATA

!write (6,*) __FILE__, __LINE__; call flush(6)

IF (LHOOK) CALL DR_HOOK('ELEDIR_MOD:ELEDIR',1,ZHOOK_HANDLE)

END SUBROUTINE ELEDIR
END MODULE ELEDIR_MOD