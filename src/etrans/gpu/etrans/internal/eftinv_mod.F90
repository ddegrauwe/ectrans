MODULE EFTINV_MOD
CONTAINS
SUBROUTINE EFTINV(ALLOCATOR,PREEL,KF_FS)

USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DISTR       ,ONLY : D, D_NUMP
USE TPM_DIM         ,ONLY : R
USE TPMALD_DIM      ,ONLY : RALD
USE TPMALD_FFT      ,ONLY : TALD

USE TPM_HICFFT      ,ONLY : EXECUTE_INV_FFT

USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

USE ISO_C_BINDING
USE BUFFERED_ALLOCATOR_MOD

!

IMPLICIT NONE

TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_FS
REAL(KIND=JPRB),    INTENT(INOUT)  :: PREEL(:)   ! (IRLEN+2)*KF_FS*NDGL_FS

INTEGER(KIND=JPIM) :: JLOT, IRLEN, JJ
INTEGER(KIND=JPIM) :: OFFSETS(2)
INTEGER(KIND=JPIM) :: LOENS(1)
integer :: istat
character(len=32) :: cfrmt
REAL(KIND=JPRB), POINTER :: PREEL3D(:,:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

IF (LHOOK) CALL DR_HOOK('EFTINV_MOD:EFTINV',0,ZHOOK_HANDLE)

IRLEN=R%NDLON+R%NNOEXTZG

! write (6,*) __FILE__, __LINE__; call flush(6)

LOENS(1)=IRLEN
JLOT=SIZE(PREEL)/(IRLEN+2)
DO JJ=1,SIZE(OFFSETS)
  OFFSETS(JJ)=(JJ-1)*(IRLEN+2)
ENDDO
! write (*,*) __FILE__, __LINE__; call flush(6)

IF (JLOT==0) THEN
  IF (LHOOK) CALL DR_HOOK('EFTINV_MOD:EFTINV',1,ZHOOK_HANDLE)
  RETURN
ENDIF

! write (6,*) __FILE__, __LINE__; call flush(6)


CALL C_F_POINTER(C_LOC(PREEL), PREEL3D, (/ SIZE(PREEL)/(KF_FS*D%NDGL_FS),KF_FS,D%NDGL_FS /))


!$ACC DATA PRESENT(PREEL) COPYIN(LOENS,OFFSETS)

#ifdef gnarls
write (6,*) __FILE__, __LINE__; call flush(6)
!$acc update host(preel)
write (6,*) 'performing FFT with batch size ',JLOT,' on data with shape ',IRLEN+2,SIZE(PREEL)/REAL(IRLEN+2)
write (6,*) 'input:'
write (cfrmt,*) '(4X,',IRLEN+2,'F10.5)'
write (6,cfrmt) PREEL
call flush(6)
#endif

#ifndef gnarls
CALL EXECUTE_INV_FFT(PREEL(:),PREEL(:),JLOT, &
    & LOENS=LOENS, &
    & OFFSETS=OFFSETS,ALLOC=ALLOCATOR%PTR)
#endif

#ifdef gnarls
write (6,*) __FILE__, __LINE__; call flush(6)
!$acc update host(preel)
write (*,*) 'output:'
write (cfrmt,*) '(4X,',IRLEN+2,'F10.5)'
write (*,cfrmt) PREEL
call flush(6)
#endif

!$ACC END DATA

! write (6,*) __FILE__, __LINE__; call flush(6)

IF (LHOOK) CALL DR_HOOK('EFTINV_MOD:EFTINV',1,ZHOOK_HANDLE)

END SUBROUTINE EFTINV
END MODULE EFTINV_MOD