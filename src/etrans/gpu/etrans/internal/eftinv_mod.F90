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
REAL(KIND=JPRB),    INTENT(INOUT)  :: PREEL(:)   ! (IRLEN+3)*NDGLG*KF_FS

INTEGER(KIND=JPIM) :: JLOT, IRLEN, JJ
INTEGER(KIND=JPIM), ALLOCATABLE :: OFFSETS(:)
INTEGER(KIND=JPIM) :: LOENS(1)
integer :: istat
character(len=32) :: cfrmt

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

IF (LHOOK) CALL DR_HOOK('ELEINV_MOD:ELEINV',0,ZHOOK_HANDLE)

IRLEN=R%NDLON+R%NNOEXTZG

! write (6,*) __FILE__, __LINE__; call flush(6)

LOENS(1)=IRLEN
JLOT=SIZE(PREEL)/(IRLEN+2)   ! daand: why is it +3, instead of expected +2 ?
ALLOCATE(OFFSETS(JLOT))
! compute offsets; TODO: avoid recomputing/putting on device every time.
DO JJ=1,JLOT
  OFFSETS(JJ)=(JJ-1)*(IRLEN+2)
ENDDO
! write (*,*) __FILE__, __LINE__; call flush(6)

IF (JLOT==0) THEN
  IF (LHOOK) CALL DR_HOOK('ELEINV_MOD:ELEINV',1,ZHOOK_HANDLE)
  RETURN
ENDIF

! write (6,*) __FILE__, __LINE__; call flush(6)



!$ACC DATA PRESENT(PREEL) COPYIN(LOENS,OFFSETS)

#ifdef gnarls
write (6,*) __FILE__, __LINE__; call flush(6)
!$acc update host(preel)
write (*,*) 'performing FFT with batch size ',JLOT,' on data with shape ',IRLEN+2,SIZE(PREEL)/REAL(IRLEN+2)
write (*,*) 'input:'
write (cfrmt,*) '(4X,',IRLEN+2,'F8.2)'
write (*,cfrmt) PREEL
call flush(6)
#endif

CALL EXECUTE_INV_FFT(PREEL(:),PREEL(:),JLOT, &
    & LOENS=LOENS, &
    & OFFSETS=OFFSETS,ALLOC=ALLOCATOR%PTR)


#ifdef gnarls
!$acc update host(preel)
write (*,*) 'output:'
write (cfrmt,*) '(4X,',IRLEN+2,'F8.2)'
write (*,cfrmt) PREEL
call flush(6)
#endif

DEALLOCATE(OFFSETS)

!$ACC END DATA

! write (6,*) __FILE__, __LINE__; call flush(6)

IF (LHOOK) CALL DR_HOOK('ELEINV_MOD:ELEINV',1,ZHOOK_HANDLE)

END SUBROUTINE EFTINV
END MODULE EFTINV_MOD