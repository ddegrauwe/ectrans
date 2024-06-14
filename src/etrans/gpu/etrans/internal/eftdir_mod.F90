MODULE EFTDIR_MOD
CONTAINS
SUBROUTINE EFTDIR(ALLOCATOR,PREEL,KF_FS,AUX_PROC)

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
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_FS
REAL(KIND=JPRB),    INTENT(INOUT)  :: PREEL(:)   ! (IRLEN+2)*NDGLG*KF_FS
EXTERNAL AUX_PROC
OPTIONAL AUX_PROC

INTEGER(KIND=JPIM) :: JLOT, IRLEN, JJ
INTEGER(KIND=JPIM), ALLOCATABLE :: OFFSETS(:)
INTEGER(KIND=JPIM) :: LOENS(1)
integer :: istat
character(len=32) :: cfrmt
REAL(KIND=JPRB) :: ZDUM
INTEGER(KIND=JPIM) :: INUL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

IF (LHOOK) CALL DR_HOOK('ELEINV_MOD:ELEINV',0,ZHOOK_HANDLE)

IRLEN=R%NDLON+R%NNOEXTZG


! Periodization of auxiliary fields in x direction
IF(R%NNOEXTZL>0) THEN
  !!! FIXME !!! CALL EXTPER(PREEL,R%NDLON+R%NNOEXTZL,1,R%NDLON,KF_FS,D%NDGL_FS,D%NSTAGTF,0)
  CALL ABORT('FIXME')
ENDIF
IF (PRESENT(AUX_PROC)) THEN
  !!! FIXME !!! CALL AUX_PROC(PREEL,ZDUM,KF_FS,D%NLENGTF,1,D%NDGL_FS,0,.TRUE.,&
  !!! FIXME !!!  & D%NSTAGTF,INUL,INUL,INUL)
  CALL ABORT('FIXME')
ENDIF


write (*,*) __FILE__, __LINE__; call flush(6)
write (*,*) 'KF_FS = ',KF_FS
write (*,*) 'shape(PREEL) = ',shape(PREEL)
write (*,*) __FILE__, __LINE__; call flush(6)


write (*,*) __FILE__, __LINE__; call flush(6)
LOENS(1)=IRLEN
JLOT=SIZE(PREEL)/(IRLEN+2)
ALLOCATE(OFFSETS(JLOT))
! compute offsets; TODO: avoid recomputing/putting on device every time.
DO JJ=1,JLOT
  OFFSETS(JJ)=(JJ-1)*(IRLEN+2)
ENDDO
write (*,*) __FILE__, __LINE__; call flush(6)

IF (JLOT==0) THEN
  IF (LHOOK) CALL DR_HOOK('ELEINV_MOD:ELEINV',1,ZHOOK_HANDLE)
  RETURN
ENDIF

write (6,*) __FILE__, __LINE__; call flush(6)



!$ACC DATA PRESENT(PREEL) COPYIN(LOENS,OFFSETS)

#ifdef gnarls
!$acc update host(preel)
write (*,*) 'performing FFT with batch size ',JLOT,' on data with shape ',IRLEN+2,SIZE(PREEL)/REAL(IRLEN+2)
write (*,*) 'input:'
write (cfrmt,*) '(4X,',IRLEN+2,'F8.2)'
write (*,cfrmt) PREEL
call flush(6)
#endif

write (6,*) __FILE__, __LINE__; call flush(6)
CALL EXECUTE_DIR_FFT(PREEL(:),PREEL(:),JLOT, &
    & LOENS=LOENS, &
    & OFFSETS=OFFSETS,ALLOC=ALLOCATOR%PTR)
write (6,*) __FILE__, __LINE__; call flush(6)


#ifdef gnarls
!$acc update host(preel)
write (*,*) 'output:'
write (cfrmt,*) '(4X,',IRLEN+2,'F8.2)'
write (*,cfrmt) PREEL
call flush(6)
#endif

DEALLOCATE(OFFSETS)

!$ACC END DATA

write (6,*) __FILE__, __LINE__; call flush(6)

IF (LHOOK) CALL DR_HOOK('EFTDIR_MOD:EFTDIR',1,ZHOOK_HANDLE)

END SUBROUTINE EFTDIR
END MODULE EFTDIR_MOD