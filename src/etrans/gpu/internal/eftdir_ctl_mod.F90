MODULE EFTDIR_CTL_MOD
CONTAINS
SUBROUTINE EFTDIR_CTL(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KF_GPB, &
 & KVSETUV,KVSETSC,KPTRGP,&
 & KVSETSC3A,KVSETSC3B,KVSETSC2,&
 & PGP,PGPUV,PGP3A,PGP3B,PGP2,AUX_PROC)

!**** *EFTDIR_CTL - Direct Fourier transform control

!     Purpose. Control routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!     CALL FTDIR_CTL(..)

!     Explicit arguments :
!     --------------------
!     KF_UV_G      - global number of spectral u-v fields
!     KF_SCALARS_G - global number of scalar spectral fields
!     KF_GP        - total number of output gridpoint fields
!     KF_FS        - total number of fields in fourier space
!     KF_GPB       - total global number of output gridpoint fields
!     PGP     -  gridpoint array
!     KVSETUV - "B" set in spectral/fourier space for
!                u and v variables
!     KVSETSC - "B" set in spectral/fourier space for
!                scalar variables
!     KPTRGP  -  pointer array to fields in gridpoint space

!     Method.
!     -------

!     Externals.  TRGTOL      - transposition routine
!     ----------  FOURIER_OUT - copy fourier data to Fourier buffer
!                 FTDIR       - fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        G. Radnoti 01-03-13 adaptation to aladin (coupling)
!     01-08-28 : G. Radnoti & R. El Khatib Fix for NPROMATR /= 0
!     19-11-01 : G. Radnoti    bug corection by introducing cpl_int interface
!     02-09-30 : P. Smolikova  AUX_PROC for d4 in NH
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM          ,ONLY : R
USE TPM_TRANS        ,ONLY : FOUBUF_IN
USE TPM_DISTR        ,ONLY : D
!USE TPM_GEN          ,ONLY : LALLOPERM2
USE EFTDATA_MOD      ,ONLY : ZGTF_PERM

USE TRGTOL_MOD       ,ONLY : TRGTOL, TRGTOL_CUDAAWARE
USE EFOURIER_OUT_MOD ,ONLY : EFOURIER_OUT
USE EFTDIR_MOD       ,ONLY : EFTDIR
USE EXTPER_MOD       ,ONLY : EXTPER
!use cudafor

!

IMPLICIT NONE

! Dummy arguments

INTEGER(KIND=JPIM),INTENT(IN) :: KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KF_GPB
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPTRGP(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3A(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3B(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC2(:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP2(:,:,:)
EXTERNAL AUX_PROC
OPTIONAL AUX_PROC

! Local variables
REAL(KIND=JPRB) :: ZDUM
REAL(KIND=JPRB), POINTER :: ZGTF (:,:)
INTEGER(KIND=JPIM) :: IST,INUL,JGL,IGL,IBLEN
INTEGER(KIND=JPIM) :: IVSETUV(KF_UV_G)
INTEGER(KIND=JPIM) :: IVSETSC(KF_SCALARS_G)
INTEGER(KIND=JPIM) :: IVSET(KF_GP)
INTEGER(KIND=JPIM) :: IFGP2,IFGP3A,IFGP3B,IOFF,J3
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

! Field distribution in Spectral/Fourier space

IF (LHOOK) CALL DR_HOOK('EFTDIR_CTL_MOD:EFTDIR_CTL',0,ZHOOK_HANDLE)

IF (ALLOCATED (ZGTF_PERM)) THEN
  IF ((UBOUND (ZGTF_PERM, 1) /= D%NLENGTF) .OR. (UBOUND (ZGTF_PERM, 2) < KF_FS)) THEN
    !$acc exit data delete (ZGTF_PERM)
    DEALLOCATE (ZGTF_PERM)
  ENDIF
ENDIF

IF (.NOT. ALLOCATED (ZGTF_PERM)) THEN
  ALLOCATE (ZGTF_PERM (D%NLENGTF,KF_FS))
  !$acc enter data create (ZGTF_PERM)
ENDIF

ZGTF => ZGTF_PERM (:,1:KF_FS)

IF(PRESENT(KVSETUV)) THEN
  IVSETUV(:) = KVSETUV(:)
ELSE
  IVSETUV(:) = -1
ENDIF
IVSETSC(:) = -1
IF(PRESENT(KVSETSC)) THEN
  IVSETSC(:) = KVSETSC(:)
ELSE
  IOFF=0
  IF(PRESENT(KVSETSC2)) THEN
    IFGP2=UBOUND(KVSETSC2,1)
    IVSETSC(1:IFGP2)=KVSETSC2(:)
    IOFF=IOFF+IFGP2
  ENDIF
  IF(PRESENT(KVSETSC3A)) THEN
    IFGP3A=UBOUND(KVSETSC3A,1)
    DO J3=1,UBOUND(PGP3A,3)
      IVSETSC(IOFF+1:IOFF+IFGP3A)=KVSETSC3A(:)
      IOFF=IOFF+IFGP3A
    ENDDO
  ENDIF
  IF(PRESENT(KVSETSC3B)) THEN
    IFGP3B=UBOUND(KVSETSC3B,1)
    DO J3=1,UBOUND(PGP3B,3)
      IVSETSC(IOFF+1:IOFF+IFGP3B)=KVSETSC3B(:)
      IOFF=IOFF+IFGP3B
    ENDDO
  ENDIF
ENDIF

IST = 1
IF(KF_UV_G > 0) THEN
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
ENDIF
IF(KF_SCALARS_G > 0) THEN
  IVSET(IST:IST+KF_SCALARS_G-1) = IVSETSC(:)
  IST = IST+KF_SCALARS_G
ENDIF

! Transposition

CALL GSTATS(158,0)

!#ifdef USE_CUDA_AWARE_MPI_FT
!CALL TRGTOL_CUDAAWARE(ZGTF,KF_FS,KF_GP,KF_SCALARS_G,IVSET,KPTRGP,&
! &PGP,PGPUV,PGP3A,PGP3B,PGP2,LDTRANSPOSED=.TRUE.)
!#else
CALL TRGTOL(ZGTF,KF_FS,KF_GP,KF_SCALARS_G,IVSET,KPTRGP,&
 &PGP,PGPUV,PGP3A,PGP3B,PGP2,LDTRANSPOSED=.TRUE.)
!$acc update device(ZGTF)
!#endif

CALL GSTATS(158,1)
CALL GSTATS(106,0)

! Periodization of auxiliary fields in x direction
IF(R%NNOEXTZL>0) THEN
  CALL ABOR1 ('EFTDIR_CTL: BIPERIODICIZATION NOT SUPPORTED')
ELSE
  IF (PRESENT(AUX_PROC)) THEN
    CALL ABOR1 ('EFTDIR_CTL: AUX_PROC not supported yet; check ZGTF dimensions')
    CALL AUX_PROC(ZGTF,ZDUM,KF_FS,D%NLENGTF,1,D%NDGL_FS,0,.TRUE.,&
     & D%NSTAGTF,INUL,INUL,INUL)
  ENDIF
ENDIF


! Fourier transform

IBLEN=D%NLENGT0B*2*KF_FS
IF (ALLOCATED(FOUBUF_IN)) THEN
  IF (MAX(1,IBLEN) > SIZE(FOUBUF_IN)) THEN
!$acc exit data delete (FOUBUF_IN)
    DEALLOCATE(FOUBUF_IN)
    ALLOCATE(FOUBUF_IN(MAX(1,IBLEN)))
!$acc enter data create (FOUBUF_IN)
  ENDIF
ELSE
  ALLOCATE(FOUBUF_IN(MAX(1,IBLEN)))
!$acc enter data create (FOUBUF_IN)
ENDIF

CALL GSTATS(1640,0)

IF(KF_FS>0) THEN
  CALL EFTDIR (ZGTF, KF_FS)
ENDIF

! Save Fourier data in FOUBUF_IN

CALL EFOURIER_OUT (ZGTF, KF_FS)

CALL GSTATS(1640,1)
CALL GSTATS(106,1)

!IF (.NOT. LALLOPERM2) THEN
  !$acc exit data delete (ZGTF_PERM)
  DEALLOCATE (ZGTF_PERM)
!ENDIF

IF (LHOOK) CALL DR_HOOK('EFTDIR_CTL_MOD:EFTDIR_CTL',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EFTDIR_CTL
END MODULE EFTDIR_CTL_MOD
