MODULE EINV_TRANS_CTL_MOD
CONTAINS
SUBROUTINE EINV_TRANS_CTL(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KF_OUT_LT,&
 & KF_UV,KF_SCALARS,KF_SCDERS,&
 & PSPVOR,PSPDIV,PSPSCALAR,KVSETUV,KVSETSC,PGP,FSPGL_PROC,&
 & PSPSC3A,PSPSC3B,PSPSC2,KVSETSC3A,KVSETSC3B,KVSETSC2,PGPUV,PGP3A,PGP3B,PGP2,&
 & PSPMEANU,PSPMEANV)

!**** *EINV_TRANS_CTL* - Control routine for inverse spectral transform.

!     Purpose.
!     --------
!        Control routine for the inverse spectral transform

!**   Interface.
!     ----------
!     CALL EINV_TRANS_CTL(...)

!     Explicit arguments :
!     --------------------
!     KF_UV_G      - global number of spectral u-v fields
!     KF_SCALARS_G - global number of scalar spectral fields
!     KF_GP        - total number of output gridpoint fields
!     KF_FS        - total number of fields in fourier space
!     KF_OUT_LT    - total number of fields coming out from inverse LT
!     KF_UV        - local number of spectral u-v fields
!     KF_SCALARS   - local number of scalar spectral fields
!     KF_SCDERS    - local number of derivatives of scalar spectral fields
!     PSPVOR(:,:)  - spectral vorticity (input)
!     PSPDIV(:,:)  - spectral divergence (input)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
!     KVSETUV(:)  - indicating which 'b-set' in spectral space owns a
!                   vor/div field. Equivalant to NBSETLEV in the IFS.
!                   The length of KVSETUV should be the GLOBAL number
!                   of u/v fields which is the dimension of u and v releated
!                   fields in grid-point space.
!     KVESETSC(:) - indicating which 'b-set' in spectral space owns a
!                   scalar field. As for KVSETUV this argument is required
!                   if the total number of processors is greater than
!                   the number of processors used for distribution in
!                   spectral wave space.
!     FSPGL_PROC  - external procedure to be executed in fourier space
!                   before transposition
!     PGP(:,:,:)  - gridpoint fields (output)

!                  The ordering of the output fields is as follows (all
!                  parts are optional depending on the input switches):

!       vorticity     : KF_UV_G fields
!       divergence    : KF_UV_G fields
!       u             : KF_UV_G fields
!       v             : KF_UV_G fields
!       scalar fields : KF_SCALARS_G fields
!       N-S derivative of scalar fields : KF_SCALARS_G fields
!       E-W derivative of u : KF_UV_G fields
!       E-W derivative of v : KF_UV_G fields
!       E-W derivative of scalar fields : KF_SCALARS_G fields

!     Method.
!     -------

!     Externals.  SHUFFLE     - reshuffle fields for load balancing
!     ----------  FIELD_SPLIT - split fields in NPROMATR packets
!                 LTINV_CTL   - control of Legendre transform
!                 FTINV_CTL   - control of Fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 01-01-03
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_GEN         ,ONLY : NPROMATR
USE TPM_TRANS       ,ONLY : LDIVGP, LSCDERS, LUVDER, LVORGP, GROWING_ALLOCATION
USE TPM_DISTR, ONLY : D

! USE SHUFFLE_MOD     ,ONLY : SHUFFLE
! USE FIELD_SPLIT_MOD ,ONLY : FIELD_SPLIT
! USE ELTINV_CTL_MOD  ,ONLY : ELTINV_CTL
! USE EFTINV_CTL_MOD  ,ONLY : EFTINV_CTL
!

USE ELTINV_MOD
USE TRMTOL_PACK_UNPACK, ONLY : TRMTOL_UNPACK, TRMTOL_UNPACK_HANDLE, PREPARE_TRMTOL_UNPACK
USE TRMTOL_MOD
USE EFSC_MOD
USE EFTINV_MOD
USE TRLTOG_MOD
USE BUFFERED_ALLOCATOR_MOD


IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV_G
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM), INTENT(IN) :: KF_GP
INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_OUT_LT
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCDERS
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPVOR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPDIV(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPSC3B(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPSC2(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3A(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3B(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC2(:)
REAL(KIND=JPRB) ,OPTIONAL   ,INTENT(OUT) :: PGP(:,:,:)
EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PGP2(:,:,:)
REAL(KIND=JPRB)      ,OPTIONAL, INTENT(IN) :: PSPMEANU(:)
REAL(KIND=JPRB)      ,OPTIONAL, INTENT(IN) :: PSPMEANV(:)

! Local variables

INTEGER(KIND=JPIM) :: IPTRGP(KF_GP),IPTRSPUV(NPROMATR),IPTRSPSC(NPROMATR)
INTEGER(KIND=JPIM) :: ISHFUV_G(KF_GP),ISHFSC_G(KF_GP)
INTEGER(KIND=JPIM) :: IVSETUV(KF_GP),IVSETSC(KF_GP)
INTEGER(KIND=JPIM) :: IBLKS,JBLK,ISTUV_G,IENUV_G
INTEGER(KIND=JPIM) :: IF_UV_G,IF_UV,ISTUV,IF_SCALARS,IF_SCALARS_G,IF_FS,IF_GP
INTEGER(KIND=JPIM) :: JFLD,ISTSC_G,IENSC_G,ISTSC,IENSC,IENUV,IF_SCDERS,IF_OUT_LT
INTEGER(KIND=JPIM) :: IOFFD,IOFFU,IOFFV,IOFFUVD,IOFFSC,IOFFSCNS,IOFFSCEW,IOFF,IF_GPB
INTEGER(KIND=JPIM) :: KUV_OFFSET, KSCALARS_OFFSET, KSCALARS_NSDER_OFFSET, &
    & KUV_EWDER_OFFSET, KSCALARS_EWDER_OFFSET
! REAL(KIND=JPRB), POINTER :: FOUBUF(:), FOUBUF_IN(:), PREEL(:)
REAL(KIND=JPRB), POINTER :: PREEL(:)
REAL(KIND=JPRB), POINTER, SAVE :: FOUBUF(:), FOUBUF_IN(:)

INTEGER(KIND=JPIM) :: ILEI2, IDIM1, IBLEN
TYPE(BUFFERED_ALLOCATOR) :: ALLOCATOR
TYPE(ELTINV_HANDLE) :: HELTINV
TYPE(TRMTOL_HANDLE) :: HTRMTOL
TYPE(TRMTOL_UNPACK_HANDLE) :: HTRMTOL_UNPACK
TYPE(TRLTOG_HANDLE) :: HTRLTOG

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

! Perform transform

IF (LHOOK) CALL DR_HOOK('EINV_TRANS_CTL_MOD:EINV_TRANS_CTL',0,ZHOOK_HANDLE)

IF(NPROMATR > 0) THEN
  print *, "This is currently not supported and/or tested (NPROMATR > 0)"
  stop 24
ENDIF


ILEI2 = 8*KF_UV + 2*KF_SCALARS + 2*KF_SCDERS
IDIM1 = 2*KF_OUT_LT


write (6,*) __FILE__, __LINE__
write( 6,*) 'KF_FS = ',KF_FS
write (6,*) 'KF_OUT_LT = ',KF_OUT_LT
write (6,*) 'D%NUMP = ',D%NUMP
write (6,*) 'D%NDGL_FS = ',D%NDGL_FS
call flush(6)

#ifndef gnarls
! daand: manual allocation of foubuf and foubuf_in, because using allocator leads to errors in TRLTOM mpi communications !?!?
IBLEN=MAX(1,D%NLENGT1B*2*KF_OUT_LT)
IF (ALLOCATED(FOUBUF_IN)) THEN
  IF (IBLEN > SIZE(FOUBUF_IN)) THEN
!$acc exit data delete (FOUBUF_IN)
    DEALLOCATE(FOUBUF_IN)
    ALLOCATE(FOUBUF_IN(IBLEN))
!$acc enter data create (FOUBUF_IN)
  ENDIF
ELSE
  ALLOCATE(FOUBUF_IN(IBLEN))
!$acc enter data create (FOUBUF_IN)
ENDIF
#endif

#ifdef gnarls
IBLEN=MAX(1,D%NLENGT0B*2*KF_OUT_LT)
IF (ALLOCATED(FOUBUF)) THEN
  IF (IBLEN > SIZE(FOUBUF)) THEN
!$acc exit data delete (FOUBUF)
    DEALLOCATE(FOUBUF)
    ALLOCATE(FOUBUF(IBLEN))
!$acc enter data create (FOUBUF)
  ENDIF
ELSE
  ALLOCATE(FOUBUF(IBLEN))
!$acc enter data create (FOUBUF)
ENDIF
#endif

ALLOCATOR = MAKE_BUFFERED_ALLOCATOR()
HELTINV = PREPARE_ELTINV(ALLOCATOR,KF_UV,KF_SCALARS,KF_SCDERS,KF_OUT_LT)   ! ZFFT, FOUBUF_IN
! HTRMTOL_PACK = PREPARE_TRMTOL_PACK(ALLOCATOR,IF_LEG)  ! FOUBUF_IN, but incompatible with LAM
HTRMTOL = PREPARE_TRMTOL(ALLOCATOR,KF_OUT_LT)  ! FOUBUF
HTRMTOL_UNPACK = PREPARE_TRMTOL_UNPACK(ALLOCATOR,KF_FS)   ! HREEL
! HFSC = PREPARE_FSC(ALLOCATOR)
! HFTINV = PREPARE_FTINV(ALLOCATOR)
HTRLTOG = PREPARE_TRLTOG(ALLOCATOR,KF_FS,KF_GP)  ! COMBUFR and COMBUFS
CALL INSTANTIATE_ALLOCATOR(ALLOCATOR, GROWING_ALLOCATION)

IF(KF_OUT_LT > 0) THEN
  CALL GSTATS(1647,0)
! write (6,*) __FILE__, __LINE__; call flush(6)
  CALL ELTINV(ALLOCATOR, HELTINV, KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,ILEI2,IDIM1,FOUBUF_IN,&
     & PSPVOR,PSPDIV,PSPSCALAR ,&
     & PSPSC3A,PSPSC3B,PSPSC2 , &
     & FSPGL_PROC=FSPGL_PROC,PSPMEANU=PSPMEANU,PSPMEANV=PSPMEANV)
! write (6,*) __FILE__, __LINE__; call flush(6)
  CALL GSTATS(1647,1)

CALL GSTATS(152,0)
! write (6,*) __FILE__, __LINE__; call flush(6)
CALL TRMTOL_CUDAAWARE(ALLOCATOR,HTRMTOL,FOUBUF_IN,FOUBUF,KF_OUT_LT)
! write (6,*) __FILE__, __LINE__; call flush(6)
CALL TRMTOL_UNPACK(ALLOCATOR,HTRMTOL_UNPACK,FOUBUF,PREEL,KF_OUT_LT,KF_FS)   ! Formerly known as fourier_in routine
! write (6,*) __FILE__, __LINE__; call flush(6)
CALL GSTATS(152,1)
ENDIF

!!$acc update host(preel)
!write (*,*) 'PREEL = ',PREEL

IF(KF_UV > 0 .OR. KF_SCDERS > 0) THEN
  CALL EFSC(PREEL, KF_UV, KF_SCALARS, KF_SCDERS,KF_FS)
ENDIF

!!$acc update host(preel)
!write (*,*) 'after EFSC:'
!write (*,*) 'PREEL = ',PREEL

IF ( KF_FS > 0 ) THEN
  CALL EFTINV(ALLOCATOR,PREEL,KF_FS)

  CALL TRLTOG_CUDAAWARE(ALLOCATOR,HTRLTOG,PREEL,KF_FS,KF_GP,KF_UV_G,KF_SCALARS_G,&
   & KVSETUV=KVSETUV,KVSETSC=KVSETSC,&
   & KVSETSC3A=KVSETSC3A,KVSETSC3B=KVSETSC3B,KVSETSC2=KVSETSC2,&
   & PGP=PGP,PGPUV=PGPUV,PGP3A=PGP3A,PGP3B=PGP3B,PGP2=PGP2)

ENDIF

  ! from global
      ! Legendre transformations
  ! CALL GSTATS(102,0)
  ! CALL LTINV(ALLOCATOR,HLTINV,KF_UV,KF_SCALARS,&
      ! & PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2, &
      ! & ZOUTS,ZOUTA,ZOUTS0,ZOUTA0)
  ! CALL GSTATS(102,1)
  ! ! Packing into send buffer, to fourier space and unpack
  ! WRITE(NOUT,*) 'inv_trans_ctl_mod:TRMTOL_CUDAAWARE'
  ! CALL GSTATS(152,0)
  ! CALL TRMTOL_PACK(ALLOCATOR,HTRMTOL_PACK,ZOUTS,ZOUTA,ZOUTS0,ZOUTA0,FOUBUF_IN,IF_LEG)
  ! CALL TRMTOL_CUDAAWARE(ALLOCATOR,HTRMTOL,FOUBUF_IN,FOUBUF,IF_LEG)
  ! CALL TRMTOL_UNPACK(ALLOCATOR,HTRMTOL_UNPACK,FOUBUF,PREEL_COMPLEX,IF_LEG,IF_FOURIER)
  ! CALL GSTATS(152,1)

  ! CALL GSTATS(107,0)
  ! ! compute NS derivatives
  ! CALL FSC(ALLOCATOR,HFSC,PREEL_COMPLEX, IF_FOURIER, KF_UV, KF_SCALARS, KUV_OFFSET, KSCALARS_OFFSET, &
      ! & KSCALARS_NSDER_OFFSET, KUV_EWDER_OFFSET, KSCALARS_EWDER_OFFSET)
  ! !Legendre transformations
  ! CALL FTINV(ALLOCATOR, HFTINV, PREEL_COMPLEX,PREEL_REAL,IF_FOURIER)
  ! CALL GSTATS(107,1)

  ! ! Transposition into grid-point space
  ! CALL GSTATS(157,0)
  ! WRITE(NOUT,*) 'inv_trans_ctl_mod:TRLTOG_CUDAAWARE'
  ! CALL TRLTOG_CUDAAWARE(ALLOCATOR,HTRLTOG,PREEL_REAL,IF_FOURIER,KF_GP,KF_UV_G,KF_SCALARS_G,&
   ! & KVSETUV=KVSETUV,KVSETSC=KVSETSC,&
   ! & KVSETSC3A=KVSETSC3A,KVSETSC3B=KVSETSC3B,KVSETSC2=KVSETSC2,&
   ! & PGP=PGP,PGPUV=PGPUV,PGP3A=PGP3A,PGP3B=PGP3B,PGP2=PGP2)
  ! CALL GSTATS(157,1)



! original code

! CALL ELTINV_CTL(KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS, &
 ! & PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPSCALAR=PSPSCALAR,&
 ! & PSPSC3A=PSPSC3A,PSPSC3B=PSPSC3B,PSPSC2=PSPSC2,&
 ! & PSPMEANU=PSPMEANU,PSPMEANV=PSPMEANV,FSPGL_PROC=FSPGL_PROC)

! CALL EFTINV_CTL(KF_UV_G,KF_SCALARS_G,&
 ! & KF_UV,KF_SCALARS,KF_SCDERS,KF_GP,KF_FS,KF_OUT_LT,&
 ! & KVSETUV=KVSETUV,KVSETSC=KVSETSC,&
 ! & KVSETSC3A=KVSETSC3A,KVSETSC3B=KVSETSC3B,KVSETSC2=KVSETSC2,&
 ! & PGP=PGP,PGPUV=PGPUV,PGP3A=PGP3A,PGP3B=PGP3B,PGP2=PGP2)

#ifdef gnarls
write (6,*) __FILE__, __LINE__
write (6,*) 'PGP2 = '
write (6,'(4F10.5)')  PGP2
write (6,*) 'PGP3A = '
write (6,'(4F10.5)')  PGP3A
call flush(6)
#endif

write (6,*) __FILE__, __LINE__
call flush(6)

IF (LHOOK) CALL DR_HOOK('EINV_TRANS_CTL_MOD:EINV_TRANS_CTL',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EINV_TRANS_CTL
END MODULE EINV_TRANS_CTL_MOD
