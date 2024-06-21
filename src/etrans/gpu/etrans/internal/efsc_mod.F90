MODULE EFSC_MOD
CONTAINS
SUBROUTINE EFSC(PREEL,KF_UV,KF_SCALARS,KF_SCDERS,KF_FS)
!SUBROUTINE EFSC(KF_UV,KF_SCALARS,KF_SCDERS,&
! & PUV,PSCALAR,PNSDERS,PEWDERS,PUVDERS)

!**** *FSC - Division by a*cos(theta), east-west derivatives

!     Purpose.
!     --------
!        In Fourier space divide u and v and all north-south
!        derivatives by a*cos(theta). Also compute east-west derivatives
!        of u,v,thermodynamic, passiv scalar variables and surface
!        pressure.

!**   Interface.
!     ----------
!        CALL FSC(..)
!        Explicit arguments :  PUV     - u and v
!        --------------------  PSCALAR - scalar valued varaibles
!                              PNSDERS - N-S derivative of S.V.V.
!                              PEWDERS - E-W derivative of S.V.V.
!                              PUVDERS - E-W derivative of u and v
!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03 (From SC2FSC)
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_TRANS       ,ONLY : LUVDER, LVORGP, LDIVGP
USE TPM_DISTR       ,ONLY : D, MYSETW, D_NPTRLS, D_NSTAGTF
USE TPM_GEOMETRY    ,ONLY : G, G_NMEN, G_NLOEN_MAX
USE TPMALD_GEO      ,ONLY : GALD
!

IMPLICIT NONE

REAL(KIND=JPRB) , INTENT(INOUT) :: PREEL(:)
INTEGER(KIND=JPIM) , INTENT(IN) :: KF_UV,KF_SCALARS,KF_SCDERS, KF_FS

INTEGER(KIND=JPIM) :: IMEN,ISTAGTF

INTEGER(KIND=JPIM) :: JF,IGLG,II,IR,JM,JGL
REAL(KIND=JPRB) :: ZIM
INTEGER(KIND=JPIM) :: I_UV_OFFSET, I_SC_OFFSET, I_SCDERS_OFFSET, I_UVDERS_OFFSET, IST
INTEGER(KIND=JPIM) :: IOFF_LAT,IOFF_UV,IOFF_UV_EWDER, IOFF_SCALARS, IOFF_SCALARS_EWDER
REAL(KIND=JPRB)    :: RET_REAL, RET_COMPLEX
character(len=32) :: frmt
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EFSC_MOD:EFSC',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------


IF(KF_UV > 0 .OR. KF_SCDERS > 0) THEN
  IST = 0
  IF(LVORGP) THEN
    IST = IST+KF_UV
  ENDIF
  IF(LDIVGP) THEN
    IST = IST+KF_UV
  ENDIF
  I_UV_OFFSET=IST

  IST = IST+2*KF_UV
  I_SC_OFFSET=IST

  IST = IST+KF_SCALARS
  !I_NSDERS_OFFSET=IST
  
  IST = IST+KF_SCDERS
  IF(LUVDER) THEN
    I_UVDERS_OFFSET=IST
    IST = IST+2*KF_UV
  ENDIF

  IF(KF_SCDERS > 0) THEN
    I_SCDERS_OFFSET=IST
  ENDIF
ENDIF

#ifdef ACCGPU
!$ACC DATA &
!$ACC& PRESENT(D_NPTRLS,D_NSTAGTF,PREEL,G_NMEN, G_NLOEN_MAX, D)
#endif


!     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.1      U AND V.

IF (LUVDER) THEN
#ifdef OMPGPU
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(NONE) SHARED(KF_UVPREEL)
#endif
#ifdef ACCGPU
  !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(IGLG,IOFF_LAT,IOFF_UV,IOFF_UV_EWDER,ZIM,RET_REAL,RET_COMPLEX,JM,JF,JGL) &
  !$ACC& FIRSTPRIVATE(KF_UV,I_UVDERS_OFFSET,I_UV_OFFSET,KF_FS) ASYNC(1)
#endif
  DO JGL=1,D%NDGL_FS
    DO JF=1,2*KF_UV
      DO JM=0,G_NLOEN_MAX/2
        IGLG = D_NPTRLS(MYSETW)+JGL-1
        IOFF_LAT = KF_FS*D_NSTAGTF(JGL)
        IOFF_UV = IOFF_LAT+(I_UV_OFFSET+JF-1)*(D_NSTAGTF(JGL+1)-D_NSTAGTF(JGL))
        IOFF_UV_EWDER = IOFF_LAT+(I_UVDERS_OFFSET+JF-1)*(D_NSTAGTF(JGL+1)-D_NSTAGTF(JGL))

        RET_REAL = 0.0_JPRBT
        RET_COMPLEX = 0.0_JPRBT

        IF (JM <= G_NMEN(IGLG)) THEN
          ZIM  = REAL(JM,JPRB)*GALD%EXWN

          RET_REAL = &
              & -PREEL(IOFF_UV+2*JM+2)*ZIM
          RET_COMPLEX =  &
              &  PREEL(IOFF_UV+2*JM+1)*ZIM
        ENDIF
        PREEL(IOFF_UV_EWDER+2*JM+1) = RET_REAL
        PREEL(IOFF_UV_EWDER+2*JM+2) = RET_COMPLEX
      ENDDO
    ENDDO
  ENDDO
ENDIF

!*       2.2     SCALAR VARIABLES

IF (KF_SCDERS > 0) THEN
#ifdef OMPGPU
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(NONE) SHARED(KF_SCALARS,PEWDERS,PSCALAR)
#endif
#ifdef ACCGPU
  !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(IGLG,IOFF_LAT,IOFF_SCALARS_EWDER,IOFF_SCALARS,ZIM,RET_REAL,RET_COMPLEX) &
  !$ACC& FIRSTPRIVATE(KF_SCALARS,I_SCDERS_OFFSET,I_SC_OFFSET,KF_FS) ASYNC(1)
#endif
  DO JGL=1,D%NDGL_FS
    DO JF=1,KF_SCALARS
      DO JM=0,G_NLOEN_MAX/2
        IGLG = D_NPTRLS(MYSETW)+JGL-1
        IOFF_LAT = KF_FS*D_NSTAGTF(JGL)
        IOFF_SCALARS_EWDER = IOFF_LAT+(I_SCDERS_OFFSET+JF-1)*(D_NSTAGTF(JGL+1)-D_NSTAGTF(JGL))
        IOFF_SCALARS = IOFF_LAT+(I_SC_OFFSET+JF-1)*(D_NSTAGTF(JGL+1)-D_NSTAGTF(JGL))

        RET_REAL = 0.0_JPRBT
        RET_COMPLEX = 0.0_JPRBT

        IF (JM <= G_NMEN(IGLG)) THEN
          ZIM  = REAL(JM,JPRB)*GALD%EXWN

          RET_REAL = &
              & -PREEL(IOFF_SCALARS+2*JM+2)*ZIM
          RET_COMPLEX = &
              &  PREEL(IOFF_SCALARS+2*JM+1)*ZIM
        ENDIF
        ! The rest from G_NMEN(IGLG+1)...MAX is zero truncated
        PREEL(IOFF_SCALARS_EWDER+2*JM+1) = RET_REAL
        PREEL(IOFF_SCALARS_EWDER+2*JM+2) = RET_COMPLEX
      ENDDO
    ENDDO
  ENDDO
ENDIF

#ifdef ACCGPU
!$ACC END DATA
#endif

#ifdef gnarls

IF (LUVDER) THEN
  DO JGL=1,D%NDGL_FS
    DO JF=1,2*KF_UV
      DO JM=0,G_NLOEN_MAX/2
        IGLG = D_NPTRLS(MYSETW)+JGL-1
        IOFF_LAT = KF_FS*D_NSTAGTF(JGL)
        IOFF_UV = IOFF_LAT+(I_UV_OFFSET+JF-1)*(D_NSTAGTF(JGL+1)-D_NSTAGTF(JGL))
        IOFF_UV_EWDER = IOFF_LAT+(I_UVDERS_OFFSET+JF-1)*(D_NSTAGTF(JGL+1)-D_NSTAGTF(JGL))

        IF (JM <= G_NMEN(IGLG)) THEN
          write (6,*) 'PREEL(',IOFF_UV_EWDER+2*JM+1,') <- PREEL(',IOFF_UV+2*JM+2,')'
          write (6,*) 'PREEL(',IOFF_UV_EWDER+2*JM+2,') <- PREEL(',IOFF_UV+2*JM+1,')'
        ELSE
          write (6,*) 'PREEL(',IOFF_UV_EWDER+2*JM+1,') <- 0.'
          write (6,*) 'PREEL(',IOFF_UV_EWDER+2*JM+2,') <- 0.'
        ENDIF
      ENDDO
    ENDDO
  ENDDO
ENDIF

!*       2.2     SCALAR VARIABLES

IF (KF_SCDERS > 0) THEN
  DO JGL=1,D%NDGL_FS
    DO JF=1,KF_SCALARS
      DO JM=0,G_NLOEN_MAX/2
        IGLG = D_NPTRLS(MYSETW)+JGL-1
        IOFF_LAT = KF_FS*D_NSTAGTF(JGL)
        IOFF_SCALARS_EWDER = IOFF_LAT+(I_SCDERS_OFFSET+JF-1)*(D_NSTAGTF(JGL+1)-D_NSTAGTF(JGL))
        IOFF_SCALARS = IOFF_LAT+(I_SC_OFFSET+JF-1)*(D_NSTAGTF(JGL+1)-D_NSTAGTF(JGL))

        IF (JM <= G_NMEN(IGLG)) THEN
          write (6,*) 'PREEL(',IOFF_SCALARS_EWDER+2*JM+1,') <- PREEL(',IOFF_SCALARS+2*JM+2,')'
          write (6,*) 'PREEL(',IOFF_SCALARS_EWDER+2*JM+2,') <- PREEL(',IOFF_SCALARS+2*JM+1,')'
        ELSE
          write (6,*) 'PREEL(',IOFF_SCALARS_EWDER+2*JM+1,') <- 0.'
          write (6,*) 'PREEL(',IOFF_SCALARS_EWDER+2*JM+2,') <- 0.'
        ENDIF
      ENDDO
    ENDDO
  ENDDO
ENDIF

call abort()


#endif










! ==== BELOW IS COPIED FROM ORIGINAL ETRANS-GPU, NOT COMPATIBLE WITH REDGREEN-OPTIMIZED ====
#ifdef gnarls
#ifndef gnarls
!$acc update host(preel)
IF(LUVDER)THEN

write (6,*) 'shape(PREEL) = ',shape(PREEL)
write (6,*) 'G_NMEN = ',G_NMEN
write (6,*) 'D_NSTAGTF = ',D_NSTAGTF
write (6,*) 'D_NPTRLS(MYSETW) = ',D_NPTRLS(MYSETW)
write (6,*) 'D%NDGL_FS = ',D%NDGL_FS

write (6,*) 'beginning of efsc: PREEL = '
write (frmt,*) '(',D%NLENGTF,'E14.4)'
write (*,frmt) PREEL

write (6,*) 'I_UV_OFFSET = ',I_UV_OFFSET
write (6,*) 'I_UVDERS_OFFSET = ',I_UVDERS_OFFSET

  DO JF=1,2*KF_UV
    DO JGL = 1, D%NDGL_FS
      IGLG    = D_NPTRLS(MYSETW)+JGL-1
      IMEN    = G_NMEN(IGLG)
      ISTAGTF = D_NSTAGTF(JGL)
      DO JM=0,IMEN
        ZIM=REAL(JM,JPRB)*GALD%EXWN
        IR = ISTAGTF+2*JM+1
        II = IR+1
        !PUVDERS(IR,JF) = -PUV(II,JF)*ZIM
        !PUVDERS(II,JF) =  PUV(IR,JF)*ZIM
        write (6,*) 'JF = ',JF,'; JGL = ',JGL,'; JM = ',JM
        write (6,*) '  PREEL(',D%NLENGTF*(I_UVDERS_OFFSET+JF-1)+IR,') <- -PREEL(',D%NLENGTF*(I_UV_OFFSET+JF-1)+II,') = ',PREEL(D%NLENGTF*(I_UV_OFFSET+JF-1)+II)
        write (6,*) '  PREEL(',D%NLENGTF*(I_UVDERS_OFFSET+JF-1)+II,') <-  PREEL(',D%NLENGTF*(I_UV_OFFSET+JF-1)+IR,') = ',PREEL(D%NLENGTF*(I_UV_OFFSET+JF-1)+IR)
      ENDDO
    ENDDO
  ENDDO
ENDIF

call abort()

!*       2.2     SCALAR VARIABLES
  
write (6,*) 'I_SC_OFFSET = ',I_SC_OFFSET
write (6,*) 'I_SCDERS_OFFSET = ',I_SCDERS_OFFSET

IF(KF_SCDERS > 0)THEN
  DO JF=1,KF_SCALARS
    DO JGL = 1, D%NDGL_FS
      IGLG    = D_NPTRLS(MYSETW)+JGL-1
      IMEN    = G_NMEN(IGLG)
      ISTAGTF = D_NSTAGTF(JGL)
      DO JM=0,IMEN
        ZIM=REAL(JM,JPRB)*GALD%EXWN
        IR = ISTAGTF+2*JM+1
        II = IR+1
        !PEWDERS(IR,JF) = -PSCALAR(II,JF)*ZIM
        !PEWDERS(II,JF) =  PSCALAR(IR,JF)*ZIM
        write (6,*) 'PREEL(',D%NLENGTF*(I_SCDERS_OFFSET+JF-1)+IR,') <- -PREEL(',D%NLENGTF*(I_SC_OFFSET+JF-1)+II,') = ',PREEL(D%NLENGTF*(I_SC_OFFSET+JF-1)+II)
        write (6,*) 'PREEL(',D%NLENGTF*(I_SCDERS_OFFSET+JF-1)+II,') <-  PREEL(',D%NLENGTF*(I_SC_OFFSET+JF-1)+IR,') = ',PREEL(D%NLENGTF*(I_SC_OFFSET+JF-1)+IR)
      ENDDO
    ENDDO
  ENDDO
ENDIF

write (6,*) 'efsc on cpu: PREEL = '
write (*,'(6E14.4)') PREEL

#endif


!*           EAST-WEST DERIVATIVES
!              ---------------------
  
!*       2.1      U AND V.
D_NLENGTF=D%NLENGTF
GALD_EXWN=GALD%EXWN
D_NDGL_FS=D%NDGL_FS

IF(LUVDER)THEN
!$acc parallel loop collapse (2) private (JF, JGL, IGLG, IMEN, ISTAGTF, JM, ZIM, IR, II) &
!$acc & present (D_NPTRLS, G_NMEN, D_NSTAGTF, PREEL) &
!$acc & copyin(I_UV_OFFSET, I_UVDERS_OFFSET, D_NLENGTF, GALD_EXWN,MYSETW,D_NDGL_FS)
  DO JF=1,2*KF_UV
    DO JGL = 1, D_NDGL_FS
      IGLG    = D_NPTRLS(MYSETW)+JGL-1
      IMEN    = G_NMEN(IGLG)
      ISTAGTF = D_NSTAGTF(JGL)
      DO JM=0,IMEN
        ZIM=REAL(JM,JPRB)*GALD_EXWN
        IR = ISTAGTF+2*JM+1
        II = IR+1
        !PUVDERS(IR,JF) = -PUV(II,JF)*ZIM
        !PUVDERS(II,JF) =  PUV(IR,JF)*ZIM
        PREEL(D_NLENGTF*(I_UVDERS_OFFSET+JF-1)+IR) = -PREEL(D_NLENGTF*(I_UV_OFFSET+JF-1)+II)*ZIM
        PREEL(D_NLENGTF*(I_UVDERS_OFFSET+JF-1)+II) =  PREEL(D_NLENGTF*(I_UV_OFFSET+JF-1)+IR)*ZIM
      ENDDO
    ENDDO
  ENDDO
!$acc end parallel loop
ENDIF
  
!*       2.2     SCALAR VARIABLES
  
IF(KF_SCDERS > 0)THEN
!$acc parallel loop collapse (2) private (JF, JGL, IGLG, IMEN, ISTAGTF, JM, ZIM, IR, II) &
!$acc & present (D_NPTRLS, G_NMEN, D_NSTAGTF, PREEL) &
!$acc & copyin(I_SC_OFFSET, I_SCDERS_OFFSET, D_NLENGTF, GALD_EXWN,MYSETW,D_NDGL_FS)
  DO JF=1,KF_SCALARS
    DO JGL = 1, D_NDGL_FS
      IGLG    = D_NPTRLS(MYSETW)+JGL-1
      IMEN    = G_NMEN(IGLG)
      ISTAGTF = D_NSTAGTF(JGL)
      DO JM=0,IMEN
        ZIM=REAL(JM,JPRB)*GALD_EXWN
        IR = ISTAGTF+2*JM+1
        II = IR+1
        !PEWDERS(IR,JF) = -PSCALAR(II,JF)*ZIM
        !PEWDERS(II,JF) =  PSCALAR(IR,JF)*ZIM
        PREEL(D_NLENGTF*(I_SCDERS_OFFSET+JF-1)+IR) = -PREEL(D_NLENGTF*(I_SC_OFFSET+JF-1)+II)*ZIM
        PREEL(D_NLENGTF*(I_SCDERS_OFFSET+JF-1)+II) =  PREEL(D_NLENGTF*(I_SC_OFFSET+JF-1)+IR)*ZIM
      ENDDO
    ENDDO
  ENDDO
!$acc end parallel loop
ENDIF


!$acc update host(preel)
write (6,*) 'efsc from gpu: PREEL = '
write (*,'(6E14.4)') PREEL
call flush(6)

#endif


IF (LHOOK) CALL DR_HOOK('EFSC_MOD:EFSC',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EFSC
END MODULE EFSC_MOD