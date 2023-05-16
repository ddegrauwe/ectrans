MODULE EFOURIER_IN_MOD
CONTAINS
SUBROUTINE EFOURIER_IN(PREEL,KFIELDS)

!**** *FOURIER_IN* - Copy fourier data from buffer to local array

!     Purpose.
!     --------
!        Routine for copying fourier data from buffer to local array

!**   Interface.
!     ----------
!     CALL FOURIER_IN(...)

!     Explicit arguments :  PREEL - local fourier/GP array
!     --------------------  KFIELDS - number of fields
!
!     Externals.  None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM
USE PARKIND_ECTRANS, ONLY : JPRBT

USE TPM_DISTR       ,ONLY : D, MYSETW, MYPROC, NPROC, &
             & D_NSTAGTF, D_MSTABF, D_NSTAGT0B, D_NPNTGTB0, D_NPROCM, D_NPTRLS
USE TPM_TRANS       ,ONLY : FOUBUF
USE TPM_GEOMETRY    ,ONLY : G, G_NMEN,G_NMEN_MAX
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
REAL(KIND=JPRBT), INTENT(OUT) :: PREEL(:,:)

INTEGER(KIND=JPIM) :: JM, JF, IGLG, IPROC, IR, II, ISTA
INTEGER(KIND=JPIM) :: IOFF, JGL
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC
INTEGER(KIND=JPIM) :: JF_MAX,JGL_MAX
INTEGER(KIND=JPIM) :: IOFF_LAT,OFFSET_VAR
INTEGER(KIND=JPIM) :: OFFSET_ISTA

IF(MYPROC > NPROC/2)THEN
   IBEG=1
   IEND=D%NDGL_FS
   IINC=1
 ELSE
   IBEG=D%NDGL_FS
   IEND=1
   IINC=-1
 ENDIF

JF_MAX=SIZE(PREEL,2)
JGL_MAX=SIZE(PREEL,1)

!$acc data &
!$acc& PRESENT(D_NPTRLS,D_NSTAGTF,D_MSTABF,D_NSTAGT0B,D_NPNTGTB0,G_NMEN,G_NMEN_MAX,D_NPROCM) &
!$acc& present(PREEL,FOUBUF)

!$acc kernels
DO JF=1,JF_MAX
  DO JGL=1,JGL_MAX
    PREEL(JGL,JF)=0._JPRBT
  ENDDO
ENDDO
!$acc end kernels
 
!     ------------------------------------------------------------------

OFFSET_VAR=D%NPTRLS(MYSETW)
IPROC = D_NPROCM(0)
OFFSET_ISTA = (D_NSTAGT0B(D_MSTABF(IPROC)))

!$acc parallel loop collapse(3) private(IGLG,IPROC,ISTA,IOFF) 
DO JGL=IBEG,IEND,IINC
   DO JF=1,KFIELDS     
     DO JM=0,G_NMEN_MAX
         IGLG    = OFFSET_VAR+JGL-1         
         IF ( JM .LE. G_NMEN(IGLG)) THEN            
            IOFF  = 1+D_NSTAGTF(JGL)
            ISTA  = (OFFSET_ISTA+D_NPNTGTB0(JM,JGL))*2*KFIELDS            
            PREEL(IOFF+2*JM+0,JF) = FOUBUF(ISTA+2*JF-1) 
            PREEL(IOFF+2*JM+1,JF) = FOUBUF(ISTA+2*JF  ) 
         END IF
      ENDDO
   ENDDO
END DO

!$acc end data

END SUBROUTINE EFOURIER_IN
END MODULE EFOURIER_IN_MOD

