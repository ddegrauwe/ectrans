MODULE EFOURIER_OUT_MOD
CONTAINS
SUBROUTINE EFOURIER_OUT(PREEL,KFIELDS)

!**** *FOURIER_OUT* - Copy fourier data from local array to buffer

!     Purpose.
!     --------
!        Routine for copying fourier data from local array to buffer

!**   Interface.
!     ----------
!     CALL FOURIER_OUT(...)

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
            &               D_NPTRLS, D_NSTAGTF, D_MSTABF, D_NSTAGT1B, &
            &               D_NPNTGTB0, D_NPROCM, D_NPROCL
USE TPM_TRANS       ,ONLY : FOUBUF_IN
USE TPM_GEOMETRY    ,ONLY : G, G_NMEN,G_NMEN_MAX
!

IMPLICIT NONE

REAL(KIND=JPRBT),  INTENT(IN) :: PREEL(:,:)
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS

INTEGER(KIND=JPIM) :: JGL
INTEGER(KIND=JPIM) :: JM, JF, IGLG, IPROC, IR, II, ISTA, JMMAX
INTEGER(KIND=JPIM) :: IOFF

!     ------------------------------------------------------------------

!$acc data &
!$acc& copy(D_NPTRLS,D_NSTAGTF,D_MSTABF,D_NSTAGT1B,D_NPNTGTB0,G_NMEN,G_NMEN_MAX,D_NPROCM) &
!$acc& present(PREEL,FOUBUF_IN)

!$acc parallel loop collapse(3) private(IGLG,JMMAX,IPROC,ISTA,IOFF)
DO JGL = 1, D%NDGL_FS
   DO JM=0,G_NMEN_MAX      
      DO JF=1,KFIELDS
         IGLG = D_NPTRLS(MYSETW)+JGL-1
         JMMAX = G_NMEN(IGLG)
         IF (JM .LE. JMMAX) THEN
            IPROC = D_NPROCM(JM)
            ISTA  = (D_NSTAGT1B(D_MSTABF(IPROC))+D_NPNTGTB0(JM,JGL))*2*KFIELDS
            IOFF  = 1+D_NSTAGTF(JGL)
            FOUBUF_IN(ISTA+2*JF-1) = PREEL(IOFF+2*JM+0,JF)
            FOUBUF_IN(ISTA+2*JF  ) = PREEL(IOFF+2*JM+1,JF)
         END IF         
      ENDDO
   ENDDO   
END DO
!$acc end data

END SUBROUTINE EFOURIER_OUT
END MODULE EFOURIER_OUT_MOD

