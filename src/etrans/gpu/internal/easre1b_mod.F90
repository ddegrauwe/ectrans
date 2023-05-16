MODULE EASRE1B_MOD
CONTAINS
SUBROUTINE EASRE1B(KFIELD,PFFT)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R, R_NDGL
USE TPM_TRANS       ,ONLY : FOUBUF_IN
USE TPM_DISTR       ,ONLY : D, D_NUMP, D_NSTAGT0B, D_NPNTGTB1, D_NPROCL

!**** *ASRE1B* - Recombine antisymmetric and symmetric parts

!     Purpose.
!     --------
!        To recombine the antisymmetric and symmetric parts of the
!        Fourier arrays and update the correct parts of the state
!        variables.

!**   Interface.
!     ----------
!        *CALL* *ASRE1B(..)

!        Explicit arguments :
!        -------------------   KFIELD - number of fields (input-c)
!                              KM - zonal wavenumber(input-c)
!                              KMLOC - local version of KM (input-c)
!                              PAOA - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM (input)
!                              PSOA - symmetric part of Fourier
!                              fields for zonal wavenumber KM (input)

!        Implicit arguments :  FOUBUF_IN - output buffer (output)
!        --------------------

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From ASRE1B in IFS CY22R1
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KFIELD
REAL(KIND=JPRB), INTENT(IN)    :: PFFT(:,:,:)

INTEGER(KIND=JPIM) :: JFLD, JGL ,IPROC
INTEGER(KIND=JPIM) :: IISTAN
INTEGER(KIND=JPIM) :: JM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: OFFSET_IISTAN
INTEGER(KIND=JPIM) :: KFIELD_2
!     ------------------------------------------------------------------

!*       1.    RECOMBINATION  OF SYMMETRIC AND ANTSYMMETRIC PARTS.
!              ---------------------------------------------------


IF (LHOOK) CALL DR_HOOK('EASRE1B_MOD:EASRE1B',0,ZHOOK_HANDLE)

IPROC=D_NPROCL(1)
OFFSET_IISTAN=D_NSTAGT0B(IPROC)
KFIELD_2=2*KFIELD;
!$acc parallel loop collapse(3) ASYNC(1) &
!$acc& present (FOUBUF_IN, PFFT, D_NSTAGT0B, D_NPNTGTB1, D_NPROCL, D_NUMP, R_NDGL) &
!$acc& private (JM, JGL, JFLD, IPROC, IISTAN) 
DO JGL=1,R_NDGL  !400             
  DO JM = 1, D_NUMP  !100      
      DO JFLD  =1,KFIELD_2!500
        IISTAN=(OFFSET_IISTAN + D_NPNTGTB1(JM,JGL))*KFIELD_2
        FOUBUF_IN(IISTAN+JFLD)=PFFT(JGL,JM,JFLD)
      ENDDO
  ENDDO
ENDDO
!$acc end parallel loop
!$ACC WAIT(1)
IF (LHOOK) CALL DR_HOOK('EASRE1B_MOD:EASRE1B',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE EASRE1B
END MODULE EASRE1B_MOD
