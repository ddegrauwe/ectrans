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
!     ------------------------------------------------------------------

!*       1.    RECOMBINATION  OF SYMMETRIC AND ANTSYMMETRIC PARTS.
!              ---------------------------------------------------


IF (LHOOK) CALL DR_HOOK('EASRE1B_MOD:EASRE1B',0,ZHOOK_HANDLE)


!$acc parallel loop collapse(3) private (JM, JGL, JFLD, IPROC, IISTAN) &
!$acc& present (FOUBUF_IN, PFFT, D_NSTAGT0B, D_NPNTGTB1, D_NPROCL, D_NUMP, R_NDGL) &
!$acc& copyin(KFIELD) default(none)
DO JM = 1, D_NUMP  !100
  DO JGL=1,R_NDGL  !400
    DO JFLD  =1,2*KFIELD !500
      IPROC=D_NPROCL(JGL)
      IISTAN=(D_NSTAGT0B(IPROC) + D_NPNTGTB1(JM,JGL))*2*KFIELD
      FOUBUF_IN(IISTAN+JFLD)=PFFT(JGL,JM,JFLD)
    ENDDO
  ENDDO
ENDDO
!$acc end parallel loop





#ifdef gnarls
!$acc data present(FOUBUF_IN)
!$acc update host(FOUBUF_IN)
write (20,*) __FILE__,__LINE__
write (20,*) 'FOUBUF_IN = '
write (20,'(6E18.8)') FOUBUF_IN
call flush(20)
!$acc end data
#endif

IF (LHOOK) CALL DR_HOOK('EASRE1B_MOD:EASRE1B',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE EASRE1B
END MODULE EASRE1B_MOD
