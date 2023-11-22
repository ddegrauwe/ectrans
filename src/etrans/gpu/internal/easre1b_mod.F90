MODULE EASRE1B_MOD
CONTAINS

SUBROUTINE EASRE1B_HIPH_WRAPPER(PFFT, FOUBUF_IN, D_NPNTGTB1, K_NPNTGTB1, R_NDGL, D_NUMP,  KJGL, KJM, KFIELD,LDS) 
            USE PARKIND_ECTRANS, ONLY : JPIM, JPRBT,JPRB
            USE, INTRINSIC:: ISO_C_BINDING

            REAL(KIND = JPRBT), TARGET:: PFFT 
            REAL(KIND = JPRBT), TARGET:: FOUBUF_IN
            INTEGER(KIND = C_INT), TARGET:: D_NPNTGTB1
            INTEGER(KIND = C_INT), INTENT(IN):: K_NPNTGTB1
            INTEGER(KIND = C_INT), INTENT(IN):: R_NDGL
            INTEGER(KIND = C_INT), INTENT(IN):: D_NUMP
            INTEGER(KIND = C_INT), INTENT(IN):: KJGL
            INTEGER(KIND = C_INT), INTENT(IN):: KJM
            INTEGER(KIND = C_INT), INTENT(IN):: KFIELD
            LOGICAL, INTENT(IN), VALUE:: LDS
     INTERFACE
         SUBROUTINE EASRE1B_HIPH(PFFT, FOUBUF_IN, D_NPNTGTB1, K_NPNTGTB1, R_NDGL, D_NUMP,  KJGL, KJM, KFIELD, LDS) bind(c, NAME="easre1b_hiph_")

             USE, INTRINSIC:: ISO_C_BINDING
             type(c_ptr), value :: PFFT 
             type(c_ptr), value :: FOUBUF_IN
             type(c_ptr), value :: D_NPNTGTB1
             INTEGER(KIND = C_INT), INTENT(IN),VALUE:: K_NPNTGTB1
             INTEGER(KIND = C_INT), INTENT(IN),VALUE:: R_NDGL
             INTEGER(KIND = C_INT), INTENT(IN),VALUE:: D_NUMP
             INTEGER(KIND = C_INT), INTENT(IN),VALUE:: KJGL
             INTEGER(KIND = C_INT), INTENT(IN),VALUE:: KJM
             INTEGER(KIND = C_INT), INTENT(IN),VALUE:: KFIELD
             LOGICAL, INTENT(IN), VALUE:: LDS
         END SUBROUTINE EASRE1B_HIPH

     END INTERFACE

     !ACC DATA PRESENT(PFFT, FOUBUF_IN, D_NPNTGTB1, R_NDGL, D_NUMP) copyin(KJGL, KJM, KFIELD, K_NPNTGTB1)

     !$ACC HOST_DATA USE_DEVICE(PFFT, FOUBUF_IN, D_NPNTGTB1) 
   
     CALL EASRE1B_HIPH(c_loc(PFFT), c_loc(FOUBUF_IN), c_loc(D_NPNTGTB1), K_NPNTGTB1, R_NDGL, D_NUMP, KJGL, KJM, KFIELD, LDS)
     
     !$ACC END HOST_DATA
    

 END SUBROUTINE EASRE1B_HIPH_WRAPPER

 
 SUBROUTINE EASRE1B_HIP(PFFT, FOUBUF_IN, D_NPNTGTB1, R_NDGL, D_NUMP, KFIELD, LDS) 
    USE PARKIND1  ,ONLY : JPIM     ,JPRB

    REAL(KIND = JPRB),INTENT(IN):: PFFT(:,:,:) 
    REAL(KIND = JPRB),INTENT(IN):: FOUBUF_IN(:)
    INTEGER(KIND = JPIM),INTENT(IN):: D_NPNTGTB1(:,:)
    INTEGER(KIND = JPIM):: KJGL
    INTEGER(KIND = JPIM):: KJM
    INTEGER(KIND = JPIM):: K_NPNTGTB1
    INTEGER(KIND = JPIM), INTENT(IN):: R_NDGL
    INTEGER(KIND = JPIM), INTENT(IN):: D_NUMP
    INTEGER(KIND = JPIM), INTENT(IN):: KFIELD
    LOGICAL, INTENT(IN), VALUE:: LDS

    KJGL = SIZE(PFFT, 1)
    KJM = SIZE(PFFT, 2)
    K_NPNTGTB1 = SIZE(D_NPNTGTB1, 1)

    CALL EASRE1B_HIPH_WRAPPER(PFFT(1,1,1), FOUBUF_IN(1), D_NPNTGTB1(1,1), K_NPNTGTB1, R_NDGL, D_NUMP, KJGL, KJM, KFIELD, LDS)

END SUBROUTINE EASRE1B_HIP


SUBROUTINE EASRE1B(KFIELD, PFFT)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R, R_NDGL
USE TPM_TRANS       ,ONLY : FOUBUF_IN
USE TPM_DISTR       ,ONLY : D, D_NUMP, D_NSTAGT0B, D_NPNTGTB1, D_NPROCL
USE TPM_GEN         ,ONLY : NERR, NOUT
!**** *ASRE1B*- Recombine antisymmetric and symmetric parts

!     Purpose.
!     --------
!        To recombine the antisymmetric/D°Nand symmetric parts of the
!        Fourier arrays and update the correct parts of the state
!        variables.

!**   Interface.
!     ----------
!        *CALL* *ASRE1B(..)

!        Explicit arguments :
!        -------------------   KFIELD-number of fields (input-c)
!                              KM-zonal wavenumber(input-c)
!                              KMLOC-local version of KM (input-c)
!                              PAOA-antisymmetric part of Fourier
!                              fields for zonal wavenumber KM (input)
!                              PSOA-symmetric part of Fourier
!                              fields for zonal wavenumber KM (input)

!        Implicit arguments :  FOUBUF_IN-output buffer (output)
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

INTEGER(KIND = JPIM), INTENT(IN):: KFIELD
REAL(KIND = JPRB), INTENT(IN)    :: PFFT(:,:,:)

INTEGER(KIND=JPIM) :: JFLD, JGL ,IPROC
INTEGER(KIND=JPIM) :: IISTAN
INTEGER(KIND=JPIM) :: JM
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE
LOGICAL            :: LLDS
INTEGER(KIND=JPIM) :: IGPU_TYPE ! 0 = OPENACC, 1 = OPENACC_LOOPS_REVERSED, 2 = HIP, 3 = HIP_LDS

PARAMETER(IGPU_TYPE=0)

!     ------------------------------------------------------------------

!*       1.    RECOMBINATION  OF SYMMETRIC AND ANTSYMMETRIC PARTS.
!              ---------------------------------------------------


IF (LHOOK) CALL DR_HOOK('EASRE1B_MOD:EASRE1B',0, ZHOOK_HANDLE)

IPROC = D_NPROCL(1)

IF (IGPU_TYPE .EQ. 0 )  THEN

WRITE(NOUT,*) "Running EASRE1B OpenACC - original loop order (IGPU_TYPE:", IGPU_TYPE,")"

!$acc parallel loop collapse(3) private (JM, JGL, JFLD, IPROC, IISTAN) &
!$acc& present (FOUBUF_IN, PFFT, D_NSTAGT0B, D_NPNTGTB1, D_NPROCL, D_NUMP, R_NDGL)

    DO JM = 1, D_NUMP  !100
       DO JFLD  =1,2*KFIELD !500        
        DO JGL=1,R_NDGL  !400
          IPROC=D_NPROCL(JGL)
          IISTAN=(D_NSTAGT0B(IPROC) + D_NPNTGTB1(JM,JGL))*2*KFIELD
          FOUBUF_IN(IISTAN+JFLD)=PFFT(JGL,JM,JFLD)
        ENDDO
    ENDDO
ENDDO
!$acc end parallel loop

ELSE IF (IGPU_TYPE .EQ. 1 )  THEN

WRITE(NOUT,*) "Running EASRE1B OpenACC - optimized loop order (IGPU_TYPE:", IGPU_TYPE,")"

!$acc parallel loop collapse(3) private (JM, JGL, JFLD, IPROC, IISTAN) &
!$acc& present (FOUBUF_IN, PFFT, D_NSTAGT0B, D_NPNTGTB1, D_NPROCL, D_NUMP, R_NDGL)

DO JGL=1,R_NDGL  !400
    DO JM = 1, D_NUMP  !100
        DO JFLD  =1,2*KFIELD !500       
          IPROC=D_NPROCL(JGL)
          IISTAN=(D_NSTAGT0B(IPROC) + D_NPNTGTB1(JM,JGL))*2*KFIELD
          FOUBUF_IN(IISTAN+JFLD)=PFFT(JGL,JM,JFLD)
        ENDDO
      ENDDO
    ENDDO
!$acc end parallel loop

ELSE

LLDS = IGPU_TYPE .EQ. 3

    WRITE(NOUT,*) "Running EASRE1B HIP (IGPU_TYPE:",IGPU_TYPE," LLDS:", LLDS,")"

    CALL EASRE1B_HIP(PFFT, FOUBUF_IN, D_NPNTGTB1, R_NDGL, D_NUMP, KFIELD, LLDS)
ENDIF


IF (LHOOK) CALL DR_HOOK('EASRE1B_MOD:EASRE1B',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE EASRE1B
END MODULE EASRE1B_MOD
