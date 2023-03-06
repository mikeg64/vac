C
*DECK VISSETARGS
      SUBROUTINE VISSETARGS(DIM1, DIM2, DIM3, DIM4, DIM5,
     A                      RESO, CORD, SLIC, DIMS, SERI,
     A                      NFLD, FORM, GEOM)
C     ******************************************************************
C     * SUBROUTINE TO SET VARIOUS ARGUMENTS WHICH ARE PASSED THROUGH   *
C     * COMMONS 'DIMS' AND 'INFO'.                                     *
C     * DIM1    CONTAINS THE SIZE OF THE VECTLEN DIMENSION.            *
C     * DIM2    CONTAINS THE SIZE OF THE FIRST COMPUTATIONAL           *
C     *         DIMENSION.                                             *
C     * DIM3    CONTAINS THE SIZE OF THE SECOND COMPUTATIONAL          *
C     *         DIMENSION.                                             *
C     * DIM4    CONTAINS THE SIZE OF THE THIRD COMPUTATIONAL           *
C     *         DIMENSION.                                             *
C     * DIM5    CONTAINS THE SIZE OF THE FOURTH COMPUTATIONAL          *
C     *         DIMENSION.                                             *
C     * RESO    CONTAINS THE RESOLUTION OF THE DATA SET                *
C     *         (VISRESO(VISCOORRAD)<=VISDIM2),                        *
C     *         (VISRESO(VISCOORPOL)<=VISDIM3),                        *
C     *         (VISRESO(VISCOORTOR)<=VISDIM4),                        *
C     *         (VISRESO(VISCOORTIM)<=VISDIM5).                        *
C     * DIMS    CONTAINS THE NUMBER OF DIMENSION TO BE PRINTED.        *
C     * CORD    CONTAINS THE DIMENSIONS WHICH WILL BE PRINTED.         *
C     * SLIC    CONTAINS THE PLANE NUMBERS FOR THE SLICED DIMENSIONS.  *
C     * SERI    DETERMINES IF THE TIME DEPENDENCE IS PRINTED.          *
C     *         VISSERIES=.TRUE. -> TIME DEPENDENCE USED AS SERIES     *
C     *                             ELEMENTS.                          *
C     * NFLD    CONTAINS THE NUMBER OF FIELDS PER NODE.                *
C     * FORM    DETERMINES WHETHER OR NOT THE DATA WILL BE PRINTED     *
C     *         IN ASCII OR IN BINARY FORMAT.                          *
C     *         VISFORMATS=VISFORMASC -> ASCII                         *
C     *         VISFORMATS=VISFORMBIN -> BINARY                        *
C     *         VISFORMATS=VISFORMXDR -> BINARY XDR                    *
C     * GEOM    DETERMINES THE TYPE OF GEOMETRY:                       *
C     *         VISGEOM=VISGEOMCYL -> CYLINDRICAL GEOMETRY,            *
C     *         VISGEOM=VISGEOMTOK -> TOKAMAK GEOMETRY                 *
C     *                               (CIRCULAR CONCENTRIC),           *
C     *         VISGEOM=VISGEOMLOO -> CORONAL LOOP GEOMETRY,           *
C     *         VISGEOM=VISGEOMAUX -> AUXILIAR GEOMETRY.               *
C     ******************************************************************
C
C#include "comform"
	INCLUDE 'INCLUDE/comform'
C
      LOGICAL       SERI
      INTEGER       DIM1,    DIM2,    DIM3,    DIM4,    DIM5,
     A              RESO(4), CORD(4), SLIC(4), DIMS,
     B              NFLD,    FORM,    GEOM
C
      COMMON /DIMS/ VISDIM1, VISDIM2, VISDIM3, VISDIM4, VISDIM5
      INTEGER       VISDIM1, VISDIM2, VISDIM3, VISDIM4, VISDIM5
      COMMON /INFO/ VISRESO,    VISCORD,    VISSLIC,    VISDIMS,
     A              VISSERI,    VISVLEN,    VISNFLD,    VISFORM,
     B              VISGEOM,    VISXDR
      LOGICAL       VISSERI,    VISXDR
      INTEGER       VISRESO(4), VISCORD(4), VISSLIC(4), VISDIMS,
     A              VISVLEN,    VISNFLD,    VISFORM,    VISGEOM
      COMMON /CTIM/ VISTIME1,   VISTRES,    VISNTST,
     A              VISFIRS,    VISLAST,    VIS1TST
      LOGICAL       VISFIRS,    VISLAST,    VIS1TST
      INTEGER       VISTRES,    VISNTST
      REAL          VISTIME1
C
C     * COUNTERS
      INTEGER I
C
      VISDIM1 = DIM1
      VISDIM2 = DIM2
      VISDIM3 = DIM3
      VISDIM4 = DIM4
      VISDIM5 = DIM5
C
      DO I=1,4
         VISRESO(I) = RESO(I)
         VISCORD(I) = CORD(I)
         VISSLIC(I) = SLIC(I)
      ENDDO
C
      VISDIMS = DIMS
      VISSERI = SERI
      VISNFLD = NFLD
      IF (FORM.EQ.VISFORMXDR) THEN
         VISFORM = VISFORMBIN
         VISXDR  = .TRUE.
      ELSE
         VISFORM = FORM
         VISXDR  = .FALSE.
      ENDIF
C
      VISGEOM = GEOM
C
C     * EVERY CALL TO VISSETARGS RESETS THESE VALUES TO THEIR DEFAULTS.
      VIS1TST=.FALSE.
      VISFIRS=.TRUE.
      VISLAST=.TRUE.
C
      RETURN
      END
