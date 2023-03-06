C
*DECK VISTRANS
      SUBROUTINE VISTRANS(VISGEOMINIT, VISGEOMFINA)
C     ******************************************************************
C     * THIS SUBROUTINE TRANSFORMS THE PHYSICAL QUANTITIES AND THE     *
C     * COORDINATES. CURRENTLY, ONLY THE FOLLOWING TRANSFORMATIONS ARE *
C     * IMPLEMENTED.                                                   *
C     *              CYLINDER-> CORONAL LOOP,                          *
C     *              CYLINDER-> CIRCULAR CONCENTRIC TOKAMAK.           *
C     ******************************************************************
C
C#include "comgeom"
C#include "comerrn"

	INCLUDE 'INCLUDE/comgeom'
	INCLUDE 'INCLUDE/comerrn'

C
      INTEGER VISGEOMINIT, VISGEOMFINA
C
      IF (VISGEOMINIT.EQ.VISGEOMCYL) THEN
         IF (VISGEOMFINA.EQ.VISGEOMTOK) THEN
            VISGEOMINIT=VISGEOMFINA
         ELSE IF (VISGEOMFINA.EQ.VISGEOMLOO) THEN
            VISGEOMINIT=VISGEOMFINA
         ELSE
            CALL VISPRERR(VISERRTRAN,'VISTRANS')
         ENDIF
      ELSE
         CALL VISPRERR(VISERRTRAN,'VISTRANS')
      ENDIF
      END
