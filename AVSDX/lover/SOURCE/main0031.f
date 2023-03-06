C
*DECK VISSIZEFORMAT
      FUNCTION VISSIZEFORMAT(TYPE)
C     ************************************************************************
C     * IT RETURNS THE SIZE IN BYTES OF THE LENGTH OF THE FORMATS IN WHICH   *
C     * THE VARIOUS TYPES ARE WRITTEN TO THE OUTPUT FILES.                   *
C     *                                                                      *
C     * ON INPUT:                                                            *
C     * TYPE     THE PRIMITIVE DATA TYPE                                     *
C     ************************************************************************
C
C#include "comtype"

	INCLUDE 'INCLUDE/comtype'

C
      INTEGER       VISSIZEFORMAT
      CHARACTER*(*) TYPE
C
      IF (TYPE.EQ.VISTYPEDOU) THEN
         VISSIZEFORMAT=17
      ELSE IF (TYPE.EQ.VISTYPEFLO) THEN
         VISSIZEFORMAT=17
      ELSE IF (TYPE.EQ.VISTYPEINT) THEN
         VISSIZEFORMAT=5
      ELSE IF (TYPE.EQ.VISTYPESHO) THEN
         VISSIZEFORMAT=5
      ELSE IF (TYPE.EQ.VISTYPECHA) THEN
         VISSIZEFORMAT=2
      ENDIF
      END
