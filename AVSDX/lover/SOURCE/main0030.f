C
*DECK VISSIZEOF
      FUNCTION VISSIZEOF(TYPE)
C     ******************************************************************
C     * THIS FUNCTION OUTPUTS HARDWARE DEPENDENT INFORMATION.          *
C     * IT RETURNS THE SIZE IN BYTES OF TE VARIOUS PRIMITVE DATA TYPES.*
C     *                                                                *
C     * ON INPUT:                                                      *
C     * TYPE     THE PRIMITIVE DATA TYPE                               *
C     ******************************************************************
C
C#include "comtype"

	INCLUDE 'INCLUDE/comtype'

C
      INTEGER       VISSIZEOF
      CHARACTER*(*) TYPE
C
      IF (TYPE.EQ.VISTYPEDOU) THEN
         VISSIZEOF=8
      ELSE IF (TYPE.EQ.VISTYPEFLO) THEN
         VISSIZEOF=4
      ELSE IF (TYPE.EQ.VISTYPEINT) THEN
         VISSIZEOF=4
      ELSE IF (TYPE.EQ.VISTYPESHO) THEN
         VISSIZEOF=2
      ELSE IF (TYPE.EQ.VISTYPECHA) THEN
         VISSIZEOF=1
      ENDIF
      END
