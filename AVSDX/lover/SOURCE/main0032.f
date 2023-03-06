C
C
*DECK VISSTRING
      FUNCTION VISSTRING(FIRST, LAST, VISSTRUCTU, VISSTRUTYP)
C     ******************************************************************
C     * THIS FUNCTION OUTPUTS A FORMAT STRING BASED ON THE NUMBER OF   *
C     * FIELDS PER NODE AND ON THE PRIMITIVE DATA TYPE OF THE FIELDS.  *
C     *                                                                *
C     * ON INPUT:                                                      *
C     * FIRST      THE NUMBER OF THE FIRST FIELD THAT NEEDS FORMATTING *
C     * LAST       THE NUMBER OF THE LAST  FIELD THAT NEEDS FORMATTING *
C     * VISSTRUCTU DETERMINES THE TYPE OF THE FIELDS                   *
C     *            VISSTRUCTU()=VISRANKSCA -> SCALAR FIELD             *
C     *            VISSTRUCTU()=VISRANKVC1 -> 1-VECTOR FIELD           *
C     *            VISSTRUCTU()=VISRANKVC2 -> 2-VECTOR FIELD           *
C     *            VISSTRUCTU()=VISRANKVC3 -> 3-VECTOR FIELD           *
C     *            VISSTRUCTU()=VISRANKVEC -> 3-VECTOR FIELD           *
C     * VISSTRUTYP CONTAINS THE PRIMITIVE TYPE OF THE FIELD            *
C     *            VISTYPEDOU -> DOUBLE                                *
C     *            VISTYPEFLO -> FLOAT                                 *
C     *            VISTYPEINT -> INTEGER                               *
C     *            VISTYPESHO -> INTEGER*2                             *
C     *            VISTYPEBYT -> CHARACTER*1                           *
C     *            VISTYPECHA -> CHARACTER*1                           *
C     ******************************************************************
C
C#include "comerrn"
C#include "comrank"
C#include "comtype"

	INCLUDE 'INCLUDE/comerrn'
	INCLUDE 'INCLUDE/comrank'
	INCLUDE 'INCLUDE/comtype'

C
      INTEGER VISSTRUCTU(*), FIRST, LAST
      CHARACTER*(*) VISSTRUTYP(*)
      CHARACTER*100  VISSTRING
C
C     * COUNTERS
      INTEGER I,J
C
C     * LOCAL VARIABLES
      CHARACTER*20 STRING2
C
      IF (FIRST.LT.1) THEN
         CALL VISPRWAR(VISWARFIRS,'VISSTRING')
         FIRST=1
      ENDIF
C
      J=0
      VISSTRING(1:1)='('
      DO 20 I=FIRST,LAST-1
         J=J+1
         IF (    VISSTRUTYP(I).EQ.VISTYPEFLO
     A       .OR.VISSTRUTYP(I).EQ.VISTYPEDOU) THEN
            STRING2(1:11)='(E16.8,1X),'
         ELSE IF (    VISSTRUTYP(I).EQ.VISTYPEINT
     A            .OR.VISSTRUTYP(I).EQ.VISTYPESHO) THEN
            STRING2(1:11)='(I4,1X),   '
         ELSE IF (    VISSTRUTYP(I).EQ.VISTYPEBYT) THEN
            STRING2(1:11)='(A1,1X),   '
         ENDIF
         IF (VISSTRUCTU(I).EQ.VISRANKSCA) THEN
            VISSTRING((J-1)*12+2:(J-1)*12+2)='1'
            VISSTRING((J-1)*12+3:J*12+1)=STRING2(1:11)
         ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC1) THEN
            VISSTRING((J-1)*12+2:(J-1)*12+2)='1'
            VISSTRING((J-1)*12+3:J*12+1)=STRING2(1:11)
         ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC2) THEN
            VISSTRING((J-1)*12+2:(J-1)*12+2)='2'
            VISSTRING((J-1)*12+3:J*12+1)=STRING2(1:11)
         ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC3) THEN
            VISSTRING((J-1)*12+2:(J-1)*12+2)='3'
            VISSTRING((J-1)*12+3:J*12+1)=STRING2(1:11)
         ELSE IF (VISSTRUCTU(I).EQ.VISRANKVEC) THEN
            VISSTRING((J-1)*12+2:(J-1)*12+2)='3'
            VISSTRING((J-1)*12+3:J*12+1)=STRING2(1:11)
         ENDIF
   20 CONTINUE
C
C     * LAST FIELD
      J=J+1
      IF (    VISSTRUTYP(I).EQ.VISTYPEFLO
     A    .OR.VISSTRUTYP(I).EQ.VISTYPEDOU) THEN
         STRING2(1:11)='(E16.8,1X) '
      ELSE IF (    VISSTRUTYP(I).EQ.VISTYPEINT
     A         .OR.VISSTRUTYP(I).EQ.VISTYPESHO) THEN
         STRING2(1:11)='(I4,1X)    '
      ELSE IF (    VISSTRUTYP(I).EQ.VISTYPEBYT) THEN
         STRING2(1:11)='(A1,1X)    '
      ENDIF
      IF (VISSTRUCTU(I).EQ.VISRANKSCA) THEN
         VISSTRING((J-1)*12+2:(J-1)*12+2)='1'
         VISSTRING((J-1)*12+3:J*12+1)=STRING2(1:11)
      ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC1) THEN
         VISSTRING((J-1)*12+2:(J-1)*12+2)='1'
         VISSTRING((J-1)*12+3:J*12+1)=STRING2(1:11)
      ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC2) THEN
         VISSTRING((J-1)*12+2:(J-1)*12+2)='2'
         VISSTRING((J-1)*12+3:J*12+1)=STRING2(1:11)
      ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC3) THEN
         VISSTRING((J-1)*12+2:(J-1)*12+2)='3'
         VISSTRING((J-1)*12+3:J*12+1)=STRING2(1:11)
      ELSE IF (VISSTRUCTU(I).EQ.VISRANKVEC) THEN
         VISSTRING((J-1)*12+2:(J-1)*12+2)='3'
         VISSTRING((J-1)*12+3:J*12+1)=STRING2(1:11)
      ENDIF
      VISSTRING(J*12+2:J*12+4)=')  '
      END
