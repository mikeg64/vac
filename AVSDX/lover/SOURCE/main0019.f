
C
*DECK VISDXINTHEAD
      SUBROUTINE VISDXINTHEAD(VISPOSITIM, VISNT,
     A                        VISFILNAME, VISRES,     VISNDIM,  VISCOOR,
     B                        VISSERIES,  VISNFIELDS, VISFIENAME,
     C                        VISSTRUCTU, VISSTRUTYP, VISFORMATS,
     D                        VISWORKINT, VISGEOM)
C     ******************************************************************
C     * THIS ROUTINE WRITES THE FIRST PART OF DX'S INTERNAL FILE       *
C     * FORMAT (OBJECT DEFINITIONS ETC.). THE SECOND PART CONTAINING   *
C     * THE DATA (BOTH COORDINATES AND DATA VALUES) HAS TO BE WRITTEN  *
C     * WITH CALLS TO VISDXWRCOO AND VISDXWRDAT. THE DATA CAN BE       *
C     * IMMEDIATELY READ INTO DX USING THE 'IMPORT' MODULE.            *
C     * THE GENERAL FILE FORMAT AS WELL AS THE INTERNAL DX FORMAT ARE  *
C     * SUPPORTED IN THIS LIBRARY. HOWEVER, WE DO RECOMMEND TO USE THE *
C     * INTERNAL FORMAT ONLY SINCE IT IS MUCH MORE EFFICIENT WITH      *
C     * RESPECT TO I/O PERFORMANCE. THE GENERAL FILE FORMAT IS ONLY    *
C     * SUITED FOR PORTABILITY REASONS. IT ALLOWS ONE TO OUTPUT        *
C     * EVERYTHING IN ASCII AND IT SEPERATES THE HEADER FROM THE DATA  *
C     * PART, I.E. DIFFERENT FILES ARE USED. THE OUTPUT IN THE DATA    *
C     * FILE IS COLUMN BASED AND CAN BE READ BY OTHER APPLICATIONS     *
C     * EASILY.                                                        *
C     *                                                                *
C     * ON INPUT:                                                      *
C     * VISPOSITIM CONTAINS THE TIME COORDINATES (LABELS).             *
C     * VISNT      CONTAINS THE SIZE OF THE FOURTH COMPUTATIONAL       *
C     *            DIMENSION.                                          *
C     * VISFILNAME CONTAINS THE FILENAME OF THE OUTPUT FILE. IF IT     *
C     *            DOES NOT ALLREADY CONTAIN THE .dx POSTFIX IT WILL   *
C     *            BE ADDED.                                           *
C     * VISRES     CONTAINS THE RESOLUTION OF THE DATA SET             *
C     *            (VISRESINIT(VISCOORRAD)<=VISNI),                    *
C     *            (VISRESINIT(VISCOORPOL)<=VISNJ),                    *
C     *            (VISRESINIT(VISCOORTOR)<=VISNK),                    *
C     *            (VISRESINIT(VISCOORTIM)<=VISNT).                    *
C     * VISNDIM    CONTAINS THE NUMBER OF DIMENSION TO BE PRINTED.     *
C     * VISCOOR    CONTAINS THE DIMENSIONS WHICH WILL BE PRINTED.      *
C     * VISSERIES  DETERMINES IF THE TIME DEPENDENCE IS PRINTED.       *
C     *            VISSERIES=.TRUE. -> TIME DEPENDENCE USED AS SERIES  *
C     *                                ELEMENTS.                       *
C     * VISNFIELDS CONTAINS THE NUMBER OF FIELDS PER NODE.             *
C     * VISFIENAME CONTAINS THE NAMES OF THE DIFFERENT FIELDS.         *
C     * VISSTRUCTU DETERMINES THE TYPE OF THE FIELDS                   *
C     *            VISSTRUCTU()=VISRANKSCA -> SCALAR FIELD,            *
C     *            VISSTRUCTU()=VISRANKVC1 -> 1-VECTOR FIELD.          *
C     *            VISSTRUCTU()=VISRANKVC2 -> 2-VECTOR FIELD.          *
C     *            VISSTRUCTU()=VISRANKVC3 -> 3-VECTOR FIELD.          *
C     *            VISSTRUCTU()=VISRANKVEC -> 3-VECTOR FIELD.          *
C     * VISFORMATS DETERMINES WHETHER OR NOT THE DATA WILL BE PRINTED  *
C     *            IN ASCII OR IN BINARY FORMAT.                       *
C     *            VISFORMATS=VISFORMASC -> ASCII                      *
C     *            VISFORMATS=VISFORMBIN -> BINARY                     *
C     * VISSTRUTYP CONTAINS THE PRIMITIVE TYPE OF THE FIELD            *
C     *            VISTYPEDOU -> DOUBLE                                *
C     *            VISTYPEFLO -> FLOAT                                 *
C     * VISWORKINT IS A INTEGER WORKING ARRAY SUPPLIED BY THE USER.    *
C     *            ITS DIMENSION SHOULD EQUAL VISNFIELDS*VISNT.        *
C     * VISGEOM    DETERMINES THE TYPE OF GEOMETRY:                    *
C     *            VISGEOM=VISGEOMCYL -> CYLINDRICAL GEOMETRY,         *
C     *            VISGEOM=VISGEOMTOK -> TOKAMAK GEOMETRY              *
C     *                                  (CIRCULAR CONCENTRIC),        *
C     *            VISGEOM=VISGEOMLOO -> CORONAL LOOP GEOMETRY,        *
C     *            VISGEOM=VISGEOMAUX -> AUXILIAR GEOMETRY.            *
C     *                                                                *
C     * VISDXINTHEAD IS CALLED BY THE USER PROGRAM OR BY VISDX.        *
C     ******************************************************************
C
C#include "comcoor"
C#include "comerrn"
C#include "cominou"
C#include "comcate"
C#include "comrank"
C#include "comform"
C#include "comtype"
C#include "comorde"
C#include "comcarr"

	INCLUDE 'INCLUDE/comcoor'
	INCLUDE 'INCLUDE/comerrn'
	INCLUDE 'INCLUDE/cominou'
	INCLUDE 'INCLUDE/comcate'
	INCLUDE 'INCLUDE/comrank'
	INCLUDE 'INCLUDE/comform'
	INCLUDE 'INCLUDE/comtype'
	INCLUDE 'INCLUDE/comorde'
	INCLUDE 'INCLUDE/comcarr'


C
      LOGICAL       VISSERIES
      INTEGER       VISNT,      VISRES(*),      VISNDIM,    VISCOOR(*),
     A              VISNFIELDS, VISSTRUCTU(*),  VISFORMATS, VISGEOM,
     B              VISWORKINT(*)
      CHARACTER*(*) VISFILNAME, VISFIENAME(*),  VISSTRUTYP(*)
      REAL          VISPOSITIM(*)
C
C     *COUNTERS
      INTEGER I,J
      DATA J/1/
C
C     * LOCAL VARIABLES
      INTEGER   RANK, SHAPE, OFFSET,
     A          ITEMS,
     B          OBJNUMBER, MEMNUMBER,
     C          IND, NSCA, NVEC
      CHARACTER FILENAME*80, CATEGORY*7, TYPE*6, FORMATS*8,
     A          VISTMPCHAR*22
      LOGICAL   NONAMES
C
C     * FUNCTIONS
      INTEGER   VISSIZEOF, VISSIZEFORMAT
      LOGICAL   VISNAMEOK
C
C     * CHECKING DIMENSIONS
      IF (VISNDIM.LT.1.OR.VISNDIM.GT.3) THEN
         CALL VISPRERR(VISERRNDIM,'VISDXINTHEAD')
      ENDIF
C
      DO 1 I=1,VISNDIM
         IF (VISCOOR(I).EQ.VISCOORTIM) THEN
            CALL VISPRERR(VISERRSERI,'VISDXINTHEAD')
         ENDIF
    1 CONTINUE
C
C     * SETTING THE POSITIONS AS THE LAST FIELD
C     * DX DEMANDS THE COORDINATES TO BE FLOATS
      VISSTRUCTU(VISNFIELDS+1)=VISRANKVC3
      VISSTRUTYP(VISNFIELDS+1)=VISTYPEFLO
      VISFIENAME(VISNFIELDS+1)='locations'
C
      IF (VISFORMATS.EQ.VISFORMBIN) THEN
         FORMATS=BYTEORDER//'ieee'
      ELSE IF (VISFORMATS.EQ.VISFORMASC) THEN
         FORMATS='ascii'
      ELSE
         CALL VISPRWAR(VISWARFORM,'VISDXINTHEAD')
         FORMATS='ascii'
      ENDIF
C
C     * FIELD NAMES
      NONAMES=.FALSE.
      DO 10 I=1,VISNFIELDS
         IF (.NOT.VISNAMEOK(VISFIENAME(I))) THEN
            NONAMES=.TRUE.
         ENDIF
   10 CONTINUE
C
      IF (NONAMES) THEN
         NSCA=0
         NVEC=0
         DO 20 I=1,VISNFIELDS
            IF (VISSTRUCTU(I).EQ.VISRANKSCA) THEN
               NSCA=NSCA+1
               WRITE(VISFIENAME(I),'(A,I1)') 'scalar', NSCA
            ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC1) THEN
               NVEC=NVEC+1
               WRITE(VISFIENAME(I),'(A,I1)') 'vector', NVEC
            ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC2) THEN
               NVEC=NVEC+1
               WRITE(VISFIENAME(I),'(A,I1)') 'vector', NVEC
            ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC3) THEN
               NVEC=NVEC+1
               WRITE(VISFIENAME(I),'(A,I1)') 'vector', NVEC
            ELSE IF (VISSTRUCTU(I).EQ.VISRANKVEC) THEN
               NVEC=NVEC+1
               WRITE(VISFIENAME(I),'(A,I1)') 'vector', NVEC
            ENDIF
   20    CONTINUE
      ENDIF
C
C     * NUMBER OF DATA ELEMENTS PER TIMESTEP
      ITEMS=1
      DO 30 I=1,VISNDIM
            ITEMS=ITEMS*VISRES(VISCOOR(I))
   30 CONTINUE
C
C     * TYPE CHECKING
      DO 40 I=1,VISNFIELDS+1
         IF (     VISSTRUTYP(I).NE.VISTYPEFLO
     A       .AND.VISSTRUTYP(I).NE.VISTYPEDOU) THEN
            CALL VISPRERR(VISERRTYPE,'VISDXINTHEAD')
         ENDIF
   40 CONTINUE
C
C     * FILENAME OK?
      IF (.NOT.VISNAMEOK(VISFILNAME)) THEN
         CALL VISPRWAR(VISWARNAME,'VISDXINTHEAD')
         VISFILNAME='dxfile.dx'
      ENDIF
C
C     * OPEN FILE FOR DX FILE FORMAT
      FILENAME             = VISFILNAME
      IND=INDEX(VISFILNAME,'.dx')
      IF (IND.EQ.0) THEN
         IND                 = MIN(78,INDEX(VISFILNAME,' '))
         IF (IND.EQ.0)   IND = MIN(78,LEN(VISFILNAME)+1)
         FILENAME(IND:IND+2) = '.dx'
      ENDIF
      OPEN(UNIT=VISDXOUTPU,FILE=FILENAME(1:IND+2),STATUS='UNKNOWN')
C
C     * IRREGULAR GRID
      OBJNUMBER= 1
      TYPE     = VISTYPEFLO
      CATEGORY = VISCATREAL
      RANK     = 1
      SHAPE    = 3
      OFFSET   = 0
C
      WRITE(VISDXOUTPU,310) OBJNUMBER, TYPE, CATEGORY, RANK,
     A                      SHAPE, ITEMS, FORMATS, OFFSET
      WRITE(VISDXOUTPU,300)
C
      IF (VISFORMATS.EQ.VISFORMBIN) THEN
         OFFSET   = OFFSET+ITEMS*VISSIZEOF(TYPE)*SHAPE
      ELSE IF (VISFORMATS.EQ.VISFORMASC) THEN
         OFFSET   = OFFSET+ITEMS*(VISSIZEFORMAT(TYPE)*SHAPE+CR)
      ENDIF
C
C     * REGULAR GRID CONNECTIONS CLASS
      OBJNUMBER= OBJNUMBER+1
      WRITE(VISDXOUTPU,340) OBJNUMBER,'gridconnections counts',
     A                      (VISRES(VISCOOR(I)), I=1,VISNDIM)
      IF (VISNDIM.EQ.3) THEN
         WRITE(VISDXOUTPU,350) '"element type"', '"cubes"'
      ELSE IF (VISNDIM.EQ.2) THEN
         WRITE(VISDXOUTPU,350) '"element type"', '"quads"'
      ELSE
         WRITE(VISDXOUTPU,350) '"element type"', '"lines"'
      ENDIF
      WRITE(VISDXOUTPU,350) '"ref"', '"positions"'
      WRITE(VISDXOUTPU,300)
C
C     * THE DATA CLASS
      IF (.NOT.VISSERIES) THEN
         DO 60 I=1,VISNFIELDS
            OBJNUMBER=OBJNUMBER+1
            TYPE     =VISSTRUTYP(I)
            CATEGORY =VISCATREAL
            IF (VISSTRUCTU(I).EQ.VISRANKSCA) THEN
               RANK  =0
               WRITE(VISDXOUTPU,360) OBJNUMBER, TYPE, CATEGORY,
     A                               RANK, ITEMS, FORMATS, OFFSET
               IF (VISFORMATS.EQ.VISFORMBIN) THEN
                  OFFSET   = OFFSET+ITEMS*VISSIZEOF(TYPE)
               ELSE IF (VISFORMATS.EQ.VISFORMASC) THEN
                  OFFSET   = OFFSET+(ITEMS*VISSIZEFORMAT(TYPE)+CR)
               ENDIF
            ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC1) THEN
               RANK  =1
               SHAPE =1
               WRITE(VISDXOUTPU,310) OBJNUMBER, TYPE, CATEGORY,
     A                        RANK, SHAPE, ITEMS, FORMATS, OFFSET
               IF (VISFORMATS.EQ.VISFORMBIN) THEN
                  OFFSET   = OFFSET+ITEMS*VISSIZEOF(TYPE)*SHAPE
               ELSE IF (VISFORMATS.EQ.VISFORMASC) THEN
                  OFFSET   = OFFSET+ITEMS*(VISSIZEFORMAT(TYPE)*SHAPE+CR)
               ENDIF
            ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC2) THEN
               RANK  =1
               SHAPE =2
               WRITE(VISDXOUTPU,310) OBJNUMBER, TYPE, CATEGORY,
     A                        RANK, SHAPE, ITEMS, FORMATS, OFFSET
               IF (VISFORMATS.EQ.VISFORMBIN) THEN
                  OFFSET   = OFFSET+ITEMS*VISSIZEOF(TYPE)*SHAPE
               ELSE IF (VISFORMATS.EQ.VISFORMASC) THEN
                  OFFSET   = OFFSET+ITEMS*(VISSIZEFORMAT(TYPE)*SHAPE+CR)
               ENDIF
            ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC3) THEN
               RANK  =1
               SHAPE =3
               WRITE(VISDXOUTPU,310) OBJNUMBER, TYPE, CATEGORY,
     A                        RANK, SHAPE, ITEMS, FORMATS, OFFSET
               IF (VISFORMATS.EQ.VISFORMBIN) THEN
                  OFFSET   = OFFSET+ITEMS*VISSIZEOF(TYPE)*SHAPE
               ELSE IF (VISFORMATS.EQ.VISFORMASC) THEN
                  OFFSET   = OFFSET+ITEMS*(VISSIZEFORMAT(TYPE)*SHAPE+CR)
               ENDIF
            ELSE IF (VISSTRUCTU(I).EQ.VISRANKVEC) THEN
               RANK  =1
               SHAPE =3
               WRITE(VISDXOUTPU,310) OBJNUMBER, TYPE, CATEGORY,
     A                        RANK, SHAPE, ITEMS, FORMATS, OFFSET
               IF (VISFORMATS.EQ.VISFORMBIN) THEN
                  OFFSET   = OFFSET+ITEMS*VISSIZEOF(TYPE)*SHAPE
               ELSE IF (VISFORMATS.EQ.VISFORMASC) THEN
                  OFFSET   = OFFSET+ITEMS*(VISSIZEFORMAT(TYPE)*SHAPE+CR)
               ENDIF
            ENDIF
            WRITE(VISDXOUTPU,350) '"dep"', '"positions"'
            WRITE(VISDXOUTPU,300)
C
C           * CLASS FIELD
            VISTMPCHAR=VISFIENAME(I)
            WRITE(VISDXOUTPU,380)
     A            '"'//VISTMPCHAR(1:INDEX(VISFIENAME(I),' ')-1)//'"',
     B            'field'
            WRITE(VISDXOUTPU,400) '"positions"', 1
            WRITE(VISDXOUTPU,400) '"connections"', 2
            WRITE(VISDXOUTPU,400) '"data"',OBJNUMBER
            WRITE(VISDXOUTPU,300)
   60    CONTINUE
C
C        * GROUP CLAUSE
         OBJNUMBER=OBJNUMBER+1
         WRITE(VISDXOUTPU,340) OBJNUMBER, 'group'
C
C        ** THE MEMBERS OF THE GROUP
         DO 70 I=1,VISNFIELDS
            VISTMPCHAR=VISFIENAME(I)
            WRITE(VISDXOUTPU,430)
     A            '"'//VISTMPCHAR(1:INDEX(VISFIENAME(I),' ')-1)//'"',
     B            '"'//VISTMPCHAR(1:INDEX(VISFIENAME(I),' ')-1)//'"'
   70    CONTINUE
C
C        ** GEOMETRY ATTRIBUTE
         WRITE(VISDXOUTPU,450) '"geom"', VISGEOM
C
      ELSE
         DO 80 J=1,VISRES(VISCOORTIM)
            DO 90 I=1,VISNFIELDS
               OBJNUMBER=OBJNUMBER+1
               TYPE     =VISTYPEDOU
               CATEGORY =VISCATREAL
               RANK     =0
C
C              * SERIES POSITION
               WRITE(VISDXOUTPU,360) OBJNUMBER, TYPE, CATEGORY,
     A                               RANK, 1, FORMATS, OFFSET
C
               IF (VISFORMATS.EQ.VISFORMBIN) THEN
                  OFFSET   = OFFSET+VISSIZEOF(TYPE)
               ELSE IF (VISFORMATS.EQ.VISFORMASC) THEN
                  OFFSET   = OFFSET+(VISSIZEFORMAT(TYPE)+CR)
               ENDIF
               WRITE(VISDXOUTPU,300)
C
C              * DATA ARRAY
               OBJNUMBER=OBJNUMBER+1
               TYPE     =VISSTRUTYP(I)
               CATEGORY =VISCATREAL
               IF (VISSTRUCTU(I).EQ.VISRANKSCA) THEN
                  RANK  =0
C
                  WRITE(VISDXOUTPU,360) OBJNUMBER, TYPE, CATEGORY,
     A                                  RANK, ITEMS, FORMATS, OFFSET
C
                  IF (VISFORMATS.EQ.VISFORMBIN) THEN
                     OFFSET   = OFFSET+ITEMS*VISSIZEOF(TYPE)
                  ELSE IF (VISFORMATS.EQ.VISFORMASC) THEN
                     OFFSET   = OFFSET+ITEMS*(VISSIZEFORMAT(TYPE)+CR)
                  ENDIF
               ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC1) THEN
                  RANK  =1
                  SHAPE =1
C
                  WRITE(VISDXOUTPU,310) OBJNUMBER, TYPE, CATEGORY,
     A                                  RANK, SHAPE, ITEMS, FORMATS,
     B                                  OFFSET
C
                  IF (VISFORMATS.EQ.VISFORMBIN) THEN
                     OFFSET = OFFSET+ITEMS*VISSIZEOF(TYPE)*SHAPE
                  ELSE IF (VISFORMATS.EQ.VISFORMASC) THEN
                    OFFSET = OFFSET+ITEMS*(VISSIZEFORMAT(TYPE)*SHAPE+CR)
                  ENDIF
               ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC2) THEN
                  RANK  =1
                  SHAPE =2
C
                  WRITE(VISDXOUTPU,310) OBJNUMBER, TYPE, CATEGORY,
     A                                  RANK, SHAPE, ITEMS, FORMATS,
     B                                  OFFSET
C
                  IF (VISFORMATS.EQ.VISFORMBIN) THEN
                     OFFSET = OFFSET+ITEMS*VISSIZEOF(TYPE)*SHAPE
                  ELSE IF (VISFORMATS.EQ.VISFORMASC) THEN
                    OFFSET = OFFSET+ITEMS*(VISSIZEFORMAT(TYPE)*SHAPE+CR)
                  ENDIF
               ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC3) THEN
                  RANK  =1
                  SHAPE =3
C
                  WRITE(VISDXOUTPU,310) OBJNUMBER, TYPE, CATEGORY,
     A                                  RANK, SHAPE, ITEMS, FORMATS,
     B                                  OFFSET
C
                  IF (VISFORMATS.EQ.VISFORMBIN) THEN
                     OFFSET = OFFSET+ITEMS*VISSIZEOF(TYPE)*SHAPE
                  ELSE IF (VISFORMATS.EQ.VISFORMASC) THEN
                    OFFSET = OFFSET+ITEMS*(VISSIZEFORMAT(TYPE)*SHAPE+CR)
                  ENDIF
               ELSE IF (VISSTRUCTU(I).EQ.VISRANKVEC) THEN
                  RANK  =1
                  SHAPE =3
C
                  WRITE(VISDXOUTPU,310) OBJNUMBER, TYPE, CATEGORY,
     A                                  RANK, SHAPE, ITEMS, FORMATS,
     B                                  OFFSET
C
                  IF (VISFORMATS.EQ.VISFORMBIN) THEN
                     OFFSET = OFFSET+ITEMS*VISSIZEOF(TYPE)*SHAPE
                  ELSE IF (VISFORMATS.EQ.VISFORMASC) THEN
                    OFFSET = OFFSET+ITEMS*(VISSIZEFORMAT(TYPE)*SHAPE+CR)
                  ENDIF
               ENDIF
               WRITE(VISDXOUTPU,350) '"dep"', '"positions"'
               WRITE(VISDXOUTPU,300)
C
C              * CLASS FIELD
               OBJNUMBER=OBJNUMBER+1
               VISWORKINT(I+(J-1)*VISNFIELDS)=OBJNUMBER
               WRITE(VISDXOUTPU,340) OBJNUMBER, 'field'
               WRITE(VISDXOUTPU,400) '"positions"', 1
               WRITE(VISDXOUTPU,400) '"connections"', 2
               WRITE(VISDXOUTPU,400) '"data"',OBJNUMBER-1
               WRITE(VISDXOUTPU,390) '"series position"', OBJNUMBER-2
               WRITE(VISDXOUTPU,300)
   90       CONTINUE
   80    CONTINUE
C
C        * SERIES CLAUSE
         DO 100 I=1,VISNFIELDS
            VISTMPCHAR=VISFIENAME(I)
            MEMNUMBER=-1
            WRITE(VISDXOUTPU,380)
     A               '"'//VISTMPCHAR(1:INDEX(VISFIENAME(I),' ')-1)//'"',
     B               'series'
            DO 110 J=1,VISRES(VISCOORTIM)
               MEMNUMBER=MEMNUMBER+1
               WRITE(VISDXOUTPU,410) MEMNUMBER, VISPOSITIM(J),
     A                               VISWORKINT(I+(J-1)*VISNFIELDS)
  110       CONTINUE
            WRITE(VISDXOUTPU,300)
  100    CONTINUE
C
C
C        * GROUP CLAUSE
         OBJNUMBER=OBJNUMBER+1
         WRITE(VISDXOUTPU,340) OBJNUMBER, 'group'
C
C        ** THE MEMBERS OF THE GROUP
         DO 120 I=1,VISNFIELDS
            VISTMPCHAR=VISFIENAME(I)
            WRITE(VISDXOUTPU,430)
     A               '"'//VISTMPCHAR(1:INDEX(VISFIENAME(I),' ')-1)//'"',
     B               '"'//VISTMPCHAR(1:INDEX(VISFIENAME(I),' ')-1)//'"'
  120    CONTINUE
C
C        ** GEOMETRY ATTRIBUTE
         WRITE(VISDXOUTPU,450) '"geom"', VISGEOM
      ENDIF
C
C     * END CLAUSE
      WRITE(VISDXOUTPU,440) 'end'
C
C     * CLOSE FILE (ONLY FOR BINARY FORMAT)
      IF (VISFORMATS.EQ.VISFORMBIN) THEN
         CLOSE(UNIT=VISDXOUTPU)
      ENDIF
C
C     * FORMATS
  300 FORMAT('#')
  310 FORMAT('object',1X,I4,1X,'class array type',1X,A,1X,'category',1X,
     A       A,1X,'rank',1X,I2,1X,'shape',1X,I4,1X,'items',1X,I7,1X,A,
     B       1X,'data',1X,I10)
  330 FORMAT(3F12.4)
  340 FORMAT('object',1X,I5,1X,'class',1X,A,1X,I3,3(1X,I3))
  350 FORMAT('attribute',1X,A,1X,'string',1X,A)
  360 FORMAT('object',1X,I5,1X,'class array type',1X,A,1X,'category',1X,
     A       A,1X,'rank',1X,I2,1X,'items',1X,I7,1X,A,1X,'data',1X,I10)
  380 FORMAT('object',1X,A,1X,'class',1X,A,1X,I4,3(1X,I4))
  390 FORMAT('attribute',1X,A,1X,'value',1X,I5)
  400 FORMAT('component',1X,A,1X,'value',1X,I5)
  410 FORMAT('member',1X,I5,1X,'position',1X,F12.4,1X,'value',1X,I5)
  420 FORMAT('member',1X,A,1X,'value',1X,I5)
  430 FORMAT('member',1X,A,1X,'value',1X,A)
  440 FORMAT(A)
  450 FORMAT('attribute',1X,A,1X,'number',1X,I1)
      END
