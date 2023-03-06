C
C
C     ==================================================================
C     =                   LOW-LEVEL ROUTINES BELOW.                    =
C     =                              |                                 =
C     =                              |                                 =
C     =                              |                                 =
C     =                            \   /                               =
C     =                             \ /                                =
C     =                              .                                 =
C     ==================================================================
*DECK VISDXGENHEAD
      SUBROUTINE VISDXGENHEAD(VISPOSITIM, VISNT,      VISFILNAME,
     A                        VISRES,     VISNDIM,    VISCOOR,
     B                        VISSERIES,  VISNFIELDS, VISFIENAME,
     C                        VISSTRUCTU, VISSTRUTYP, VISFORMATS)
C     ******************************************************************
C     * THIS ROUTINE WRITES THE HEADER FILE FOR DX'S GENERAL ARRAY FILE*
C     * FORMAT. THE CORRESPONDING DATA FILE HAS TO BE WRITTEN WITH A   *
C     * CALL TO VISDATAOUT. THE DATA CAN BE IMMEDIATELY READ INTO DX   *
C     * USING THE 'IMPORT' MODULE (FILENAME=<VISFILNAME>.gai).         *
C     * A DESCRIPTION OF THE GENERAL ARRAY FILE FORMAT CAN BE FOUND IN *
C     * THE USER'S MANUAL OF DX.                                       *
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
C     * VISFILNAME CONTAINS THE BASENAME OF THE OUTPUT FILES.          *
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
C     *            IN ASCII OR IN BINARY FORMAT                        *
C     *            VISFORMATS=VISFORMASC -> ASCII,                     *
C     *            VISFORMATS=VISFORMBIN -> BINARY.                    *
C     * VISSTRUTYP CONTAINS THE PRIMITIVE TYPE OF THE FIELD            *
C     *            VISTYPEDOU -> DOUBLE                                *
C     *            VISTYPEFLO -> FLOAT                                 *
C     *                                                                *
C     * VISDXGENHEAD IS CALLED BY THE USER PROGRAM.                    *
C     ******************************************************************
C
C#include "comcoor"
C#include "comerrn"
C#include "comrank"
C#include "cominou"
C#include "comform"
C#include "comtype"

	INCLUDE 'INCLUDE/comcoor'
	INCLUDE 'INCLUDE/comerrn'
	INCLUDE 'INCLUDE/comrank'
	INCLUDE 'INCLUDE/cominou'
	INCLUDE 'INCLUDE/comform'
	INCLUDE 'INCLUDE/comtype'



C
      LOGICAL       VISSERIES
      INTEGER       VISRES(*),  VISCOOR(*),     VISNDIM,
     A              VISNFIELDS, VISSTRUCTU(*),  VISFORMATS,
     B              VISNT
      CHARACTER*(*) VISFILNAME, VISFIENAME(*),  VISSTRUTYP(*)
      REAL          VISPOSITIM(*)
C
C     *COUNTERS
      INTEGER I
C
C     * LOCAL VARIABLES
      REAL          START, DELTA
      INTEGER       IND,
     A              NSCA, NVEC
      CHARACTER     HEADERNAME*80
      LOGICAL       NONAMES
C
C     * TO AVOID ANOTHER WORKING ARRAY THIS VARIABLE IS DEFINED FOR 50
C       FIELDS WHICH SHOULD SUFFICE.
      CHARACTER     VISTMPCHAR*450
C
C     * FUNCTIONS
      LOGICAL       VISNAMEOK
C
C     * CHECKING DIMENSIONS
      IF (VISNDIM.LT.1.OR.VISNDIM.GT.3) THEN
         CALL VISPRERR(VISERRNDIM,'VISDXGENHEAD')
      ENDIF
C
      DO 10 I=1,VISNDIM
         IF (VISCOOR(I).EQ.VISCOORTIM) THEN
            CALL VISPRERR(VISERRSERI,'VISDXGENHEAD')
         ENDIF
   10 CONTINUE
C
C     * SETTING THE POSITIONS AS THE LAST FIELD
C     * DX DEMANDS THE COORDINATES TO BE FLOATS
      VISSTRUCTU(VISNFIELDS+1)=VISRANKVC3
      VISSTRUTYP(VISNFIELDS+1)=VISTYPEFLO
      VISFIENAME(VISNFIELDS+1)='locations'
C
C     * FILE NAME OK?
      IF (.NOT.VISNAMEOK(VISFILNAME)) THEN
         CALL VISPRWAR(VISWARNAME,'VISDXGEN')
         VISFILNAME='dxfile'
      ENDIF
C
C     * OPEN HEADER FILE FOR GENERAL ARRAY IMPORTER FORMAT
      HEADERNAME           = VISFILNAME
      IND                  = INDEX(VISFILNAME,'.gai')
      IF (IND.EQ.0) THEN
         IND                   = MIN(77,INDEX(VISFILNAME,' '))
         IF (IND.EQ.0)     IND = MIN(77,LEN(VISFILNAME)+1)
         HEADERNAME(IND:IND+3) = '.gai'
      ENDIF
      OPEN(UNIT=VISDXOUTPU,FILE=HEADERNAME(1:IND+3),STATUS='UNKNOWN')
C
C     * NAME OF DATAFILE
      WRITE(VISDXOUTPU,200) 'file', VISFILNAME(1:IND-1)
C
C     * GRID
      WRITE(VISDXOUTPU,210) 'grid', VISRES(VISCOOR(1)),
     A                      ('x', VISRES(VISCOOR(I)), I=2,VISNDIM)
C
C     * SERIES (VISCOORTIM IS CONSIDERED AS SERIES ELEMENT ONLY.)
      IF (VISSERIES) THEN
         IF (VISRES(VISCOORTIM).NE.1) THEN
            START= VISPOSITIM(1)
            DELTA=(VISPOSITIM(VISRES(VISCOORTIM))-VISPOSITIM(1))
     A            /REAL(VISRES(VISCOORTIM)-1)
            WRITE(VISDXOUTPU,220) 'series', VISRES(VISCOORTIM),
     A                             START,   DELTA
         ELSE
            WRITE(VISDXOUTPU,220) 'series', 1, VISPOSITIM(1), 0.0
         ENDIF
      ENDIF
C
C     * FIELD NAMES
      NONAMES=.FALSE.
      DO 40 I=1,VISNFIELDS
         IF (.NOT.VISNAMEOK(VISFIENAME(I))) THEN
            NONAMES=.TRUE.
         ENDIF
   40 CONTINUE
C
      IF (NONAMES) THEN
         NSCA=0
         NVEC=0
         DO 50 I=1,VISNFIELDS
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
   50    CONTINUE
      ENDIF
C
      WRITE(VISDXOUTPU,230) 'field',
     A           VISFIENAME(1)(1:INDEX(VISFIENAME(1),' ')-1),
     B      (',',VISFIENAME(I)(1:INDEX(VISFIENAME(I),' ')-1),
     C                                       I=2,VISNFIELDS+1)
C
C     * STRUCTURE
      INDSN=1
      IF (VISSTRUCTU(1).EQ.VISRANKSCA) THEN
         WRITE(VISTMPCHAR(INDSN:INDSN+5),'(A6)') 'scalar'
         INDSN = INDSN+6
      ELSE IF (VISSTRUCTU(1).EQ.VISRANKVC1) THEN
         WRITE(VISTMPCHAR(INDSN:INDSN+7),'(A8)') '1-vector'
         INDSN = INDSN+8
      ELSE IF (VISSTRUCTU(1).EQ.VISRANKVC2) THEN
         WRITE(VISTMPCHAR(INDSN:INDSN+7),'(A8)') '2-vector'
         INDSN = INDSN+8
      ELSE IF (VISSTRUCTU(1).EQ.VISRANKVC3) THEN
         WRITE(VISTMPCHAR(INDSN:INDSN+7),'(A8)') '3-vector'
         INDSN = INDSN+8
      ELSE IF (VISSTRUCTU(1).EQ.VISRANKVEC) THEN
         WRITE(VISTMPCHAR(INDSN:INDSN+7),'(A8)') '3-vector'
         INDSN = INDSN+8
      ENDIF
      DO 70 I=2,VISNFIELDS
         IF (VISSTRUCTU(I).EQ.VISRANKSCA) THEN
            WRITE(VISTMPCHAR(INDSN:INDSN+7),'('','',1X,A6)') 'scalar'
            INDSN = INDSN+8
         ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC1) THEN
            WRITE(VISTMPCHAR(INDSN:INDSN+7),'('','',1X,A8)') '1-vector'
            INDSN = INDSN+10
         ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC2) THEN
            WRITE(VISTMPCHAR(INDSN:INDSN+8),'('','',1X,A8)') '2-vector'
            INDSN = INDSN+10
         ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC3) THEN
            WRITE(VISTMPCHAR(INDSN:INDSN+9),'('','',1X,A8)') '3-vector'
            INDSN = INDSN+10
         ELSE IF (VISSTRUCTU(I).EQ.VISRANKVEC) THEN
            WRITE(VISTMPCHAR(INDSN:INDSN+9),'('','',1X,A8)') '3-vector'
            INDSN = INDSN+10
         ENDIF
   70 CONTINUE
C
C     ** LOCATIONS STRUCTURE
      WRITE(VISTMPCHAR(INDSN:INDSN+9),'('','',1X,A8)') '3-vector'
      INDSN = INDSN+9
C
      WRITE(VISDXOUTPU,200) 'structure', VISTMPCHAR(1:INDSN)
C
C     * FORMAT
      IF (VISFORMATS.EQ.VISFORMASC) THEN
         WRITE(VISDXOUTPU,200) 'format', 'ascii'
      ELSE IF (VISFORMATS.EQ.VISFORMBIN) THEN
         WRITE(VISDXOUTPU,200) 'format', 'binary'
      ELSE
         CALL VISPRWAR(VISWARFORM,'VISDXGENHEAD')
         VISFORMATS=VISFORMASC
         WRITE(VISDXOUTPU,200) 'format', 'ascii'
      ENDIF
C
C     * TYPE
      DO 80 I=1,VISNFIELDS+1
         IF ( VISSTRUTYP(I)(1:INDEX(VISSTRUTYP(I),' ')-1).NE.VISTYPEFLO
     A   .AND.VISSTRUTYP(I)(1:INDEX(VISSTRUTYP(I),' ')-1).NE.VISTYPEDOU)
     B      THEN
            CALL VISPRERR(VISERRTYPE,'VISDXGENHEAD')
         ENDIF
   80 CONTINUE
C
      WRITE(VISDXOUTPU,230) 'type',
     A            VISSTRUTYP(1)(1:INDEX(VISSTRUTYP(1),' ')-1),
     B            (',',VISSTRUTYP(I)(1:INDEX(VISSTRUTYP(I),' ')-1),
     C            I=2, VISNFIELDS+1)
C
C     * MAJORITY
      WRITE(VISDXOUTPU,200) 'majority', 'column'
C
C     * HEADER
      IF (VISFORMATS.EQ.VISFORMASC) THEN
         WRITE(VISDXOUTPU,240) 'header','lines', VISNFIELDS+1
      ELSE IF (VISFORMATS.EQ.VISFORMBIN) THEN
         WRITE(VISDXOUTPU,250) 'header','bytes', 0
      ENDIF
C
C     * INTERLEAVING
      WRITE(VISDXOUTPU,200) 'interleaving', 'field'
C
C     * END
      WRITE(VISDXOUTPU,260) 'end'
C
C     * CLOSE FILE
      CLOSE(UNIT=VISDXOUTPU)
C
C     * FORMATS
  200 FORMAT(A,' = ',A)
  210 FORMAT(A,' = ',I3,3(1X,A,1X,I3))
  220 FORMAT(A,' = ',(I4,2(',',F12.4)))
  230 FORMAT(A,' = ',A,8(A,1X,A))
  240 FORMAT(A,' = ',A,1X,I5)
  250 FORMAT(A,' = ',A,1X,I5)
  260 FORMAT(A)
      END
