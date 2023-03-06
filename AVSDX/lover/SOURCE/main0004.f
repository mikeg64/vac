C
*DECK VISGOFORIT
      SUBROUTINE VISGOFORIT(VISPACKAGE, VISHOW,     VISMATDATA,
     A                      VISPOSX,    VISPOSY,    VISPOSZ,
     B                      VISPOSITIM, VISFILNAME, VISFIENAME,
     C                      VISSTRUCTU, VISSTRUTYP, VISWORK,
     D                      VISWORKINT)
C     ******************************************************************
C     * THIS IS THE HIGH LEVEL DRIVER ROUTINE WHICH CALLS THE VAROUS   *
C     * LOW LEVEL ROUTINES. FOR MORE SPECIFIC INFORMATION SEE THE      *
C     * SPECIFIC SUBROUTINES.                                          *
C     *                                                                *
C     * ON INPUT:                                                      *
C     * VISPACKAGE DETERMINES FOR WHICH VISUALIZATION PACKAGE TO       *
C     *            FORMAT: 'A' --> AVS, 'D' --> DX                     *
C     * VISHOW     DETERMINES THE FILE FORMAT                          *
C     *            VISHOW=VISGENERAL -> GENERAL  FILE FORMAT,          *
C     *            VISHOW=VISINTERNA -> INTERNAL FILE FORMAT.          *
C     * VISMATDATA CONTAINS THE DATA.                                  *
C     * VISPOSX    CONTAINS THE X-COORDINATES.                         *
C     * VISPOSY    CONTAINS THE Y-COORDINATES.                         *
C     * VISPOSZ    CONTAINS THE Z-COORDINATES.                         *
C     * VISPOSITIM CONTAINS THE TIME COORDINATES (LABELS).             *
C     * VISFILNAME CONTAINS THE FILENAME OF THE OUTPUT FILE. IF IT     *
C     *            DOES NOT ALLREADY CONTAIN THE .fld (AVS) or .dx     *
C     *            (DX INTERNAL) POSTFIX IT WILL BE ADDED.             *
C     * VISFIENAME CONTAINS THE NAMES OF THE DIFFERENT FIELDS.         *
C     * VISSTRUCTU DETERMINES THE TYPE OF THE FIELDS                   *
C     *            VISSTRUCTU()=VISRANKSCA -> SCALAR FIELD,            *
C     *            VISSTRUCTU()=VISRANKVC1 -> 1-VECTOR FIELD.          *
C     *            VISSTRUCTU()=VISRANKVC2 -> 2-VECTOR FIELD.          *
C     *            VISSTRUCTU()=VISRANKVC3 -> 3-VECTOR FIELD.          *
C     *            VISSTRUCTU()=VISRANKVEC -> 3-VECTOR FIELD.          *
C     * VISSTRUTYP CONTAINS THE PRIMITIVE TYPE OF THE FIELD            *
C     *            VISTYPEDOU -> DOUBLE                                *
C     *            VISTYPEFLO -> FLOAT                                 *
C     * VISWORK    IS A REAL WORKING ARRAY DECLARED BY THE USER AND    *
C     *            WITH A DIMENSION >=  VISDIM1*                       *
C     *                               4*VISDIM2*VISDIM3*VISDIM4        *
C     * VISWORKINT IS A INTEGER WORK ARRAY DECLARED BY THE USER AND    *
C     *            WITH A DIMENSION >=  VISDIM5*VISNFLD                *
C     *                                                                *
C     * CALLS:                                                         *
C     * VISDXGENHEAD, VISDXINTHEAD, VISAVSHEAD, VISWRDXCOO, VISWRDXDAT,*
C     * VISDATAOUT,   VISPRERR,     VISLOVER                           *
C     *                                                                *
C     * FUNCTIONS:                                                     *
C     * VISOPENDX, VISOPENAVS                                          *
C     ******************************************************************
C#include "comrank"
C#include "comform"
C#include "comerrn"
	INCLUDE 'INCLUDE/comrank'
	INCLUDE 'INCLUDE/comform'
	INCLUDE 'INCLUDE/comerrn'
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
      CHARACTER*(*) VISPACKAGE, VISFILNAME, VISFIENAME(*), VISSTRUTYP(*)
      INTEGER       VISHOW, VISSTRUCTU(*), VISWORKINT(*)
      REAL          VISMATDATA(VISDIM1,VISDIM2,VISDIM3,VISDIM4,VISDIM5),
     A              VISPOSX(VISDIM2,VISDIM3,VISDIM4),
     B              VISPOSY(VISDIM2,VISDIM3,VISDIM4),
     C              VISPOSZ(VISDIM2,VISDIM3,VISDIM4),
     D              VISPOSITIM(VISDIM5),
     E              VISWORK(*)
C
C     * COUNTERS
      INTEGER I
C
C     * LOCAL VARIABLES
      INTEGER IOUNIT, OFFSET
C
C     * SAVING LOCAL VARIABLES
      SAVE IOUNIT
C
C     * EXTERNAL SUBROUTINES
      EXTERNAL VISDXGENHEAD, VISDXINTHEAD, VISDXWRCOO, VISDXWRDAT,
     A         VISAVSHEAD,   VISDATAOUT,   VISPRERR,   VISLOVER
C
C     * FUNCTIONS
      INTEGER VISOPENDX, VISOPENAVS
C
C     * COPYRIGHT NOTICE
      IF (VISFIRS) CALL VISLOVER
C
C     * DETERMINING THE VECTOR LENGTH
      VISVLEN = 0
      DO 10 I=1, VISNFLD
         IF (VISSTRUCTU(I).EQ.VISRANKSCA) THEN
            VISVLEN = VISVLEN + 1
         ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC1) THEN
            VISVLEN = VISVLEN + 1
         ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC2) THEN
            VISVLEN = VISVLEN + 2
         ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC3) THEN
            VISVLEN = VISVLEN + 3
         ELSE IF (VISSTRUCTU(I).EQ.VISRANKVEC) THEN
            VISVLEN = VISVLEN + 3
         ELSE
            CALL VISPRERR(VISERRWRST,'VISGOFORIT')
         ENDIF
   10 CONTINUE
C
C     * DETERMINING OFFSET FOR VISWORK USED BY VISDXWRDAT
      OFFSET = 3*VISDIM2*VISDIM3*VISDIM4+1
C
      IF (VISPACKAGE.EQ.'A') GO TO 30
      IF (VISPACKAGE.NE.'D') THEN
         CALL VISPRERR(VISERRPACK,'VISGOFORIT')
      ENDIF
C
      IF (VISHOW.EQ.VISGENERAL) GO TO 20
      IF (VISHOW.NE.VISINTERNA) THEN
         CALL VISPRERR(VISERRFIFO,'VISGOFORIT')
      ENDIF
C
C     * DRIVER PART FOR DX INTERNAL FILE FORMAT
C
C     ** first call ?
      IF (.NOT.VISFIRS) GOTO 15
C
      CALL VISDXINTHEAD(VISPOSITIM, VISDIM5,    VISFILNAME, VISRESO,
     A                  VISDIMS,    VISCORD,    VISSERI,    VISNFLD,
     B                  VISFIENAME, VISSTRUCTU, VISSTRUTYP, VISFORM,
     C                  VISWORKINT, VISGEOM)
      CALL VISDXWRCOO(VISPOSX,       VISPOSY,    VISPOSZ,    VISDIM2,
     A                VISDIM3,       VISDIM4, VISFILNAME,    VISRESO,
     B                VISDIMS,       VISCORD,    VISSLIC,    VISNFLD,
     C                VISSTRUCTU, VISSTRUTYP,    VISFORM,    VISWORK)
      VISFIRS=.FALSE.
      IF (VIS1TST) THEN
         VISNTST    = 0
         VISTIME1   = VISPOSITIM(1)
         VISTRES    = VISRESO(4)
         VISRESO(4) = 1
         RETURN
      ENDIF
C
C     ** all calls
   15 CONTINUE
C
C     * RECONFIGURING FOR ONE TIMESTEP A TIME
      IF (VIS1TST) THEN
         VISNTST       = VISNTST+1
         IF (VISNTST.EQ.VISTRES) VISLAST=.TRUE.
         VISPOSITIM(1) = VISPOSITIM(VISNTST)
      ENDIF
C
      CALL VISDXWRDAT(VISMATDATA, VISPOSITIM,    VISDIM1,    VISDIM2,
     A                VISDIM3,       VISDIM4,    VISDIM5, VISFILNAME,
     B                VISRESO,       VISDIMS,    VISCORD,    VISSERI,
     C                VISSLIC,       VISVLEN,    VISNFLD, VISSTRUCTU,
     D                VISSTRUTYP,    VISFORM, VISWORK(OFFSET))
C
C     * RESTORING VARIABLES ON LAST CALL
      IF (VIS1TST.AND.VISLAST) THEN
         VISRESO(4)    = VISTRES
         VISPOSITIM(1) = VISTIME1
         VISFIRS       = .TRUE.
         VISLAST       = .FALSE.
      ENDIF
      RETURN
C
   20 CONTINUE
C
C     * DRIVER PART FOR DX GENERAL FILE FORMAT
C
C     ** first call ?
      IF (.NOT.VISFIRS) GOTO 25
C
      CALL VISDXGENHEAD(VISPOSITIM, VISDIM5,    VISFILNAME, VISRESO,
     A                  VISDIMS,    VISCORD,    VISSERI,    VISNFLD,
     B                  VISFIENAME, VISSTRUCTU, VISSTRUTYP, VISFORM)
C
      IOUNIT=VISOPENDX(VISFILNAME,VISFORM)
C
      VISFIRS=.FALSE.
      IF (VIS1TST) THEN
         VISNTST    = 0
         VISTIME1   = VISPOSITIM(1)
         VISTRES    = VISRESO(4)
         VISRESO(4) = 1
         RETURN
      ENDIF
C
C     ** all calls
   25 CONTINUE
C
C     * RECONFIGURING FOR ONE TIMESTEP A TIME
      IF (VIS1TST) THEN
         VISNTST       = VISNTST+1
         IF (VISNTST.EQ.VISTRES) VISLAST=.TRUE.
         VISPOSITIM(1) = VISPOSITIM(VISNTST)
      ENDIF
C
C     * last call ?
      IF (VISLAST)THEN
         IOUNIT=-ABS(IOUNIT)
      ENDIF
      CALL VISDATAOUT(VISMATDATA,    VISPOSX,    VISPOSY,    VISPOSZ,
     A                VISPOSITIM,    VISDIM1,    VISDIM2,    VISDIM3,
     B                VISDIM4,       VISDIM5, VISFILNAME,    VISRESO,
     C                VISDIMS,       VISCORD,    VISSERI,    VISSLIC,
     D                VISVLEN,       VISNFLD, VISFIENAME, VISSTRUCTU,
     E                VISSTRUTYP,    VISFORM,     VISXDR,     IOUNIT,
     F                VISWORK(OFFSET),           VISWORK)
C
C     * RESTORING VARIABLES ON LAST CALL
      IF (VIS1TST.AND.VISLAST) THEN
         VISRESO(4)    = VISTRES
         VISPOSITIM(1) = VISTIME1
         VISFIRS       = .TRUE.
         VISLAST       = .FALSE.
      ENDIF
      RETURN
C
   30 CONTINUE
C
      IF (VISHOW.EQ.VISGENERAL) GO TO 40
      IF (VISHOW.NE.VISINTERNA) THEN
         CALL VISPRERR(VISERRFIFO,'VISGOFORIT')
      ENDIF
C     * DRIVER PART FOR AVS INTERNAL FILE FORMAT
C
C     ** first call ?
      IF (.NOT.VISFIRS) GOTO 35
C
      CALL VISAVSHEAD(VISFILNAME,    VISRESO,    VISDIMS,    VISCORD,
     A                VISVLEN,       VISNFLD, VISFIENAME, VISSTRUCTU,
     B                VISSTRUTYP,    VISFORM,     VISXDR, VISINTERNA)
C
      VISFIRS=.FALSE.
      IF (VIS1TST) THEN
         VISNTST    = 0
         VISTIME1   = VISPOSITIM(1)
         VISTRES    = VISRESO(4)
         VISRESO(4) = 1
         RETURN
      ENDIF
C
C     ** all calls
   35 CONTINUE
C
C     * RECONFIGURING FOR ONE TIMESTEP A TIME
      IF (VIS1TST) THEN
         VISNTST       = VISNTST+1
         IF (VISNTST.EQ.VISTRES) VISLAST=.TRUE.
         VISPOSITIM(1) = VISPOSITIM(VISNTST)
         VISSLIC(4)    = 1
      ENDIF
C
      IOUNIT=99
      CALL VISDATAOUT(VISMATDATA,    VISPOSX,    VISPOSY,    VISPOSZ,
     A                VISPOSITIM,    VISDIM1,    VISDIM2,    VISDIM3,
     B                VISDIM4,       VISDIM5, VISFILNAME,    VISRESO,
     C                VISDIMS,       VISCORD,    VISSERI,    VISSLIC,
     D                VISVLEN,       VISNFLD, VISFIENAME, VISSTRUCTU,
     E                VISSTRUTYP,    VISFORM,     VISXDR,     IOUNIT,
     F                VISWORK(OFFSET),           VISWORK)
C
C     * RESTORING VARIABLES ON LAST CALL
      IF (VIS1TST.AND.VISLAST) THEN
         VISRESO(4)    = VISTRES
         VISPOSITIM(1) = VISTIME1
         VISFIRS       = .TRUE.
         VISLAST       = .FALSE.
      ENDIF
      RETURN
C
   40 CONTINUE
C
C     * DRIVER PART FOR AVS GENERAL FILE FORMAT
C
C     ** first call ?
      IF (.NOT.VISFIRS) GOTO 45
C
      CALL VISAVSHEAD(VISFILNAME,    VISRESO,    VISDIMS,    VISCORD,
     A                VISVLEN,       VISNFLD, VISFIENAME, VISSTRUCTU,
     B                VISSTRUTYP,    VISFORM,     VISXDR, VISGENERAL)
C
      IOUNIT=VISOPENAVS(VISFILNAME,VISFORM)
C
      VISFIRS=.FALSE.
      IF (VIS1TST) THEN
         VISNTST    = 0
         VISTIME1   = VISPOSITIM(1)
         VISTRES    = VISRESO(4)
         VISRESO(4) = 1
         RETURN
      ENDIF
C
C     ** all calls
   45 CONTINUE
C
C     * RECONFIGURING FOR ONE TIMESTEP A TIME
      IF (VIS1TST) THEN
         VISNTST       = VISNTST+1
         IF (VISNTST.EQ.VISTRES) VISLAST=.TRUE.
         VISPOSITIM(1) = VISPOSITIM(VISNTST)
         VISSLIC(4)    = 1
      ENDIF
C
C     * last call ?
      IF (VISLAST)THEN
         IOUNIT=-ABS(IOUNIT)
      ENDIF
C
      CALL VISDATAOUT(VISMATDATA,    VISPOSX,    VISPOSY,    VISPOSZ,
     A                VISPOSITIM,    VISDIM1,    VISDIM2,    VISDIM3,
     B                VISDIM4,       VISDIM5, VISFILNAME,    VISRESO,
     C                VISDIMS,       VISCORD,    VISSERI,    VISSLIC,
     D                VISVLEN,       VISNFLD, VISFIENAME, VISSTRUCTU,
     E                VISSTRUTYP,    VISFORM,     VISXDR,     IOUNIT,
     F                VISWORK(OFFSET),           VISWORK)
C
C     * RESTORING VARIABLES ON LAST CALL
      IF (VIS1TST.AND.VISLAST) THEN
         VISRESO(4)    = VISTRES
         VISPOSITIM(1) = VISTIME1
         VISFIRS       = .TRUE.
         VISLAST       = .FALSE.
      ENDIF
      RETURN
      END
