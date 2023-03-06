*DECK COPYRIGHT
************************************************************************
* MAIN program: @(#)liblover.F	2.1
*                                                                      *
* Author: A.J.C. Belien                                                *
*         FOM-Instituut voor Plasmafysica `Rijnhuizen'                 *
*         P.O. Box 1207                                                *
*         3430 BE Nieuwegein                                           *
*                                                                      *
* (C) 1994-1996: FOM-Institute for Plasma Physics `Rijnhuizen',        *
*                Academic Computing Services Amsterdam (SARA)          *
* (C) 1997-2000: FOM-Institute for Plasma Physics `Rijnhuizen'         *
************************************************************************
*DECK VISLOVER
      SUBROUTINE VISLOVER
C     ******************************************************************
C     ******************************************************************
C
#include "cominou"
C
      LOGICAL VISPRINTCOPYRIGHT
      SAVE VISPRINTCOPYRIGHT
      DATA VISPRINTCOPYRIGHT/.TRUE./

      IF (VISPRINTCOPYRIGHT) THEN
         WRITE(VISSTDOUTP,100)
         VISPRINTCOPYRIGHT=.FALSE.
      END IF
C
C     * FORMATS
  100 FORMAT(/
     C/'**************************************************************',
     C/'*                                                            *',
     C/'*         Library Of Visualization Exporting Routines        *',
     C/'*                                                            *',
     C/'*                             BY                             *',
     C/'*                                                            *',
     C/'*                       SANDER BELIEN                        *',
     C/'*                                                            *',
     C/'*  (C) 1994-1996:                                            *',
     C/'*            FOM-Institute for Plasma Physics Rijnhuizen,    *',
     C/'*            Academic Computing Services Amsterdam (SARA)    *',
     C/'*  (C) 1997-2000:                                            *',
     C/'*            FOM-Institute for Plasma Physics Rijnhuizen     *',
     C/'**************************************************************')
      END
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
#include "comform"
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
C
*DECK VISONESTEP
      SUBROUTINE VISONESTEP
C     ******************************************************************
C     * THIS ROUTIME IS CALLED TO SETUP ITEMS IF ONE WANTS TO OUTPUT   *
C     * ONE TIMESTEP A TIME.                                           *
C     ******************************************************************
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
      VIS1TST=.TRUE.
      VISFIRS=.TRUE.
      VISLAST=.FALSE.
C
      RETURN
      END
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
#include "comrank"
#include "comform"
#include "comerrn"
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
C
*DECK VISCROP7
      SUBROUTINE VISCROP7(VISMATDATA, VISPOSIRAD, VISPOSIPOL, 
     A                    VISPOSITOR, VISPOSITIM, VISCROPSIZ, 
     B                    VISSTRUCTU)
C     ******************************************************************
C     * ARGUMENT DESCRIPTION: SEE VISCROP                              *
C     ******************************************************************
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
      INTEGER       VISCROPSIZ, VISSTRUCTU(*)
      REAL          VISMATDATA(VISDIM1,VISDIM2,VISDIM3,VISDIM4,VISDIM5),
     A              VISPOSIRAD(VISDIM2),
     B              VISPOSIPOL(VISDIM3),
     C              VISPOSITOR(VISDIM4),
     D              VISPOSITIM(VISDIM5)
C
#include "comerrn"
#include "comrank"
C
C     * LOCAL VARIABLE
      INTEGER VISTMP
C
C     * COUNTERS
      INTEGER I
C
C     * EXTERNAL SUBROUTINES
      EXTERNAL VISPRERR, VISCROP
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
            CALL VISPRERR(VISERRWRST,'VISCROP7')
         ENDIF
   10 CONTINUE
C
      IF (VIS1TST) THEN
         VISTMP     = VISRESO(4)
         VISRESO(4) = 1
      ENDIF
C
      CALL VISCROP(VISMATDATA, 
     A             VISPOSIRAD, VISPOSIPOL, VISPOSITOR, VISPOSITIM,
     B             VISDIM1,       VISDIM2,    VISDIM3,    VISDIM4,
     C             VISDIM5,    VISCROPSIZ,    VISRESO,    VISVLEN)
C
      IF (VIS1TST) THEN
         VISRESO(4) = VISTMP
      ENDIF
C
      RETURN
      END
C
*DECK VISCROP7I
      SUBROUTINE VISCROP7I(VISMATDATA, VISPOSX,    VISPOSY,     VISPOSZ, 
     A                     VISPOSITIM, VISCROPSIZ, VISSTRUCTU)
C     ******************************************************************
C     * ARGUMENT DESCRIPTION: SEE VISCROPI                             *
C     ******************************************************************
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
      INTEGER       VISCROPSIZ, VISSTRUCTU(*)
      REAL          VISMATDATA(VISDIM1,VISDIM2,VISDIM3,VISDIM4,VISDIM5),
     A              VISPOSX(VISDIM2,VISDIM3,VISDIM4),
     B              VISPOSY(VISDIM2,VISDIM3,VISDIM4),
     C              VISPOSZ(VISDIM2,VISDIM3,VISDIM4),
     D              VISPOSITIM(VISDIM5)
C
#include "comerrn"
#include "comrank"
C
C     * LOCAL VARIABLE
      INTEGER VISTMP
C
C     * COUNTERS
      INTEGER I
C
C     * EXTERNAL SUBROUTINES
      EXTERNAL VISPRERR, VISCROPI
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
            CALL VISPRERR(VISERRWRST,'VISCROP7I')
         ENDIF
   10 CONTINUE
C
      IF (VIS1TST) THEN
         VISTMP     = VISRESO(4)
         VISRESO(4) = 1
      ENDIF
C
      CALL VISCROPI(VISMATDATA,
     A              VISPOSX,       VISPOSX,    VISPOSZ,    VISPOSITIM,
     B              VISDIM1,       VISDIM2,    VISDIM3,       VISDIM4,
     C              VISDIM5,    VISCROPSIZ,    VISRESO,       VISVLEN)
C
      IF (VIS1TST) THEN
         VISRESO(4) = VISTMP
      ENDIF
C
      RETURN
      END
C
*DECK VISINTERPOL7
      SUBROUTINE VISINTERPOL7(VISMATDATA, VISPOSIRAD, VISPOSIPOL, 
     A                        VISPOSITOR, VISPOSITIM, VISRESINIT, 
     B                        VISSTRUCTU)
C     ******************************************************************
C     * ARGUMENT DESCRIPTION: SEE VISINTERPOL                          *
C     ******************************************************************
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
      INTEGER       VISRESINIT(*), VISSTRUCTU(*)
      REAL          VISMATDATA(VISDIM1,VISDIM2,VISDIM3,VISDIM4,VISDIM5),
     A              VISPOSIRAD(VISDIM2),
     B              VISPOSIPOL(VISDIM3),
     C              VISPOSITOR(VISDIM4),
     D              VISPOSITIM(VISDIM5)
C
#include "comerrn"
#include "comrank"
C
C     * LOCAL VARIABLE
      INTEGER VISTMP
C
C     * COUNTERS
      INTEGER I
C
C     * EXTERNAL SUBROUTINES
      EXTERNAL VISPRERR, VISINTERPOL
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
            CALL VISPRERR(VISERRWRST,'VISINTERPOL7')
         ENDIF
   10 CONTINUE
C
      IF (VIS1TST) THEN
         VISTMP     = VISRESO(4)
         VISRESO(4) = 1
      ENDIF
C
      CALL VISINTERPOL(VISMATDATA, 
     A                 VISPOSIRAD, VISPOSIPOL, VISPOSITOR, VISPOSITIM,
     B                 VISDIM1,       VISDIM2,    VISDIM3,    VISDIM4,
     C                 VISDIM5,    VISRESINIT,    VISRESO,    VISVLEN)
C
      IF (VIS1TST) THEN
         VISRESO(4) = VISTMP
      ENDIF
C
      RETURN
      END
C
*DECK VISADDIGNO10
      SUBROUTINE VISADDIGNO10(VISMATDATA, VISPOSIRAD, VISPOSIPOL, 
     A                       VISPOSITOR, VISPOSITIM, VISRESINIT, 
     B                       VISSTRUCTU, MODES,      FINTIM, VISCOS)
C     ******************************************************************
C     * ARGUMENT DESCRIPTION: SEE VISADDIGNO                           *
C     ******************************************************************
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
      INTEGER       VISRESINIT(*), VISSTRUCTU(*)
      REAL          VISMATDATA(VISDIM1,VISDIM2,VISDIM3,VISDIM4,VISDIM5),
     A              VISPOSIRAD(VISDIM2),
     B              VISPOSIPOL(VISDIM3),
     C              VISPOSITOR(VISDIM4),
     D              VISPOSITIM(VISDIM5),
     E              MODES(*), FINTIM
      LOGICAL       VISCOS
C
#include "comerrn"
#include "comrank"
C
C     * LOCAL VARIABLE
      INTEGER VISTMP
C
C     * COUNTERS
      INTEGER I
C
C     * EXTERNAL SUBROUTINES
      EXTERNAL VISPRERR, VISTOCART
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
            CALL VISPRERR(VISERRWRST,'VISADDIGNO9')
         ENDIF
   10 CONTINUE
C
      IF (VIS1TST) THEN
         VISTMP     = VISRESO(4)
         VISRESO(4) = 1
      ENDIF
C
      CALL VISADDIGNO(VISMATDATA, 
     A               VISPOSIRAD, VISPOSIPOL, VISPOSITOR, 
     B               VISPOSITIM,
     C               VISDIM1,       VISDIM2,    VISDIM3,    VISDIM4,
     D               VISDIM5,    VISRESINIT,    VISRESO,    VISVLEN,
     E               MODES,          FINTIM,    VISCOS)
C
      IF (VIS1TST) THEN
         VISRESO(4) = VISTMP
      ENDIF
C
      RETURN
      END
C
*DECK VISADDIGNO10I
      SUBROUTINE VISADDIGNO10I(VISMATDATA,    VISPOSX,    VISPOSY, 
     A                         VISPOSZ,    VISPOSITIM, VISRESINIT, 
     B                         VISSTRUCTU, MODES, FINTIM, IGNO, VISCOS)
C     ******************************************************************
C     * ARGUMENT DESCRIPTION: SEE VISADDIGNOI                          *
C     ******************************************************************
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
      INTEGER       VISRESINIT(*), VISSTRUCTU(*)
      REAL          VISMATDATA(VISDIM1,VISDIM2,VISDIM3,VISDIM4,VISDIM5),
     A              VISPOSX(VISDIM2,VISDIM3,VISDIM4),
     B              VISPOSY(VISDIM2,VISDIM3,VISDIM4),
     C              VISPOSZ(VISDIM2,VISDIM3,VISDIM4),
     D              VISPOSITIM(VISDIM5),
     E              MODES(*), FINTIM
      LOGICAL       VISCOS
C
#include "comerrn"
#include "comrank"
C
C     * LOCAL VARIABLE
      INTEGER VISTMP
C
C     * COUNTERS
      INTEGER I
C
C     * EXTERNAL SUBROUTINES
      EXTERNAL VISPRERR, VISTOCART, IGNO
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
            CALL VISPRERR(VISERRWRST,'VISADDIGNO9I')
         ENDIF
   10 CONTINUE
C
      IF (VIS1TST) THEN
         VISTMP     = VISRESO(4)
         VISRESO(4) = 1
      ENDIF
C
      CALL VISADDIGNOI(VISMATDATA, 
     A                 VISPOSX,       VISPOSY,    VISPOSZ, 
     B                 VISPOSITIM,
     C                 VISDIM1,       VISDIM2,    VISDIM3,    VISDIM4,
     D                 VISDIM5,    VISRESINIT,    VISRESO,    VISVLEN,
     E                 MODES,          FINTIM,       IGNO,    VISCOS)
C
      IF (VIS1TST) THEN
         VISRESO(4) = VISTMP
      ENDIF
C
      RETURN
      END
C
*DECK VISTOCART9
      SUBROUTINE VISTOCART9(VISMATDATA,    VISPOSX,    VISPOSY,VISPOSZ, 
     A                      VISPOSIRAD, VISPOSIPOL, VISPOSITOR, 
     B                      VISSTRUCTU, INVASPECT)
C     ******************************************************************
C     * ARGUMENT DESCRIPTION: SEE VISTOCART                            *
C     ******************************************************************
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
      INTEGER       VISSTRUCTU(*)
      REAL          VISMATDATA(VISDIM1,VISDIM2,VISDIM3,VISDIM4,VISDIM5),
     A              VISPOSX(VISDIM2,VISDIM3,VISDIM4),
     B              VISPOSY(VISDIM2,VISDIM3,VISDIM4),
     C              VISPOSZ(VISDIM2,VISDIM3,VISDIM4),
     D              VISPOSIRAD(VISDIM2),
     E              VISPOSIPOL(VISDIM3),
     F              VISPOSITOR(VISDIM4),
     G              INVASPECT
C
#include "comerrn"
#include "comrank"
#include "comgeom"
C
C     * LOCAL VARIABLE
      INTEGER VISTMP
C
C     * COUNTERS
      INTEGER I
C
C     * EXTERNAL SUBROUTINES
      EXTERNAL VISPRERR, VISPRWAR, VISTOCART
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
            CALL VISPRERR(VISERRWRST,'VISTOCART10')
         ENDIF
   10 CONTINUE
C
      IF (     VISGEOM.NE.VISGEOMCYL.AND.VISGEOM.NE.VISGEOMTOK
     A    .AND.VISGEOM.NE.VISGEOMLOO.AND.VISGEOM.NE.VISGEOMAUX) THEN
          CALL VISPRWAR(VISWARGEOM,'VISTOCART10')
      ENDIF
C
      IF (VIS1TST) THEN
         VISTMP     = VISRESO(4)
         VISRESO(4) = 1
      ENDIF
C
      CALL VISTOCART(VISMATDATA, 
     A               VISPOSX,    VISPOSY,       VISPOSZ,
     A               VISPOSIRAD, VISPOSIPOL, VISPOSITOR, 
     B               VISDIM1,       VISDIM2,    VISDIM3,    VISDIM4,
     D               VISDIM5,       VISRESO,    VISGEOM,    VISVLEN,
     E               VISNFLD,    VISSTRUCTU,  INVASPECT)
C
      IF (VIS1TST) THEN
         VISRESO(4) = VISTMP
      ENDIF
C
      RETURN
      END
C
C     ==================================================================
C     =                   UTILITY SECTION BELOW.                       =
C     =                              |                                 =
C     =                              |                                 =
C     =                              |                                 =
C     =                            \   /                               =
C     =                             \ /                                =
C     =                              .                                 =
C     ==================================================================
*DECK VISCROP
      SUBROUTINE VISCROP(VISMATDATA, 
     A                   VISPOSIRAD, VISPOSIPOL, VISPOSITOR, VISPOSITIM,
     B                   VISNV,      VISNI,      VISNJ,      VISNK,      
     C                   VISNT,      VISCROPSIZ, VISRESINIT, VISVECTLEN)
C     ******************************************************************
C     * THIS ROUTINE REDUCES THE NUMBER OF ELEMENTS WITH A FACTOR      *
C     * VISCROPSIZ IN EACH DIMENSION. IT IS A PREREQUISIT THAT THE     *
C     * MAPPING FROM COMPUTATIONAL (I,J,K) SPACE TO PHYSICAL (X,Y,Z)   *
C     * SPACE IS WRITTEN AS A PRODUCT OF ARRAYS. I.E.                  *
C     * X(I,J,K)=R(I)*THETA(J)*PSI(K),                                 *
C     * Y(I,J,K)=R(I)*THETA(J)*PSI(K),                                 *
C     * Z(I,J,K)=R(I)*THETA(J)*PSI(K),                                 *
C     * WHERE WE HAVE TAKEN R, THETA, AND PSI AS AN EXAMPLE OF A       *
C     * CYLINDRICAL GEOMETRY. IN THIS LIBRARY WE HAVE ADOPTED A NAMING *
C     * CONVENTION RELATED TO THIS GEOMETRY. THE FIRST COMPUTATIONAL   *
C     * DIMENSION WILL ALWAYS BE DENOTED BY VISCOORRAD, THE SECOND BY  *
C     * VISCOORPOL (POLOIDAL), THE THIRD BY VISCOORTOR (TOROIDAL), AND *
C     * THE FOURTH BY VISCOORTIM (TIME).                               *
C     *                                                                *
C     * ON INPUT:                                                      *
C     * VISMATDATA CONTAINS THE UNCROPPED DATA.                        *
C     * VISNV      CONTAINS THE SIZE OF THE VECTLEN DIMENSION.         *
C     * VISNI      CONTAINS THE SIZE OF THE FIRST COMPUTATIONAL        *
C     *            DIMENSION.                                          *
C     * VISNJ      CONTAINS THE SIZE OF THE SECOND COMPUTATIONAL       *
C     *            DIMENSION.                                          *
C     * VISNK      CONTAINS THE SIZE OF THE THIRD COMPUTATIONAL        *
C     *            DIMENSION.                                          *
C     * VISNT      CONTAINS THE SIZE OF THE FOURTH COMPUTATIONAL       *
C     *            DIMENSION.                                          *
C     * VISCROPSIZ CONTAINS THE FACTOR WITH WHICH THE DATA SET IS      *
C     *            CROPPED (>=1).                                      *
C     * VISRESINIT CONTAINS THE INITIAL RESOLUTION OF THE DATA SET     *
C     *            (VISRESINIT(VISCOORRAD)<=VISNI),                    *
C     *            (VISRESINIT(VISCOORPOL)<=VISNJ),                    *
C     *            (VISRESINIT(VISCOORTOR)<=VISNK),                    *
C     *            (VISRESINIT(VISCOORTIM)<=VISNT).                    *
C     * VISVECTLEN CONTAINS THE NUMBER OF DATA ELEMENTS AT A NODE      *
C     *            (VISVECTLEN<=VISNV).                                *
C     * VISPOSIRAD CONTAINS THE ORIGINAL RADIAL   COORDINATES.         *
C     * VISPOSIPOL CONTAINS THE ORIGINAL POLOIDAL COORDINATES.         *
C     * VISPOSITOR CONTAINS THE ORIGINAL TOROIDAL COORDINATES.         *
C     * VISPOSITIM CONTAINS THE ORIGINAL TIME     COORDINATES.         *
C     *                                                                *
C     * ON OUTPUT:                                                     *
C     * VISRESINIT CONTAINS THE NEW RESOLUTION.                        *
C     * VISMATDATA CONTAINS THE CROPPED DATA SET.                      *
C     * VISPOSIRAD CONTAINS THE CROPPED RADIAL        COORDINATES.     *
C     * VISPOSIPOL CONTAINS THE CROPPED POLOIDAL      COORDINATES.     *
C     * VISPOSITOR CONTAINS THE CROPPED TOROIDAL      COORDINATES.     *
C     * VISPOSITIM CONTAINS THE CROPPED TIME          COORDINATES.     *
C     *                                                                *
C     * VISCROP IS CALLED BY THE USER PROGRAM.                         *
C     ******************************************************************
C
#include "comcoor"
#include "comerrn"
C
      INTEGER       VISCROPSIZ, VISRESINIT(*), VISVECTLEN, 
     A              VISNV,  VISNI,  VISNJ,  VISNK,  VISNT
      REAL          VISMATDATA(VISNV,VISNI,VISNJ,VISNK,VISNT),
     A              VISPOSIRAD(VISNI),
     B              VISPOSIPOL(VISNJ),
     C              VISPOSITOR(VISNK),
     D              VISPOSITIM(VISNT)
C
C     * COUNTERS
      INTEGER  I,  J,  K,  L,  M,
     A        II, JJ, KK, LL
C
C     * INITIAL VALUES
      DATA JJ/0/
C
C     * CHECKING DIMENSION VARIABLES
      CALL VISDIMVAR(VISNV, VISNI, VISNJ, VISNK, VISNT, VISVECTLEN, 
     A               VISRESINIT, VISRESINIT, 'VISCROP')
C
    1 IF (VISCROPSIZ.LT.1) THEN
         CALL VISPRERR(VISWARCROP,'VISCROP')
         VISCROPSIZ=1
      ELSE
         II=0
         DO 10 I=1, VISRESINIT(VISCOORTIM), VISCROPSIZ
            II=II+1
            JJ=0
            KK=0
            LL=0
            DO 20 J=1, VISRESINIT(VISCOORTOR), VISCROPSIZ
               JJ=JJ+1
               KK=0
               LL=0
               DO 30 K=1, VISRESINIT(VISCOORPOL), VISCROPSIZ
                  KK=KK+1
                  LL=0
                  DO 40 L=1, VISRESINIT(VISCOORRAD), VISCROPSIZ
                     LL=LL+1
                     DO 50 M=1, VISVECTLEN
                        VISMATDATA(M,LL,KK,JJ,II)=VISMATDATA(M,L,K,J,I)
   50                CONTINUE
   40             CONTINUE
   30          CONTINUE
   20       CONTINUE
   10    CONTINUE
C
         KK=0
         DO 60 K=1, VISRESINIT(VISCOORRAD), VISCROPSIZ
            KK=KK+1
            VISPOSIRAD(KK)=VISPOSIRAD(K)  
   60    CONTINUE
         VISRESINIT(VISCOORRAD)=KK
C
         JJ=0
         DO 70 J=1, VISRESINIT(VISCOORPOL), VISCROPSIZ
            JJ=JJ+1
            VISPOSIPOL(JJ)=VISPOSIPOL(J)  
   70    CONTINUE
         VISRESINIT(VISCOORPOL)=JJ
C
         II=0
         DO 80 I=1, VISRESINIT(VISCOORTOR), VISCROPSIZ
            II=II+1
            VISPOSITOR(II)=VISPOSITOR(I)  
   80    CONTINUE
         VISRESINIT(VISCOORTOR)=II
C
         LL=0
         DO 120 L=1, VISRESINIT(VISCOORTIM), VISCROPSIZ
            LL=LL+1
            VISPOSITIM(LL)=VISPOSITIM(L)  
  120    CONTINUE
         VISRESINIT(VISCOORTIM)=LL
      ENDIF
      RETURN
C
      END
C
*DECK VISCROPI
      SUBROUTINE VISCROPI(VISMATDATA, 
     A                   VISPOSX, VISPOSY, VISPOSZ, VISPOSITIM,
     B                   VISNV,      VISNI,      VISNJ,      VISNK,      
     C                   VISNT,      VISCROPSIZ, VISRESINIT, VISVECTLEN)
C     ******************************************************************
C     * THIS ROUTINE REDUCES THE NUMBER OF ELEMENTS WITH A FACTOR      *
C     * VISCROPSIZ IN EACH DIMENSION. IT IS A PREREQUISIT THAT THE     *
C     * THE MAPPING FROM COMPUTATIONAL (I,J,K) SPACE TO PHYSICAL       *
C     * (X,Y,Z) IS EXPLICITLY GIVEN BY THREE MATRICES VISPOSX(I,J,K),  *
C     * VISPOSY(I,J,K), AND VISPOSZ(I,J,K).                            * 
C     *                                                                *
C     * ON INPUT:                                                      *
C     * VISMATDATA CONTAINS THE UNCROPPED DATA.                        *
C     * VISNV      CONTAINS THE SIZE OF THE VECTLEN DIMENSION.         *
C     * VISNI      CONTAINS THE SIZE OF THE FIRST COMPUTATIONAL        *
C     *            DIMENSION.                                          *
C     * VISNJ      CONTAINS THE SIZE OF THE SECOND COMPUTATIONAL       *
C     *            DIMENSION.                                          *
C     * VISNK      CONTAINS THE SIZE OF THE THIRD COMPUTATIONAL        *
C     *            DIMENSION.                                          *
C     * VISNT      CONTAINS THE SIZE OF THE FOURTH COMPUTATIONAL       *
C     *            DIMENSION.                                          *
C     * VISCROPSIZ CONTAINS THE FACTOR WITH WHICH THE DATA SET IS      *
C     *            CROPPED (>=1).                                      *
C     * VISRESINIT CONTAINS THE INITIAL RESOLUTION OF THE DATA SET     *
C     *            (VISRESINIT(VISCOORRAD)<=VISNI),                    *
C     *            (VISRESINIT(VISCOORPOL)<=VISNJ),                    *
C     *            (VISRESINIT(VISCOORTOR)<=VISNK),                    *
C     *            (VISRESINIT(VISCOORTIM)<=VISNT).                    *
C     * VISVECTLEN CONTAINS THE NUMBER OF DATA ELEMENTS AT A NODE      *
C     *            (VISVECTLEN<=VISNV).                                *
C     * VISPOSX CONTAINS THE ORIGINAL CARTESIAN X COORDINATES.         *
C     * VISPOSY CONTAINS THE ORIGINAL CARTESIAN Y COORDINATES.         *
C     * VISPOSZ CONTAINS THE ORIGINAL CARTESIAN Z COORDINATES.         *
C     * VISPOSITIM CONTAINS THE ORIGINAL TIME     COORDINATES.         *
C     *                                                                *
C     * ON OUTPUT:                                                     *
C     * VISRESINIT CONTAINS THE NEW RESOLUTION.                        *
C     * VISMATDATA CONTAINS THE CROPPED DATA SET.                      *
C     * VISPOSX    CONTAINS THE CROPPED CARTESIAN X COORDINATES.       *
C     * VISPOSY    CONTAINS THE CROPPED CARTESIAN Z COORDINATES.       *
C     * VISPOSZ    CONTAINS THE CROPPED CARTESIAN Z COORDINATES.       *
C     * VISPOSITIM CONTAINS THE CROPPED TIME        COORDINATES.       *
C     *                                                                *
C     * VISCROPI IS CALLED BY THE USER PROGRAM.                        *
C     ******************************************************************
C
#include "comcoor"
#include "comerrn"
C
      INTEGER       VISCROPSIZ, VISRESINIT(*), VISVECTLEN, 
     A              VISNV,  VISNI,  VISNJ,  VISNK,  VISNT
      REAL          VISMATDATA(VISNV,VISNI,VISNJ,VISNK,VISNT),
     D              VISPOSITIM(VISNT),
     E              VISPOSX(VISNI,VISNJ,VISNK),
     F              VISPOSY(VISNI,VISNJ,VISNK),
     G              VISPOSZ(VISNI,VISNJ,VISNK)
C
C     * COUNTERS
      INTEGER  I,  J,  K,  L,  M,
     A        II, JJ, KK, LL
C
C     * INITIAL VALUES
      DATA JJ/0/
C
C     * CHECKING DIMENSION VARIABLES
      CALL VISDIMVAR(VISNV, VISNI, VISNJ, VISNK, VISNT, VISVECTLEN, 
     A               VISRESINIT, VISRESINIT, 'VISCROP')
C
    1 IF (VISCROPSIZ.LT.1) THEN
         CALL VISPRWAR(VISWARCROP,'VISCROP')
         VISCROPSIZ=1
      ELSE
         II=0
         DO 10 I=1, VISRESINIT(VISCOORTIM), VISCROPSIZ
            II=II+1
            JJ=0
            KK=0
            LL=0
            DO 20 J=1, VISRESINIT(VISCOORTOR), VISCROPSIZ
               JJ=JJ+1
               KK=0
               LL=0
               DO 30 K=1, VISRESINIT(VISCOORPOL), VISCROPSIZ
                  KK=KK+1
                  LL=0
                  DO 40 L=1, VISRESINIT(VISCOORRAD), VISCROPSIZ
                     LL=LL+1
                     DO 50 M=1, VISVECTLEN
                        VISMATDATA(M,LL,KK,JJ,II)=VISMATDATA(M,L,K,J,I)
   50                CONTINUE
   40             CONTINUE
   30          CONTINUE
   20       CONTINUE
   10    CONTINUE
C
         KK=0
         DO 90 K=1,VISRESINIT(VISCOORTOR),VISCROPSIZ       
            KK=KK+1
            JJ=0
            II=0
            DO 100 J=1,VISRESINIT(VISCOORPOL),VISCROPSIZ
               JJ=JJ+1
               II=0
               DO 110 I=1,VISRESINIT(VISCOORRAD),VISCROPSIZ
                  II=II+1
                  VISPOSX(II,JJ,KK)=VISPOSX(I,J,K)
                  VISPOSY(II,JJ,KK)=VISPOSY(I,J,K)
                  VISPOSZ(II,JJ,KK)=VISPOSZ(I,J,K)
  110          CONTINUE
  100       CONTINUE
   90    CONTINUE
         VISRESINIT(VISCOORRAD)=II
         VISRESINIT(VISCOORPOL)=JJ
         VISRESINIT(VISCOORTOR)=KK
C
         LL=0
         DO 120 L=1, VISRESINIT(VISCOORTIM), VISCROPSIZ
            LL=LL+1
            VISPOSITIM(LL)=VISPOSITIM(L)  
  120    CONTINUE
         VISRESINIT(VISCOORTIM)=LL
      ENDIF
      RETURN
C
      END
C
*DECK VISINTERPOL
      SUBROUTINE VISINTERPOL(VISMATDATA, 
     A                       VISPOSIRAD, VISPOSIPOL, VISPOSITOR, 
     B                       VISPOSITIM,
     C                       VISNV,      VISNI,      VISNJ,      VISNK,      
     D                       VISNT,      VISRESINIT, VISRESFINA, 
     E                       VISVECTLEN)
C     ******************************************************************
C     * THIS ROUTINE LINEARLY INTERPOLATES THE DATA SET TO GET THE     * 
C     * FINAL RESOLUTION. INTERPOLATION IS DONE FOR EACH DIMENSION.    * 
C     *                                                                * 
C     * ON INPUT:                                                      *
C     * VISMATDATA CONTAINS THE ORIGINAL DATA.                         *
C     * VISNV      CONTAINS THE SIZE OF THE VECTLEN DIMENSION.         *
C     * VISNI      CONTAINS THE SIZE OF THE FIRST COMPUTATIONAL        *
C     *            DIMENSION.                                          *
C     * VISNJ      CONTAINS THE SIZE OF THE SECOND COMPUTATIONAL       *
C     *            DIMENSION.                                          *
C     * VISNK      CONTAINS THE SIZE OF THE THIRD COMPUTATIONAL        *
C     *            DIMENSION.                                          *
C     * VISNT      CONTAINS THE SIZE OF THE FOURTH COMPUTATIONAL       *
C     *            DIMENSION.                                          *
C     * VISRESINIT CONTAINS THE INITIAL RESOLUTION OF THE DATA SET     *
C     *            (VISRESINIT(VISCOORRAD)<=VISNI),                    *
C     *            (VISRESINIT(VISCOORPOL)<=VISNJ),                    *
C     *            (VISRESINIT(VISCOORTOR)<=VISNK),                    *
C     *            (VISRESINIT(VISCOORTIM)<=VISNT).                    *
C     * VISRESFINA CONTAINS THE FINAL   RESOLUTION OF THE DATA SET     *
C     *            (VISRESFINA(VISCOORRAD)<=VISNI),                    *
C     *            (VISRESFINA(VISCOORPOL)<=VISNJ),                    *
C     *            (VISRESFINA(VISCOORTOR)<=VISNK),                    *
C     *            (VISRESFINA(VISCOORTIM)<=VISNT).                    *
C     * VISVECTLEN CONTAINS THE NUMBER OF DATA ELEMENTS AT A NODE.     *
C     *            (VISVECTLEN<=VISNV).                                *
C     * VISPOSIRAD CONTAINS THE ORIGINAL RADIAL   COORDINATES.         *
C     * VISPOSIPOL CONTAINS THE ORIGINAL POLOIDAL COORDINATES.         *
C     * VISPOSITOR CONTAINS THE ORIGINAL TOROIDAL COORDINATES.         *
C     * VISPOSITIM CONTAINS THE ORIGINAL TIME     COORDINATES.         *
C     *                                                                *
C     * ON OUTPUT:                                                     *
C     * VISRESINIT CONTAINS THE NEW RESOLUTION.                        *
C     * VISMATDATA CONTAINS THE INTERPOLATED DATA SET.                 *
C     * VISPOSITIM CONTAINS THE INTERPOLATED TIME     COORDINATES.     *
C     *                                                                *
C     * VISPOSIRAD CONTAINS THE INTERPOLATED RADIAL   COORDINATES.     *
C     * VISPOSIPOL CONTAINS THE INTERPOLATED POLOIDAL COORDINATES.     *
C     * VISPOSITOR CONTAINS THE INTERPOLATED TOROIDAL COORDINATES.     *
C     *                                                                *
C     * VISINTERPOL ARE CALLED BY THE USER PROGRAM.                    *
C     ******************************************************************
C
#include "comcoor"
#include "comerrn"
C
      INTEGER VISNV,     VISNI,      VISNJ,      VISNK,      VISNT,
     A        VISRESINIT(*), VISRESFINA(*),
     B        VISVECTLEN 
      REAL    VISMATDATA(VISNV,VISNI,VISNJ,VISNK,VISNT),
     A        VISPOSIRAD(*),
     B        VISPOSIPOL(VISNJ),
     C        VISPOSITOR(VISNK),
     D        VISPOSITIM(VISNT)
C
C     * COUNTERS
      INTEGER I, J, K, L, M,
     A       II,JJ,KK,LL
C
C     * INTERPOLATION (TO AVOID ANOTHER WORKING ARRAY R IS DECLARED
C                      TO HAVE 1000 ELEMENTS WHICH SHOULD SUFFICE.)
      REAL    RANGE,
     A        VU, VL,
     B        RU, RL, R(1000)
      INTEGER KMAX, LMAX, MMAX
      DATA KMAX/1/
     A     LMAX/1/
     B     MMAX/1/
C
C     * SUPPRESSING UNDERFLOW ON IBM
C
C     * CHECKING DIMENSION VARIABLES
  7   CALL VISDIMVAR(VISNV, VISNI, VISNJ, VISNK, VISNT, VISVECTLEN,
     A               VISRESINIT, VISRESFINA, 'VISINTERPOL')
C
      DO 10 I=VISCOORRAD, VISCOORTIM
C
         IF (    VISRESINIT(I).EQ.VISRESFINA(I)
     A       .OR.VISRESINIT(I).EQ.1) GO TO 10
C
         RANGE=0.0
         IF (I.EQ.VISCOORRAD) THEN
            RANGE=VISPOSIRAD(VISRESINIT(VISCOORRAD))-VISPOSIRAD(1)
         ELSE IF (I.EQ.VISCOORPOL) THEN 
            RANGE=VISPOSIPOL(VISRESINIT(VISCOORPOL))-VISPOSIPOL(1) 
         ELSE IF (I.EQ.VISCOORTOR) THEN
            RANGE=VISPOSITOR(VISRESINIT(VISCOORTOR))-VISPOSITOR(1) 
         ELSE IF (I.EQ.VISCOORTIM) THEN
            RANGE=VISPOSITIM(VISRESINIT(VISCOORTIM))-VISPOSITIM(1) 
         ENDIF
C
         IF (I.EQ.VISCOORRAD) THEN
            R(1)=VISPOSIRAD(1)
            RL  =VISPOSIRAD(1)
            RU  =VISPOSIRAD(2)
         ELSE IF (I.EQ.VISCOORPOL) THEN
            R(1)=VISPOSIPOL(1)
            RL  =VISPOSIPOL(1)
            RU  =VISPOSIPOL(2)
         ELSE IF (I.EQ.VISCOORTOR) THEN
            R(1)=VISPOSITOR(1)
            RL  =VISPOSITOR(1)
            RU  =VISPOSITOR(2)
         ELSE IF (I.EQ.VISCOORTIM) THEN
            R(1)=VISPOSITIM(1)
            RL  =VISPOSITIM(1)
            RU  =VISPOSITIM(2)
         ENDIF
C
         DO 20 J=2, VISRESFINA(I)
            IF (I.EQ.VISCOORRAD) THEN
               IF (VISRESFINA(VISCOORRAD).NE.1) THEN
                  R(J)=VISPOSIRAD(1)+REAL(J-1)
     A              /REAL(VISRESFINA(VISCOORRAD)-1)*RANGE
               ELSE
                   GOTO 10
               ENDIF
            ELSE IF (I.EQ.VISCOORPOL) THEN 
               IF (VISRESFINA(VISCOORPOL).NE.1) THEN
                  R(J)=VISPOSIPOL(1)+REAL(J-1)
     A              /REAL(VISRESFINA(VISCOORPOL)-1)*RANGE
               ELSE
                   GOTO 10
               ENDIF
            ELSE IF (I.EQ.VISCOORTOR) THEN
               IF (VISRESFINA(VISCOORTOR).NE.1) THEN
                  R(J)=VISPOSITOR(1)+REAL(J-1)
     A              /REAL(VISRESFINA(VISCOORTOR)-1)*RANGE
               ELSE
                   GOTO 10
               ENDIF
            ELSE IF (I.EQ.VISCOORTIM) THEN
               IF (VISRESFINA(VISCOORTIM).NE.1) THEN
                  R(J)=VISPOSITIM(1)+REAL(J-1)
     A              /REAL(VISRESFINA(VISCOORTIM)-1)*RANGE
               ELSE
                   GOTO 10
               ENDIF
            ENDIF
C
            L=1
            DO 30 K = 1, VISRESINIT(I)-1
               IF (I.EQ.VISCOORRAD) THEN
                  IF (VISPOSIRAD(K).LT.R(J).AND.VISPOSIRAD(K+1).GE.R(J)) 
     A            THEN
                     RL = VISPOSIRAD(K)
                     RU = VISPOSIRAD(K+1)
                     L  = K
                  ENDIF
               ELSE IF(I.EQ.VISCOORPOL) THEN
                  IF (VISPOSIPOL(K).LT.R(J).AND.VISPOSIPOL(K+1).GE.R(J)) 
     A            THEN
                     RL = VISPOSIPOL(K)
                     RU = VISPOSIPOL(K+1)
                     L  = K
                  ENDIF
               ELSE IF(I.EQ.VISCOORTOR) THEN
                  IF (VISPOSITOR(K).LT.R(J).AND.VISPOSITOR(K+1).GE.R(J)) 
     A            THEN
                     RL = VISPOSITOR(K)
                     RU = VISPOSITOR(K+1)
                     L  = K
                  ENDIF
               ELSE IF(I.EQ.VISCOORTIM) THEN
                  IF (VISPOSITIM(K).LT.R(J).AND.VISPOSITIM(K+1).GE.R(J)) 
     A            THEN
                     RL = VISPOSITIM(K)
                     RU = VISPOSITIM(K+1)
                     L  = K
                  ENDIF
               ENDIF
   30       CONTINUE
C
            DO 40 II=1, VISVECTLEN
               IF (I.EQ.VISCOORRAD) THEN
                  DO 50 JJ=1,VISRESINIT(VISCOORPOL)
                     DO 60 KK=1,VISRESINIT(VISCOORTOR)
                        DO 70 LL=1,VISRESINIT(VISCOORTIM)
                           VU=VISMATDATA(II,L+1,JJ,KK,LL)
                           VL=VISMATDATA(II,L,JJ,KK,LL)
                           VISMATDATA(II,J,JJ,KK,LL)=VL+(VU-VL)
     A                                              *(R(J)-RL)/(RU-RL)
   70                   CONTINUE
   60                CONTINUE
   50             CONTINUE
               ELSEIF (I.EQ.VISCOORPOL) THEN
                  DO 80 JJ=1,VISRESFINA(VISCOORRAD)
                     DO 90 KK=1,VISRESINIT(VISCOORTOR)
                        DO 100 LL=1,VISRESINIT(VISCOORTIM)
                           VU=VISMATDATA(II,JJ,L+1,KK,LL)
                           VL=VISMATDATA(II,JJ,L,KK,LL)
                           VISMATDATA(II,JJ,J,KK,LL)=VL+(VU-VL)
     A                                              *(R(J)-RL)/(RU-RL)
  100                   CONTINUE
   90                CONTINUE
   80             CONTINUE
               ELSEIF (I.EQ.VISCOORTOR) THEN
                  DO 110 JJ=1,VISRESFINA(VISCOORRAD)
                     DO 120 KK=1,VISRESFINA(VISCOORPOL)
                        DO 130 LL=1,VISRESINIT(VISCOORTIM)
                           VU=VISMATDATA(II,JJ,KK,L+1,LL)
                           VL=VISMATDATA(II,JJ,KK,L,LL)
                           VISMATDATA(II,JJ,KK,J,LL)=VL+(VU-VL)
     A                                              *(R(J)-RL)/(RU-RL)
  130                   CONTINUE
  120                CONTINUE
  110             CONTINUE
               ELSEIF (I.EQ.VISCOORTIM) THEN
                  DO 140 JJ=1,VISRESFINA(VISCOORRAD)
                     DO 150 KK=1,VISRESFINA(VISCOORPOL)
                        DO 160 LL=1,VISRESFINA(VISCOORTOR)
                           VU=VISMATDATA(II,JJ,KK,LL,L+1)
                           VL=VISMATDATA(II,JJ,KK,LL,L)
                           VISMATDATA(II,JJ,KK,LL,J)=VL+(VU-VL)
     A                                              *(R(J)-RL)/(RU-RL)
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
               ENDIF
   40       CONTINUE
C
   20    CONTINUE
         DO 170 M=1, VISRESFINA(I)
            IF (I.EQ.VISCOORRAD) THEN
               VISPOSIRAD(M)=R(M)
            ELSE IF (I.EQ.VISCOORPOL) THEN
               VISPOSIPOL(M)=R(M)
            ELSE IF (I.EQ.VISCOORTOR) THEN
               VISPOSITOR(M)=R(M)
            ELSE IF (I.EQ.VISCOORTIM) THEN
               VISPOSITIM(M)=R(M)
            ENDIF
  170    CONTINUE
   10 CONTINUE
C
  270 DO 280 I=VISCOORRAD, VISCOORTIM
         VISRESINIT(I)=VISRESFINA(I)
  280 CONTINUE
      RETURN
C
      END
C
*DECK VISADDIGNO
      SUBROUTINE VISADDIGNO(VISMATDATA,
     A                      VISPOSIRAD, VISPOSIPOL, VISPOSITOR,
     B                      VISPOSITIM,
     C                      VISNV,      VISNI,      VISNJ,      VISNK,      
     D                      VISNT,
     E                      VISRESINIT, VISRESFINA, VISVECTLEN, MODES,
     F                      FINTIM, VISCOS)
C     ******************************************************************
C     * THIS ROUTINE ADDS THE DEPENDENCE OF AN IGNORABLE COORDINATE TO *
C     * THE PHYSICAL QUANTITIES. A COORDINATE IS IGNORABLE IF ITS      *
C     * RESOLUTION EQUALS 1. THE DEPENDENCE OF THE PHYSICAL QUANTITIES *
C     * ON THE IGNORABLE COORDINATES IS ASSUMED TO BE COSINUS-LIKE.    *
C     *                                                                *
C     * ON INPUT:                                                      * 
C     * VISMATDATA CONTAINS THE ORIGINAL DATA.                         * 
C     * VISNV      CONTAINS THE SIZE OF THE VECTLEN DIMENSION.         *
C     * VISNI      CONTAINS THE SIZE OF THE FIRST COMPUTATIONAL        *
C     *            DIMENSION.                                          *
C     * VISNJ      CONTAINS THE SIZE OF THE SECOND COMPUTATIONAL       *
C     *            DIMENSION.                                          *
C     * VISNK      CONTAINS THE SIZE OF THE THIRD COMPUTATIONAL        *
C     *            DIMENSION.                                          *
C     * VISNT      CONTAINS THE SIZE OF THE FOURTH COMPUTATIONAL       *
C     *            DIMENSION.                                          *
C     * VISRESINIT CONTAINS THE INITIAL RESOLUTION OF THE DATA SET     *
C     *            (VISRESINIT(VISCOORRAD)<=VISNI),                    *
C     *            (VISRESINIT(VISCOORPOL)<=VISNJ),                    *
C     *            (VISRESINIT(VISCOORTOR)<=VISNK),                    *
C     *            (VISRESINIT(VISCOORTIM)<=VISNT).                    *
C     * VISRESFINA CONTAINS THE FINAL   RESOLUTION OF THE DATA SET     *
C     *            (VISRESFINA(VISCOORRAD)<=VISNI),                    *
C     *            (VISRESFINA(VISCOORPOL)<=VISNJ),                    *
C     *            (VISRESFINA(VISCOORTOR)<=VISNK),                    *
C     *            (VISRESFINA(VISCOORTIM)<=VISNT).                    *
C     *            IF VISRESINIT(I)=1 AND VISRESFINA(I)>1,             *
C     *               COORDINATE I IS CONSIDERED TO BE IGNORABLE       *
C     *               AND ITS DEPENDENCE WILL BE ADDED.                *
C     * VISVECTLEN CONTAINS THE NUMBER OF DATA ELEMENTS AT A NODE      *
C     *            (VISVECTLEN<=VISNV).                                *
C     * VISPOSIRAD CONTAINS THE ORIGINAL RADIAL   COORDINATES,         *
C     * VISPOSIPOL CONTAINS THE ORIGINAL POLOIDAL COORDINATES,         *
C     * VISPOSITOR CONTAINS THE ORIGINAL TOROIDAL COORDINATES.         *
C     * VISPOSITIM CONTAINS THE ORIGINAL TIME     COORDINATES.         *
C     * MODES      CONTAINS THE FOURIER MODE NUMBERS REPRESENTING THE  *
C     *            DEPENDENCE OF THE CORRESPONDING IGNORABLE           *
C     *            COORDINATE.                                         *
C     * FINTIM     CONTAINS THE TIME RANGE NEEDED IN CASE T IS THE     *
C     *            IGNORABLE COORDINATE.                               *
C     * VISCOS     IF .TRUE. IGNORABLE DEPENDENCE IS COSINUS-LIKE ELSE * 
C     *            SINUS-LIKE                                          *
C     *                                                                *
C     * ON OUTPUT:                                                     *
C     * VISPOSIRAD CONTAINS THE NEW RADIAL   COORDINATES.              *
C     * VISPOSIPOL CONTAINS THE NEW POLOIDAL COORDINATES.              *
C     * VISPOSITOR CONTAINS THE NEW TOROIDAL COORDINATES.              *
C     * VISPOSITIM CONTAINS THE NEW TIME COORDINATES.                  *
C     *                                                                *
C     * VISADDIGNO CALLS VISIGNOCOS.                                   *       
C     * VISADDIGNO IS CALLED BY THE USER PROGRAM.                      *
C     ******************************************************************
C
#include "comgeom"
#include "comcoor"
#include "comerrn"
#include "comcons"
#include "cominou"
C
      REAL    MODES(*)
      INTEGER VISNV,      VISNI,      VISNJ,      VISNK,     VISNT,
     A        VISRESINIT(*), VISRESFINA(*), VISVECTLEN
      REAL    VISMATDATA(VISNV,VISNI,VISNJ,VISNK,VISNT),
     D        VISPOSIRAD(VISNI),
     E        VISPOSIPOL(VISNJ),
     F        VISPOSITOR(VISNK),
     G        VISPOSITIM(VISNT),
     H        FINTIM
      LOGICAL VISCOS
C
C     *COUNTERS
      INTEGER  I, J, K, L, M,
     A        II,JJ,KK,LL
C
C     * LOCAL VARIABLES
      REAL    IGNORDEP, VISIGNOCOS
      INTEGER IGNINDICATION, IGNORABLE, NIGNO, IGNOCOOR(4)
C
C     * CHECKING DIMENSION VARIABLES
      CALL VISDIMVAR(VISNV, VISNI, VISNJ, VISNK, VISNT, VISVECTLEN,
     A               VISRESINIT, VISRESFINA, 'VISADDIGNO')
C
    1 IGNORDEP     =1.0
      IGNORABLE    =1
      IGNINDICATION=0
C
      DO 10 I=VISCOORRAD, VISCOORTIM
         IF (VISRESINIT(I).EQ.IGNORABLE
     A       .AND.VISRESFINA(I).GT.IGNORABLE) THEN
            IF (I.EQ.VISCOORRAD) THEN
               IGNINDICATION=IGNINDICATION +1
            ELSE IF (I.EQ.VISCOORPOL) THEN
               IGNINDICATION=IGNINDICATION +10
            ELSE IF (I.EQ.VISCOORTOR) THEN
               IGNINDICATION=IGNINDICATION +100
            ELSE IF (I.EQ.VISCOORTIM) THEN
               IGNINDICATION=IGNINDICATION +1000
            ENDIF
         ENDIF
   10 CONTINUE
C
      NIGNO=0
      DO 20 I=VISCOORRAD,VISCOORTIM
         IF (VISRESINIT(I).EQ.IGNORABLE) THEN
            NIGNO          =NIGNO+1
            IGNOCOOR(NIGNO)=I
            DO 30 J=1,VISRESFINA(I)
               IF (I.EQ.VISCOORRAD.AND.VISRESFINA(I).NE.1) THEN
                  VISPOSIRAD(J)=REAL(J-1)/REAL(VISRESFINA(I)-1)
               ELSE IF (I.EQ.VISCOORPOL.AND.VISRESFINA(I).NE.1) THEN
                  VISPOSIPOL(J)=2*PI*REAL(J-1)/REAL(VISRESFINA(I)-1)
               ELSE IF (I.EQ.VISCOORTOR.AND.VISRESFINA(I).NE.1) THEN
                  VISPOSITOR(J)=2*PI*REAL(J-1)/REAL(VISRESFINA(I)-1)
               ELSE IF (I.EQ.VISCOORTIM.AND.VISRESFINA(I).NE.1) THEN
                  VISPOSITIM(J)=FINTIM*REAL(J-1)/REAL(VISRESFINA(I)-1)
               ENDIF
   30       CONTINUE
         ENDIF
   20 CONTINUE
C
      DO 40 I=1,VISRESFINA(VISCOORTIM)
         IF (VISRESINIT(VISCOORTIM).EQ.IGNORABLE) THEN
            II=IGNORABLE
         ELSE
            II=I
         ENDIF
         DO 50 J=1,VISRESFINA(VISCOORTOR)
            IF (VISRESINIT(VISCOORTOR).EQ.IGNORABLE) THEN
               JJ=IGNORABLE
            ELSE
               JJ=J
            ENDIF
            DO 60 K=1,VISRESFINA(VISCOORPOL)
               IF (VISRESINIT(VISCOORPOL).EQ.IGNORABLE) THEN
                  KK=IGNORABLE
               ELSE
                  KK=K
               ENDIF
               DO 70 L=1,VISRESFINA(VISCOORRAD)
                  IF (VISRESINIT(VISCOORRAD).EQ.IGNORABLE) THEN
                     LL=IGNORABLE
                  ELSE
                     LL=L
                  ENDIF
                  IF (VISCOS) THEN
                     IGNORDEP=COS(VISIGNOCOS(VISPOSIRAD,VISPOSIPOL,
     A                                       VISPOSITOR,VISPOSITIM,
     B                                       VISNI,VISNJ,VISNK,VISNT,
     C                                       IGNINDICATION,MODES,
     D                                       L,K,J,I))
                  ELSE
                     IGNORDEP=SIN(VISIGNOCOS(VISPOSIRAD,VISPOSIPOL,
     A                                       VISPOSITOR,VISPOSITIM,
     B                                       VISNI,VISNJ,VISNK,VISNT,
     C                                       IGNINDICATION,MODES,
     D                                       L,K,J,I))
                  ENDIF
                  DO 80 M=1,VISVECTLEN
                     VISMATDATA(M,L,K,J,I)=VISMATDATA(M,LL,KK,JJ,II)
     A                                     *IGNORDEP
   80             CONTINUE
   70          CONTINUE
   60       CONTINUE
   50    CONTINUE
   40 CONTINUE
C
      RETURN
      END
C
*DECK VISADDIGNOI
      SUBROUTINE VISADDIGNOI(VISMATDATA,
     A                      VISPOSX, VISPOSY, VISPOSZ,
     B                      VISPOSITIM,
     C                      VISNV,      VISNI,      VISNJ,      VISNK,      
     D                      VISNT,
     E                      VISRESINIT, VISRESFINA, VISVECTLEN, MODES,
     F                      FINTIM,     IGNO, VISCOS)
C     ******************************************************************
C     * THIS ROUTINE ADDS THE DEPENDENCE OF AN IGNORABLE COORDINATE TO *
C     * THE PHYSICAL QUANTITIES. A COORDINATE IS IGNORABLE IF ITS      *
C     * RESOLUTION EQUALS 1. THE DEPENDENCE OF THE PHYSICAL QUANTITIES *
C     * ON THE IGNORABLE COORDINATES IS ASSUMED TO BE FOURIER-LIKE.    *
C     *                                                                *
C     * THE VARIABLE IGNO IS PROVIDED BY THE                           *
C     * USER PROGRAM AS THE NAME OF A ROUTINE WHICH ADDS THE           *
C     * DEPENDENCE ON IGNORABLE COORDINATES TO THE CARTESIAN ONES.     *
C     * THIS IS NECESSARY AS THE DEPENDENCE BETWEEN THE CARTESIAN      *
C     * COORDINATES AND THE COMPUTATIONAL COORDINATES IS NOT KNOWN A   *
C     * PRIORI BY THE LIBRARY.                                         *
C     *                                                                *
C     * ON INPUT:                                                      * 
C     * VISMATDATA CONTAINS THE ORIGINAL DATA.                         * 
C     * VISNV      CONTAINS THE SIZE OF THE VECTLEN DIMENSION.         *
C     * VISNI      CONTAINS THE SIZE OF THE FIRST COMPUTATIONAL        *
C     *            DIMENSION.                                          *
C     * VISNJ      CONTAINS THE SIZE OF THE SECOND COMPUTATIONAL       *
C     *            DIMENSION.                                          *
C     * VISNK      CONTAINS THE SIZE OF THE THIRD COMPUTATIONAL        *
C     *            DIMENSION.                                          *
C     * VISNT      CONTAINS THE SIZE OF THE FOURTH COMPUTATIONAL       *
C     *            DIMENSION.                                          *
C     * VISRESINIT CONTAINS THE INITIAL RESOLUTION OF THE DATA SET     *
C     *            (VISRESINIT(VISCOORRAD)<=VISNI),                    *
C     *            (VISRESINIT(VISCOORPOL)<=VISNJ),                    *
C     *            (VISRESINIT(VISCOORTOR)<=VISNK),                    *
C     *            (VISRESINIT(VISCOORTIM)<=VISNT).                    *
C     * VISRESFINA CONTAINS THE FINAL   RESOLUTION OF THE DATA SET     *
C     *            (VISRESFINA(VISCOORRAD)<=VISNI),                    *
C     *            (VISRESFINA(VISCOORPOL)<=VISNJ),                    *
C     *            (VISRESFINA(VISCOORTOR)<=VISNK),                    *
C     *            (VISRESFINA(VISCOORTIM)<=VISNT).                    *
C     *            IF VISRESINIT(I)=1 AND VISRESFINA(I)>1,             *
C     *               COORDINATE I IS CONSIDERED TO BE IGNORABLE       *
C     *               AND ITS DEPENDENCE WILL BE ADDED.                *
C     * VISVECTLEN CONTAINS THE NUMBER OF DATA ELEMENTS AT A NODE      *
C     *            (VISVECTLEN<=VISNV).                                *
C     *                                                                *
C     * VISPOSX    CONTAINS THE ORIGINAL CARTESIAN X COORDINATES,      *
C     * VISPOSY    CONTAINS THE ORIGINAL CARTESIAN Y COORDINATES,      *
C     * VISPOSZ    CONTAINS THE ORIGINAL CARTESIAN Z COORDINATES.      *
C     * VISPOSITIM CONTAINS THE ORIGINAL TIME     COORDINATES.         *
C     * MODES      CONTAINS THE FOURIER MODE NUMBERS REPRESENTING THE  *
C     *            DEPENDENCE OF THE CORRESPONDING IGNORABLE           *
C     *            COORDINATE.                                         *
C     * FINTIM     CONTAINS THE TIME RANGE NEEDED IN CASE T IS THE     *
C     *            IGNORABLE COORDINATE.                               *
C     * IGNO       IS THE USER-SUPPLIED ROUTINE THAT IS USED TO PRODUCE*
C     *            TO CORRECT DEPENDENCE OF THE CARTESIAN X,Y,Z        *
C     *            COORDINATES ON THE IGNORABLE ONES.                  *
C     * VISCOS     IF .TRUE. IGNORABLE DEPENDENCE IS COSINUS-LIKE ELSE * 
C     *            SINUS-LIKE                                          *
C     *                                                                *
C     * ON OUTPUT:                                                     *
C     * VISPOSX    CONTAINS THE CARTESIAN X COORDINATES.               *
C     * VISPOSY    CONTAINS THE CARTESIAN Z COORDINATES.               *
C     * VISPOSZ    CONTAINS THE CARTESIAN Z COORDINATES.               *
C     * VISPOSITIM CONTAINS THE NEW TIME COORDINATES.                  *
C     *                                                                *
C     * VISADDIGNOI CALLS VISIGNOCOSI AND IGNI.                        *       
C     * VISADDIGNOI ARE CALLED BY THE USER PROGRAM.                    *
C     ******************************************************************
C
#include "comgeom"
#include "comcoor"
#include "comerrn"
#include "comcons"
#include "cominou"
C
      REAL    MODES(*)
      INTEGER VISNV,      VISNI,      VISNJ,      VISNK,     VISNT,
     A        VISRESINIT(*), VISRESFINA(*), VISVECTLEN
      REAL    VISMATDATA(VISNV,VISNI,VISNJ,VISNK,VISNT),
     A        VISPOSX(VISNI,VISNJ,VISNK),
     B        VISPOSY(VISNI,VISNJ,VISNK),
     C        VISPOSZ(VISNI,VISNJ,VISNK),
     G        VISPOSITIM(VISNT),
     H        FINTIM
      LOGICAL VISCOS
C
C     *COUNTERS
      INTEGER  I, J, K, L, M,
     A        II,JJ,KK,LL
C
C     * LOCAL VARIABLES
      REAL    IGNORDEP, VISIGNOCOSI
      INTEGER IGNINDICATION, IGNORABLE, NIGNO, IGNOCOOR(4)
C
C     * CHECKING DIMENSION VARIABLES
      CALL VISDIMVAR(VISNV, VISNI, VISNJ, VISNK, VISNT, VISVECTLEN,
     A               VISRESINIT, VISRESFINA, 'VISADDIGNO')
C
    1 IGNORDEP     =1.0
      IGNORABLE    =1
      IGNINDICATION=0
C
      DO 10 I=VISCOORRAD, VISCOORTIM
         IF (VISRESINIT(I).EQ.IGNORABLE
     A       .AND.VISRESFINA(I).GT.IGNORABLE) THEN
            IF (I.EQ.VISCOORRAD) THEN
               IGNINDICATION=IGNINDICATION +1
            ELSE IF (I.EQ.VISCOORPOL) THEN
               IGNINDICATION=IGNINDICATION +10
            ELSE IF (I.EQ.VISCOORTOR) THEN
               IGNINDICATION=IGNINDICATION +100
            ELSE IF (I.EQ.VISCOORTIM) THEN
               IGNINDICATION=IGNINDICATION +1000
            ENDIF
         ENDIF
   10 CONTINUE
C
      NIGNO=0
      DO 20 I=VISCOORRAD,VISCOORTIM
         IF (VISRESINIT(I).EQ.IGNORABLE) THEN
            NIGNO          =NIGNO+1
            IGNOCOOR(NIGNO)=I
            DO 30 J=1,VISRESFINA(I)
               IF (I.EQ.VISCOORTIM.AND.VISRESFINA(I).NE.1) THEN
                  VISPOSITIM(J)=FINTIM*REAL(J-1)/REAL(VISRESFINA(I)-1)
               ENDIF
   30       CONTINUE
         ENDIF
   20 CONTINUE
C
C     ******************************************************************
C     * IN CASE VISCOORTYP=VISIRRCOOR THE COORDINATES SUPPLIED ARE     *
C     * CARTESIAN. HOWEVER, THEY DO NOT CORRESPOND, IN GENERAL, TO THE *
C     * COORDINATES USED IN THE ORIGINAL PROBLEM. THE CARTESIAN        *
C     * COORDINATES DO DEPEND ON THESE ONES, E.G. X(I,J,K), Y(I,J,K),  *
C     * Z(I,J,K), WHERE I,J, AND K INDICATE THE ORIGINAL COORDINATES.  *
C     * IF ONE OF THE THREE I,J, AND, K IS AN IGNORABLE COORDINATE,    *
C     * THE USER-SUPPLIED ROUTINE 'IGNO' IS USED TO ADD THE IGNORABLE  *
C     * COORDINATE DEPENDENCE TO THE CARTESIAN COORDINATES.            *
C     * FOR EXAMPLE, IF X(I,J,1), Y(I,J,1), AND Z(I,J,1) (HENCE K IS   *
C     * IGNORABLE) THE ROUTINE 'IGNO' IS USED TO GET X(I,J,K), Y(I,J,K)*
C     * AND Z(I,J,K).                                                  *
C     ******************************************************************
C
         CALL IGNO(VISPOSX,    VISPOSY,    VISPOSZ, VISNI, VISNJ, VISNK,
     A             VISRESINIT, VISRESFINA)

C
      DO 40 I=1,VISRESFINA(VISCOORTIM)
         IF (VISRESINIT(VISCOORTIM).EQ.IGNORABLE) THEN
            II=IGNORABLE
         ELSE
            II=I
         ENDIF
         DO 50 J=1,VISRESFINA(VISCOORTOR)
            IF (VISRESINIT(VISCOORTOR).EQ.IGNORABLE) THEN
               JJ=IGNORABLE
            ELSE
               JJ=J
            ENDIF
            DO 60 K=1,VISRESFINA(VISCOORPOL)
               IF (VISRESINIT(VISCOORPOL).EQ.IGNORABLE) THEN
                  KK=IGNORABLE
               ELSE
                  KK=K
               ENDIF
               DO 70 L=1,VISRESFINA(VISCOORRAD)
                  IF (VISRESINIT(VISCOORRAD).EQ.IGNORABLE) THEN
                     LL=IGNORABLE
                  ELSE
                     LL=L
                  ENDIF
                  IF (VISCOS) THEN
                     IGNORDEP=COS(VISIGNOCOSI(VISPOSX,VISPOSY,VISPOSZ,
     A                                        VISPOSITIM,
     B                                        VISNI,VISNJ,VISNK,VISNT,
     C                                        IGNINDICATION,MODES,
     D                                        L,K,J,I))
                  ELSE
                     IGNORDEP=SIN(VISIGNOCOSI(VISPOSX,VISPOSY,VISPOSZ,
     A                                        VISPOSITIM,
     B                                        VISNI,VISNJ,VISNK,VISNT,
     C                                        IGNINDICATION,MODES,
     D                                        L,K,J,I))
                  ENDIF
                  DO 80 M=1,VISVECTLEN
                     VISMATDATA(M,L,K,J,I)=VISMATDATA(M,LL,KK,JJ,II)
     A                                     *IGNORDEP
   80             CONTINUE
   70          CONTINUE
   60       CONTINUE
   50    CONTINUE
   40 CONTINUE
      RETURN
C
      END
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
#include "comgeom"
#include "comerrn"
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
C
*DECK VISTOCART 
      SUBROUTINE VISTOCART(VISMATDATA,
     A                     VISPOSX,    VISPOSY,    VISPOSZ,  
     B                     VISPOSIRAD, VISPOSIPOL, VISPOSITOR,
     C                     VISNV,      VISNI,      VISNJ,      VISNK,      
     D                     VISNT,      VISRES,     VISGEOM,    
     E                     VISVECTLEN, VISNFIELDS, VISSTRUCTU, 
     F                     INVASPECT)
C     ******************************************************************
C     * THIS ROUTINE CONVERTS COORDINATES ADAPTED TO SPECIAL GEOMETRIES*
C     * TO CARTESIAN COORDINATES WHICH ARE USED BY THE VISUALIZATION   *
C     * PACKAGES. AUXILIAR GEOMETRIES SHOULD ALWAYS BE DEFINED WITH    *
C     * CARTESIAN COORDINATES. HENCE, THIS ROUTINE DOES NOT HAVE TO BE *
C     * USED FOR THESE GEOMETRIES.                                     *
C     *                                                                *
C     * ON INPUT:                                                      *
C     * VISMATDATA CONTAINS THE ORIGINAL  DATA.                        *
C     * VISNV      CONTAINS THE SIZE OF THE VECTLEN DIMENSION.         *
C     * VISNI      CONTAINS THE SIZE OF THE FIRST COMPUTATIONAL        *
C     *            DIMENSION.                                          *
C     * VISNJ      CONTAINS THE SIZE OF THE SECOND COMPUTATIONAL       *
C     *            DIMENSION.                                          *
C     * VISNK      CONTAINS THE SIZE OF THE THIRD COMPUTATIONAL        *
C     *            DIMENSION.                                          *
C     * VISNT      CONTAINS THE SIZE OF THE FOURTH COMPUTATIONAL       *
C     *            DIMENSION.                                          *
C     * VISGEOM    DETERMINES THE TYPE OF GEOMETRY:                    *
C     *            VISGEOM=VISGEOMCYL -> CYLINDRICAL GEOMETRY,         *
C     *            VISGEOM=VISGEOMTOK -> TOKAMAK GEOMETRY              *
C     *                                  (CIRCULAR CONCENTRIC),        *
C     *            VISGEOM=VISGEOMLOO -> CORONAL LOOP GEOMETRY,        *
C     *            VISGEOM=VISGEOMAUX -> AUXILIAR GEOMETRY.            *
C     * VISNFIELDS CONTAINS THE NUMBER OF FIELDS PER NODE.             *
C     * VISSTRUCTU DETERMINES THE TYPE OF THE FIELDS.                  *
C     *            VISSTRUCTU()=VISRANKSCA -> SCALAR FIEL.             *
C     *            VISSTRUCTU()=VISRANKVC1 -> 1-VECTOR FIELD.          *
C     *            VISSTRUCTU()=VISRANKVC2 -> 2-VECTOR FIELD.          *
C     *            VISSTRUCTU()=VISRANKVC3 -> 3-VECTOR FIELD.          *
C     *            VISSTRUCTU()=VISRANKVEC -> 3-VECTOR FIELD.          *
C     * VISPOSIRAD CONTAINS THE RADIAL   COORDINATES.                  *
C     * VISPOSIPOL CONTAINS THE POLOIDAL COORDINATES.                  *
C     * VISPOSITOR CONTAINS THE TOROIDAL COORDINATES.                  *
C     * INVASPECT  CONTAINS THE INVERSE ASPECT RATIO WHICH IS NEEDED   *
C     *            FOR THE THREE SPECIAL GEOMETRIES.                   *
C     *                                                                *
C     * ON OUTPUT:                                                     *
C     * VISPOSX CONTAINS THE X-COORDINATES.                            *
C     * VISPOSY CONTAINS THE Y-COORDINATES.                            *
C     * VISPOSZ CONTAINS THE Z-COORDINATES.                            *
C     *                                                                *
C     * VISTOCART IS CALLED BY THE USER PROGRAM.                       *
C     ******************************************************************
C
#include "comgeom"
#include "comcoor"
#include "comcons"
#include "comrank"
#include "comerrn"
C 
      INTEGER VISRES(*), VISSTRUCTU(*), VISGEOM, VISNFIELDS, VISVECTLEN,
     A        VISNV,     VISNI,         VISNJ,   VISNK,      VISNT
      REAL    VISMATDATA(VISNV,VISNI,VISNJ,VISNK,VISNT),
     A        VISPOSX(VISNI,VISNJ,VISNK),
     B        VISPOSY(VISNI,VISNJ,VISNK),
     C        VISPOSZ(VISNI,VISNJ,VISNK),
     D        VISPOSIRAD(VISNI),
     E        VISPOSIPOL(VISNJ),
     F        VISPOSITOR(VISNK),
     G        INVASPECT
C
C     * CYLINDRICAL COORDINATES
      REAL    R,  THETA,  ZCYL
C     * CYLINDRICAL VECTOR COMPONENTS
      REAL   VR, VTHETA, VZCYL
C
C     * TOKAMAK AND LOOP COORDINATES
      REAL    PHI
C     * TOKAMAK AND LOOP COORDINATES
      REAL   VPHI
C
C     * CARTESIAN VECTOR COMPONENTS
      REAL   VX, VY, VZ
C
C     * COUNTERS
      INTEGER I,J,K,L,M
C
C     * TEMP VARIABELES
      INTEGER  OFFSET
C
C     * CHECKING DIMENSION VARIABLES
      CALL VISDIMVAR(VISNV, VISNI, VISNJ, VISNK, VISNT, VISVECTLEN,
     A               VISRES, VISRES, 'VISTOCART')
C
C     * THE CYLINDER
      IF (VISGEOM.EQ.VISGEOMCYL) THEN
         DO 10 K=1,VISRES(VISCOORTOR)
            ZCYL=VISPOSITOR(K)/INVASPECT
            DO 20 J=1,VISRES(VISCOORPOL)
               THETA=VISPOSIPOL(J)
               DO 30 I=1,VISRES(VISCOORRAD)
                  R=VISPOSIRAD(I)
                  VISPOSX(I,J,K)=R*COS(THETA)
                  VISPOSY(I,J,K)=R*SIN(THETA)
                  VISPOSZ(I,J,K)=ZCYL
C
                  OFFSET=1
                  DO 40 M=1,VISNFIELDS
                     IF (VISSTRUCTU(M).EQ.VISRANKVEC  .OR.
     A                   VISSTRUCTU(M).EQ.VISRANKVC3) THEN
C                       * 3-VECTORS DO HAVE TO BE TRANSFORMED
                        DO 50 L=1,VISRES(VISCOORTIM)
                           VR    =VISMATDATA(OFFSET,I,J,K,L)
                           VTHETA=VISMATDATA(OFFSET+1,I,J,K,L)
                           VZCYL =VISMATDATA(OFFSET+2,I,J,K,L)
                           VX    =VR*COS(THETA)-VTHETA*SIN(THETA)
                           VY    =VR*SIN(THETA)+VTHETA*COS(THETA)
                           VZ    =VZCYL
                           VISMATDATA(OFFSET  ,I,J,K,L)=VX
                           VISMATDATA(OFFSET+1,I,J,K,L)=VY
                           VISMATDATA(OFFSET+2,I,J,K,L)=VZ
   50                   CONTINUE
                        OFFSET=OFFSET+VISRANKVC3
                     ELSE 
C                       * SCALARS AND 1- AND 2-VECTORS DO NOT HAVE TO 
C                       * BE TRANSFORMED
                        OFFSET=OFFSET+VISSTRUCTU(M)
                     ENDIF
   40             CONTINUE
   30          CONTINUE
   20       CONTINUE
   10    CONTINUE         
      ENDIF
C
C     * THE CIRCULAR CONCENTRIC TOKAMAK
      IF (VISGEOM.EQ.VISGEOMTOK) THEN
         DO 60 K=1,VISRES(VISCOORTOR)
            PHI=VISPOSITOR(K)
            DO 70 J=1,VISRES(VISCOORPOL)
               THETA=VISPOSIPOL(J)
               DO 80 I=1,VISRES(VISCOORRAD)
                  R=VISPOSIRAD(I)
                  VISPOSX(I,J,K)=(1.0/INVASPECT+R*COS(THETA))*COS(PHI)
                  VISPOSY(I,J,K)=(1.0/INVASPECT+R*COS(THETA))*SIN(PHI)
                  VISPOSZ(I,J,K)=R*SIN(THETA)
C
                  OFFSET=1
                  DO 90 M=1,VISNFIELDS
                     IF (VISSTRUCTU(M).EQ.VISRANKVEC  .OR.
     A                   VISSTRUCTU(M).EQ.VISRANKVC3) THEN
C                       * VECTORS DO HAVE TO BE TRANSFORMED
                        DO 100 L=1,VISRES(VISCOORTIM)
                           VR    =VISMATDATA(OFFSET,I,J,K,L)
                           VTHETA=VISMATDATA(OFFSET+1,I,J,K,L)
                           VPHI  =VISMATDATA(OFFSET+2,I,J,K,L)
                           VX    =(VR*COS(THETA)-VTHETA*SIN(THETA))
     A                            *COS(PHI)
                           VY    =(VR*SIN(THETA)+VTHETA*COS(THETA))
     A                            *SIN(PHI)
                           VZ    =VR*SIN(THETA)+VTHETA*COS(THETA)
                           VISMATDATA(OFFSET  ,I,J,K,L)=VX
                           VISMATDATA(OFFSET+1,I,J,K,L)=VY
                           VISMATDATA(OFFSET+2,I,J,K,L)=VZ
  100                   CONTINUE
                        OFFSET=OFFSET+VISRANKVC3
                     ELSE
C                       * SCALARS AND 1- AND 2-VECTORS DO NOT HAVE TO 
C                       * BE TRANSFORMED
                        OFFSET=OFFSET+VISSTRUCTU(M)
                     ENDIF
   90             CONTINUE
   80          CONTINUE
   70       CONTINUE
   60    CONTINUE         
      ENDIF
C
C     * THE CORONAL LOOP
      IF (VISGEOM.EQ.VISGEOMLOO) THEN
         DO 110 K=1,VISRES(VISCOORTOR)
            PHI=VISPOSITOR(K)/VISPOSITOR(VISRES(VISCOORTOR))
     A          *4.0*ATAN(1.0)
            DO 120 J=1,VISRES(VISCOORPOL)
               THETA=VISPOSIPOL(J)
               DO 130 I=1,VISRES(VISCOORRAD)
                  R=VISPOSIRAD(I)
                  VISPOSX(I,J,K)=(1./INVASPECT+R*SIN(THETA))*COS(PHI)
                  VISPOSY(I,J,K)=(1./INVASPECT+R*SIN(THETA))*SIN(PHI)
                  VISPOSZ(I,J,K)=R*COS(THETA)
C
                  OFFSET=1
                  DO 140 M=1,VISNFIELDS
                     IF (VISSTRUCTU(M).EQ.VISRANKVEC  .OR.
     A                   VISSTRUCTU(M).EQ.VISRANKVC3) THEN
C                       * VECTORS DO HAVE TO BE TRANSFORMED
                        DO 150 L=1,VISRES(VISCOORTIM)
                           VR    =VISMATDATA(OFFSET,I,J,K,L)
                           VTHETA=VISMATDATA(OFFSET+1,I,J,K,L)
                           VPHI  =VISMATDATA(OFFSET+2,I,J,K,L)
                           VX=(VR*SIN(THETA)+VTHETA*COS(THETA))*SIN(PHI)
     A                        -VPHI*SIN(PHI)
                           VY=(VR*SIN(THETA)+VTHETA*COS(THETA))*COS(PHI)
     A                        +VPHI*COS(PHI)
                           VZ=VR*COS(THETA)-VTHETA*SIN(THETA)
                           VISMATDATA(OFFSET  ,I,J,K,L)=VX
                           VISMATDATA(OFFSET+1,I,J,K,L)=VY
                           VISMATDATA(OFFSET+2,I,J,K,L)=VZ
  150                   CONTINUE
                        OFFSET=OFFSET+VISRANKVC3
                     ELSE
C                       * SCALARS AND 1- AND 2-VECTORS DO NOT HAVE TO 
C                       * BE TRANSFORMED
                        OFFSET=OFFSET+VISSTRUCTU(M)
                     ENDIF
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE         
      ENDIF
C
C     * AUXILIARY GEOMETRY
      IF (VISGEOM.EQ.VISGEOMAUX) THEN
         CALL VISPRERR(VISERRAUXI,'VISTOCART')
      ENDIF
      END
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
#include "comcoor"
#include "comerrn"
#include "comrank"
#include "cominou"
#include "comform"
#include "comtype"
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
#include "comcoor"
#include "comerrn"
#include "cominou"
#include "comcate"
#include "comrank"
#include "comform"
#include "comtype"
#include "comorde"
#include "comcarr"
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
C
*DECK VISAVSHEAD
      SUBROUTINE VISAVSHEAD(VISFILNAME, VISRES,     VISNDIM,    VISCOOR,
     A                      VISVECTLEN, VISNFIELDS, VISFIENAME, 
     B                      VISSTRUCTU, VISSTRUTYP, VISFORMATS,  VISXDR,
     C                      VISHOW)
C     ******************************************************************
C     * THIS ROUTINE TOGETHER WITH VISOPENAVS (GENERAL FORMAT ONLY)    *
C     * AND VISDATAOUT                                                 *
C     * OUTPUT THE DATA IN SUCH A FORMAT THAT IT CAN BE IMMEDIATELY    *
C     * READ INTO AVS USING THE 'READ FIELD' MODULE. THE GENERAL FILE  *
C     * FORMAT AS WELL AS THE INTERNAL AVS FIELD FORMAT ARE SUPPORTED. *
C     * HOWEVER WE DO RECOMMEND TO USE THE INTERNAL FORMAT ONLY SINCE  *
C     * IT IS MUCH MORE EFFICIENT WITH RESPECT TO I/O PERFORMANCE.     *
C     * THE GENERAL FILE FORMAT IS ONLY SUITED FOR PORTABILITY REASONS.*
C     * IT ALLOWS ONE TO OUTPUT EVERYTHING IN ASCII AND IT SEPERATES   *
C     * THE HEADER FROM THE DATA PART, I.E. DIFFERENT FILES ARE USED.  *
C     * THE OUTPUT IN THE DATA FILE IS COLUMN BASED AND CAN BE READ    *
C     * BY OTHER APPLICATIONS EASILY.                                  *
C     *                                                                *
C     * ON INPUT:                                                      *
C     * VISFILNAME CONTAINS THE FILENAME OF THE OUTPUT FILE. IF IT     *
C     *            DOES NOT ALLREADY CONTAIN THE .fld POSTFIX IT WILL  *
C     *            BE ADDED.                                           *
C     * VISRES     CONTAINS THE RESOLUTION OF THE DATA SET             *
C     *            (VISRESINIT(VISCOORRAD)<=VISNI),                    *
C     *            (VISRESINIT(VISCOORPOL)<=VISNJ),                    *
C     *            (VISRESINIT(VISCOORTOR)<=VISNK),                    *
C     *            (VISRESINIT(VISCOORTIM)<=VISNT).                    *
C     * VISNDIM    CONTAINS THE NUMBER OF DIMENSION TO BE PRINTED.     *
C     * VISCOOR    CONTAINS THE DIMENSIONS WHICH WILL BE PRINTED.      *
C     * VISVECTLEN CONTAINS THE NUMBER OF DATA ELEMENTS PER NODE.      *
C     * VISNFIELDS CONTAINS THE NUMBER OF FIELDS PER NODE.             *
C     * VISFIENAME CONTAINS THE NAMES OF THE DIFFERENT FIELDS.         *
C     * VISHOW     DETERMINES THE FILE FORMAT                          *
C     *            VISHOW=VISGENERAL -> GENERAL  FILE FORMAT,          *
C     *            VISHOW=VISINTERNA -> INTERNAL FILE FORMAT.          *
C     * VISSTRUCTU DETERMINES THE TYPE OF THE FIELDS                   *
C     *            VISSTRUCTU()=VISRANKSCA -> SCALAR FIELD,            *
C     *            VISSTRUCTU()=VISRANKVC1 -> 1-VECTOR FIELD.          *
C     *            VISSTRUCTU()=VISRANKVC2 -> 2-VECTOR FIELD.          *
C     *            VISSTRUCTU()=VISRANKVC3 -> 3-VECTOR FIELD.          *
C     *            VISSTRUCTU()=VISRANKVEC -> 3-VECTOR FIELD.          *
C     * VISSTRUTYP CONTAINS THE PRIMITIVE TYPE OF THE FIELD            *
C     *            VISTYPEDOU -> DOUBLE                                *
C     *            VISTYPEFLO -> FLOAT                                 *
C     * VISFORMATS DETERMINES WHETHER OR NOT THE DATA WILL BE PRINTED  *
C     *            IN ASCII OR IN BINARY FORMAT.                       *
C     *            VISFORMATS=VISFORMASC -> ASCII                      *
C     *            VISFORMATS=VISFORMBIN -> BINARY                     *
C     * VISXDR     DETERMINES IF BINARY IS IN XDR FORMAT OR NOT        *
C     *            (AVS ONLY)                                          *
C     ******************************************************************
C
#include "comcoor"
#include "comerrn"
#include "cominou"
#include "comrank"
#include "comtype"
#include "comform"
C
      LOGICAL       VISXDR
      INTEGER       VISRES(*),  VISNDIM,    VISCOOR(*),     
     A              VISVECTLEN, VISNFIELDS, VISSTRUCTU(*),  VISFORMATS, 
     B              VISHOW
      CHARACTER*(*) VISFILNAME, VISFIENAME(*),  VISSTRUTYP(*)
C
C     *COUNTERS
      INTEGER I,J
C
C     * LOCAL VARIABLES
      INTEGER NSPACE, SKIP, OFFSET, STRIDE,
     A        IND, IND2
      CHARACTER VISCOORTYP*(*), HEADERNAME*80, TYPE*9,
     A          FILETYPE*11
      PARAMETER (
     1    NSPACE=3,
     2    VISCOORTYP='irregular'
     3)
C
C     * INITIAL VALUES
      DATA SKIP/0/
C
C     * FUNCTIONS
      INTEGER       VISSIZEOF
      LOGICAL       VISNAMEOK
C
C     * CHECKING DIMENSIONS
      IF (VISNDIM.LT.1.OR.VISNDIM.GT.4) THEN
         CALL VISPRERR(VISERRNDIM,'VISAVSHEAD')
      ENDIF
C     
C     * FILE NAME OK?
      IF (.NOT.VISNAMEOK(VISFILNAME)) THEN
         CALL VISPRWAR(VISWARNAME,'VISAVSHEAD')
         VISFILNAME='avsfile'
      ENDIF
C
C     * OPEN HEADER FILE FOR GENERAL FILE FORMAT
C     * FOR AVS INTERNAL FIELD FORMAT THIS HEADER FILE ALSO CONTAINS THE
C     * BINARY DATA PART
      HEADERNAME = VISFILNAME
      IND = INDEX(HEADERNAME,'.fld')
      IF (IND.EQ.0) THEN
         IND = MIN(76,INDEX(HEADERNAME,' '))
         IF (IND.EQ.0) IND = MAX(1,LEN(HEADERNAME)-3)
      ENDIF
      IND = MIN(76,IND)
      HEADERNAME(IND:IND+3)   = '.fld'
      HEADERNAME(IND+4:IND+4) = CHAR(0)
      OPEN(UNIT=VISAVSOUTP,FILE=HEADERNAME(1:IND+3),STATUS='UNKNOWN')
C
C     * AVS INDICATION
      WRITE(VISAVSOUTP,100)
C
C     * NUMBER OF DIMENSIONS
      WRITE(VISAVSOUTP,110) VISNDIM
C   
C     * DIMENSION SIZES
      WRITE(VISAVSOUTP,120) (I, VISRES(VISCOOR(I)), I=1,VISNDIM)
C
C     * DIMENSION OF PHYSICAL SPACE
      WRITE(VISAVSOUTP,130) NSPACE
C
C     * NUMBER OF DATA VALUES PER DATA NODE
      WRITE(VISAVSOUTP,140) VISVECTLEN
C
C     * SETTING THE POSITIONS AS THE LAST FIELD
      VISSTRUCTU(VISNFIELDS+1)=VISRANKVC3
      VISSTRUTYP(VISNFIELDS+1)=VISSTRUTYP(1)
      VISFIENAME(VISNFIELDS+1)='locations'
C
C     * PRIMITIVE DATA TYPE (ALL STRUCTURES MUST HAVE THE SAME TYPE IN 
C       AVS)
      TYPE=VISSTRUTYP(1)
      DO 10 I=2,VISNFIELDS+1
         IF (VISSTRUTYP(I).NE.TYPE) THEN
            CALL VISPRERR(VISERRSTRU,'VISAVSHEAD')
         ENDIF
   10 CONTINUE
      IF (VISXDR) THEN
         WRITE(VISAVSOUTP,150) 'xdr_'//TYPE
      ELSE
         WRITE(VISAVSOUTP,150) TYPE
      ENDIF
C
C     * COORDINATE TYPE (IN THIS VERSION ONLY IRREGULAR IS SUPPORTED)
      WRITE(VISAVSOUTP,160) VISCOORTYP
C
C     * STRUCTURE NAME LABELS
      DO 20 I=1,VISNFIELDS
         DO 30 J=1,VISSTRUCTU(I)
            IND2=INDEX(VISFIENAME(I),' ')
            IF (IND2.EQ.0) IND2 = LEN(VISFIENAME(I))+1
            IF (VISSTRUCTU(I).EQ.VISRANKSCA) THEN
               WRITE(VISAVSOUTP,170) VISFIENAME(I)(1:IND2-1)
            ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC1) THEN
               WRITE(VISAVSOUTP,170) VISFIENAME(I)(1:IND2-1),J
            ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC2) THEN
               WRITE(VISAVSOUTP,170) VISFIENAME(I)(1:IND2-1),J
            ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC3) THEN
               WRITE(VISAVSOUTP,170) VISFIENAME(I)(1:IND2-1),J
            ELSE IF (VISSTRUCTU(I).EQ.VISRANKVEC) THEN
               WRITE(VISAVSOUTP,170) VISFIENAME(I)(1:IND2-1),J
            ENDIF
   30    CONTINUE
   20 CONTINUE
C
      IF (.NOT.(VISFORMATS.EQ.VISFORMASC.OR.
     A    VISFORMATS.EQ.VISFORMBIN.OR.
     B    VISFORMATS.EQ.VISFORMUNF)) THEN
          CALL VISPRWAR(VISWARFORM,'VISAVSHEAD')
          VISFORMATS=VISFORMASC
      ENDIF
C
      IF (VISFORMATS.EQ.VISFORMBIN.AND.VISHOW.EQ.VISINTERNA) THEN
C        * ONLY FOR INTERNAL AVS FIELD
         FILETYPE='binary'
      ENDIF
C
C     * DATA VARIABLES
      DO 40 I=1,VISVECTLEN
         IF (VISFORMATS.EQ.VISFORMASC) THEN
            FILETYPE='ascii'
            SKIP    =VISNFIELDS+1
            OFFSET=I-1
            STRIDE=VISNFIELDS+NSPACE
            WRITE(VISAVSOUTP,190) I, VISFILNAME(1:IND-1), FILETYPE(1:5), 
     A                            SKIP, OFFSET, STRIDE
         ELSE IF (VISFORMATS.EQ.VISFORMUNF) THEN
            FILETYPE='unformatted'
            SKIP  =4+(I-1)*VISSIZEOF(TYPE)
            STRIDE=NSPACE*VISSIZEOF(VISSTRUTYP(VISNFIELDS+1))
     A            +VISVECTLEN*VISSIZEOF(TYPE)
            WRITE(VISAVSOUTP,195) I, VISFILNAME(1:IND-1),FILETYPE(1:11),
     A                            SKIP, STRIDE
         ENDIF
   40 CONTINUE
C
C     * COORDINATES
      DO 50 I=1,NSPACE
         IF (VISFORMATS.EQ.VISFORMASC) THEN
            OFFSET=VISNFIELDS+(I-1)
            STRIDE=VISNFIELDS+NSPACE
            WRITE(VISAVSOUTP,180) I, VISFILNAME(1:IND-1), FILETYPE(1:5),
     A                            SKIP, OFFSET, STRIDE
         ELSE IF (VISFORMATS.EQ.VISFORMUNF) THEN
            SKIP  =4+(VISVECTLEN+I-1)*VISSIZEOF(TYPE)
            STRIDE=NSPACE*VISSIZEOF(VISSTRUTYP(VISNFIELDS+1))
     A            +VISVECTLEN*VISSIZEOF(VISSTRUTYP(1))
            WRITE(VISAVSOUTP,185) I, VISFILNAME(1:IND-1),FILETYPE(1:11),
     A                            SKIP, STRIDE
         ENDIF
   50 CONTINUE
C
C     * CLOSE FILE 
      CLOSE(UNIT=VISAVSOUTP)
C
C     * ADDING TWO FORM FEEDS IN CASE OF BINARY FORMAT
      IF (VISFORMATS.EQ.VISFORMBIN) THEN
         CALL writeff(HEADERNAME)
      ENDIF
C
C     * FORMATS
  100 FORMAT('# AVS field file',/,'#')
  110 FORMAT('ndim=',I1) 
  120 FORMAT('dim',I1,'=',I4) 
  130 FORMAT('nspace=',I1) 
  140 FORMAT('veclen=',I1) 
  150 FORMAT('data=',A) 
  160 FORMAT('field=',A) 
  170 FORMAT('label=',A,I1)
  180 FORMAT('coord',1X,I1,1X,'file=',A,1X,'filetype=',A,1X,
     A       'skip=',1X,I7,1X,'offset=',I7,1X,'stride=',I7)
  185 FORMAT('coord',1X,I1,1X,'file=',A,1X,'filetype=',A,1X,
     A       'skip=',1X,I7,1X,'stride=',I7)
  190 FORMAT('variable',1X,I1,1X,'file=',A,1X,'filetype=',A,1X,
     A       'skip=',1X,I7,1X,'offset=',I7,1X,'stride=',I7)
  195 FORMAT('variable',1X,I1,1X,'file=',A,1X,'filetype=',A,1X,
     A       'skip=',1X,I7,1X,'stride=',I7)
      END
C
*DECK VISOPENDX
      FUNCTION VISOPENDX(VISFILNAME, VISFORMATS)
C     ******************************************************************
C     * THIS FUNCTION IS USED TO OPEN THE DATA FILE FOR DX'S GENERAL   *
C     * FILE FORMAT. IT IS CALLED BY THE USER AFTER THE HEADER HAS     *
C     * BEEN WRITTEN. THE FUNCTIONS RETURNS THE UNIT NUMBER ON WHICH   *
C     * THE DATA FILE HAS BEEN OPENED FOR WRITING.                     *
C     *                                                                *
C     * ON INPUT:                                                      *
C     * VISFLNAME CONTAINS THE NAME OF THE DATA FILE TO BE OPENED.     *
C     * VISFORMATS CONTAINS THE FORM IN WHICH DATA WILL BE WRITTEN     *
C     *    VISFORMATS=VISFORMASC -> ASCII,                             *
C     *    VISFORMATS=VISFORMBIN -> BINARY,                            *
C     *    VISFORMATS=VISFORMUNF -> NOT SUPPORTED IN DX.               *
C     *                                                                *
C     * CALLED BY: USER PROGRAM                                        *
C     * TYPICAL USE:                                                   *
C     *         ...                                                    *
C     *         call visdxheadgen(...)                                 *
C     *         visoutput=VISOPENDX(...)                               *
C     *         ...                                                    *
C     *         do i=1,number_of_timesteps                             *
C     *            call visdataout(...)                                *
C     *         end do                                                 *
C     ******************************************************************
C
#include "cominou"
#include "comerrn"
#include "comform"
C
      INTEGER       VISOPENDX, VISFORMATS
      CHARACTER*(*) VISFILNAME
C
C     * LOCAL VARIABLES
      INTEGER IND
C
C     * FUNCTIONS
      LOGICAL VISNAMEOK
C     
C     * FILE NAME OK?
      IF (.NOT.VISNAMEOK(VISFILNAME)) THEN
         CALL VISPRWAR(VISWARNAME,'VISDXINT')
         VISFILNAME='dxfile'
      ENDIF
C
      IND=INDEX(VISFILNAME,'.gai')
      IF (IND.EQ.0) THEN
         IND               = INDEX(VISFILNAME,' ')
         IF (IND.EQ.0) IND = LEN(VISFILNAME)+1
      ENDIF
      OPEN(UNIT=VISDXOUTPU, FILE=VISFILNAME(1:IND-1),STATUS='UNKNOWN')
C
      VISOPENDX=VISDXOUTPU
C
      END
C
*DECK VISOPENAVS
      FUNCTION VISOPENAVS(VISFILNAME, VISFORMATS)
C     ******************************************************************
C     * THIS FUNCTION IS USED TO OPEN THE DATA FILE FOR AVS'S GENERAL  *
C     * FILE FORMAT. IT IS CALLED BY THE USER AFTER THE HEADER HAS     *
C     * BEEN WRITTEN. THE FUNCTIONS RETURNS THE UNIT NUMBER ON WHICH   *
C     * THE DATA FILE HAS BEEN OPENED FOR WRITING.                     *
C     *                                                                *
C     * ON INPUT:                                                      *
C     * VISFILNAME CONTAINS THE NAME OF THE DATA FILE TO BE OPENED.    *
C     * VISFORMATS CONTAINS THE FORM IN WHICH DATA WILL BE WRITTEN     *
C     *    VISFORMATS=VISFORMASC -> ASCII,                             *
C     *    VISFORMATS=VISFORMBIN -> BINARY (NOT IMPLEMENTED FOR        *
C     *                                     GENERAL FILE FORMAT),      *
C     *    VISFORMATS=VISFORMUNF -> FORTRAN UNFORMATTED.               *
C     *                                                                *
C     * CALLED BY: USER PROGRAM                                        *
C     * TYPICAL USE:                                                   *
C     *         ...                                                    *
C     *         call visavshead(...)                                   *
C     *         visoutput=VISOPENAVS(...)                              *
C     *         ...                                                    *
C     *         do i=1,number_of_timesteps                             *
C     *            call visdataout(...)                                *
C     *         end do                                                 *
C     ******************************************************************
C
#include "cominou"
#include "comerrn"
#include "comform"
C
      INTEGER       VISOPENAVS, VISFORMATS
      CHARACTER*(*) VISFILNAME
C
C     * LOCAL VARIABLES
      INTEGER IND
C
C     * FUNCTIONS
      LOGICAL VISNAMEOK
C     
C     * FILE NAME OK?
      IF (.NOT.VISNAMEOK(VISFILNAME)) THEN
         CALL VISPRWAR(VISWARNAME,'VISDXINT')
         VISFILNAME='avsfile'
      ENDIF
C
      IF (VISFORMATS.NE.VISFORMBIN) THEN
         IND=INDEX(VISFILNAME,'.fld')
         IF (IND.EQ.0) THEN
            IND               = INDEX(VISFILNAME,' ')
            IF (IND.EQ.0) IND = LEN(VISFILNAME)+1
         ENDIF
         IF (VISFORMATS.EQ.VISFORMASC) THEN
            OPEN(UNIT=VISAVSOUTP, FILE=VISFILNAME(1:IND-1),
     A           STATUS='UNKNOWN')
         ELSE IF (VISFORMATS.EQ.VISFORMUNF) THEN
            OPEN(UNIT=VISAVSOUTP, FORM='UNFORMATTED',
     A           FILE=VISFILNAME(1:IND-1),STATUS='UNKNOWN')
         ENDIF
      ENDIF
C
      VISOPENAVS=VISAVSOUTP
C
      END
C
*DECK VISDXWRDAT
      SUBROUTINE VISDXWRDAT(VISMATDATA, VISPOSITIM,
     A                    VISNV,      VISNI,      VISNJ,     VISNK,     
     B                    VISNT,      
     C                    VISFILNAME, VISRES,     VISNDIM,   VISCOOR,
     D                    VISSERIES,  VISSLICE,   VISVECTLEN,VISNFIELDS,
     E                    VISSTRUCTU, VISSTRUTYP, VISFORMATS,VISWORK)
C     ******************************************************************
C     * THIS SUBROUTINE PRINTS THE MATRIX VISMATDATA. THE LAST FORTRAN *
C     * FORTRAN INDEX VARIES FASTEST.                                  *
C     *                                                                *
C     * ON INPUT:                                                      *
C     * VISMATDATA CONTAINS THE UNCROPPED DATA.                        *
C     * VISPOSITIM CONTAINS THE TIME COORDINATES (LABELS).             *
C     * VISFIENAME CONTAINS THE NAMES OF THE DIFFERENT FIELDS.         *
C     * VISNV      CONTAINS THE SIZE OF THE VECTLEN DIMENSION.         *
C     * VISNI      CONTAINS THE SIZE OF THE FIRST COMPUTATIONAL        *
C     *            DIMENSION.                                          *
C     * VISNJ      CONTAINS THE SIZE OF THE SECOND COMPUTATIONAL       *
C     *            DIMENSION.                                          *
C     * VISNK      CONTAINS THE SIZE OF THE THIRD COMPUTATIONAL        *
C     *            DIMENSION.                                          *
C     * VISNT      CONTAINS THE SIZE OF THE FOURTH COMPUTATIONAL       *
C     *            DIMENSION.                                          *
C     * VISFILNAME CONTAINS THE FILENAME OF THE OUTPUT FILE.           *
C     * VISRES     CONTAINS THE RESOLUTION OF THE DATA SET             *
C     *            (VISRESINIT(VISCOORRAD)<=VISNI),                    *
C     *            (VISRESINIT(VISCOORPOL)<=VISNJ),                    *
C     *            (VISRESINIT(VISCOORTOR)<=VISNK),                    *
C     *            (VISRESINIT(VISCOORTIM)<=VISNT).                    *
C     * VISNDIM    CONTAINS THE NUMBER OF DIMENSIONS TO BE PRINTED.    *
C     * VISCOOR    CONTAINS THE DIMENSIONS WHICH WILL BE PRINTED.      *
C     * VISSERIES  DETERMINES IF THE TIME DEPENDENCE IS PRINTED.       *
C     *            VISSERIES=.TRUE. -> TIME DEPENDENCE USED AS SERIES  *
C     *                                ELEMENTS.                       *
C     * VISSLICE   CONTAINS THE POSITION (SLICE) NUMBERS OF THE        *
C     *            DIMENSIONS WHICH ARE SLICED AND HENCE WILL NOT BE   *
C     *            VISUALIZED.                                         *
C     * VISVECTLEN CONTAINS THE NUMBER OF DATA ELEMENTS AT A NODE      *
C     *            (VISVECTLEN<=VISNV).                                *
C     * VISFORMATS DETERMINES WHETHER OR NOT THE DATA WILL BE PRINTED  *
C     *            IN ASCII OR IN BINARY FORMAT.                       *
C     *            VISFORMATS=VISFORMASC -> ASCII                      *
C     *            VISFORMATS=VISFORMBIN -> BINARY                     *
C     * VISNFIELDS CONTAINS THE NUMBER OF FIELDS PER NODE.             *
C     * VISSTRUCTU DETERMINES THE TYPE OF THE FIELDS                   *
C     *            VISSTRUCTU()=VISRANKSCA -> SCALAR FIELD,            *
C     *            VISSTRUCTU()=VISRANKVC1 -> 1-VECTOR FIELD.          *
C     *            VISSTRUCTU()=VISRANKVC2 -> 2-VECTOR FIELD.          *
C     *            VISSTRUCTU()=VISRANKVC3 -> 3-VECTOR FIELD.          *
C     *            VISSTRUCTU()=VISRANKVEC -> 3-VECTOR FIELD.          *
C     * VISSTRUTYP CONTAINS THE PRIMITIVE TYPE OF THE FIELD            *
C     *            VISTYPEDOU -> DOUBLE                                *
C     *            VISTYPEFLO -> FLOAT                                 *
C     * VISWORK    IS A WORKING ARRAY SUPPLIED BY THE USER. ITS        *
C     *            DIMENSION SHOULD EQUAL VISNV*VISNI*VISNJ*VISNK*4    *
C     ******************************************************************
C
#include "comcoor"
#include "cominou"
#include "comform"
C
      LOGICAL VISSERIES
      INTEGER VISRES(*),  VISCOOR(*), VISSLICE(*),   VISNDIM,   
     A        VISVECTLEN, VISFORMATS, VISNFIELDS, VISSTRUCTU(*), 
     B        VISNV,      VISNI,      VISNJ,
     C        VISNK,      VISNT
      CHARACTER*(*) VISFILNAME, VISSTRUTYP(*)
      REAL          VISMATDATA(VISNV,VISNI,VISNJ,VISNK,VISNT),
     A              VISPOSITIM(VISNT),
     B              VISWORK(VISNV*VISNI*VISNJ*VISNK*VISNT)
C
C     * COUNTERS
      INTEGER  I, J, K, L, M,
     A        II, IFIELD,
     B        I0, I1, J0, J1, K0, K1, L0, L1, M0, M1
C
C     *LOCAL VARIABLES ASCII
      INTEGER      FIRST, LAST, DELTA
#ifdef CRAY
      INTEGER      IERR
#endif
      CHARACTER    STRING*100, FILENAME*80
C
C     *LOCAL VARIABLES ASCII
      INTEGER      NOBJ, IND
C
C     * FUNCTIONS
      INTEGER VISGETRES
#ifdef CRAY
      INTEGER CRAY2IEG
#endif
      CHARACTER*100 VISSTRING
C
C     * IF VISGETRES=1 LOOK UP SLICE NUMBER
      J1=VISGETRES(VISRES, VISNDIM, VISCOOR, VISCOORRAD)
      K1=VISGETRES(VISRES, VISNDIM, VISCOOR, VISCOORPOL)
      L1=VISGETRES(VISRES, VISNDIM, VISCOOR, VISCOORTOR)
      IF (VISSERIES) THEN
         M1=VISRES(VISCOORTIM)
      ELSE
         M1=1
      ENDIF
      J0=1
      K0=1
      L0=1
      M0=1
      IF (J1.EQ.1.AND.J1.LT.VISNI) THEN
          J0=VISSLICE(VISCOORRAD)
          J1=VISSLICE(VISCOORRAD)
      ENDIF
      IF (K1.EQ.1.AND.K1.LT.VISNJ) THEN
          K0=VISSLICE(VISCOORPOL)
          K1=VISSLICE(VISCOORPOL)
      ENDIF
      IF (L1.EQ.1.AND.L1.LT.VISNK) THEN
          L0=VISSLICE(VISCOORTOR)
          L1=VISSLICE(VISCOORTOR)
      ENDIF
      IF (M1.EQ.1.AND.M1.LT.VISNT) THEN
          M0=VISSLICE(VISCOORTIM)
          M1=VISSLICE(VISCOORTIM)
      ENDIF
C
      II = 0
      DO 70 M=M0,M1
         I0 = 1
         DO 60 IFIELD=1,VISNFIELDS
            I1 = I0 + VISSTRUCTU(IFIELD) - 1
C
C           * ASCII
            IF (VISFORMATS.EQ.VISFORMASC) THEN
C              FIRST=I0
C              LAST =I1
               FIRST=IFIELD
               LAST =IFIELD
               DELTA=LAST-FIRST
               STRING=VISSTRING(FIRST,LAST,VISSTRUCTU, VISSTRUTYP)
               IF (VISSERIES) THEN
                  IF (VISRES(VISCOORTIM).GT.1) THEN
                     WRITE(VISDXOUTPU,'(F16.8,1X)') VISPOSITIM(M)
                  ELSE
                     WRITE(VISDXOUTPU,'(F16.8,1X)') VISPOSITIM(1)
                  ENDIF
               ENDIF
               IF (VISRES(VISCOORTIM).GT.1) THEN
                  WRITE(VISDXOUTPU,STRING(1:(DELTA+1)*12+4))
     A            ((((VISMATDATA(I,J,K,L,M),I=I0,I1), 
     B                           L=L0,L1), K=K0,K1), J=J0,J1)
               ELSE
                  WRITE(VISDXOUTPU,STRING(1:(DELTA+1)*12+4))
     A            ((((VISMATDATA(I,J,K,L,1),I=I0,I1), 
     B                        L=L0,L1), K=K0,K1), J=J0,J1)
               ENDIF
C
C           * BINARY
            ELSE IF (VISFORMATS.EQ.VISFORMBIN) THEN
               IF (VISSERIES) THEN
                  II=II+1
                  VISWORK(II)=VISPOSITIM(M)
               ENDIF
               DO 20 J=J0,J1
                  DO 30 K=K0,K1
                     DO 40 L=L0,L1
                        DO 50 I=I0,I1
                           II=II+1
                           IF (VISRES(VISCOORTIM).EQ.1) THEN
                              VISWORK(II)=VISMATDATA(I,J,K,L,1)
                           ELSE
                              VISWORK(II)=VISMATDATA(I,J,K,L,M)
                           ENDIF
  50                    CONTINUE
  40                 CONTINUE
  30              CONTINUE
  20           CONTINUE
C
            ENDIF
C
            I0 = I0 + VISSTRUCTU(IFIELD)
  60     CONTINUE
  70  CONTINUE 
C
      IF (VISFORMATS.EQ.VISFORMBIN) THEN
         NOBJ=II
C
C        ** WRITING BINARY FILE USING 
C
C        * CHECKING FILENAME, ADDING NULL CHARACTER
         FILENAME = VISFILNAME
         IND = INDEX(VISFILNAME,'.dx')
         IF (IND.EQ.0) THEN
            IND = MIN(77,INDEX(VISFILNAME,' '))
            IF (IND.EQ.0) IND = MAX(1,LEN(VISFILNAME)-2)
            FILENAME(IND:IND+2) = '.dx'
         ENDIF
         IND = MIN(77,IND)
         FILENAME(IND+3:IND+3)=CHAR(0)
#ifdef CRAY
         IERR=CRAY2IEG(8,NOBJ,VISWORK,0,VISWORK)
         IF (IERR.GT.0) CALL VISPRERR(VISERRIEEE,'VISDXWRDAT')
         CALL DXBINWRITEFLO(FILENAME,VISWORK,NOBJ,8)
#else
         CALL DXBINWRITEDOU(FILENAME,VISWORK,NOBJ,8)
#endif
      ENDIF
      END
C
*DECK VISDXWRCOO
      SUBROUTINE VISDXWRCOO(VISPOSX,    VISPOSY,    VISPOSZ,  
     A                    VISNI,      VISNJ,      VISNK,     
     B                    VISFILNAME, VISRES,     VISNDIM,   VISCOOR,
     C                    VISSLICE,   VISNFIELDS, VISSTRUCTU,VISSTRUTYP,
     D                    VISFORMATS, VISWORK)
C     ******************************************************************
C     * THIS SUBROUTINE PRINTS THE CARTESIAN COORDINATES. THE LAST     *
C     * FORTRAN INDEX VARIES FASTEST.                                  *
C     *                                                                *
C     * ON INPUT:                                                      *
C     * VISPOSX    CONTAINS THE X-COORDINATES.                         *
C     * VISPOSY    CONTAINS THE Y-COORDINATES.                         *
C     * VISPOSZ    CONTAINS THE Z-COORDINATES.                         *
C     * VISNI      CONTAINS THE SIZE OF THE FIRST COMPUTATIONAL        *
C     *            DIMENSION.                                          *
C     * VISNJ      CONTAINS THE SIZE OF THE SECOND COMPUTATIONAL       *
C     *            DIMENSION.                                          *
C     * VISNK      CONTAINS THE SIZE OF THE THIRD COMPUTATIONAL        *
C     *            DIMENSION.                                          *
C     * VISFILNAME CONTAINS THE FILENAME OF THE OUTPUT FILE.           *
C     * VISRES     CONTAINS THE RESOLUTION OF THE DATA SET             *
C     *            (VISRESINIT(VISCOORRAD)<=VISNI),                    *
C     *            (VISRESINIT(VISCOORPOL)<=VISNJ),                    *
C     *            (VISRESINIT(VISCOORTOR)<=VISNK),                    *
C     *            (VISRESINIT(VISCOORTIM)<=VISNT).                    *
C     * VISNDIM    CONTAINS THE NUMBER OF DIMENSIONS TO BE PRINTED.    *
C     * VISCOOR    CONTAINS THE DIMENSIONS WHICH WILL BE PRINTED.      *
C     * VISSLICE   CONTAINS THE POSITION (SLICE) NUMBERS OF THE        *
C     *            DIMENSIONS WHICH ARE SLICED AND HENCE WILL NOT BE   *
C     *            VISUALIZED.                                         *
C     * VISFORMATS DETERMINES WHETHER OR NOT THE DATA WILL BE PRINTED  *
C     *            IN ASCII OR IN BINARY FORMAT.                       *
C     *            VISFORMATS=VISFORMASC -> ASCII                      *
C     *            VISFORMATS=VISFORMBIN -> BINARY                     *
C     * VISNFIELDS CONTAINS THE NUMBER OF FIELDS PER NODE.             *
C     * VISSTRUCTU DETERMINES THE TYPE OF THE FIELDS                   *
C     *            VISSTRUCTU()=VISRANKSCA -> SCALAR FIELD,            *
C     *            VISSTRUCTU()=VISRANKVC1 -> 1-VECTOR FIELD.          *
C     *            VISSTRUCTU()=VISRANKVC2 -> 2-VECTOR FIELD.          *
C     *            VISSTRUCTU()=VISRANKVC3 -> 3-VECTOR FIELD.          *
C     *            VISSTRUCTU()=VISRANKVEC -> 3-VECTOR FIELD.          *
C     * VISSTRUTYP CONTAINS THE PRIMITIVE TYPE OF THE FIELD            *
C     *            VISTYPEDOU -> DOUBLE                                *
C     *            VISTYPEFLO -> FLOAT                                 *
C     * VISWORK    IS A WORKING ARRAY SUPPLIED BY THE USER. ITS        *
C     *            DIMENSION SHOULD EQUAL 3*VISNI*VISNJ*VISNK.         *
C     ******************************************************************
C
#include "comcoor"
#include "cominou"
#include "comform"
C
      INTEGER VISRES(*),  VISCOOR(*), VISSLICE(*),   VISNDIM, 
     A        VISFORMATS, VISNFIELDS, VISSTRUCTU(*), 
     B        VISNI,      VISNJ,      VISNK
      CHARACTER*(*) VISFILNAME, VISSTRUTYP(*)
      REAL  VISPOSX(VISNI,VISNJ,VISNK),
     A      VISPOSY(VISNI,VISNJ,VISNK),
     B      VISPOSZ(VISNI,VISNJ,VISNK),
     C      VISWORK(3*VISNI*VISNJ*VISNK)
C
C     * COUNTERS
      INTEGER  I, J, K,
     A        II,
     B        I0, I1, J0, J1, K0, K1
C
C     *LOCAL VARIABLES ASCII
      INTEGER      FIRST, LAST, DELTA
      CHARACTER    STRING*100, FILENAME*80
C
C     *LOCAL VARIABLES BINARY
#ifdef CRAY
      REAL    TMP
      INTEGER IERR
#endif
      INTEGER NOBJ, IND
C
C     * FUNCTIONS
      INTEGER VISGETRES
#ifdef CRAY
      INTEGER CRAY2IEG
#endif
      CHARACTER*100 VISSTRING
C
C     * IF VISGETRES=1 LOOK UP SLICE NUMBER
      I1=VISGETRES(VISRES, VISNDIM, VISCOOR, VISCOORRAD)
      J1=VISGETRES(VISRES, VISNDIM, VISCOOR, VISCOORPOL)
      K1=VISGETRES(VISRES, VISNDIM, VISCOOR, VISCOORTOR)
      I0=1
      J0=1
      K0=1
      IF (I1.EQ.1.AND.I1.LT.VISNI) THEN
          I0=VISSLICE(VISCOORRAD)
          I1=VISSLICE(VISCOORRAD)
      ENDIF
      IF (J1.EQ.1.AND.J1.LT.VISNJ) THEN
          J0=VISSLICE(VISCOORPOL)
          J1=VISSLICE(VISCOORPOL)
      ENDIF
      IF (K1.EQ.1.AND.K1.LT.VISNK) THEN
          K0=VISSLICE(VISCOORTOR)
          K1=VISSLICE(VISCOORTOR)
      ENDIF
C
C     * ASCII
      IF (VISFORMATS.EQ.VISFORMASC) THEN
         FIRST=VISNFIELDS+1
         LAST =VISNFIELDS+1
         DELTA=LAST-FIRST
         STRING=VISSTRING(FIRST, LAST, VISSTRUCTU, VISSTRUTYP)
         WRITE(VISDXOUTPU,STRING(1:(DELTA+1)*12+4))
     A         (((VISPOSX(I,J,K), VISPOSY(I,J,K), VISPOSZ(I,J,K),          
     B                    K=K0,K1), J=J0,J1), I=I0,I1)
C
C     * BINARY
      ELSE IF (VISFORMATS.EQ.VISFORMBIN) THEN
         II=1  
         DO 10 I=I0,I1
            DO 20 J=J0,J1
               DO 30 K=K0,K1
                  VISWORK(II  )=VISPOSX(I,J,K)
                  VISWORK(II+1)=VISPOSY(I,J,K)
                  VISWORK(II+2)=VISPOSZ(I,J,K) 
                  II           =II+3 
   30          CONTINUE 
   20       CONTINUE
   10    CONTINUE
C
         NOBJ= VISGETRES(VISRES, VISNDIM, VISCOOR, VISCOORRAD)
     A        *VISGETRES(VISRES, VISNDIM, VISCOOR, VISCOORPOL)
     B        *VISGETRES(VISRES, VISNDIM, VISCOOR, VISCOORTOR)
     C        *3
C
C        ** WRITING BINARY FILE 
C
C        * CHECKKING FILENAME, ADDING NULL CHARACTER
         FILENAME = VISFILNAME
         IND = INDEX(VISFILNAME,'.dx')
         IF (IND.EQ.0) THEN
            IND = MIN(77,INDEX(VISFILNAME,' '))
            IF (IND.EQ.0) IND = MAX(1,LEN(VISFILNAME)-2)
            FILENAME(IND:IND+2) = '.dx'
         ENDIF
         IND = MIN(77,IND)
         FILENAME(IND+3:IND+3)=CHAR(0)
#ifdef CRAY
C        * SOMETHING GOES WRONG IF WE WRITE ALL NOBJ ELEMENTS 
C          AT ONCE. IT IS CURRENTLY FIXED BY WRITING TWICE A HALF
C          OF THE DATA SET.
         IERR = CRAY2IEG(2,NOBJ/2,VISWORK,0,VISWORK)
         IF (IERR.GT.0) CALL VISPRERR(VISERRIEEE,'VISDXWRCOO')
         CALL DXBINWRITEFLO(FILENAME,VISWORK,NOBJ/2,4)
         IERR = CRAY2IEG(2,NOBJ-NOBJ/2,VISWORK,0,VISWORK(NOBJ/2+1))
         IF (IERR.GT.0) CALL VISPRERR(VISERRIEEE,'VISDXWRCOO')
         CALL DXBINWRITEFLO(FILENAME,VISWORK,NOBJ-NOBJ/2,4)
#else
         CALL DXBINWRITEDOU(FILENAME,VISWORK,NOBJ,4)
#endif
      ENDIF
      END
C
*DECK VISDATAOUT
      SUBROUTINE VISDATAOUT(VISMATDATA,
     A                      VISPOSX,    VISPOSY,    VISPOSZ, VISPOSITIM,
     B                      VISNV,      VISNI,      VISNJ,      VISNK,   
     C                      VISNT,
     D                      VISFILNAME, VISRES,     VISNDIM,    VISCOOR,
     E                      VISSERIES,  VISSLICE,   VISVECTLEN,
     F                      VISNFIELDS, VISFIENAME, VISSTRUCTU, 
     G                      VISSTRUTYP, VISFORMATS, VISXDR,   VISOUTPUT,
     H                      VISWORK,    VISWORK2)
C     ******************************************************************
C     * THIS ROUTINE IS USED TO WRITE THE DATA FILE FOR THE GENERAL    *
C     * FILE FORMATS AND AVS'S INTERNAL FIELD FORMAT.                  *
C     *                                                                *
C     * ON INPUT:                                                      *
C     * VISMATDATA CONTAINS THE DATA.                                  *
C     * VISPOSX    CONTAINS THE X-COORDINATES.                         *
C     * VISPOSY    CONTAINS THE Y-COORDINATES.                         *
C     * VISPOSZ    CONTAINS THE Z-COORDINATES.                         *
C     * VISPOSITIM CONTAINS THE TIME COORDINATES (LABELS).             *
C     * VISNV      CONTAINS THE SIZE OF THE VECTLEN DIMENSION.         *
C     * VISNI      CONTAINS THE SIZE OF THE FIRST COMPUTATIONAL        *
C     *            DIMENSION.                                          *
C     * VISNJ      CONTAINS THE SIZE OF THE SECOND COMPUTATIONAL       *
C     *            DIMENSION.                                          *
C     * VISNK      CONTAINS THE SIZE OF THE THIRD COMPUTATIONAL        *
C     *            DIMENSION.                                          *
C     * VISNT      CONTAINS THE SIZE OF THE FOURTH COMPUTATIONAL       *
C     *            DIMENSION.                                          *
C     * VISFILNAME CONTAINS THE FILENAME OF THE OUTPUT FILE.           *
C     * VISRES     CONTAINS THE RESOLUTION OF THE DATA SET             *
C     *            (VISRESINIT(VISCOORRAD)<=VISNI),                    *
C     *            (VISRESINIT(VISCOORPOL)<=VISNJ),                    *
C     *            (VISRESINIT(VISCOORTOR)<=VISNK),                    *
C     *            (VISRESINIT(VISCOORTIM)<=VISNT).                    *
C     * VISNDIM    CONTAINS THE NUMBER OF DIMENSIONS TO BE PRINTED.    *
C     * VISCOOR    CONTAINS THE DIMENSIONS WHICH WILL BE PRINTED.      *
C     * VISSERIES  DETERMINES IF THE TIME DEPENDENCE IS PRINTED.       *
C     *            VISSERIES=.TRUE. -> TIME DEPENDENCE USED AS SERIES  *
C     *                                ELEMENTS.                       *
C     *            THIS VARIABLE IS ONLY USED AFTER A CALL TO          *
C     *            VISDXGENHEAD. WHEN USED IN CONJUNCTION WITH         *
C     *            VISAVSHEAD IT SHOULD BE SET EQUAL TO .FALSE. .      *
C     * VISSLICE   CONTAINS THE POSITION (SLICE) NUMBERS OF THE        *
C     *            DIMENSIONS WHICH ARE SLICED AND HENCE WILL NOT BE   *
C     *            VISUALIZED.                                         *
C     * VISVECTLEN CONTAINS THE NUMBER OF DATA ELEMENTS AT A NODE      *
C     *            (VISVECTLEN<=VISNV).                                *
C     * VISNFIELDS CONTAINS THE NUMBER OF FIELDS PER NODE.             *
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
C     * VISFORMATS DETERMINES WHETHER OR NOT THE DATA WILL BE PRINTED  *
C     *            IN ASCII OR IN BINARY FORMAT.                       *
C     *            VISFORMATS=VISFORMASC -> ASCII                      *
C     *            VISFORMATS=VISFORMBIN -> BINARY                     *
C     * VISXDR     DETERMINES IF BINARY IS IN XDR FORMAT OR NOT        *
C     *            (AVS ONLY)                                          *
C     * VISOUTPUT  CONTAINS THE UNIT NUMBER ON WHICH TO WRITE.         *
C     *            IF VISOUTPUT<0 IT INDICATES THE LAST CALL TO THIS   *
C     *            ROUTINE. THIS INDICATION IS NECESSARY TO WRITE THE  *
C     *            INTERNAL AVS FIELD FILE AT THE RIGHT TIME WHEN ONE  *
C     *            TIMESTEP A TIME IS BEING WRITTEN.                   *
C     * VISWORK    IS A WORKING ARRAY SUPPLIED BY THE USER. ITS        *
C     *            DIMENSION SHOULD EQUAL VISNV*VISNI*VISNJ*VISNK*4    *
C     *            THE USER MUST ALSO DECLARE IT AS A SAVED ARRAY,     *
C     *            I.E. SAVE VISWORK.                                  *
C     * VISWORK2   IS A WORKING ARRAY SUPPLIED BY THE USER. ITS        *
C     *            DIMENSION SHOULD EQUAL 3*VISNI*VISNJ*VISNK.         *
C     *            THE USER MUST ALSO DECLARE IT AS A SAVED ARRAY,     *
C     *            I.E. SAVE VISWORK2.                                 *
C     *                                                                *
C     * VISDATAOUT IS CALLED BY THE USER PROGRAM. IT IS USED IN        *
C     * CONJUNCTION WITH VISDXGENHEAD AND VISAVS(GENERAL FILE FORMAT). *
C     ******************************************************************
C
#include "comcoor"
#include "cominou"
#include "comrank"
#include "comform"
#include "comtype"
C
      LOGICAL       VISSERIES,  VISXDR
      INTEGER       VISRES(*),  VISNDIM,    VISCOOR(*),     VISSLICE(*),
     A              VISVECTLEN, VISNFIELDS, VISSTRUCTU(*),  VISFORMATS,
     B              VISOUTPUT,  VISNV,      VISNI,          VISNJ,          
     C              VISNK,      VISNT
      CHARACTER*(*) VISFILNAME, VISFIENAME(*),  VISSTRUTYP(*)
      REAL          VISMATDATA(VISNV,VISNI,VISNJ,VISNK,VISNT),
     A              VISPOSX(VISNI,VISNJ,VISNK),
     B              VISPOSY(VISNI,VISNJ,VISNK),
     C              VISPOSZ(VISNI,VISNJ,VISNK),
     D              VISPOSITIM(VISNT),
     E              VISWORK(VISNV*VISNI*VISNJ*VISNK*4),
     F              VISWORK2(VISNI*VISNJ*VISNK,3)
      COMMON /CTIM/ VISTIME1,   VISTRES,    VISNTST,
     A              VISFIRS,    VISLAST,    VIS1TST
      LOGICAL       VISFIRS,    VISLAST,    VIS1TST
      INTEGER       VISTRES,    VISNTST
      REAL          VISTIME1
C
C     * COUNTERS
      INTEGER I,J,K,L,M,
     A        II, JJ, SUMII, SUMJJ,
     B        J0, J1, K0, K1, L0, L1, M0, M1
C
C     * LOCAL VARIABLES
      INTEGER   IND, FIRSTCOLUMN, LASTCOLUMN, 
     A          OFFSET, DOFFSET, XOFFSET, YOFFSET, ZOFFSET,
     B          DLENGTH, XLENGTH, YLENGTH, ZLENGTH
#ifdef CRAY
      INTEGER      IERR
#endif
      CHARACTER STRING*40, FILENAME*80
      LOGICAL   FIRSTCALL, LASTCALL
C
C     * INITIAL VALUES
      DATA  FIRSTCALL/.TRUE./
     A      LASTCALL /.FALSE./
     B      SUMII   /0/
     C      SUMJJ   /0/
     D      DOFFSET /0/
     E      XOFFSET /0/
     F      YOFFSET /0/
     G      DLENGTH /0/
     H      XLENGTH /0/
     I      YLENGTH /0/
C
C     * FUNCTIONS
      INTEGER VISGETRES, VISSIZEOF, findoffset, INDEX
#ifdef CRAY
      INTEGER CRAY2IEG
#endif
      CHARACTER*100 VISSTRING
C
C     * SAVING VARIABLES FOR SUBSEQUENT CALLS (CALLS USED TO OUTPUT 
C     * SINGLE TIMESTEPS)
      SAVE  SUMII, SUMJJ, FIRSTCALL, 
     A      DOFFSET, XOFFSET, YOFFSET, ZOFFSET,
     B      DLENGTH, XLENGTH, YLENGTH
C
C     * IF VISGETRES=1 LOOK UP SLICE NUMBER
      J1=VISGETRES(VISRES, VISNDIM, VISCOOR, VISCOORRAD)
      K1=VISGETRES(VISRES, VISNDIM, VISCOOR, VISCOORPOL)
      L1=VISGETRES(VISRES, VISNDIM, VISCOOR, VISCOORTOR)
      IF (VISSERIES) THEN
         M1=VISRES(VISCOORTIM)
      ELSE
         M1=VISGETRES(VISRES, VISNDIM, VISCOOR, VISCOORTIM)
      ENDIF
      J0=1
      K0=1
      L0=1
      M0=1
      IF (J1.EQ.1.AND.J1.LT.VISNI) THEN
          J0=VISSLICE(VISCOORRAD)
          J1=VISSLICE(VISCOORRAD)
      ENDIF
      IF (K1.EQ.1.AND.K1.LT.VISNJ) THEN
          K0=VISSLICE(VISCOORPOL)
          K1=VISSLICE(VISCOORPOL)
      ENDIF
      IF (L1.EQ.1.AND.L1.LT.VISNK) THEN
          L0=VISSLICE(VISCOORTOR)
          L1=VISSLICE(VISCOORTOR)
      ENDIF
      IF (M1.EQ.1.AND.M1.LT.VISNT) THEN
          M0=VISSLICE(VISCOORTIM)
          M1=VISSLICE(VISCOORTIM)
      ENDIF
C
C     * LAST CALL ?
      IF (VISOUTPUT.LT.0.OR.VISLAST) THEN
         VISOUTPUT=IABS(VISOUTPUT)
         LASTCALL=.TRUE.
      ENDIF
C
C     * WRITE INFO IF FORMAT=ASCII (ONLY ON FIRST CALL)
      IF (FIRSTCALL) THEN
         IF (VISFORMATS.EQ.VISFORMASC) THEN
            FIRSTCOLUMN=1
            DO 10 I=1,VISNFIELDS
               LASTCOLUMN =FIRSTCOLUMN+VISSTRUCTU(I)-1
               WRITE(VISOUTPUT,200) FIRSTCOLUMN,LASTCOLUMN,VISFIENAME(I)
               FIRSTCOLUMN=FIRSTCOLUMN+VISSTRUCTU(I)
   10       CONTINUE
            WRITE(VISOUTPUT,200) FIRSTCOLUMN, FIRSTCOLUMN+2,
     A                           VISFIENAME(VISNFIELDS+1)
         ELSE IF (VISFORMATS.EQ.VISFORMBIN) THEN
C
C           * CHECKING FILENAME, ADDING NULL CHARACTER
            FILENAME = VISFILNAME
            IND = INDEX(FILENAME,'.fld')
            IF (IND.EQ.0) THEN
               IND = MIN(76,INDEX(FILENAME,' '))
               IF (IND.EQ.0) IND = MAX(1,LEN(FILENAME)-3)
            ENDIF
            IND = MIN(76,IND)
            FILENAME(IND:IND+3)   = '.fld'
            FILENAME(IND+4:IND+4) = CHAR(0)
C
C           * DETERMINING THE OFFSETS OF THE VARIOUS BINARY PARTS
            DOFFSET=findoffset(FILENAME)
C
C           ** IF ONE TIMESTEP A TIME, NUMBER OF TIMESTEPS TO BE 
C              WRITTEN OUT IS STORED IN VISTRES.
            IF (VIS1TST) THEN
               DLENGTH=VISVECTLEN*(J1-J0+1)*(K1-K0+1)*(L1-L0+1)
     A                 *VISTRES*VISSIZEOF(VISSTRUTYP(1))
            ELSE
               DLENGTH=VISVECTLEN*(J1-J0+1)*(K1-K0+1)*(L1-L0+1)
     A                 *(M1-M0+1)*VISSIZEOF(VISSTRUTYP(1))
            ENDIF
            IF (VIS1TST) THEN
               XLENGTH=(J1-J0+1)*(K1-K0+1)*(L1-L0+1)
     A                 *VISTRES*VISSIZEOF(VISSTRUTYP(VISNFIELDS+1))
            ELSE
               XLENGTH=(J1-J0+1)*(K1-K0+1)*(L1-L0+1)
     A                 *(M1-M0+1)*VISSIZEOF(VISSTRUTYP(VISNFIELDS+1))
            ENDIF
            YLENGTH=XLENGTH
            ZLENGTH=XLENGTH
            XOFFSET=DOFFSET+DLENGTH
            YOFFSET=DOFFSET+DLENGTH+XLENGTH
            ZOFFSET=DOFFSET+DLENGTH+XLENGTH+YLENGTH
C
         ENDIF
         FIRSTCALL=.FALSE.
      ENDIF
C
      IF (VISFORMATS.EQ.VISFORMASC) THEN
C        * WRITE DATA AND POSITIONS IN FIELD FORMAT
         STRING=VISSTRING(1,VISNFIELDS+1, VISSTRUCTU, VISSTRUTYP)
         DO 30 M=M0,M1
            DO 40 L=L0,L1
               DO 50 K=K0,K1
                  DO 60 J=J0,J1
                     WRITE(VISOUTPUT,STRING(1:(VISNFIELDS+1)*12+4)) 
     A                    (VISMATDATA(I,J,K,L,M), I=1,VISVECTLEN), 
     B                     VISPOSX(J,K,L),
     C                     VISPOSY(J,K,L),
     D                     VISPOSZ(J,K,L)
   60             CONTINUE
   50          CONTINUE
   40       CONTINUE
   30    CONTINUE
C
C     * FORTRAN UNFORMATTED DATA FOR AVS GENERAL FIELD FORMAT
      ELSE IF (VISFORMATS.EQ.VISFORMUNF) THEN
          WRITE(VISOUTPUT) (((((VISMATDATA(I,J,K,L,M), I=1,VISVECTLEN), 
     B                          VISPOSX(J,K,L),
     C                          VISPOSY(J,K,L),
     D                          VISPOSZ(J,K,L),
     E                          J=J0,J1), K=K0,K1), L=L0,L1), M=M0,M1)
C
C     * BINARY DATA FOR AVS INTERNAL FIELD FORMAT
      ELSE IF (VISFORMATS.EQ.VISFORMBIN) THEN
C
C        * CHECKING FILENAME, ADDING NULL CHARACTER
         FILENAME = VISFILNAME
         IND = INDEX(FILENAME,'.fld')
         IF (IND.EQ.0) THEN
            IND = MIN(76,INDEX(FILENAME,' '))
            IF (IND.EQ.0) IND = MAX(1,LEN(FILENAME)-3)
         ENDIF
         IND = MIN(76,IND)
         FILENAME(IND:IND+3)   = '.fld'
         FILENAME(IND+4:IND+4) = CHAR(0)
C
         II=0
         JJ=0        
         DO 70 M=M0,M1
            DO 80 L=L0,L1
               DO 90 K=K0,K1
                  DO 100 J=J0,J1
                     DO 110 I=1,VISVECTLEN
                        II=II+1
                        VISWORK(II)=VISMATDATA(I,J,K,L,M)
  110                CONTINUE
                     JJ=JJ+1
                     VISWORK2(JJ,1)=VISPOSX(J,K,L)
                     VISWORK2(JJ,2)=VISPOSY(J,K,L)
                     VISWORK2(JJ,3)=VISPOSZ(J,K,L)
  100             CONTINUE
   90          CONTINUE
   80       CONTINUE
   70    CONTINUE
C
         IF (VISSTRUTYP(1).EQ.VISTYPEFLO) THEN
#ifdef CRAY
            CALL VISPRERR(VISERRCFLO,'VISDATAOUT')
#else
            IF (VISXDR) CALL VISPRERR(VISERRCXDR,'VISDATAOUT')
C
            OFFSET=DOFFSET+SUMII
            CALL AVSBINWRITEFLO(FILENAME,VISWORK,II,OFFSET,VISXDR)
            OFFSET=XOFFSET+SUMJJ
            CALL AVSBINWRITEFLO(FILENAME,VISWORK2(1,1),JJ,OFFSET,
     A                          VISXDR)
            OFFSET=YOFFSET+SUMJJ
            CALL AVSBINWRITEFLO(FILENAME,VISWORK2(1,2),JJ,OFFSET,
     A                          VISXDR)
            OFFSET=ZOFFSET+SUMJJ
            CALL AVSBINWRITEFLO(FILENAME,VISWORK2(1,3),JJ,OFFSET,
     A                          VISXDR)
#endif
         ELSE IF (VISSTRUTYP(1).EQ.VISTYPEDOU) THEN
            IF (VISXDR) THEN
#ifdef CRAY
C
C              * WRITING RECORD CONTAINING DATA VALUES
               IERR=CRAY2IEG(8,II,VISWORK,0,VISWORK)
               IF (IERR.GT.0) CALL VISPRERR(VISERRIEEE,'VISDATAOUT')
               OFFSET=DOFFSET+SUMII
               CALL AVSBINWRITEDOU(FILENAME,VISWORK,II,OFFSET,VISXDR)
C
C              * WRITING RECORD CONTAINING X-COORDINATES
               IERR=CRAY2IEG(2,JJ/2,VISWORK2,0,VISWORK2)
               IF (IERR.GT.0) CALL VISPRERR(VISERRIEEE,'VISDATAOUT')
               OFFSET=XOFFSET+SUMJJ
               CALL AVSBINWRITEFLO(FILENAME,VISWORK2(1,1),JJ/2,OFFSET,
     A                             VISXDR)
               IERR=CRAY2IEG(2,JJ-JJ/2,VISWORK2(1,1),0,
     A                                 VISWORK2(JJ/2+1,1))
               IF (IERR.GT.0) CALL VISPRERR(VISERRIEEE,'VISDATAOUT')
               OFFSET=XOFFSET+SUMJJ+(JJ/2)*VISSIZEOF(VISSTRUTYP(1))
               CALL AVSBINWRITEFLO(FILENAME,VISWORK2(1,1),
     A                             JJ-JJ/2,OFFSET,VISXDR)
C
C              * WRITING RECORD CONTAINING Y-COORDINATES
               IERR=CRAY2IEG(2,JJ/2,VISWORK2(1,2),0,VISWORK2(1,2))
               IF (IERR.GT.0) CALL VISPRERR(VISERRIEEE,'VISDATAOUT')
               OFFSET=YOFFSET+SUMJJ
               CALL AVSBINWRITEFLO(FILENAME,VISWORK2(1,2),JJ/2,OFFSET,
     A                             VISXDR)
               IERR=CRAY2IEG(2,JJ-JJ/2,VISWORK2(1,2),0,
     A                                 VISWORK2(JJ/2+1,2))
               IF (IERR.GT.0) CALL VISPRERR(VISERRIEEE,'VISDATAOUT')
               OFFSET=YOFFSET+SUMJJ+(JJ/2)*VISSIZEOF(VISSTRUTYP(1))
               CALL AVSBINWRITEFLO(FILENAME,VISWORK2(1,2),
     A                             JJ-JJ/2,OFFSET,VISXDR)
C
C              * WRITING RECORD CONTAINING Z-COORDINATES
               IERR=CRAY2IEG(2,JJ/2,VISWORK2(1,3),0,VISWORK2(1,3))
               IF (IERR.GT.0) CALL VISPRERR(VISERRIEEE,'VISDATAOUT')
               OFFSET=ZOFFSET+SUMJJ
               CALL AVSBINWRITEFLO(FILENAME,VISWORK2(1,3),JJ/2,OFFSET,
     A                             VISXDR)
               IERR=CRAY2IEG(2,JJ-JJ/2,VISWORK2(1,3),0,
     A                                 VISWORK2(JJ/2+1,3))
               IF (IERR.GT.0) CALL VISPRERR(VISERRIEEE,'VISDATAOUT')
               OFFSET=ZOFFSET+SUMJJ+(JJ/2)*VISSIZEOF(VISSTRUTYP(1))
               CALL AVSBINWRITEFLO(FILENAME,VISWORK2(1,3),
     A                             JJ-JJ/2,OFFSET,VISXDR)
#else
               CALL VISPRERR(VISERRCXDR,'VISDATAOUT')
#endif
            ELSE
#ifdef CRAY
               OFFSET=DOFFSET+SUMII
               CALL AVSBINWRITEFLO(FILENAME,VISWORK,II,OFFSET,VISXDR)
               OFFSET=XOFFSET+SUMJJ
               CALL AVSBINWRITEFLO(FILENAME,VISWORK2(1,1),JJ,OFFSET,
     A                             VISXDR)
               OFFSET=YOFFSET+SUMJJ
               CALL AVSBINWRITEFLO(FILENAME,VISWORK2(1,2),JJ,OFFSET,
     A                             VISXDR)
               OFFSET=ZOFFSET+SUMJJ
               CALL AVSBINWRITEFLO(FILENAME,VISWORK2(1,3),JJ,OFFSET,
     A                             VISXDR)
#else
               OFFSET=DOFFSET+SUMII
               CALL AVSBINWRITEDOU(FILENAME,VISWORK,II,OFFSET,VISXDR)
               OFFSET=XOFFSET+SUMJJ
               CALL AVSBINWRITEDOU(FILENAME,VISWORK2(1,1),JJ,OFFSET,
     A                             VISXDR)
               OFFSET=YOFFSET+SUMJJ
               CALL AVSBINWRITEDOU(FILENAME,VISWORK2(1,2),JJ,OFFSET,
     A                             VISXDR)
               OFFSET=ZOFFSET+SUMJJ
               CALL AVSBINWRITEDOU(FILENAME,VISWORK2(1,3),JJ,OFFSET,
     A                             VISXDR)
#endif
            ENDIF
         ENDIF
      ENDIF
C
C     * SUB TOTALS OF WRITTEN ITEMS
      SUMII=SUMII+II*VISSIZEOF(VISSTRUTYP(1))
      SUMJJ=SUMJJ+JJ*VISSIZEOF(VISSTRUTYP(1))
C
      IF (LASTCALL) THEN
         FIRSTCALL=.TRUE.
         LASTCALL =.FALSE.
      ENDIF
C
C     * FORMATS
  200 FORMAT('COLUMNS ',I2,' -',I2,': ',A)
      END
C
*DECK VISGETRES
C     ==================================================================
C     = FUNCTIONS AND SUBROUTINES WHICH ARE CALLED BY THE ABOVE        =
C     = SUBROUTINES. THEY SHOULD NOT BE CALLED BY THE USER. THEY MAY   =
C     = FURTHERMORE CONTAIN HARDWARE DEPENDENT INFORMATION!!!!!!       =
C     =                              |                                 =
C     =                              |                                 =
C     =                              |                                 =
C     =                            \   /                               =
C     =                             \ /                                =
C     =                              .                                 =
C     ==================================================================
      FUNCTION VISGETRES(VISRES, VISNDIM, VISCOOR, WHICHCOOR)
C     ******************************************************************
C     * THIS FUNCTION RETURNS THE RESOLUTION BELONGING TO WHICHCOOR    *
C     * IF WHICHCOOR HAS TO BE OUTPUTTED. OTHERWISE A VALUE OF 1 IS    *
C     * RETURNED.                                                      *
C     *                                                                *
C     * ON INPUT:                                                      *
C     * VISRES     CONTAINS THE RESOLUTION OF THE DATA                 *
C     * VISNDIM    CONTAINS THE NUMBER OF DIMENSION TO BE PRINTED      *
C     * VISCOOR    CONTAINS THE DIMENSIONS WHICH WILL BE PRINTED       *
C     ******************************************************************
C
      INTEGER VISGETRES, VISRES(*), VISCOOR(*), VISNDIM, WHICHCOOR
C
C     * COUNTERS
      INTEGER I
C
      DO 10 I=1, VISNDIM
         IF (VISCOOR(I).EQ.WHICHCOOR) THEN
            VISGETRES = VISRES(WHICHCOOR)
            RETURN
         ENDIF
   10 CONTINUE
      VISGETRES = 1
      END
C
*DECK VISNAMEOK
      FUNCTION VISNAMEOK(NAME)
C     ******************************************************************
C     * FUNCTION TO CHECK IF A CHARACTER STRING STARTS WITH A LETTER.  *
C     * IF THIS IS TRUE THEN .TRUE. WILL BE RETURNED, OTHERWISE .FALSE.*
C     *                                                                *
C     * ON INPUT:                                                      *
C     * NAME     THE CHARACTER STRING                                  *
C     ******************************************************************
C
      CHARACTER*(*) NAME
      LOGICAL       VISNAMEOK
C
      VISNAMEOK=.TRUE.
      IF (.NOT.((NAME(1:1).GE.'A'.AND.NAME(1:1).LE.'Z')
     A    .OR.(NAME(1:1).GE.'a'.AND.NAME(1:1).LE.'z'))) THEN
         VISNAMEOK=.FALSE.
      ENDIF
      END   
C
*DECK VISPRWAR
      SUBROUTINE VISPRWAR(VISERRORNO,VISSUBNAME)
C     ******************************************************************
C     * THIS ROUTINE PRINTS WARNING MESSAGES. VISERRORNO CONTAINS THE  *
C     * NUMBER OF THE MESSAGE. VISSUBNAME CONTAINS THE NAME OF THE     *
C     * CALLING SUBROUTINE.                                            *
C     ******************************************************************
C
#include "comerrn"
#include "cominou"
C
      CHARACTER*(*) VISSUBNAME
      INTEGER       VISERRORNO
C
      IF (VISERRORNO.EQ.VISWARCROP) WRITE(VISSTDERRO,1)  1,VISSUBNAME
      IF (VISERRORNO.EQ.VISWARNAME) WRITE(VISSTDERRO,9)  9,VISSUBNAME
      IF (VISERRORNO.EQ.VISWARFIRS) WRITE(VISSTDERRO,11)  11,VISSUBNAME
      IF (VISERRORNO.EQ.VISWARLAST) WRITE(VISSTDERRO,12)  12,VISSUBNAME
      IF (VISERRORNO.EQ.VISWARFORM) WRITE(VISSTDERRO,13)  13,VISSUBNAME
C
C     *FORMATS
    1 FORMAT(//I2,1X,'***',A,'*WARNING: CROPSIZ LESS THAN 1. NO', 
     A                   ' ACTIONS PERFORMED.')
    9 FORMAT(//I2,1X,'***',A,'*WARNING: FILE NAME DOES NOT BEGIN WITH', 
     A                   ' A LETTER. STANDARD NAME HAS BEEN TAKEN.')
   11 FORMAT(//I2,1X,'***',A,'*WARNING: FIRST MUST SATISFY FIRST>=1.')
   12 FORMAT(//I2,1X,'***',A,'*WARNING: LAST  MUST SATISFY LAST',
     A                   ' <=51')
   13 FORMAT(//I2,1X,'***',A,'*WARNING: WRONG FORMAT, ACSII ASSUMED.')
      END
C
*DECK VISPRERR
      SUBROUTINE VISPRERR(VISERRORNO,VISSUBNAME)
C     ******************************************************************
C     * THIS ROUTINE PRINTS ERROR MESSAGES. VISERRORNO CONTAINS THE    *
C     * NUMBER OF THE MESSAGE. VISSUBNAME CONTAINS THE NAME OF THE     *
C     * CALLING SUBROUTINE.                                            *
C     ******************************************************************
C
#include "comerrn"
#include "cominou"
C
      CHARACTER*(*) VISSUBNAME
      INTEGER       VISERRORNO
C
      IF (VISERRORNO.EQ.VISERRAUXI) WRITE(VISSTDERRO,4)  4,VISSUBNAME
      IF (VISERRORNO.EQ.VISERRTRAN) WRITE(VISSTDERRO,5)  5,VISSUBNAME
      IF (VISERRORNO.EQ.VISERRNDIM) WRITE(VISSTDERRO,7)  7,VISSUBNAME
      IF (VISERRORNO.EQ.VISERRSERI) WRITE(VISSTDERRO,8)  8,VISSUBNAME
      IF (VISERRORNO.EQ.VISERRFIFO) WRITE(VISSTDERRO,10)  10,VISSUBNAME
      IF (VISERRORNO.EQ.VISERRTYPE) WRITE(VISSTDERRO,14)  14,VISSUBNAME
      IF (VISERRORNO.EQ.VISERRSTRU) WRITE(VISSTDERRO,15)  15,VISSUBNAME
      IF (VISERRORNO.EQ.VISERRWRST) WRITE(VISSTDERRO,21)  21,VISSUBNAME
      IF (VISERRORNO.EQ.VISERRPACK) WRITE(VISSTDERRO,22)  22,VISSUBNAME
      IF (VISERRORNO.EQ.VISERRCFLO) WRITE(VISSTDERRO,23)  23,VISSUBNAME
      IF (VISERRORNO.EQ.VISERRCXDR) WRITE(VISSTDERRO,24)  24,VISSUBNAME
      IF (VISERRORNO.EQ.VISERRIEEE) WRITE(VISSTDERRO,25)  25,VISSUBNAME
C
      WRITE(VISSTDERRO,*)
      WRITE(VISSTDERRO,'(A)') '***VISPRERR: TERMINATING DUE TO ERROR.'
C
      STOP
C
C     *FORMATS
    4 FORMAT(//I2,1X,'***',A,'*ERROR: COORDINATES FOR AUXILIARY',
     A                   ' GEOMETRY SHOULD ALWAYS BE DEFINED WITH',
     B                   ' CARTESIAN COORDINATES.')
    5 FORMAT(//I2,1X,'***',A,'*ERROR: TRANSFORMATION ONLY POSSIBLE FOR',
     A                   '  CYLINDER -> TOKAMAK AND',
     B                   '  CYLINDER -> CORONAL LOOP.')
    7 FORMAT(//I2,1X,'***',A,'*ERROR: NDIM SHOULD SATISFY:',
     A                   ' 1 <= NDIM <= 3.')
    8 FORMAT(//I2,1X,'***',A,'*ERROR: THE TIME DEPENDENCE CAN ONLY BE',
     A                   ' VISUALIZED BY MAKING VISSERIES .TRUE.',
     B                   ' IT MAY NOT APPEAR IN VISCOOR.')
   10 FORMAT(//I2,1X,'***',A,'*ERROR: WRONG FILE FORMAT INDICATION. IT',
     A                   ' SHOULD EQUAL VISGENERAL OR VISINTERNA.')
   14 FORMAT(//I2,1X,'***',A,'*ERROR: WRONG PRIMITIVE TYPE')
   15 FORMAT(//I2,1X,'***',A,'*ERROR: ALL THE DATA STRUCTURES AT A',
     A                   ' NODE SHOULD HAVE THE SAME PRIMITIVE DATA',
     B                   ' TYPE.')
   21 FORMAT(//I2,1X,'***',A,'*ERROR: LOVER ONLY SUPPORTS THE SCALAR',
     A                   ' FIELD STRUCTURE (VISRANKSCA) AND THE',
     B                   ' 1-VECTOR, 2-VECTOR, AND 3-VECTOR',/,
     B               11X,' FIELD STRUCTURE (VISRANKVC1,2,3).')
   22 FORMAT(//I2,1X,'***',A,'*ERROR: INPUT PARAMETER VISPACKAGE MUST',
     A                   ' EQUAL ''A'' FOR AVS OR ''D'' FOR DX.')
   23 FORMAT(//I2,1X,'***',A,'*ERROR: ONLY LOVER TYPE DOUBLE (= CRAY',
     A                   ' TYPE FLOAT) SUPPORTED.')
   24 FORMAT(//I2,1X,'***',A,'*ERROR: CURRENTLY, XDR ONLY SUPPORTED ON',
     A                   ' CRAY.')
   25 FORMAT(//I2,1X,'***',A,'*ERROR: OVERFLOWS DETECTED IN CALL TO',
     A                   ' CRAY2IEG.')
      END
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
#include "comtype"
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
#include "comtype"
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
#include "comerrn"
#include "comrank"
#include "comtype"
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
C
*DECK VISDIMVAR
      SUBROUTINE VISDIMVAR(VISNV, VISNI, VISNJ, VISNK, VISNT,VISVECTLEN, 
     A                     VISRESINIT,   VISRESFINA,   VISSUBNAME)
C     ******************************************************************
C     * CHECKING VARIABLES                                             *
C     ******************************************************************
C
#include "comcoor"
#include "cominou"
#include "comerrn"
C
      INTEGER VISNV, VISNI, VISNJ, VISNK, VISNT, VISVECTLEN,
     A        VISRESINIT(*), VISRESFINA(*)
      CHARACTER*(*) VISSUBNAME
C
      IF (VISNV.LT.VISVECTLEN) THEN
          WRITE(VISSTDERRO,1) 16,VISSUBNAME, VISNV
          STOP
      ENDIF
      IF (VISNI.LT.VISRESINIT(VISCOORRAD)) THEN
          WRITE(VISSTDERRO,2) 17,VISSUBNAME, VISNI
          STOP
      ENDIF
      IF (VISNJ.LT.VISRESINIT(VISCOORPOL)) THEN
          WRITE(VISSTDERRO,3) 18,VISSUBNAME, VISNJ
          STOP
      ENDIF
      IF (VISNK.LT.VISRESINIT(VISCOORTOR)) THEN
          WRITE(VISSTDERRO,4) 19,VISSUBNAME, VISNK
          STOP
      ENDIF
      IF (VISNT.LT.VISRESINIT(VISCOORTIM)) THEN
          WRITE(VISSTDERRO,5) 20,VISSUBNAME, VISNT
          STOP
      ENDIF
C
      IF (VISNI.LT.VISRESFINA(VISCOORRAD)) THEN
          WRITE(VISSTDERRO,6) 17,VISSUBNAME, VISNI
          STOP
      ENDIF
      IF (VISNJ.LT.VISRESFINA(VISCOORPOL)) THEN
          WRITE(VISSTDERRO,7) 18,VISSUBNAME, VISNJ
          STOP
      ENDIF
      IF (VISNK.LT.VISRESFINA(VISCOORTOR)) THEN
          WRITE(VISSTDERRO,8) 19,VISSUBNAME, VISNK
          STOP
      ENDIF
      IF (VISNT.LT.VISRESFINA(VISCOORTIM)) THEN
          WRITE(VISSTDERRO,9) 20,VISSUBNAME, VISNT
          STOP
      ENDIF
C
C     * FORMATS
    1 FORMAT(//I2,1X,'***',A,'*ERROR: NUMBER OF DATA ELEMENTS PER NODE',
     A                   ' SHOULD SATISFY VISVECTLEN <= ',I1)
    2 FORMAT(//I2,1X,'***',A,'*ERROR: VISRES(1) SHOULD BE <= ',I2)
    3 FORMAT(//I2,1X,'***',A,'*ERROR: VISRES(2) SHOULD BE <= ',I2)
    4 FORMAT(//I2,1X,'***',A,'*ERROR: VISRES(3) SHOULD BE <= ',I2)
    5 FORMAT(//I2,1X,'***',A,'*ERROR: VISRES(4) SHOULD BE <= ',I2)
    6 FORMAT(//I2,1X,'***',A,'*ERROR: VISRESFINA(1) SHOULD BE <= ',I2)
    7 FORMAT(//I2,1X,'***',A,'*ERROR: VISRESFINA(2) SHOULD BE <= ',I2)
    8 FORMAT(//I2,1X,'***',A,'*ERROR: VISRESFINA(3) SHOULD BE <= ',I2)
    9 FORMAT(//I2,1X,'***',A,'*ERROR: VISRESFINA(4) SHOULD BE <= ',I2)
      END
C
*DECK VISIGNOCOS
      FUNCTION VISIGNOCOS(VISPOSIRAD, VISPOSIPOL, VISPOSITOR, 
     A                    VISPOSITIM,
     B                    VISNI,      VISNJ,      VISNK,      VISNT,
     C                    IGNINDICATION, MODES, IRAD, IPOL, ITOR, ITIM)
C     ******************************************************************
C     * THIS FUNCTION CALCULATES THE ARGUMENT OF THE COSINUS. THIS     *
C     * COSINUS CONTAINS ALL THE INFORMATION OF THE IGNORABLE          *
C     * COORDINATES. IT IS ASSUMED THAT THE DEPENDENCE ON IGNORABLE    *
C     * COORDINATES IS FOURIER LIKE. FOR EACH OGNORABLE COORDINATE     *
C     * ONLY ONE FOURIER TERM IS TAKEN INTO ACCOUNT.                   *
C     *                                                                *
C     * E.G. COS(M*THETA+N*ZETA) IF THETA AND ZETA ARE IGNORABLE.      * 
C     *                                                                *
C     * ON INPUT:                                                      *
C     * IGNINDICATION INDICATES WHICH VARIABLES ARE IGNORABLE.         *
C     *               ITS VALUE IS OBTAINED BY SUMMATION.              *
C     *               ADD    1 IF RADIAL   COORDINATE IS IGNORABLE     *
C     *               ADD   10 IF POLOIDAL COORDINATE IS IGNORABLE     *
C     *               ADD  100 IF TOROIDAL COORDINATE IS IGNORABLE     *
C     *               ADD 1000 IF TIME     COORDINATE IS IGNORABLE     *
C     * MODES CONTAINS THE FOURIER MODE NUMBERS.                       *
C     * IRAD  CONTAINS THE RADIAL   INDEX FOR VISPOSIRAD.              *
C     * IPOL  CONTAINS THE POLOIDAL INDEX FOR VISPOSIRAD.              *
C     * ITOR  CONTAINS THE TOROIDAL INDEX FOR VISPOSIRAD.              *
C     * ITIM  CONTAINS THE TIME     INDEX FOR VISPOSIRAD.              *
C     *                                                                *
C     * VISIGNOCOS IS CALLED BY VISADDIGNO.                            *
C     ******************************************************************
C
#include "comcoor"
C
      INTEGER VISNI,      VISNJ,      VISNK,     VISNT,
     A        IGNINDICATION, 
     B        IRAD, IPOL, ITOR, ITIM
      REAL    VISPOSIRAD(VISNI),
     E        VISPOSIPOL(VISNJ),
     F        VISPOSITOR(VISNK),
     G        VISPOSITIM(VISNT)
      REAL    VISIGNOCOS, MODES(*)
C
C     * TEMP VARIABLES
      REAL    KR, M, N, OMEGA
      INTEGER I,J,K,L
C
C     * MODE NUMBERS
    1 KR   =MODES(VISCOORRAD)
      M    =MODES(VISCOORPOL)
      N    =MODES(VISCOORTOR)
      OMEGA=MODES(VISCOORTIM)
C
C     * WHICH COORDINATES ARE IGNORABLE?
      I=IGNINDICATION/1000
      J=MOD(IGNINDICATION,1000)/100
      K=MOD(MOD(IGNINDICATION,1000),100)/10
      L=MOD(MOD(MOD(IGNINDICATION,1000),100),10)
C
      VISIGNOCOS=    KR*VISPOSIRAD(IRAD)*L
     1          +    M*VISPOSIPOL(IPOL)*K
     2          +    N*VISPOSITOR(ITOR)*J
     3          -OMEGA*VISPOSITIM(ITIM)*I
      RETURN
C
      END
C
*DECK VISIGNOCOSI	
      FUNCTION VISIGNOCOSI(VISPOSX, VISPOSY, VISPOSZ, 
     A                    VISPOSITIM,
     B                    VISNI,      VISNJ,      VISNK,      VISNT,
     C                    IGNINDICATION, MODES, IRAD, IPOL, ITOR, ITIM)
C     ******************************************************************
C     * THIS FUNCTION CALCULATES THE ARGUMENT OF THE COSINUS. THIS     *
C     * COSINUS CONTAINS ALL THE INFORMATION OF THE IGNORABLE          *
C     * COORDINATES. IT IS ASSUMED THAT THE DEPENDENCE ON IGNORABLE    *
C     * COORDINATES IS FOURIER LIKE. FOR EACH OGNORABLE COORDINATE     *
C     * ONLY ONE FOURIER TERM IS TAKEN INTO ACCOUNT.                   *
C     *                                                                *
C     * E.G. COS(M*THETA+N*ZETA) IF THETA AND ZETA ARE IGNORABLE.      * 
C     *                                                                *
C     * ON INPUT:                                                      *
C     * IGNINDICATION INDICATES WHICH VARIABLES ARE IGNORABLE.         *
C     *               ITS VALUE IS OBTAINED BY SUMMATION.              *
C     *               ADD    1 IF RADIAL   COORDINATE IS IGNORABLE     *
C     *               ADD   10 IF POLOIDAL COORDINATE IS IGNORABLE     *
C     *               ADD  100 IF TOROIDAL COORDINATE IS IGNORABLE     *
C     *               ADD 1000 IF TIME     COORDINATE IS IGNORABLE     *
C     * MODES CONTAINS THE FOURIER MODE NUMBERS.                       *
C     * IRAD  CONTAINS THE RADIAL   INDEX FOR VISPOSIRAD.              *
C     * IPOL  CONTAINS THE POLOIDAL INDEX FOR VISPOSIRAD.              *
C     * ITOR  CONTAINS THE TOROIDAL INDEX FOR VISPOSIRAD.              *
C     * ITIM  CONTAINS THE TIME     INDEX FOR VISPOSIRAD.              *
C     *                                                                *
C     * VISIGNOCOSI IS CALLED BY VISADDIGNOI.                            *
C     ******************************************************************
C
#include "comcoor"
C
      INTEGER VISNI,      VISNJ,      VISNK,     VISNT,
     A        IGNINDICATION, 
     B        IRAD, IPOL, ITOR, ITIM
      REAL    VISPOSX(VISNI,VISNJ,VISNK),
     B        VISPOSY(VISNI,VISNJ,VISNK),
     C        VISPOSZ(VISNI,VISNJ,VISNK),
     G        VISPOSITIM(VISNT)
      REAL    VISIGNOCOSI, MODES(*)
C
C     * TEMP VARIABLES
      REAL    KR, M, N, OMEGA
      INTEGER I,J,K,L
C
C     * MODE NUMBERS
    1 KR   =MODES(VISCOORRAD)
      M    =MODES(VISCOORPOL)
      N    =MODES(VISCOORTOR)
      OMEGA=MODES(VISCOORTIM)
C
C     * WHICH COORDINATES ARE IGNORABLE?
      I=IGNINDICATION/1000
      J=MOD(IGNINDICATION,1000)/100
      K=MOD(MOD(IGNINDICATION,1000),100)/10
      L=MOD(MOD(MOD(IGNINDICATION,1000),100),10)
C
      VISIGNOCOSI=    KR*VISPOSX(IRAD,IPOL,ITOR)*L
     1           +    M*VISPOSY(IRAD,IPOL,ITOR)*K
     2           +    N*VISPOSZ(IRAD,IPOL,ITOR)*J
     3           -OMEGA*VISPOSITIM(ITIM)*I
      RETURN
C
      END
