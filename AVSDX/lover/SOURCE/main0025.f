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
C#include "comcoor"
C#include "cominou"
C#include "comrank"
C#include "comform"
C#include "comtype"

	INCLUDE 'INCLUDE/comcoor'
	INCLUDE 'INCLUDE/cominou'
	INCLUDE 'INCLUDE/comrank'
	INCLUDE 'INCLUDE/comform'
	INCLUDE 'INCLUDE/comtype'


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
C#ifdef CRAY
C      INTEGER      IERR
C#endif
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
C#ifdef CRAY
C      INTEGER CRAY2IEG
C#endif
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
C#ifdef CRAY
C            CALL VISPRERR(VISERRCFLO,'VISDATAOUT')
C#else
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
C#endif
         ELSE IF (VISSTRUTYP(1).EQ.VISTYPEDOU) THEN
            IF (VISXDR) THEN

               CALL VISPRERR(VISERRCXDR,'VISDATAOUT')
            ELSE
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
