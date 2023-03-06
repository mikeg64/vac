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
C#include "comcoor"
C#include "cominou"
C#include "comform"
	INCLUDE 'INCLUDE/comcoor'
	INCLUDE 'INCLUDE/cominou'
	INCLUDE 'INCLUDE/comform'


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
C#ifdef CRAY
C      INTEGER      IERR
C#endif
      CHARACTER    STRING*100, FILENAME*80
C
C     *LOCAL VARIABLES ASCII
      INTEGER      NOBJ, IND
C
C     * FUNCTIONS
      INTEGER VISGETRES
C#ifdef CRAY
C      INTEGER CRAY2IEG
C#endif
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
C#ifdef CRAY
C         IERR=CRAY2IEG(8,NOBJ,VISWORK,0,VISWORK)
C         IF (IERR.GT.0) CALL VISPRERR(VISERRIEEE,'VISDXWRDAT')
C         CALL DXBINWRITEFLO(FILENAME,VISWORK,NOBJ,8)
C#else
         CALL DXBINWRITEDOU(FILENAME,VISWORK,NOBJ,8)
C#endif
      ENDIF
      END
