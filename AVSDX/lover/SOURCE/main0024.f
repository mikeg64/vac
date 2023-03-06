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
C#include "comcoor"
C#include "cominou"
C#include "comform"

	INCLUDE 'INCLUDE/comcoor'
	INCLUDE 'INCLUDE/cominou'
	INCLUDE 'INCLUDE/comform'


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
C#ifdef CRAY
C      REAL    TMP
C      INTEGER IERR
C#endif
      INTEGER NOBJ, IND
C
C     * FUNCTIONS
      INTEGER VISGETRES
C#ifdef CRAY
C      INTEGER CRAY2IEG
C#endif
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
C#ifdef CRAY
C        * SOMETHING GOES WRONG IF WE WRITE ALL NOBJ ELEMENTS
C          AT ONCE. IT IS CURRENTLY FIXED BY WRITING TWICE A HALF
C          OF THE DATA SET.
C         IERR = CRAY2IEG(2,NOBJ/2,VISWORK,0,VISWORK)
C         IF (IERR.GT.0) CALL VISPRERR(VISERRIEEE,'VISDXWRCOO')
C         CALL DXBINWRITEFLO(FILENAME,VISWORK,NOBJ/2,4)
C         IERR = CRAY2IEG(2,NOBJ-NOBJ/2,VISWORK,0,VISWORK(NOBJ/2+1))
C         IF (IERR.GT.0) CALL VISPRERR(VISERRIEEE,'VISDXWRCOO')
C         CALL DXBINWRITEFLO(FILENAME,VISWORK,NOBJ-NOBJ/2,4)
C#else
         CALL DXBINWRITEDOU(FILENAME,VISWORK,NOBJ,4)
C#endif
      ENDIF
      END
