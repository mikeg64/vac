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
C#include "comgeom"
C#include "comcoor"
C#include "comerrn"
C#include "comcons"
C#include "cominou"


	INCLUDE 'INCLUDE/comgeom'
	INCLUDE 'INCLUDE/comcoor'
	INCLUDE 'INCLUDE/comerrn'
	INCLUDE 'INCLUDE/comcons'
	INCLUDE 'INCLUDE/cominou'
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
