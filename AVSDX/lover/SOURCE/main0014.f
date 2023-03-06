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
