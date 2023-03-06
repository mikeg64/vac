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
C#include "comcoor"
C#include "comerrn"
	INCLUDE 'INCLUDE/comcoor'
	INCLUDE 'INCLUDE/comerrn'
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
