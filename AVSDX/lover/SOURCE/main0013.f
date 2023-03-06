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
C#include "comcoor"
C#include "comerrn"
	INCLUDE 'INCLUDE/comcoor'
	INCLUDE 'INCLUDE/comerrn'
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
