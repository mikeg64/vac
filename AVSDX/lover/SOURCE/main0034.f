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
C#include "comcoor"

	INCLUDE 'INCLUDE/comcoor'

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
