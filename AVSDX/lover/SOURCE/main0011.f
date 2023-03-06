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
C#include "comcoor"
C#include "comerrn"
	INCLUDE 'INCLUDE/comcoor'
	INCLUDE 'INCLUDE/comerrn'

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
