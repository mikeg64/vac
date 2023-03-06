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
C#include "comgeom"
C#include "comcoor"
C#include "comcons"
C#include "comrank"
C#include "comerrn"

	INCLUDE 'INCLUDE/comgeom'
	INCLUDE 'INCLUDE/comcoor'
	INCLUDE 'INCLUDE/comcons'
	INCLUDE 'INCLUDE/comrank'
	INCLUDE 'INCLUDE/comerrn'


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
