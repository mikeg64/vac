C
*DECK VISADDIGNO10
      SUBROUTINE VISADDIGNO10(VISMATDATA, VISPOSIRAD, VISPOSIPOL,
     A                       VISPOSITOR, VISPOSITIM, VISRESINIT,
     B                       VISSTRUCTU, MODES,      FINTIM, VISCOS)
C     ******************************************************************
C     * ARGUMENT DESCRIPTION: SEE VISADDIGNO                           *
C     ******************************************************************
      COMMON /DIMS/ VISDIM1, VISDIM2, VISDIM3, VISDIM4, VISDIM5
      INTEGER       VISDIM1, VISDIM2, VISDIM3, VISDIM4, VISDIM5
      COMMON /INFO/ VISRESO,    VISCORD,    VISSLIC,    VISDIMS,
     A              VISSERI,    VISVLEN,    VISNFLD,    VISFORM,
     B              VISGEOM,    VISXDR
      LOGICAL       VISSERI,    VISXDR
      INTEGER       VISRESO(4), VISCORD(4), VISSLIC(4), VISDIMS,
     A              VISVLEN,    VISNFLD,    VISFORM,    VISGEOM
      COMMON /CTIM/ VISTIME1,   VISTRES,    VISNTST,
     A              VISFIRS,    VISLAST,    VIS1TST
      LOGICAL       VISFIRS,    VISLAST,    VIS1TST
      INTEGER       VISTRES,    VISNTST
      REAL          VISTIME1
C
      INTEGER       VISRESINIT(*), VISSTRUCTU(*)
      REAL          VISMATDATA(VISDIM1,VISDIM2,VISDIM3,VISDIM4,VISDIM5),
     A              VISPOSIRAD(VISDIM2),
     B              VISPOSIPOL(VISDIM3),
     C              VISPOSITOR(VISDIM4),
     D              VISPOSITIM(VISDIM5),
     E              MODES(*), FINTIM
      LOGICAL       VISCOS
C
C#include "comerrn"
C#include "comrank"
	INCLUDE 'INCLUDE/comerrn'
	INCLUDE 'INCLUDE/comrank'
C
C     * LOCAL VARIABLE
      INTEGER VISTMP
C
C     * COUNTERS
      INTEGER I
C
C     * EXTERNAL SUBROUTINES
      EXTERNAL VISPRERR, VISTOCART
C
C     * DETERMINING THE VECTOR LENGTH
      VISVLEN = 0
      DO 10 I=1, VISNFLD
         IF (VISSTRUCTU(I).EQ.VISRANKSCA) THEN
            VISVLEN = VISVLEN + 1
         ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC1) THEN
            VISVLEN = VISVLEN + 1
         ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC2) THEN
            VISVLEN = VISVLEN + 2
         ELSE IF (VISSTRUCTU(I).EQ.VISRANKVC3) THEN
            VISVLEN = VISVLEN + 3
         ELSE IF (VISSTRUCTU(I).EQ.VISRANKVEC) THEN
            VISVLEN = VISVLEN + 3
         ELSE
            CALL VISPRERR(VISERRWRST,'VISADDIGNO9')
         ENDIF
   10 CONTINUE
C
      IF (VIS1TST) THEN
         VISTMP     = VISRESO(4)
         VISRESO(4) = 1
      ENDIF
C
      CALL VISADDIGNO(VISMATDATA,
     A               VISPOSIRAD, VISPOSIPOL, VISPOSITOR,
     B               VISPOSITIM,
     C               VISDIM1,       VISDIM2,    VISDIM3,    VISDIM4,
     D               VISDIM5,    VISRESINIT,    VISRESO,    VISVLEN,
     E               MODES,          FINTIM,    VISCOS)
C
      IF (VIS1TST) THEN
         VISRESO(4) = VISTMP
      ENDIF
C
      RETURN
      END
