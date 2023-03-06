C
*DECK VISGETRES
C     ==================================================================
C     = FUNCTIONS AND SUBROUTINES WHICH ARE CALLED BY THE ABOVE        =
C     = SUBROUTINES. THEY SHOULD NOT BE CALLED BY THE USER. THEY MAY   =
C     = FURTHERMORE CONTAIN HARDWARE DEPENDENT INFORMATION!!!!!!       =
C     =                              |                                 =
C     =                              |                                 =
C     =                              |                                 =
C     =                            \   /                               =
C     =                             \ /                                =
C     =                              .                                 =
C     ==================================================================
      FUNCTION VISGETRES(VISRES, VISNDIM, VISCOOR, WHICHCOOR)
C     ******************************************************************
C     * THIS FUNCTION RETURNS THE RESOLUTION BELONGING TO WHICHCOOR    *
C     * IF WHICHCOOR HAS TO BE OUTPUTTED. OTHERWISE A VALUE OF 1 IS    *
C     * RETURNED.                                                      *
C     *                                                                *
C     * ON INPUT:                                                      *
C     * VISRES     CONTAINS THE RESOLUTION OF THE DATA                 *
C     * VISNDIM    CONTAINS THE NUMBER OF DIMENSION TO BE PRINTED      *
C     * VISCOOR    CONTAINS THE DIMENSIONS WHICH WILL BE PRINTED       *
C     ******************************************************************
C
      INTEGER VISGETRES, VISRES(*), VISCOOR(*), VISNDIM, WHICHCOOR
C
C     * COUNTERS
      INTEGER I
C
      DO 10 I=1, VISNDIM
         IF (VISCOOR(I).EQ.WHICHCOOR) THEN
            VISGETRES = VISRES(WHICHCOOR)
            RETURN
         ENDIF
   10 CONTINUE
      VISGETRES = 1
      END
