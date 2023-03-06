C
*DECK VISNAMEOK
      FUNCTION VISNAMEOK(NAME)
C     ******************************************************************
C     * FUNCTION TO CHECK IF A CHARACTER STRING STARTS WITH A LETTER.  *
C     * IF THIS IS TRUE THEN .TRUE. WILL BE RETURNED, OTHERWISE .FALSE.*
C     *                                                                *
C     * ON INPUT:                                                      *
C     * NAME     THE CHARACTER STRING                                  *
C     ******************************************************************
C
      CHARACTER*(*) NAME
      LOGICAL       VISNAMEOK
C
      VISNAMEOK=.TRUE.
      IF (.NOT.((NAME(1:1).GE.'A'.AND.NAME(1:1).LE.'Z')
     A    .OR.(NAME(1:1).GE.'a'.AND.NAME(1:1).LE.'z'))) THEN
         VISNAMEOK=.FALSE.
      ENDIF
      END
