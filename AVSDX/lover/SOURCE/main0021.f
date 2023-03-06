C
*DECK VISOPENDX
      FUNCTION VISOPENDX(VISFILNAME, VISFORMATS)
C     ******************************************************************
C     * THIS FUNCTION IS USED TO OPEN THE DATA FILE FOR DX'S GENERAL   *
C     * FILE FORMAT. IT IS CALLED BY THE USER AFTER THE HEADER HAS     *
C     * BEEN WRITTEN. THE FUNCTIONS RETURNS THE UNIT NUMBER ON WHICH   *
C     * THE DATA FILE HAS BEEN OPENED FOR WRITING.                     *
C     *                                                                *
C     * ON INPUT:                                                      *
C     * VISFLNAME CONTAINS THE NAME OF THE DATA FILE TO BE OPENED.     *
C     * VISFORMATS CONTAINS THE FORM IN WHICH DATA WILL BE WRITTEN     *
C     *    VISFORMATS=VISFORMASC -> ASCII,                             *
C     *    VISFORMATS=VISFORMBIN -> BINARY,                            *
C     *    VISFORMATS=VISFORMUNF -> NOT SUPPORTED IN DX.               *
C     *                                                                *
C     * CALLED BY: USER PROGRAM                                        *
C     * TYPICAL USE:                                                   *
C     *         ...                                                    *
C     *         call visdxheadgen(...)                                 *
C     *         visoutput=VISOPENDX(...)                               *
C     *         ...                                                    *
C     *         do i=1,number_of_timesteps                             *
C     *            call visdataout(...)                                *
C     *         end do                                                 *
C     ******************************************************************
C
C#include "cominou"
C#include "comerrn"
C#include "comform"



	INCLUDE 'INCLUDE/cominou'
	INCLUDE 'INCLUDE/comerrn'
	INCLUDE 'INCLUDE/comform'
C
      INTEGER       VISOPENDX, VISFORMATS
      CHARACTER*(*) VISFILNAME
C
C     * LOCAL VARIABLES
      INTEGER IND
C
C     * FUNCTIONS
      LOGICAL VISNAMEOK
C
C     * FILE NAME OK?
      IF (.NOT.VISNAMEOK(VISFILNAME)) THEN
         CALL VISPRWAR(VISWARNAME,'VISDXINT')
         VISFILNAME='dxfile'
      ENDIF
C
      IND=INDEX(VISFILNAME,'.gai')
      IF (IND.EQ.0) THEN
         IND               = INDEX(VISFILNAME,' ')
         IF (IND.EQ.0) IND = LEN(VISFILNAME)+1
      ENDIF
      OPEN(UNIT=VISDXOUTPU, FILE=VISFILNAME(1:IND-1),STATUS='UNKNOWN')
C
      VISOPENDX=VISDXOUTPU
C
      END
