C
*DECK VISOPENAVS
      FUNCTION VISOPENAVS(VISFILNAME, VISFORMATS)
C     ******************************************************************
C     * THIS FUNCTION IS USED TO OPEN THE DATA FILE FOR AVS'S GENERAL  *
C     * FILE FORMAT. IT IS CALLED BY THE USER AFTER THE HEADER HAS     *
C     * BEEN WRITTEN. THE FUNCTIONS RETURNS THE UNIT NUMBER ON WHICH   *
C     * THE DATA FILE HAS BEEN OPENED FOR WRITING.                     *
C     *                                                                *
C     * ON INPUT:                                                      *
C     * VISFILNAME CONTAINS THE NAME OF THE DATA FILE TO BE OPENED.    *
C     * VISFORMATS CONTAINS THE FORM IN WHICH DATA WILL BE WRITTEN     *
C     *    VISFORMATS=VISFORMASC -> ASCII,                             *
C     *    VISFORMATS=VISFORMBIN -> BINARY (NOT IMPLEMENTED FOR        *
C     *                                     GENERAL FILE FORMAT),      *
C     *    VISFORMATS=VISFORMUNF -> FORTRAN UNFORMATTED.               *
C     *                                                                *
C     * CALLED BY: USER PROGRAM                                        *
C     * TYPICAL USE:                                                   *
C     *         ...                                                    *
C     *         call visavshead(...)                                   *
C     *         visoutput=VISOPENAVS(...)                              *
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
      INTEGER       VISOPENAVS, VISFORMATS
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
         VISFILNAME='avsfile'
      ENDIF
C
      IF (VISFORMATS.NE.VISFORMBIN) THEN
         IND=INDEX(VISFILNAME,'.fld')
         IF (IND.EQ.0) THEN
            IND               = INDEX(VISFILNAME,' ')
            IF (IND.EQ.0) IND = LEN(VISFILNAME)+1
         ENDIF
         IF (VISFORMATS.EQ.VISFORMASC) THEN
            OPEN(UNIT=VISAVSOUTP, FILE=VISFILNAME(1:IND-1),
     A           STATUS='UNKNOWN')
         ELSE IF (VISFORMATS.EQ.VISFORMUNF) THEN
            OPEN(UNIT=VISAVSOUTP, FORM='UNFORMATTED',
     A           FILE=VISFILNAME(1:IND-1),STATUS='UNKNOWN')
         ENDIF
      ENDIF
C
      VISOPENAVS=VISAVSOUTP
C
      END
