*DECK COPYRIGHT
************************************************************************
* MAIN program: @(#)liblover.F   2.1
*                                                                      *
* Author: A.J.C. Belien                                                *
*         FOM-Instituut voor Plasmafysica `Rijnhuizen'                 *
*         P.O. Box 1207                                                *
*         3430 BE Nieuwegein                                           *
*                                                                      *
* (C) 1994-1996: FOM-Institute for Plasma Physics `Rijnhuizen',        *
*                Academic Computing Services Amsterdam (SARA)          *
* (C) 1997-2000: FOM-Institute for Plasma Physics `Rijnhuizen'         *
************************************************************************
*DECK VISLOVER
      SUBROUTINE VISLOVER
C     ******************************************************************
C     ******************************************************************
C
	INCLUDE 'INCLUDE/cominou'
C#include "cominou"
C
      LOGICAL VISPRINTCOPYRIGHT
      SAVE VISPRINTCOPYRIGHT
      DATA VISPRINTCOPYRIGHT/.TRUE./

      IF (VISPRINTCOPYRIGHT) THEN
         WRITE(VISSTDOUTP,100)
         VISPRINTCOPYRIGHT=.FALSE.
      END IF
C
C     * FORMATS
  100 FORMAT(/
     C/'**************************************************************',
     C/'*                                                            *',
     C/'*         Library Of Visualization Exporting Routines        *',
     C/'*                                                            *',
     C/'*                             BY                             *',
     C/'*                                                            *',
     C/'*                       SANDER BELIEN                        *',
     C/'*                                                            *',
     C/'*  (C) 1994-1996:                                            *',
     C/'*            FOM-Institute for Plasma Physics Rijnhuizen,    *',
     C/'*            Academic Computing Services Amsterdam (SARA)    *',
     C/'*  (C) 1997-2000:                                            *',
     C/'*            FOM-Institute for Plasma Physics Rijnhuizen     *',
     C/'**************************************************************')
      END
