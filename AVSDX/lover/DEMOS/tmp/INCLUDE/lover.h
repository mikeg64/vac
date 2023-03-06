************************************************************************
* INLCUDE file: @(#)lover.h	2.1
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
      INTEGER VISGEOMCYL, VISGEOMTOK, VISGEOMLOO, VISGEOMAUX 
      PARAMETER (
     1   VISGEOMCYL = 1,
     2   VISGEOMTOK = 2,
     3   VISGEOMLOO = 3,
     4   VISGEOMAUX = 4
     5)
      INTEGER VISCOORRAD, VISCOORPOL, VISCOORTOR, VISCOORTIM
      PARAMETER (
     1   VISCOORRAD = 1,
     2   VISCOORPOL = 2,
     3   VISCOORTOR = 3,
     4   VISCOORTIM = 4
     5)
      INTEGER VISRANKSCA, VISRANKVC1, VISRANKVC2, VISRANKVC3, VISRANKVEC
      PARAMETER (
     1   VISRANKSCA = 1,
     2   VISRANKVC1 = 1,
     3   VISRANKVC2 = 2,
     4   VISRANKVC3 = 3,
     5   VISRANKVEC = 3
     6)
      REAL    PI
      PARAMETER (
     1   PI         = 3.1415926535898
     2)
      INTEGER VISFORMASC, VISFORMBIN, VISFORMUNF, VISGENERAL, 
     A        VISINTERNA, VISFORMXDR
      PARAMETER (
     1     VISFORMASC=1,
     2     VISFORMBIN=2,
     3     VISFORMUNF=3,
     4     VISFORMXDR=4,
     5     VISGENERAL=6,
     6     VISINTERNA=7
     7)
      CHARACTER VISTYPEBYT*(*), VISTYPECHA*(*), VISTYPESHO*(*), 
     A          VISTYPEINT*(*), VISTYPEFLO*(*), VISTYPEDOU*(*)
      PARAMETER (
     1     VISTYPEBYT = 'byte', 
     2     VISTYPECHA = 'char', 
     3     VISTYPESHO = 'short',
     4     VISTYPEINT = 'int',
     5     VISTYPEFLO = 'float',
     6     VISTYPEDOU = 'double'
     7)
