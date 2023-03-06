/*@(#)findoffset.c	2.1
 *                                                                      
 * Author: A.J.C. Belien                                                
 *         FOM-Instituut voor Plasmafysica `Rijnhuizen'                 
 *         P.O. Box 1207                                                
 *         3430 BE Nieuwegein                                           
 *                                                                      
 * (C) 1994-1996: FOM-Institute for Plasma Physics `Rijnhuizen',        
 *                Academic Computing Services Amsterdam (SARA)          
 * (C) 1997-2000: FOM-Institute for Plasma Physics `Rijnhuizen'         
 */
#include <stdio.h>

#ifdef CRAY
#include <fortran.h>
long FINDOFFSET(_fcd pathname)
#else
#ifdef IBMRIOS
long findoffset(char * pathname, int   pathnamelength)
#else
long findoffset_(char * pathname, int   pathnamelength)
#endif
#endif
{

FILE *fp;
long currpos;

#ifdef CRAY
   if ((fp = fopen(_fcdtocp(pathname), "ab")) != NULL)
#else
   if ((fp = fopen(pathname, "ab")) != NULL)
#endif
     {
       if ((currpos = ftell(fp)) == -1)
         fprintf(stderr,"FINDOFFSET: FTELL ERROR\n");
       
     }
   else
     
     fprintf(stderr, "FINDOFFSET: OPEN ERROR\n");

   fclose(fp);

   return currpos;
}
