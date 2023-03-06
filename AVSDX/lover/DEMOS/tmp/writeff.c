/*@(#)writeff.c	2.1
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
void WRITEFF(_fcd pathname)
#else
#ifdef IBMRIOS
void writeff(char * pathname, int   pathnamelength)
#else
void writeff_(char * pathname, int   pathnamelength)
#endif
#endif
{

char LF;
FILE *fp;

   LF = 12;

#ifdef CRAY
   if ((fp = fopen(_fcdtocp(pathname), "ab")) != NULL)
#else
   if ((fp = fopen(pathname, "ab")) != NULL)
#endif
     {
       if (fwrite(&LF, sizeof(LF), 1, fp) != 1)
         fprintf(stderr,"WRITEFF: WRITE ERROR\n");
       if (fwrite(&LF, sizeof(LF), 1, fp) != 1)
         fprintf(stderr,"WRITEFF: WRITE ERROR\n");
     }
   else

     fprintf(stderr, "WRITEFF: OPEN ERROR\n");

   fclose(fp);
}
