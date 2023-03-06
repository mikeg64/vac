/*@(#)dxbinwriteflo.c	2.1
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
void DXBINWRITEFLO(_fcd pathname, float object[], size_t *nobj, 
                   size_t *size)
#else
#ifdef IBMRIOS
void dxbinwriteflo(char *pathname, float object[], size_t *nobj, 
                   size_t *size, int pathnamelength)
#else 
#ifdef DEC
void dxbinwriteflo_(char *pathname, float object[], int *nobj, 
                   int *size, int pathnamelength)
#else
void dxbinwriteflo_(char *pathname, float object[], size_t *nobj, 
                   size_t *size, int pathnamelength)
#endif
#endif
#endif
{

FILE *fp;

#ifdef CRAY
 if ((fp = fopen(_fcdtocp(pathname), "ab")) != NULL)
#else   
 if ((fp = fopen(pathname, "ab")) != NULL)
#endif
   {
     if (fwrite(object, *size, *nobj, fp) != *nobj)
       {
         fprintf(stderr,"DXBINWRITEFLO*ERROR: WRITE ERROR\n");
         exit(1);
       }
   }
 else
   {
     fprintf(stderr, "DXBINWRITEFLO*ERROR: OPEN ERROR\n");
     exit(1);
   }

 fclose(fp);
}
