/*@(#)dxbinwritedou.c	2.1
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
void DXBINWRITEDOU(_fcd pathname, double object[], size_t *nobj, 
                   size_t *size)
#else
#ifdef IBMRIOS
void dxbinwritedou(char *pathname, double object[], size_t *nobj, 
                   size_t *size, int pathnamelength)
#else 
#ifdef DEC
void dxbinwritedou_(char *pathname, double object[], int *nobj, 
                   int *size, int pathnamelength)
#else
void dxbinwritedou_(char *pathname, double object[], size_t *nobj, 
                   size_t *size, int pathnamelength)
#endif
#endif
#endif
{

FILE *fp;
#ifndef CRAY
float *fl_object;
int   i;
#endif

#ifndef CRAY
 if (*size == sizeof(float))
   {
     fl_object = (float *) calloc(*nobj,*size);
     
     for (i=0; i<*nobj;i++)
       fl_object[i] = (float) object[i];
   }
#endif

#ifdef CRAY
 if ((fp = fopen(_fcdtocp(pathname), "ab")) != NULL)
#else
 if ((fp = fopen(pathname, "ab")) != NULL)
#endif
   {
#ifndef CRAY
     if (*size == sizeof(float))
       {
	 if (fwrite(fl_object, *size, *nobj, fp) != *nobj)
	   {
	     fprintf(stderr,"DXBINWRITEDOU*ERROR: WRITE ERROR\n");
	     exit(1);
	   }
       }
     else if (fwrite(object, *size, *nobj, fp) != *nobj)
       {
	 {
	   fprintf(stderr,"DXBINWRITEDOU*ERROR: WRITE ERROR\n");
	   exit(1);
	 }
       }
#else
     if (fwrite(object, *size, *nobj, fp) != *nobj)
       {
	 fprintf(stderr,"DXBINWRITEDOU*ERROR: WRITE ERROR\n");
	 exit(1);
       }
#endif
   }
 else
   {
     fprintf(stderr, "DXBINWRITEDOU*ERROR: OPEN ERROR\n");
     exit(1);
   }

#ifndef CRAY
 if (*size == sizeof(float))
   free(fl_object);
#endif
 fclose(fp);
}
