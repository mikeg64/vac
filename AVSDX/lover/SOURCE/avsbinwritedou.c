/*@(#)avsbinwritedou.c	2.1
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
void AVSBINWRITEDOU(_fcd pathname, double object[], size_t *nobj,
                    long *offset, int *xdr)
#else
#ifdef IBMRIOS
void avsbinwritedou(char *pathname, double object[], size_t *nobj,
                    long *offset, int *xdr, int   pathnamelength)
#else 
#ifdef DEC
void avsbinwritedou_(char *pathname, double object[], int *nobj,
                    int *offset, int *xdr, int   pathnamelength)
#else
void avsbinwritedou_(char *pathname, double object[], size_t *nobj,
                    long *offset, int *xdr, int   pathnamelength)
#endif
#endif
#endif
{

FILE *fp;

#ifdef CRAY
   if ((fp = fopen(_fcdtocp(pathname), "r+")) != NULL)
#else
   if ((fp = fopen(pathname, "r+")) != NULL)
#endif
     {
      if (fseek(fp,*offset,SEEK_SET) == -1)
	{
         fprintf(stderr,"AVSBINWRITEDOU*ERROR: FSEEK ERROR\n");
         exit(1);
	}

      if (xdr)
	{
	  if (fwrite(object, 8, *nobj, fp) != *nobj)
	    {
	      fprintf(stderr,"AVSBINWRITEDOU*ERROR: WRITE ERROR\n");
	      exit(1);
	    }
	}
      else
	{
	  if (fwrite(object, sizeof(object[0]), *nobj, fp) != *nobj)
	    {
	      fprintf(stderr,"AVSBINWRITEDOU*ERROR: WRITE ERROR\n");
	      exit(1);
	    }
	}
     }
   else
     {
       fprintf(stderr, "AVSBINWRITEDOU*ERROR: OPEN ERROR\n");
       exit(1);
     }

   fclose(fp);
}
