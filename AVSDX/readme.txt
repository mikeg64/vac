MKG mod 5/11/2009
Instructions for building lover library for convertdata utility
change directory to lover/SOURCE

use the same compiler as used with the compilation of vac and vacini to build convertdata

to build the libraries
use the command
make -f liblover1.m all

it may be necessary to alter the compilers and switches in this mk file
copy lover.h to the include folder

Modifications by MKG 26th Oct 2007
comorde  MSB was defined as the default expressing byte order here
comcarr  removed defines set the default CR=0 called by main0019.f


@(#)README	2.1
This directory contains the source programs constituing lover2.1+.

'liblover.mk ARCH ' (with ARCH one of CRAY, IBMRIOS, SGI, SUN) 
compiles the whole tree. It is automatically called by the 
installation script. To invoke this script yourself, make sure 
the environment variables are set to their correct values.

	LOVER        top directory name under which the tree can 
                     be found (defaults to /usr/local/lover)
	LIB          directory wher liblover.a is placed 
                     (defaults to /usr/local/lib)
	INCLUDE      directory in which the user include file 
                     lover.h is placed (defaults to /usr/local/include)

	F77          the name of the FORTRAN compiler 
                     (defaults to f90 on CRAY, xlf on IBMRIOS, 
                               to f77 on SGI and SUN)
	F77FLAGS     flags to be passed to the FORTRAN compiler 
                     (defaults to -O on CRAY,
                      defaults to -O -qautodbl=dbl4 on IBM,
                      defaults to -O -r8 on SGI and SUN) 
	CC           the name of the C compiler 
                     (defaults to cc on CRAY, IBMRIOS, and SGI,
                      defaults to acc on SUN)
	CFLAGS       flags to be passed to the C compiler (defaults to -O)

Notice well, on workstations the source has to be compiled with the 
autodouble option (passed as a flag with F77FLAGS). Your local 
compiler guide will tell you what the autodouble option is.

