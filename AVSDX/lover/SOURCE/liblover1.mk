#! /bin/sh
#@(#)liblover.mk	2.1
#
#Shell script to construct the LOVER library.
#
#Author:  A.J.C. Belien
#Version: 2.1
#Date:    14:04:52 00/03/14
#
#(C) 1994-1996: FOM-Institute for Plasma Physics `Rijnhuizen',
#               Academic Computing Services Amsterdam (SARA)
#(C) 1997-2000: FOM-Institute for Plasma Physics `Rijnhuizen'
#
CC= pgcc
F77= pgf77
LIB= ../../lib/.
F77FLAGS=
CFLAGS=
 


objs :
	$(F77) *.f -c $(F77FLAGS) -I./INCLUDE
	$(CC)  *.c -c $(CFLAGS)


#.f.o:
#	$(F77) $(F77FLAGS) -I./INCLUDE -o $@ -c $<



#Compile
#all : $(OBJ)
#	ar -rv liblover.a $(OBJ)

#echo =============================
#echo Creating library liblover.a        
#echo in directory $LIB          
#echo =============================

all : objs
	ar vq liblover.a *.o
	mv liblover.a $(LIB)
#

#
#clean :


	

