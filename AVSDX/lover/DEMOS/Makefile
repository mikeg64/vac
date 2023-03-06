#@(#)Makefile	2.1
#Makefile to construct the LOVER demos. 
#
#Author:  A.J.C. Belien
#Version: 2.1
#Date:    17:02:46 00/03/13
#
#(C) 1994-1996: FOM-Institute for Plasma Physics `Rijnhuizen',
#               Academic Computing Services Amsterdam (SARA)
#(C) 1997-2000: FOM-Institute for Plasma Physics `Rijnhuizen'
#
# IRIX uses inherited shell to execute commands. 
# Hence, we explicitly set it to the Bourne shell.
SHELL    =/bin/sh
#
LIB      =../../lib
INCLUDE  =../../include
INCLUDES =-I$(INCLUDE)
LIBS     =-L$(LIB) -llover
LIBOBJ   =$(LIB)/liblover.a
#LIBSRC   =SOURCE/liblover.F
#LIBMK    =SOURCE/liblover1.mk
#
DEMOS    = slice   \
           format1 \
           format2 \
           format3 \
           format4 \
           timesteps  \
           timesteps2 \
           timesteps3 
#
all:
#	make $(LIBOBJ)
	make $(DEMOS)
#
$(DEMOS): $$@.F
	g77 -o $@ $@.F $(F77FLAGS) $(INCLUDES) $(LIBS) \

#
#$(LIBOBJ): $(LIBSRC)
#	$(LIBMK) $(ARCH)
#
clean:
	-rm -f $(DEMOS)



