#! /bin/csh
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
if ($#argv != 1) then
   echo "***liblover.mk: usage: liblover.mk <CRAY|IBMRIOS|SGI|SUN>"
   exit 1
endif
#
setenv ARCH $1
#
if (!($?LOVER)) then
   setenv LOVER /usr/local/lover
endif
if (!($?LIB)) then
   setenv LIB /usr/local/lib
endif
if (!($?INCLUDE)) then
   setenv INCLUDE /usr/local/include
endif
if (!($?F77)) then
   if ($ARCH == CRAY) then
      setenv F77 f90
   else if ($ARCH == IBM) then
      setenv F77 xlf
   else 
      setenv F77 f77
   endif
endif
if (!($?F77FLAGS)) then
   if ($ARCH == CRAY) then
      setenv F77FLAGS ""
   else if ($ARCH == IBM) then
      setenv F77FLAGS "-O -qautodbl=dbl4"
   else if ($ARCH == DEC) then
      setenv F77FLAGS "-O -r8 -DLSB"
   else 
      setenv F77FLAGS "-O -r8"
   endif
endif
if (!($?CC)) then
   if ($ARCH == SUN) then
      setenv CC acc
   else
      setenv CC cc
   endif
endif
if (!($?CFLAGS)) then
   setenv CFLAGS "-O"
endif
#
if (-d . && -w .) then
   rm -rf tmp
   mkdir tmp
else
   echo ================================================
   echo No write access in working directory `pwd`
   echo ================================================
   exit 1
endif
#
echo =============================
echo Copying liblover.F
echo =============================
cp -f $LOVER/SOURCE/liblover.F tmp
cp -f $LOVER/SOURCE/dxbinwritedou.c tmp
cp -f $LOVER/SOURCE/dxbinwriteflo.c tmp
cp -f $LOVER/SOURCE/avsbinwritedou.c tmp
cp -f $LOVER/SOURCE/avsbinwriteflo.c tmp
cp -f $LOVER/SOURCE/writeff.c tmp
cp -f $LOVER/SOURCE/findoffset.c tmp
#
mkdir tmp/INCLUDE
cp -f $LOVER/SOURCE/INCLUDE/* tmp/INCLUDE
cd tmp
#
echo =============================
echo Splitting liblover.F         
echo =============================
fsplit liblover.F     
#
# fsplit on SUN splits always in .f files. 
if ($ARCH == "IBMRIOS" || $ARCH == "SUN") then
   foreach file (*.f)
      set basename = `echo $file | cut -f1 -d.`
      mv $file $basename.F
   end
endif
#
echo ===============================
echo Compiling all liblover routines    
echo ===============================
rm -f liblover.F
$F77 *.F *.f -c $F77FLAGS -D$ARCH -I./INCLUDE
$CC  *.c -c $CFLAGS -D$ARCH
#
echo =============================
echo Creating library liblover.a        
echo in directory $LIB          
echo =============================
ar vq liblover.a *.o
#
# ranlib not available under UNICOS
if ($ARCH != "CRAY" && $ARCH != "SGI") then
   ranlib liblover.a
endif
#
mv liblover.a $LIB
cd ..
rm -rf tmp

#end of shell script liblover

