C ============================================================================
C
C This is an include file for the conversion program VAC2AVSDX in src/convert.f
C
C ============================================================================

C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C THE WORK ARRAY SHOULD BE ABLE TO CONTAIN A FULL SNAPSHOT FOR ASCII-BINARY
C CONVERSION, AND ABOUT 6 SNAPSHOT SIZE FOR AVS/DX CONVERSION.
C
C THEREFORE THE SIZE "worksize" MAY HAVE TO BE ADJUSTED FOR VERY BIG FILES, 
C OR FOR MACHINES WITH VERY LITTLE MEMORY.
C
C NPICTMAX should be greater than the number of snapshots in a file.
C
C VERBOSE determines the amount of info printed. Possible values: 0, 1, 2

      INTEGER worksize,npictmax,verbose
      PARAMETER (worksize=8000000,npictmax=1000,verbose=2)

      REAL*8 work(worksize)
      COMMON /WORK/ work

C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C PARAMETERS

      INTEGER unitini, unitout
      PARAMETER (unitini=10, unitout=12)

C Maximum number of snapshots in one file to be converted
C and maximum number of scalars and vectors in a file (MHD has 4: rho,M,e,B)

      INTEGER maxntime,maxnfields
      PARAMETER (maxntime=1000,maxnfields=4)

C GLOBAL VARIABLES

      INTEGER it,ndim,signedndim,neqpar,nw,nx(3)
      INTEGER nx1,nx2,nx3,nxs,ndir,pictsize,ipict,npict
      INTEGER ios

      COMMON /INTE/ it,ndim,signedndim,neqpar,nw,nx,
     &     nx1,nx2,nx3,nxs,ndir,pictsize,ipict,npict,
     &     ios

      REAL*8 t,eqpar(100)
      COMMON /DOUB/ t,eqpar

      CHARACTER*79 filenameini,filenameout,fileheadini,varnames
      CHARACTER*10 typefileini,typefileout

      COMMON /CHAR/ filenameini,filenameout,fileheadini,varnames,
     &     typefileini,typefileout

      LOGICAL readpict(npictmax)
      COMMON /LOGI/ readpict

C ============================================================================
