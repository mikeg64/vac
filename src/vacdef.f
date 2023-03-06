!##############################################################################
! include vacdef

IMPLICIT NONE

!HPF$ PROCESSORS PP(NUMBER_OF_PROCESSORS())

! DEFINITIONS OF GLOBAL PARAMETERS AND VARIABLES
! Parameters:

! Indices for cylindrical coordinates FOR TESTS, negative value when not used:
INTEGER,PARAMETER:: r_=1, phi_=-9, z_=-8

! Indices for cylindrical coordinates FOR INDEXING, always positive
INTEGER,PARAMETER:: pphi_=1, zz_=1

include 'vacpar.f'

INTEGER,PARAMETER:: ixGlo1=1,ixGlo2=1
! The next line is edited by SETVAC
INTEGER,PARAMETER:: ixGhi1=104,ixGhi2=104,ixGhimin=104,ixGhimax=104
INTEGER,PARAMETER:: ndim=2, ndir=2

INTEGER,PARAMETER:: dixBlo=2,dixBhi=2

 !Size of work array for VACPOISSON
INTEGER,PARAMETER:: nhi=nw*ixGhi1*ixGhi2 !Maximum number of unknowns for VACIMPL
 

INTEGER,PARAMETER:: nhiB=10           ! maximum No. boundary sections

INTEGER,PARAMETER:: nsavehi=100       ! maximum No. saves into outputfiles
                                      ! defined by arrays of tsave or itsave

INTEGER,PARAMETER:: niw_=nw+1         !Indexname for size of iw index array

INTEGER,PARAMETER:: filelog_=1,fileout_=2,nfile=2 ! outputfiles

INTEGER,PARAMETER:: unitstdin=5,unitterm=6,uniterr=6,unitini=10 ! Unit names. 
                                  ! Outputfiles use unitini+1..initini+nfile
                                  ! Default parfiles uses unitini-1

INTEGER,PARAMETER:: biginteger=10000000

DOUBLE PRECISION,PARAMETER:: pi= 3.1415926535897932384626433832795
DOUBLE PRECISION,PARAMETER:: smalldouble=1.D-99, bigdouble=1.D+99
DOUBLE PRECISION,PARAMETER:: zero=0D0,one=1D0,two=2D0,half=0.5D0,quarter&
   =0.25D0

INTEGER,PARAMETER:: toosmallp_=1,toosmallr_=2,couranterr_=3,poissonerr_=4
INTEGER,PARAMETER:: nerrcode=4

include 'vacusrpar.f'



!-- Common variables:







! Unit for reading input parameters.
INTEGER:: unitpar

! Logical to set verbosity. For MPI parallel run only PE 0 is verbose
LOGICAL:: verbose

! General temporary arrays, any subroutine call may change them 
! except for subroutines which say the opposite in their header
DOUBLE PRECISION:: tmp,tmp2

! Number of errors during calculation
INTEGER:: nerror

!Kronecker delta and Levi-Civita tensors
INTEGER:: kr,lvc

!Grid parameters
INTEGER:: ixMmin1,ixMmin2,ixMmax1,ixMmax2,ixGmin1,ixGmin2,ixGmax1,ixGmax2,nx1,&
   nx2,nx
INTEGER:: dixBmin1,dixBmin2,dixBmax1,dixBmax2
! x and dx are local for HPF
DOUBLE PRECISION:: x,dx
DOUBLE PRECISION:: volume,dvolume
DOUBLE PRECISION:: area,areaC
DOUBLE PRECISION:: areadx,areaside

! Variables for generalized coordinates and polargrid
LOGICAL::          gencoord, polargrid

DOUBLE PRECISION:: surfaceC, normalC

!Boundary region parameters
DOUBLE PRECISION:: fixB1,fixB2
INTEGER:: nB,ixBmin,ixBmax,idimB,ipairB
LOGICAL:: upperB,fixedB,nofluxB,extraB
CHARACTER*10 :: typeB,typeBscalar

!Equation and method parameters
DOUBLE PRECISION:: eqpar,procpar

! Time step control parameters
DOUBLE PRECISION:: courantpar,dtpar,dtdiffpar,dtcourant,dtmrpc
LOGICAL:: dtcantgrow
INTEGER:: slowsteps

! Parameters for the implicit techniques
 

INTEGER:: nwimpl,nimpl
DOUBLE PRECISION:: implpar,impldiffpar,implerror,implrelax,impldwlimit
INTEGER:: implrestart,implrestart2,impliter,impliternr,implmrpcpar
CHARACTER*10 :: typeimplinit,typeimpliter,typeimplmat
LOGICAL:: implconserv,implnewton,implcentered,implnewmat
LOGICAL:: implpred,impl3level,impljacfast,implsource

!Method switches
INTEGER:: iw_full,iw_semi,iw_impl,iw_filter
INTEGER:: iw_vector,vectoriw
          ! The upper bound+1 in iw_vector avoids F77 warnings when nvector=0
CHARACTER*10 :: typefull1,typepred1,typeimpl1,typefilter1
CHARACTER*10 :: typelimited,typefct,typetvd,typeaxial
CHARACTER*10 :: typepoisson, typeconstrain
CHARACTER*10 :: typelimiter,typeentropy
CHARACTER*10 :: typeadvance, typedimsplit, typesourcesplit
LOGICAL:: dimsplit,sourcesplit,sourceunsplit,artcomp,useprimitive
LOGICAL:: divbfix,divbwave,divbconstrain,angmomfix,compactres,smallfix
INTEGER:: idimsplit
INTEGER:: nproc
DOUBLE PRECISION:: entropycoef,constraincoef
DOUBLE PRECISION:: smallp,smallpcoeff,smallrho,smallrhocoeff,vacuumrho
DOUBLE PRECISION:: muscleta1,muscleta2,musclomega,acmcoef,acmexpo
LOGICAL:: acmnolim, fourthorder
INTEGER:: acmwidth

!Previous time step and residuals
DOUBLE PRECISION:: wold,residual,residmin,residmax

! Flux storage for flux-CT and flux-CD methods !!! for MHD only !!! 


!Time parameters
INTEGER:: step,istep,nstep,it,itmin,itmax,nexpl,nnewton,niter,nmatvec
DOUBLE PRECISION:: t,tmax,dt,dtmin,cputimemax
LOGICAL:: tmaxexact
DOUBLE PRECISION:: tsave,tsavelast,dtsave
INTEGER:: itsave,itsavelast,ditsave
INTEGER:: isavet,isaveit

!File parameters
CHARACTER*79 :: filenameini,filenameout,filename
CHARACTER*79 :: fileheadini,fileheadout,varnames,wnames
CHARACTER*10 :: typefileini,typefileout,typefilelog
LOGICAL::             fullgridini,fullgridout
INTEGER::             snapshotini,snapshotout,isaveout

!Test parameters
CHARACTER*79 :: teststr
INTEGER:: ixtest1,ixtest2,ixtest3,iwtest,idimtest
LOGICAL:: oktest !This is a local variable for all subroutines and functions

! end include vacdef
!##############################################################################
COMMON /DOUB/ tmp(ixGlo1:ixGhi1,ixGlo2:ixGhi2),tmp2(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2),x(IXGlo1:IXGhi1,IXGlo2:IXGhi2,ndim),dx(IXGlo1:IXGhi1,&
   IXGlo2:IXGhi2,ndim),volume,dvolume(IXGlo1:IXGhi1,IXGlo2:IXGhi2),&
   area(IXGLO1:IXGHI1),areaC(IXGLO1:IXGHI1),areadx(IXGLO1:IXGHI1),&
   areaside(IXGLO1:IXGHI1),surfaceC(2,2,ndim), normalC(2,2,ndim,ndim),fixB1(&
   -dixBlo:dixBhi,ixGLO2:ixGHI2,nw),fixB2(ixGLO1:ixGHI1,-dixBlo:dixBhi,nw),&
   eqpar(neqpar+nspecialpar),procpar(nprocpar),courantpar,dtpar,dtdiffpar,&
   dtcourant(ndim),dtmrpc,implpar,impldiffpar,implerror,implrelax,impldwlimit,&
   entropycoef(nw),constraincoef,smallp,smallpcoeff,smallrho,smallrhocoeff,&
   vacuumrho,muscleta1,muscleta2,musclomega,acmcoef(nw),acmexpo,&
   wold(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),residual,residmin,residmax,t,tmax,dt,&
   dtmin,cputimemax,tsave(nsavehi,nfile),tsavelast(nfile),dtsave(nfile)
COMMON /INTE/ unitpar,nerror(nerrcode),kr(3,3),lvc(3,3,3),ixMmin1,ixMmin2,&
   ixMmax1,ixMmax2,ixGmin1,ixGmin2,ixGmax1,ixGmax2,nx1,nx2,nx(ndim),dixBmin1,&
   dixBmin2,dixBmax1,dixBmax2,nB,ixBmin(ndim,nhiB),ixBmax(ndim,nhiB),&
   idimB(nhiB),ipairB(nhiB),slowsteps,nwimpl,nimpl,implrestart,implrestart2,&
   impliter,impliternr,implmrpcpar,iw_full(niw_),iw_semi(niw_),iw_impl(niw_),&
   iw_filter(niw_),iw_vector(nvector+1),vectoriw(nw),idimsplit,nproc(nfile&
   +2),acmwidth,step,istep,nstep,it,itmin,itmax,nexpl,nnewton,niter,nmatvec,&
   itsave(nsavehi,nfile),itsavelast(nfile),ditsave(nfile),isavet(nfile),&
   isaveit(nfile),snapshotini,snapshotout,isaveout,ixtest1,ixtest2,ixtest3,&
   iwtest,idimtest
COMMON /CHAR/ typeB(nw,nhiB),typeBscalar(nhiB),typeimplinit,typeimpliter,&
   typeimplmat,typefull1,typepred1,typeimpl1,typefilter1,typelimited,typefct,&
   typetvd,typeaxial,typepoisson, typeconstrain,typelimiter(nw),&
   typeentropy(nw),typeadvance, typedimsplit, typesourcesplit,filenameini,&
   filenameout,filename(nfile),fileheadini,fileheadout,varnames,wnames,&
   typefileini,typefileout,typefilelog,teststr
COMMON /LOGI/ verbose,gencoord, polargrid,upperB(nhiB),fixedB(nw,nhiB),&
   nofluxB(nw,ndim),extraB,dtcantgrow,implconserv,implnewton,implcentered,&
   implnewmat,implpred,impl3level,impljacfast,implsource,dimsplit,sourcesplit,&
   sourceunsplit,artcomp(nw),useprimitive,divbfix,divbwave,divbconstrain,&
   angmomfix,compactres,smallfix,acmnolim, fourthorder,tmaxexact,fullgridini,&
   fullgridout
