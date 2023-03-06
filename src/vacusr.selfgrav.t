!==============================================================================
!
!    THE FOLLOWING SUBROUTINES ADD GRAVITATIONAL SOURCE TERMS.
!    GRAVITATIONAL FORCE IS CALCULATED FROM DENSITY BY SOLVING 
!    A POISSON EQUATION
!
!------------------------------------------------------------------------------
!    See vacusr.t.selfgrav and vacusrpar.t.selfgrav for an example of usage
!
!    The eqpar(grav_) coefficient is used in the Poisson equation
!    for the gravitational potential Phi:
!
!    div grad Phi = eqpar(grav_) * rho
!
!    Note that eqpar(grav_) includes the usual 1/(4 pi) factor too!
!    Set eqpar(grav_)=0 to switch the gravitational sources off.
!    The components of the gravitational force are calculated from
!
!    grav_i = - grad_i (Phi)
!
!    Gravitational force is added to the momentum equation:
!
!    d m_i/dt += rho*grav_i
!
!    Gravitational work is added to the energy equation (if present):
!
!    de/dt += Sum_i m_i*grav_i
!
!    The boundary conditions and the parameters for the Poisson solver
!    can be defined in the namelist &GRAVLIST at the end of the parameter 
!    file for VAC par/PROBLEM. For example
!
! &gravlist
!        typeBgrav=      4*'nul'
!        typestop=       'max'
!        tolerance=      0.01
!        matvecmax=      100
!        useoldphi=      T
! /
!
!    The Poisson problem is solved to an accuracy 
! 
!    max( | G*rho - div grad Phi | ) < 0.01
!
!    using at most 100 iterations, and the initial guess for Phi is taken
!    from the previous time step. The boundary conditions are Phi=0 for 
!    all boundaries since typeBgrav='nul'. Other possibilities are 
!    typestop='rel' or 'abs' and typeBgrav='periodic','cont','symm','asymm'
!
!============================================================================
subroutine addsource_grav(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)

! Add gravity source calculated from wCT to w within ixO for all variables 
! in iws. wCT is at time qtC, w is advanced from qt to qt+qdt.

include 'vacdef.f'

integer::          ixI^L,ixO^L,iws(niw_)
double precision:: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)
integer:: iiw,iw,idim
logical:: oktime
double precision:: time1,timegrav,cputime
data timegrav /0.D0/

double precision:: grav(ixG^T,ndim)
common /gravity/ grav
!-----------------------------------------------------------------------------
oktest=index(teststr,'selfgrav')>=1

oktime=index(teststr,'timegrav')>=1
if(oktest)write(*,*)'Addsource_grav old w:',w(ixtest^D,iwtest)

if(oktime)time1=cputime()

call ensurebound(1,ixI^L,ixO^L,qtC,wCT)

! Calculate gravitational force array grav
call setgrav(wCT,ixI^L,ixO^L,grav)

if(oktest)write(*,*)'grav:', (grav(ixtest^D,idim),idim=1,ndim)

! add sources from gravity
do iiw=1,iws(niw_); iw=iws(iiw)
   select case(iw)
   case(m^D_)
      ! dm_i/dt += rho*g_i
      idim=iw-m0_
      w(ixO^S,m0_+idim)=w(ixO^S,m0_+idim)+ &
          qdt*grav(ixO^S,idim)*wCT(ixO^S,rho_)

   case(e_)
      ! de/dt += Sum_i g_i*m_i
      do idim=1,ndim
         w(ixO^S,ee_)=w(ixO^S,ee_)+ &
             qdt*grav(ixO^S,idim)*wCT(ixO^S,m0_+idim)
      end do
   end select ! iw
end do        ! iiw

if(oktime)then
    time1=cputime()-time1; timegrav=timegrav+time1
    write(*,*)'Grav.time:',time1,timegrav
endif

if(oktest)write(*,*)'Addsource_grav new w:',w(ixtest^D,iwtest)

return
end
!=============================================================================
subroutine setgrav(w,ixI^L,ixO^L,grav)

! Set the gravitational acceleration within ixO based on w(ixI,rho) 

include 'vacdef.f'

double precision:: w(ixG^T,nw),grav(ixG^T,ndim)
integer:: ixI^L,ixO^L

double precision:: grho(ixG^T),phi(ixG^T),tolerance
integer:: idim,matvecmax,info,ios
character*3:: typestop
CHARACTER*^LENTYPE :: typeBgrav(nhiB)
logical:: initialized,useoldphi

! Default values for iterative solver
data typestop /'max'/
data tolerance /0.005/
data matvecmax /100/
data useoldphi /.true./

! The parameter file may contain a &GRAVLIST ... / to overwrite defaults
namelist /gravlist/ typeBgrav,typestop,tolerance,matvecmax,useoldphi

data initialized /.false./

! Gravitational potential is saved for initial guess at the next time step
! Only needed if useoldphi is true
save phi,typeBgrav
!----------------------------------------------------------------------------

oktest=index(teststr,'setgrav')>=1
if(oktest)write(*,*)'SetGrav rho:',w(ixtest^D,rho_)

if(.not.initialized)then
   phi(ixG^S)=zero

   typeBgrav(1:nhiB)='nul'

   tolerance=tolerance*eqpar(grav_)*sum(w(ixM^S,rho_))/(nx^D*)

   write(*,*)'SetGrav calculates gravitational potential using Poisson solver'

   ! Attempt reading user defined values for iteration parameters
   read(*,gravlist,iostat=ios)
   if(ios/=0)write(*,*)&
      'No &GRAVLIST in parameter file. Default values are used.'
   write(*,gravlist)
   initialized=.true.
endif

! Right hand side, boundary cond., and initial guess for the Poisson solver
grho(ixM^S)=eqpar(grav_)*w(ixM^S,rho_)
typeBscalar(1:nB)=typeBgrav(1:nB)
if(.not.useoldphi)phi(ixG^S)=zero

{^IFPOISSON
call poisson('Self grav. ',&
   grho,tolerance,typestop,matvecmax,info,useoldphi,phi)
!} call die('Error: Poisson solver is OFF! setvac -on=poisson; make vac')

if(oktest)write(*,*)'Poisson solver info:',info

! Stop if solver failed to converge
if(info<0)then
   write(*,*)'Time step, info:',it,info
   call die('Error in SelfGrav: Could not solve the Poisson equation!')
endif

! Calculate gravitational force grav=-grad(phi)

! First get the ghost cell values for the solution
call boundscalar(phi);

do idim=1,ndim
   call gradient(.true.,phi,ixO^L,idim,tmp)
   grav(ixO^S,idim)=-tmp(ixO^S)
end do

if(oktest)write(*,*)'SetGrav grav:', (grav(ixtest^D,idim),idim=1,ndim)

return
end
!=============================================================================

subroutine getdt_grav(w,ix^L)

! Time step limit is based on the condition
!
!    dt < sqrt(dx/grav)

include 'vacdef.f'

double precision:: w(ixG^T,nw)
integer:: ix^L
double precision:: dtgrav
save dtgrav

double precision:: grav(ixG^T,ndim)
common/gravity/grav
!----------------------------------------------------------------------------

oktest=index(teststr,'getdt')>=1

if(it==itmin)call setgrav(w,ixG^L,ix^L,grav)

dtgrav=one/sqrt(maxval(abs(grav(ix^S,1:ndim))/dx(ix^S,1:ndim)))

{^IFMPI call mpiallreduce(dtgrav,MPI_MIN)}

! limit the time step
dt=min(dt,dtgrav)

if(oktest)write(*,*)'Gravity limit for dt:',dtgrav

return
end

!=============================================================================
