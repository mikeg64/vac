!==============================================================================
!
!    THE FOLLOWING SUBROUTINES ADD THERMAL CONDUCTION SOURCE TERMS, 
!    SET THE THERMAL CONDUCTION COEFFICIENT KAPPA, AND CHECK DT
!
!    These subroutines were developed and tested by Rony Keppens.
!
!------------------------------------------------------------------------------
!    To include these subroutines into the VACUSR module at the top write
! 
!INCLUDE:vacusr.thermcond.t
!
!     then call addsource_tcond from the specialsource subroutine like
!
!if(abs(eqpar(kappa_))>smalldouble)&
!   call addsource_tcond(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)
!
!    where eqpar(kappa_) is the thermal conduction coefficient and kappa_ 
!    should be defined in vacusrpar.t like this
!
!INTEGER,PARAMETER:: kappa_=neqpar+1, nspecialpar=1
!
!    Positive value for eqpar(kappa_) defines a constant coefficient, 
!    use 0 for no thermal conduction.
!    The !!! comments show how a "kappa" array can be used if kappa is not 
!    constant. The "setkappa" subroutine has to be completed then.
!
!    To ensure numerical stability call the getdt_tcond subroutine from
!    getdt_special.
!    
!============================================================================
subroutine addsource_tcond(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)

! Add thermal conduction source to w within ixO based on time centered wCT

include 'vacdef.f'

integer::          ixI^L,ixO^L,iws(niw_)
double precision:: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)

integer:: ix^L,idir,iiw,iw
double precision, dimension(ixG^T):: thermcond

!!! ! Define common kappa array if the thermal conduction varies in space
!!! double precision:: kappa(ixG^T)
!!! common/tcond/kappa
!-----------------------------------------------------------------------------

oktest=index(teststr,'addsource')>=1
if(oktest)write(*,*)'addsource_tcond with kappa=',eqpar(kappa_)

! Calculating thermal conduction sources 
! involves second derivatives, two extra layers
call ensurebound(2,ixI^L,ixO^L,qtC,wCT)

ix^L=ixO^L^LADD1;

!!! ! kappa is calculated by getdt_tcond in the first time step if it varies 
!!! ! in space only. Calculate kappa in every time step if it varies in time
!!! call setkappa(wCT,ixI^L,ix^L,kappa)

! first calculate temperature from pressure, T=p/rho, i.e. k_Boltzman/mu=1.

call getpthermal(.true.,wCT,ixI^L,tmp)
tmp(ixI^S)=tmp(ixI^S)/wCT(ixI^S,rho_)
if(oktest)write(*,*)'temperature:',tmp(ixtest^D)

! Calculate source for energy from thermal conduction
do iiw=1,iws(niw_); iw=iws(iiw)
   if(iw==e_)then
      if(fourthorder.and.eqpar(kappa_)>0)then
         ! For constant kappa and uniform Cartesian grid 
         ! we can use fourth order laplace of temperature
         call laplace4(tmp,ixO^L,thermcond)
         w(ixO^S,ee_)=w(ixO^S,ee_)+thermcond(ixO^S)*eqpar(kappa_)*qdt
         return
      endif
      ! de/dt= +div(kappa*grad T), thus e=e+d_i gradt_i
      do idir=1,ndim
         ! Calculate temperature gradient within ixL: 
         ! tmp2= grad T, thus gradt_i=d_i T
         call gradient(.true.,tmp,ix^L,idir,tmp2)
         ! Multiply temperature gradient with conduction coefficient and dt
         tmp2(ix^S)=tmp2(ix^S)*eqpar(kappa_)*qdt
         !!! ! For spatially varying kappa use this instead of the line above
         !!! tmp2(ix^S)=tmp2(ix^S)*kappa(ix^S)*qdt
         ! take divergence of temperature gradient to get thermal conduction
         call gradient(.false.,tmp2,ixO^L,idir,thermcond)
         w(ixO^S,ee_)=w(ixO^S,ee_)+thermcond(ixO^S)

         if(oktest)write(*,*)'energy source:',thermcond(ixtest^D)
      enddo
   end if
end do

return
end

!=============================================================================
!!! subroutine setkappa(w,ixI^L,ixO^L,kappa)

! Set the thermal conduction coefficient kappa within ixO based on w(ixI).

!!! include 'vacdef.f'

!!! double precision:: w(ixG^T,nw),kappa(ixG^T)
!!! integer:: ixI^L,ixO^L

!----------------------------------------------------------------------------

!!! return
!!! end
!=============================================================================
subroutine getdt_tcond(w,ix^L)

! Check diffusion time limit for dt < dtdiffpar * dx**2 / ((gamma-1)*kappa/rho)

! Based on Hirsch volume 2, p.631, eq.23.2.17 and the remarks below

include 'vacdef.f'

double precision:: w(ixG^T,nw),dtdiff_tcond
integer:: ix^L,idim

!!! ! For spatially varying kappa you need a common kappa array
!!! double precision:: kappa(ixG^T)
!!! common/tcond/kappa

!-----------------------------------------------------------------------------

oktest=index(teststr,'getdt')>=1

! Calculate the thermal energy diffusivity tmp = (gamma-1)*kappa/rho
if(eqpar(kappa_)>zero)then
   ! For spatially uniform kappa:
   tmp(ix^S)=(eqpar(gamma_)-1)*eqpar(kappa_)/w(ix^S,rho_)
else
   ! For spatially varying kappa uncomment the 2 lines, and remove the 3rd
   !!! call setkappa(w,ixG^L,ix^L,kappa)
   !!! tmp(ix^S)=(eqpar(gamma_)-1)*kappa(ix^S)/w(ix^S,rho_)
   call die('For eqpar(kappa_)<0 edit vacusr.thermcond.t')
endif

do idim=1,ndim
    dtdiff_tcond=dtdiffpar/maxval(tmp(ix^S)/dx(ix^S,idim)**2)
    {^IFMPI call mpiallreduce(dtdiff_tcond,MPI_MIN)}

    ! limit the time step
    dt=min(dt,dtdiff_tcond)
    if(oktest)write(*,*)'Thermal cond. idim, dt:',idim, dtdiff_tcond
enddo 

return
end
!=============================================================================

