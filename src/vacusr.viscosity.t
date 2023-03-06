!==============================================================================
!
!    THE FOLLOWING SUBROUTINES ADD VISCOUS SOURCE TERMS, SET NU, AND CHECK DT
!
!    These subroutines were developed and tested by Rony Keppens.
!
!------------------------------------------------------------------------------
!    See vacusr.t.viscosity and vacusrpar.t.viscosity for an example of usage
!
!    Viscous forces in the momentum equations:
!
!    d m_i/dt +=  - div (eqpar(nu_) * PI)
!
!    Viscous work in the energy equation:
!
!    de/dt    += - div (v . eqpar(nu_) * PI)
!
!    where the PI stress tensor is
!
!    PI_i,j = - (dv_j/dx_i + dv_i/dx_j) + (2/3)*Sum_k dv_k/dx_k
!
!    where eqpar(nu_) is the viscosity coefficient. Positive value for nu 
!    defines a constant viscosity, and use 0 for no viscosity. 
!    The !!! comments show how a "nu" array can be used if nu is not constant.
!    The "setnu" subroutine has to be completed then.
!
!==============================================================================
subroutine addsource_visc(qdt,ixI^L,ixO^L,iws,qtC,w,qt,wnew)

! Add viscosity source to wnew within ixO 

include 'vacdef.f'

integer::          ixI^L,ixO^L,iws(niw_)
double precision:: qdt,qtC,qt,w(ixG^T,nw),wnew(ixG^T,nw)

integer:: ix,ix^L,idim,idir,jdir,iiw,iw
double precision:: lambda(ixG^T,ndir,ndir)

!!! ! Define common nu array if the viscosity varies in space
!!! double precision:: nu(ixG^T)
!!! common/visc/nu
!-----------------------------------------------------------------------------

oktest=index(teststr,'addsource')>=1
if(oktest)write(*,*)'addsource_visc with nu=',eqpar(nu_)

! Calculating viscosity sources 
! involves second derivatives, two extra layers
call ensurebound(2,ixI^L,ixO^L,qtC,w)
ix^L=ixO^L^LADD1;

!!! ! nu array is calculated in the first time step by getdt_visc if it varies
!!! ! in space only. Calculate nu in every time step if it varies in time:
!!! call setnu(w,ixI^L,ix^L,nu)

! construct lambda tensor: lambda_ij = gradv_ij + gradv_ji
! initialize

lambda(ix^S,1:ndir,1:ndir)=zero

!next construct 

do idim=1,ndim; do idir=1,ndir
! Calculate velocity gradient tensor within ixL: gradv= grad v, 
! thus gradv_ij=d_j v_i
  tmp(ixI^S)=w(ixI^S,m0_+idir)/w(ixI^S,rho_)
  call gradient(.true.,tmp,ix^L,idim,tmp2)
  if(oktest)write(*,*)'i,j,d_j v_i:',idir,idim,tmp2(ixtest^D)
  lambda(ix^S,idim,idir)= lambda(ix^S,idim,idir)+ tmp2(ix^S)
  lambda(ix^S,idir,idim)= lambda(ix^S,idir,idim)+ tmp2(ix^S)
enddo; enddo;

{^IFPHI
! special treatment of axial symmetric case 
!   added term to d_phi v_r   = -v_phi/r
!   added term to d_phi v_phi = +v_r/r
! enters lambda tensor in r_phi_,phi_r_,phi_phi_
if(typeaxial=='cylinder'.and. phi_ >1)then
   if(gencoord)then
      lambda(ix^S,phi_,r_)=lambda(ix^S,phi_,r_)  &
            - w(ix^S,phi_)/w(ix^S,rho_)/x(ix^S,r_)
      lambda(ix^S,phi_,phi_)=lambda(ix^S,phi_,phi_)  &
            + two*w(ix^S,r_)/w(ix^S,rho_)/x(ix^S,r_)
   else
     forall(ix= ix^LIM1:)lambda(ix,ix^SE,phi_,r_)=lambda(ix,ix^SE,phi_,r_) &
                  -w(ix,ix^SE,phi_)/w(ix,ix^SE,rho_)*areaside(ix)
     forall(ix= ix^LIM1:)lambda(ix,ix^SE,phi_,phi_)=lambda(ix,ix^SE,phi_,phi_)&
                  +two*w(ix,ix^SE,r_)/w(ix,ix^SE,rho_)*areaside(ix)
   endif
   lambda(ix^S,r_,phi_)=lambda(ix^S,phi_,r_)  
endif
}

! Multiply lambda with viscosity coefficient and dt

do idir=1,ndir; do jdir=1,ndir
  lambda(ix^S,idir,jdir)=lambda(ix^S,idir,jdir)*eqpar(nu_)*qdt
  !!! ! Use the following line instead of the one above if nu varies in space
  !!! lambda(ix^S,idir,jdir)=lambda(ix^S,idir,jdir)*nu(ix^S)*qdt
enddo; enddo;

if(oktest)write(*,*)'ndir,lambda(1:ndir,1:ndir):',&
          ndir,lambda(ixtest^D,1:ndir,1:ndir)

!calculate div v term through trace action separately

tmp(ix^S)=(^C&lambda(ix^S,^C,^C)+)/3

if(oktest)write(*,*)'trace:',tmp(ixtest^D)

!substract trace from diagonal elements

do idir=1,ndir
   lambda(ix^S,idir,idir)=lambda(ix^S,idir,idir)-tmp(ix^S)
enddo

! Calculate source terms from viscosity 
do iiw=1,iws(niw_); iw=iws(iiw)
   select case(iw)
   case(m^C_)
      ! dm/dt= +div(nu*[d_j v_i+d_i v_j]-(2*nu/3)* div v * kr) 
      ! hence m_j=m_j+d_i tensor_ji
      idir=iw-m0_
      do idim=1,ndim 
            tmp(ix^S)=lambda(ix^S,idir,idim)
            call gradient(.true.,tmp,ixO^L,idim,tmp2)
            wnew(ixO^S,iw)=wnew(ixO^S,iw)+tmp2(ixO^S)
            if(oktest)write(*,*)'momentum source for iw:',iw,tmp2(ixtest^D)
      enddo
      ! special treatment of axial symmetric case 
      !   added term to d m_r/dt   = (lambda_rr - lambda_phiphi )/r
      !   added term to d m_phi/dt = 2*lambda_rphi/r
      !   added term to d m_z/dt   = lambda_rz/r
      if(typeaxial=='cylinder')then
        if(oktest)write(*,*)'correcting for axial symmetry'
        select case(idir)
           case(r_)
               tmp(ixO^S)=lambda(ixO^S,r_,r_)
               if(phi_>1) tmp(ixO^S)=tmp(ixO^S)-lambda(ixO^S,pphi_,pphi_)
           case(phi_)
               tmp(ixO^S)=two*lambda(ixO^S,r_,pphi_)
           case(z_)
               tmp(ixO^S)=lambda(ixO^S,r_,zz_)
        endselect
        if(gencoord)then
            wnew(ixO^S,iw)=wnew(ixO^S,iw)+tmp(ixO^S)/x(ixO^S,r_)
            if(oktest)write(*,*)'momentum source for iw:', &
                iw,tmp(ixtest^D)/x(ixtest^D,r_)
        else
            forall(ix= ixO^LIM1:) wnew(ix,ixO^SE,iw)=wnew(ix,ixO^SE,iw)  &
                 +tmp(ix,ixO^SE)*areaside(ix)
        endif
      endif
   case(e_)
      ! de/dt= +div(v.dot.[nu*[d_j v_i+d_i v_j]-(2*nu/3)* div v *kr])
      ! thus e=e+d_i v_j tensor_ji
      do idim=1,ndim
         tmp(ix^S)=(^C&w(ix^S,m^C_)*lambda(ix^S,^C,idim)+)/w(ix^S,rho_)
         ! axial symmetric case covered using gradient(.false.)=divergence
         call gradient(.false.,tmp,ixO^L,idim,tmp2)
         wnew(ixO^S,ee_)=wnew(ixO^S,ee_)+tmp2(ixO^S)
         if(oktest)write(*,*)'energy source:',tmp2(ixtest^D)
      enddo
   end select ! iw
end do        ! iiw

return
end

!=============================================================================
!!! subroutine setnu(w,ixI^L,ixO^L,nu)

! Set the viscosity coefficient nu within ixO based on w(ixI). 

!!! include 'vacdef.f'

!!! double precision:: w(ixG^T,nw),nu(ixG^T)
!!! integer:: ixI^L,ixO^L
!----------------------------------------------------------------------------
!!! return
!!! end

!=============================================================================
subroutine getdt_visc(w,ix^L)

! Check diffusion time limit for dt < dtdiffpar * dx**2 / (nu/rho)

! Based on Hirsch volume 2, p.631, eq.23.2.17

include 'vacdef.f'

double precision:: w(ixG^T,nw),dtdiff_visc
integer:: ix^L,idim

!!! ! For spatially varying nu you need a common nu array
!!! double precision:: nu(ixG^T)
!!! common/visc/nu
!-----------------------------------------------------------------------------

oktest=index(teststr,'getdt')>=1

if(abs(eqpar(nu_))<smalldouble) return

! Calculate the dynamic viscosity tmp=nu/rho

if(eqpar(nu_)>zero)then
   ! For a constant viscosity:
   tmp(ix^S)=eqpar(nu_)/w(ix^S,rho_)
else
   ! For spatially varying nu uncomment the 2 lines below and delete the 3rd
   !!! call setnu(w,ixG^L,ixG^L,nu)
   !!! tmp(ix^S)=nu(ix^S)/w(ix^S,rho_)
   call die('For eqpar(nu_)<0 edit vacusr.viscosity.t')
endif

do idim=1,ndim
   dtdiff_visc=dtdiffpar/maxval(tmp(ix^S)/dx(ix^S,idim)**2)
   {^IFMPI call mpiallreduce(dtdiff_visc,MPI_MIN)}

   ! limit the time step
   dt=min(dt,dtdiff_visc)
   if(oktest) write(*,*)'idim, viscous dt:',idim,dtdiff_visc
enddo

return
end
!=============================================================================

