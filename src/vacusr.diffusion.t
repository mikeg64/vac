!==============================================================================
!
!    THE FOLLOWING SUBROUTINES ADD DIFFUSIVE SOURCE TERMS AND LIMIT DT
!
!------------------------------------------------------------------------------
!    See vacusr.t.diffusion and vacusrpar.t.diffusion for an example of usage
!
!    d rho/dt += div (diff_i * grad_i (rho))
!
!    The eqpar(diff1_),eqpar(diff2_),... parameters are the diffusion 
!    coefficients in each dimension. Set them to 0 for no diffusion 
!    in that direction. It is trivial to simplify this library
!    to a single isotropic diffusion coefficient.
!
!    The !!! comments show how a diff array could be used for a spatially
!    (and maybe temporally) varying diffusionn coefficient.
!    "subroutine setdiff" has to be completed then.
!
!============================================================================
subroutine addsource_diff(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)

! Add diffusion source calculated from wCT to w within ixO for all variables 
! in iws. wCT is at time qtC, w is advanced from qt to qt+qdt.

include 'vacdef.f'

integer::          ixI^L,ixO^L,iws(niw_)
double precision:: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)
integer:: ix^L,iiw,iw,idim

!!! ! Define common diff array if the diffusion coefficient varies in space
!!! double precision:: diff(ixG^T,ndim)
!!! common/diffuse/diff

!-----------------------------------------------------------------------------

oktest=index(teststr,'addsource')>=1
if(oktest)write(*,*)'addsource_diff with diff=', ^D&eqpar(diff^D_)

! Calculating diffusive sources involves second derivatives, two extra layers:
call ensurebound(2,ixI^L,ixO^L,qtC,wCT)

ix^L=ixO^L^LADD1;

!!! ! If diff needs to be calculated only once do it for the whole grid
!!! if(it==itmin)call setdiff(w,ixG^L,ixG^L,diff)
!!! ! Otherwise call setdiff in every time step
!!! call setdiff(w,ixI^L,ix^L,diff)

! add sources from diffusion
do iiw=1,iws(niw_); iw=iws(iiw)
   select case(iw)
   case(rho_)
      ! drho/dt= +div(diff*grad rho)
      do idim=1,ndim
         if(abs(eqpar(diff0_+idim))>smalldouble)then
            ! Extract density
            tmp(ixI^S)=w(ixI^S,rho_)
            ! Calculate density gradient within ixL
            call gradient(.true.,tmp,ix^L,idim,tmp2)
            ! Multiply density gradient with diffusion coefficient and dt
            tmp2(ix^S)=tmp2(ix^S)*eqpar(diff0_+idim)*qdt
            !!! For spatially varying diff use this instead of the line above
            !!! tmp2(ix^S)=tmp2(ix^S)*diff(ix^S,idim)*qdt

            ! Calculate divergence
            call gradient(.false.,tmp2,ixO^L,idim,tmp)
            ! Add source
            w(ixO^S,rho_)=w(ixO^S,rho_)+tmp(ixO^S)
         endif
      enddo
   end select ! iw
end do        ! iiw

return
end
!=============================================================================
!!! subroutine setdiff(w,ixI^L,ixO^L,diff)

! Set the diffusion coefficient array within ixO based on x(ixI,ndim) 
! and/or w(ixI,nw)

!!! include 'vacdef.f'

!!! double precision:: w(ixG^T,nw),diff(ixG^T,ndim)
!!! integer:: ixI^L,ixO^L
!----------------------------------------------------------------------------
!!! return
!!! end

!=============================================================================
subroutine getdt_diff(w,ix^L)

! Check diffusion time limit for dt

include 'vacdef.f'

double precision:: w(ixG^T,nw),dtdiff_diff
integer:: ix^L,idim
save dtdiff_diff

!!! ! For spatially varying diffusion coefficient you need a common diff array
!!! double precision:: diff(ixG^T,ndim)
!!! common/diffuse/diff
!-----------------------------------------------------------------------------

oktest=index(teststr,'getdt')>=1

if(it==itmin)then
   ! If the diffusion coefficient is defined by the equation parameters:
   dtdiff_diff=bigdouble
   do idim=1,ndim
      if(eqpar(diff0_+idim)>zero)&
      dtdiff_diff=min(dtdiff_diff,minval(dx(ix^S,idim))**2/eqpar(diff0_+idim))
   enddo
   dtdiff_diff=dtdiffpar*dtdiff_diff

   !!! ! If the diffusion coefficient varies spatially, use this instead of 
   !!! ! the lines above
   !!! call setdiff(w,ix^L,ix^L,diff)
   !!! ! If diff does not vary with time, calculate dtdiff_diff here
   !!! dtdiff_diff=dtdiffpar*minval(dx(ix^S,1:ndim)**2/diff(ix^S,1:ndim))
endif
!!! ! If diff varies with time, calculate dtdiff_diff here
!!! dtdiff_diff=dtdiffpar*minval(dx(ix^S,1:ndim)**2/diff(ix^S,1:ndim))

{^IFMPI call mpiallreduce(dtdiff,MPI_MIN)}

! limit the time step
dt=min(dt,dtdiff_diff)

return
end
!=============================================================================



