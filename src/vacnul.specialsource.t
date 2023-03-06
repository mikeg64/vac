!=============================================================================
! There are 3 subroutines in this include file:
!
! specialsource -- for sources other than resistivity
! getdt_special -- for time step conditions other than CFL or resistivity
! specialeta    -- for non-constant resistivity with eqpar(eta_)<zero
!
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)

! Calculate w(iw)=w(iw)+qdt*SOURCE(wCT,iw) within ixO for all indices
! iw=iws(iiw) where the indirect index iiw=1..iws(niw_).
! wCT is at time qCT, while w is at time qt on input and qt+qdt on output.
!
! You may want to use second order accurate time integration if 
! "sourcesplit=T" and/or the "typefull='tvd'" method is used in the par-file.
!
! If the source needs wCT values outside ixO, ensurebound should be used:
!
! call ensurebound(dix,ixI^L,ixO^L,qtC,wCT)
!
! where "dix" is the number of extra layers needed, typically 1 or 2.

include 'vacdef.f'

integer:: ixI^L,ixO^L,iws(niw_)
double precision:: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)

! These variables will most likely be needed:
!
! integer:: iw,iiw
! double precision:: s(ixG^T)
!-----------------------------------------------------------------------------

! Here are some typical structures to be used
!
! do iiw=1,iws(niw_); iw=iws(iiw)
!    select case(iw)
!    case(m1_)
!       ! The source is based on the time centered wCT
!       call getmyforce(wCT,ixO^L,s)
!       w(ixO^S,m1_)=w(ixO^S,m1_) + qdt*s(ixO^S)
!    case(e_)
!       call getmyheating(wCT,ixO^L,s)
!       w(ixO^S,e_) =w(ixO^S,e_)  + qdt*s(ixO^S)
!    end select
! end do

return
end

!=============================================================================
subroutine getdt_special(w,ix^L)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the VACPHYS 
! module have already been called. 

include 'vacdef.f'

double precision:: w(ixG^T,nw)
integer:: ix^L
!-----------------------------------------------------------------------------

! Something like this
!
! dt_mine=...
! dt=min(dt,dt_mine)

return
end

!=============================================================================
subroutine specialeta(w,ix^L,idirmin)

! Set the common "eta" array for resistive MHD based on w and the common
! "current" variable which has components between idirmin and 3.
! If resistivity should be treated implicitly, the "gradeta" array 
! has to be set. It may use the closest neighbours of "w" (jx,ix,hx), 
! and/or local values (ix) of "current" only, i.e. it has to be compact.

! REMOVE THE "!!!" COMMENT SIGNS FROM THE COMMON ARRAY DECLARATIONS!

include 'vacdef.f'

double precision:: w(ixG^T,nw)
integer:: ix^L,idirmin

!!! double precision:: current(ixG^T,7-2*ndir:3),eta(ixG^T),gradeta(ixG^T,ndim)
!!! common/resist/current,eta,gradeta
!-----------------------------------------------------------------------------

! If eta depends on local values of w only (definitely NOT on the current)
! you can set gradeta like this:
!
!   eta(ixG^S)=...
!
!   do idim=1,ndim
!      call gradient(.true.,eta,ix^L,idim,tmp)
!      gradeta(ix^S,idim)=tmp(ix^S)
!   enddo
!
! Otherwise figure out a compact formula, using e.g. grad(J)=Laplace B, etc.

call die('specialeta is not defined')
end
!=============================================================================
