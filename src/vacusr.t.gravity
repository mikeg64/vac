!##############################################################################
! module vacusr - gravity

INCLUDE:vacusr.gravity.t

INCLUDE:vacnul.specialbound.t
INCLUDE:vacnul.specialini.t
!INCLUDE:vacnul.specialsource.t
INCLUDE:vacnul.specialio.t

!=============================================================================
! There are 3 subroutines in this file:
!
! specialsource -- for sources other than resistivity
! getdt_special -- for time step conditions other than CFL or resistivity
! specialeta    -- for non-constant resistivity with eqpar(eta_)<zero
!
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)

include 'vacdef.f'

integer:: ixI^L,ixO^L,iws(niw_)
double precision:: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)
!-----------------------------------------------------------------------------

call addsource_grav(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)

return
end

!=============================================================================
subroutine getdt_special(w,ix^L)

include 'vacdef.f'

double precision:: w(ixG^T,nw)
integer:: ix^L
!-----------------------------------------------------------------------------

call getdt_grav(w,ix^L)

return
end

!=============================================================================
subroutine specialeta(w,ix^L,idirmin)

include 'vacdef.f'

double precision:: w(ixG^T,nw)
integer:: ix^L,idirmin
!-----------------------------------------------------------------------------

call die('vacusr.gravity: specialeta is not defined')
end
!=============================================================================

! end module vacusr - gravity
!##############################################################################
