!=============================================================================
subroutine specialini(ix^L,w)

! Initialize w for VACINI, user-defined

include 'vacdef.f'

integer:: ix^L
double precision:: w(ixG^T,nw)
!-----------------------------------------------------------------------------

call die('Special initial condition is not defined')
end
!=============================================================================
