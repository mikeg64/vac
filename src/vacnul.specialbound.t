!=============================================================================
subroutine specialbound(qt,ix^L,iw,iB,w)

! Calculates the boundary values in the iB-th boundary segment, user-defined

include 'vacdef.f'

integer:: ix^L,iw,iB
double precision:: qt,w(ixG^T,nw)
!-----------------------------------------------------------------------------

call die('Special boundary is not defined')
end
!=============================================================================
