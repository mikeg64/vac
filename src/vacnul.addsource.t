!=============================================================================
subroutine addsource(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)

! w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO

include 'vacdef.f'

integer::          ixI^L,ixO^L,iws(niw_)
double precision:: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)
!-----------------------------------------------------------------------------

return
end
!=============================================================================
