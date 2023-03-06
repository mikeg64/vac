!=============================================================================
subroutine process(count,idim^LIM,w)

! Process w before it is advected in directions idim^LIM, or before save
! count=1 and 2 for first and second (half step) processing during advection
! count=ifile+2 for saving results into the file indexed by ifile

include 'vacdef.f'

integer:: count,idim^LIM
double precision:: w(ixG^T,nw)

!-----------------------------------------------------------------------------

return
end
!=============================================================================
