!=============================================================================
subroutine keeppositive_rho(ix^L,w)

! Keep density positive

include 'vacdef.f'

integer::          ix^L
double precision:: w(ixG^T,nw)
logical:: toosmallr
!-----------------------------------------------------------------------------

! Overwrite smallrho if it is not yet defined
if(smallrho== -one)then
    smallrho=minval(w(ix^S,rho_))*smallrhocoeff
    {^IFMPI call mpiallreduce(smallrho,MPI_MIN)}
endif

! Check for very small density
toosmallr=any(w(ix^S,rho_)<max(zero,smallrho))

if(toosmallr)then
   nerror(toosmallr_)=nerror(toosmallr_)+1
   if(nerror(toosmallr_)==1)then
      write(*,'(a,i2,a,i7)')'Too small density (code=',toosmallr_,') at it=',it
      write(*,*)'Value < smallrho: ',minval(w(ix^S,rho_)),smallrho
!     write(*,*)'Location: ',minloc(w(ix^S,rho_)) !F77_
      {write(*,*)'Processor:',ipe ^IFMPI}
   endif
   if(smallrho>zero)w(ix^S,rho_)=max(w(ix^S,rho_),smallrho)
endif

return
end
!=============================================================================
