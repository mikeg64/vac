!#############################################################################
! module vacmc
! MacCormack scheme
!=============================================================================
subroutine maccormack(qdt,ixII^L,ixO^L,iws,idim^LIM,qt,w)

! Advance the iws flow variables from t to t+qdt within ixO^L with 
! MacCormack scheme by solving dw/dt+dF_i(w)/dx_i=S(w) type equation. 
! For dix=0 the rightside fluxes are added, for dix=1 the leftside fluxes.

include 'vacdef.f'

double precision:: qdt,qt,w(ixG^T,nw)
integer:: ixII^L,ixO^L,iws(niw_),idim^LIM

double precision :: w1(ixG^T,nw),v(ixG^T,ndim),f(ixG^T,ndim),fC(ixG^T)
integer:: ix,dix,ixI^L,ix^L,ixC^L,ijx^L,hix^L,ixS^L,iiw,iw,idim,idir,idir^LIM
logical:: transport
!-----------------------------------------------------------------------------

oktest=index(teststr,'maccormack')>=1

! determine the value of dix based on the parity of it:
dix=it-2*(it/2)

if(oktest)write(*,*)'MacCormack it,dix,w:',it,dix,w(ixtest^D,iwtest)

! Determine needed directions for velocity and flux
if(gencoord)then
   idirmin=1
   idirmax=ndim
else
   idir^LIM=idim^LIM;
endif

! An extra layer is needed in each direction for which fluxes are added.
ixI^L=ixO^L; ixC^L=ixO^L;
do idim= idir^LIM
   ixI^L=ixI^L^LADDkr(idim,^D);
   ixCmin^D=ixCmin^D-(1-dix)*kr(idim,^D);
   ixCmax^D=ixCmax^D+dix*kr(idim,^D);
enddo
if(ixII^L^LTixI^L|.or.|.or.) call die( &
   'Error in MacCormack: Non-conforming input limits')

! Calculate velocities
do idir= idir^LIM
   call getv(w,ixI^L,idir,fC)
   v(ixI^S,idir)=fC(ixI^S)
enddo

! Initialize w1
w1(ixC^S,1:nw)=w(ixC^S,1:nw)

! Calcultae w(1) with one-sided flux eq3.69a/5.13a, and part of eq.3.69b/5.13.b
do iiw=1,iws(niw_); iw=iws(iiw)
   ! Calculate all needed fluxes for iw, use fC as a temporary variable
   do idir= idir^LIM
      call getflux(w,ixI^L,iw,idir,fC,transport)
      if(transport)then
         f(ixI^S,idir)=fC(ixI^S)+v(ixI^S,idir)*w(ixI^S,iw)
      else
         f(ixI^S,idir)=fC(ixI^S)
      endif
   enddo
   ! Add gradients of fluxes for all directions
   do idim= idim^LIM
      ijx^L=ixC^L+(1-dix)*kr(idim,^D); hix^L=ixC^L-dix*kr(idim,^D);
      ixmax^D=ijxmax^D;ixmin^D=hixmin^D;
      ixSmin^D=ixCmin^D-kr(idim,^D);ixSmax^D=ixCmax^D;
      if(gencoord)then
         fC(ix^S)= ^D&normalC(ixS^S,idim,^D)*f(ix^S,^D)+
      else
         fC(ix^S)=f(ix^S,idim)
      endif

      call addflux(qdt,ixC^L,iw,idim,fC,ijx^L,fC,hix^L,w1)

   end do    !next idim
end do       !next iw

if(oktest)write(*,*)'w(1) without source:',w1(ixtest^D,iwtest)

! Add unsplit and geometrical sources to w1 as in eq.3.69a
if(typeaxial/='slab'.and.idimmin==1) &
   call addgeometry(qdt,ixC^L,iws,w,w1)
if(sourceunsplit) call addsource2(qdt*(idimmax-idimmin+one)/ndim, &
    ixII^L,ixC^L,iws,qt,w,qt,w1)

if(oktest)write(*,*)'w(1) with    source:',w1(ixtest^D,iwtest)

! Calculate velocities
do idir= idir^LIM
   call getv(w1,ixC^L,idir,fC)
   v(ixC^S,idir)=fC(ixC^S)
enddo

! Calculate w(2) from other sided flux(1)
dix=1-dix
do iiw=1,iws(niw_); iw=iws(iiw)
   ! Calculate all needed fluxes for iw
   do idir= idir^LIM
      call getflux(w1,ixC^L,iw,idir,fC,transport)
      if(transport)then
         f(ixC^S,idir)=fC(ixC^S)+v(ixC^S,idir)*w1(ixC^S,iw)
      else
         f(ixC^S,idir)=fC(ixC^S)
      endif
   enddo
   ! Add gradients of fluxes for all directions
   do idim= idim^LIM
      ijx^L=ixO^L+(1-dix)*kr(idim,^D); hix^L=ixO^L-dix*kr(idim,^D);
      ixmax^D=ijxmax^D;ixmin^D=hixmin^D;
      ixSmin^D=ixOmin^D-kr(idim,^D);ixSmax^D=ixOmax^D;

      if(gencoord)then
         fC(ix^S)= ^D&normalC(ixS^S,idim,^D)*f(ix^S,^D)+
      else
         fC(ix^S)=f(ix^S,idim)
      endif

      call addflux(qdt,ixO^L,iw,idim,fC,ijx^L,fC,hix^L,w)

   end do    !next idim
end do       !next iw

if(oktest)write(*,*)'w(2) without source:',w(ixtest^D,iwtest)

! Add unsplit and geometrical sources to w as in eq.3.69b
if(typeaxial/='slab'.and.idimmin==1) &
    call addgeometry(qdt,ixO^L,iws,w1,w)
if(sourceunsplit)call addsource2(qdt*(idimmax-idimmin+one)/ndim, &
    ixC^L,ixO^L,iws,qt+qdt,w1,qt,w)

if(oktest)write(*,*)'w(2) with    source:',w(ixtest^D,iwtest)

! Average U(1) and U(2) (which is in w), eq3.69b or eq.5.13b
do iiw=1,iws(niw_); iw=iws(iiw)
   w(ixO^S,iw)=half*(w(ixO^S,iw)+w1(ixO^S,iw))
end do

if(oktest)write(*,*)'w with averaging:',w(ixtest^D,iwtest)

return
end
!=============================================================================
! end module vacmc
!#############################################################################
