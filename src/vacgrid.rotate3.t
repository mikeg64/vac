! Subroutines used by VACGRID for 3D calculations

!=============================================================================
subroutine rotatew3(ix^L,idim,w)

! Rotate all vector variables of w into the normalC(:,idim,:) coord. system.
!
! This is the 3D version. The 3 base vectors n,m and l are
!
! n=(nx,ny,nz)
! m=(-ny,nx,0)/sqrt(1-nz**2)
! l=(-nx*nz, -ny*nz, 1-nz**2)/sqrt(1-nz**2)=((0,0,1)-nz*n)/sqrt(1-nz**2)
!
! for idim==1, and cyclicly shifted to the right for idim==2 and 3.
! Cyclic shifting preserves right handedness in 3D.
!
! This set of base vectors fails to exist for |nz|==1. We use 
!
! n=(nx,ny,nz)
! m=(0,-nz,ny)/sqrt(1-nx**2)
! l=(1-nx**2, -nx*ny, -nx*nz)/sqrt(1-nx**2)=((1,0,0)-nx*n)/sqrt(1-nx**2)
!
! for |nz|>=half, which surely exist since then nx<1.

include 'vacdef.f'

integer:: ix^L,idim
double precision:: w(ixG^T,1:nw)
integer:: ivector,iwv,jdim,kdim
double precision:: tmp3(ixG^T)
double precision:: norm(IXG^T,ndim)
common /gen3D/ norm
!-----------------------------------------------------------------------------

oktest=index(teststr,'rotatew')>=1
if(oktest)write(*,*)'Rotatew3 idim,wold:',idim,&
   w(ixtest^D,iwtest)

jdim=idim+1-ndim*(idim/ndim); kdim=jdim+1-ndim*(jdim/ndim)

do ivector=1,nvector
   iwv=iw_vector(ivector)
   if(oktest)write(*,*)'ivector,iwv:',ivector,iwv

   tmp(ix^S)=w(ix^S,iwv+1)
   tmp2(ix^S)=w(ix^S,iwv+2)
   tmp3(ix^S)=w(ix^S,iwv+3)

   w(ix^S,iwv+idim)=  normalC(ix^S,idim,1)*tmp(ix^S)&
                     +normalC(ix^S,idim,2)*tmp2(ix^S)&
                     +normalC(ix^S,idim,3)*tmp3(ix^S)

   where(abs(normalC(ix^S,idim,3))<half)
      w(ix^S,iwv+jdim)=norm(ix^S,idim)*&
         (-normalC(ix^S,idim,2)*tmp(ix^S)+normalC(ix^S,idim,1)*tmp2(ix^S))

      w(ix^S,iwv+kdim)=norm(ix^S,idim)*&
         (tmp3(ix^S)-normalC(ix^S,idim,3)*w(ix^S,iwv+idim))
   elsewhere
      w(ix^S,iwv+jdim)=norm(ix^S,idim)*&
         (-normalC(ix^S,idim,3)*tmp2(ix^S)+normalC(ix^S,idim,2)*tmp3(ix^S))

      w(ix^S,iwv+kdim)=norm(ix^S,idim)*&
         (tmp(ix^S)-normalC(ix^S,idim,1)*w(ix^S,iwv+idim))
   endwhere

   if(oktest.and.iwv<iwtest.and.iwv+ndim>=iwtest)then
       write(*,*)'vold:',tmp(ixtest^D),tmp2(ixtest^D),tmp3(ixtest^D)
       write(*,*)'vnew:',^D&w(ixtest^DD,iwv+^D)
   endif
enddo

if(oktest)write(*,*)'Rotatew3 iwtest,wnew:',iwtest,&
   w(ixtest^D,iwtest)

return
end

!=============================================================================
subroutine rotateback3(ix^L,idim,idir,idirf,rotC)

! Calculate the RotC(idir,idirf) rotation matrix element for the idim side.
! which is the idir component of the idirf-the base vector.
!
! This is the 3D version. The 3 base vectors n,m and l are defined in rotatew3.

! Can not use tmp,tmp2!!!

include 'vacdef.f'

double precision:: rotC(ixG^T)
integer:: ix^L,idim,idir,idirf,jdim
double precision:: norm(IXG^T,ndim)
common /gen3D/ norm
!-----------------------------------------------------------------------------

oktest=index(teststr,'rotateback')>0
if(oktest)write(*,*)'RotateBack idim,idir,idirf,ixL:',idim,idir,idirf,ix^L

jdim=idim+1-ndim*(idim/ndim)

if(idirf==idim)then
   rotC(ix^S)=  normalC(ix^S,idim,idir)
elseif(idirf==jdim)then
   select case(idir)
   case(1)
      where(abs(normalC(ix^S,idim,3))<half)
         rotC(ix^S)= -normalC(ix^S,idim,2)*norm(ix^S,idim)
      elsewhere
         rotC(ix^S)= zero
      endwhere
   case(2)
      where(abs(normalC(ix^S,idim,3))<half)
         rotC(ix^S)=  normalC(ix^S,idim,1)*norm(ix^S,idim)
      elsewhere
         rotC(ix^S)= -normalC(ix^S,idim,3)*norm(ix^S,idim)
      endwhere
   case(3)
      where(abs(normalC(ix^S,idim,3))<half)
         rotC(ix^S)= zero
      elsewhere
         rotC(ix^S)=  normalC(ix^S,idim,2)*norm(ix^S,idim)
      endwhere
   end select
else
   select case(idir)
   case(1)
      where(abs(normalC(ix^S,idim,3))<half)
         rotC(ix^S)= -normalC(ix^S,idim,1)*normalC(ix^S,idim,3)*norm(ix^S,idim)
      elsewhere
         rotC(ix^S)=  1/norm(ix^S,idim)
      endwhere
   case(2)
      where(abs(normalC(ix^S,idim,3))<half)
         rotC(ix^S)= -normalC(ix^S,idim,2)*normalC(ix^S,idim,3)*norm(ix^S,idim)
      elsewhere
         rotC(ix^S)= -normalC(ix^S,idim,1)*normalC(ix^S,idim,2)*norm(ix^S,idim)
      endwhere
   case(3)
      where(abs(normalC(ix^S,idim,3))<half)
         rotC(ix^S)=  1/norm(ix^S,idim)
      elsewhere
         rotC(ix^S)= -normalC(ix^S,idim,1)*normalC(ix^S,idim,3)*norm(ix^S,idim)
      endwhere
   end select
endif

if(oktest)write(*,*)'RotateBack rotC:',rotC(ixtest^D)


return
end
!=============================================================================
