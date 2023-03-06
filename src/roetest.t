!=============================================================================
program roetest

! Test the approximate Roe solver by checking the R*L*dwC=dwC identity
include 'vacdef.f'

double precision, dimension(ixG^T,nw):: w,wR,wroe
double precision, dimension(ixG^T)   :: jumpC,rjumpC,a
double precision:: smalla,dwCtest(nw),wleft(nw),wright(nw)
integer:: ix^L,ixI^L,ix^D,iw,il,idim,iws(niw_)
!-----------------------------------------------------------------------------

write(*,*)'Teststr?'
read(*,'(a)')teststr
oktest=index(teststr,'roetest')>=1

kr=0;^C&kr(^C,^C)=1;
do iw=1,nw
   iws(iw)=iw
enddo
iws(niw_)=nw

ixImin^D=1;ixImax^D=2;ix^L=1;ix^D=1;ixtest^D=ix^D;iwtest=1;idimtest=1

write(*,*)'nw,neqpar:',nw,neqpar
write(*,*)'Give equation parameters'
read(*,*)eqpar
write(*,*)'Give left  state wL:'
read(*,*)wleft
write(*,*)'Give right state wR:'
read(*,*)wright
do iw=1,nw
   w(ixI^S,iw) =wleft(iw)
   wR(ixI^S,iw)=wright(iw)
end do
if(oktest)write(*,*)'ixI:',ixI^L

do idim=1,ndim
   write(*,*)'idim=',idim
   call average(w,wR,ix^L,iws,idim,wroe)
   write(*,*)'Average: ',wroe(ix^D,1:nw)
   ! call checkmhd(wroe,ix^D)
   dwCtest=0d0
   do il=1,nw
      call geteigenjump(w,wR,wroe,ix^L,il,idim,smalla,a,jumpC)
      write(*,*)'il,L*(wR-wL),c:',il,jumpC(ix^D),a(ix^D)
      do iw=1,nw
         call rtimes(jumpC,wroe,ix^L,iw,il,idim,rjumpC)
         if(oktest)write(*,*)'r_il*L*(wR-wL):',rjumpC(ix^D)
         dwCtest(iw)=dwCtest(iw)+rjumpC(ix^D)
      end do
   end do     !il
   write(*,*)'R*L*(wR-wL)             wR-wL'
   do iw=1,nw
      write(*,*)dwCtest(iw),wR(ix^D,iw)-w(ix^D,iw)
   end do
end do        !idim

stop
end
!=============================================================================
subroutine checkmhd(wroe,ix^D)

! Check identites
include 'vacdef.f'

double precision, dimension(ixG^T,nw):: wroe
double precision, dimension(ixG^T)   :: cfast,cslow,afast,aslow
common /roe/ cfast,cslow,afast,aslow
integer:: ix^D
!-----------------------------------------------------------------------------

if(abs(afast(ix^D)**2+aslow(ix^D)**2-1d0)>1D-10)&
   write(*,*)'Problem af**2+as**2/=1; af,as:',&
   afast(ix^D),aslow(ix^D)

if(abs((afast(ix^D)*cfast(ix^D))**2+(aslow(ix^D)*cslow(ix^D))**2-&
   wroe(ix^D,ee_)**2)>1D-10)&
   write(*,*)'Problem (af*cf)**2+(as*cs)**2/=a**2; af,as,cf,cs,a:',&
   afast(ix^D),aslow(ix^D),cfast(ix^D),cslow(ix^D),wroe(ix^D,ee_)

write(*,*)'a,b2=cf**2+cs**2-a**2',&
   wroe(ix^D,ee_),cfast(ix^D)**2+cslow(ix^D)**2-wroe(ix^D,ee_)**2

return
end
!=============================================================================
subroutine ensurebound(dix,ixI^L,ixO^L,qt,w)

include 'vacdef.f'

integer:: dix,ixI^L,ixO^L
double precision:: qt,w(ixG^T,nw)
!-----------------------------------------------------------------------------

! Dummy subroutine

return
end
!=============================================================================
subroutine gradient(realgrad,q,ix^L,idir,gradq)

include 'vacdef.f'

logical:: realgrad
integer:: ix^L,idir
double precision:: q(ixG^T),gradq(ixG^T)
!-----------------------------------------------------------------------------

! Dummy subroutine

return
end
!=============================================================================
