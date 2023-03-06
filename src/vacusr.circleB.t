!#############################################################################
! module vacusr - circular boundary

! INCLUDE:vacnul.specialbound.t
INCLUDE:vacnul.specialini.t
INCLUDE:vacnul.specialsource.t
! INCLUDE:vacnul.specialio.t

!=============================================================================
subroutine specialbound(qt,ix^L,iw,iB,w)

! Reflective circular boundary for Cartesian grids.
! There are two versions, your choice depends on "teststr" in &testlist.
! Neither of these versions are perfect (that's why there are two of them).
!
! The default version conserves mass, since there is no flux through the cell
! interfaces around the circle, however this boundary is not smooth,
! which results in the formation of a boundary layer for shearing flows. 
! It is easy to change the shape of the boundary by defining the logical arrays
! "edgeminB, edgemaxB, outsideB". Note that the DEC Fortran 90 compiler does
! the where statements very inefficiently, in Fortran 77 the loops are much
! faster. This is probably a bug of the F90 compiler.
!
! "teststr='bilinear'"
! The version with bilinear interpolation uses a reflection onto a perfect
! circle, but due to the asymmetry of the grid cells relative to this imaginary
! boundary the mass conservation is only approximate. The error is
! very small for waves propagating in the radial direction, but significant
! for shearing flows. The shape of the boundary is not easy to change.
!
! To use the circular boundaries in the parfile have something like this
!
! &boundlist
!        typeB=3*'special'
!        nB=1
!        idimB=1
!        ixBmin=3,3
!        ixBmax=68,68
! &end
!
! where the full grid (including ghost cells) is assumed to be 70x70.

include 'vacdef.f'

integer:: ix^L,iw,iB
double precision:: qt,w(ixG^T,nw)
!-----------------------------------------------------------------------------

if(index(teststr,'bilinear')>=1)then
   call specialbound_bilinear(qt,ix^L,iw,iB,w)
else
   call specialbound_ragged(qt,ix^L,iw,iB,w)
endif

return
end

!=============================================================================
subroutine specialbound_ragged(qt,ix^L,iw,iB,w)

! Circular boundary implemented as reflective cell interfaces around the circle
! The size of the circle is determined by ixminB and ixmaxB in &boundlist,
! which is assumed to be symmetric around the x=0 point.
! Reflect in direction idimsplit for the cell interface closest to the circle.

include 'vacdef.f'

integer:: ix^L,iw,iB
double precision:: qt,w(ixG^T,nw)
double precision:: w0(1:nw),r(ixG^T),rcircle,rBmin,rBmax
logical:: edgeminB(ixG^T,ndim),edgemaxB(ixG^T,ndim),outsideB(ixG^T)

integer:: hx^L,jx^L,kx^L,idim

common/circle/rCircle,rBmin,rBmax
save w0,edgeminB,edgemaxB,outsideB
logical:: initialized
data initialized/.false./
!-----------------------------------------------------------------------------
if(.not.dimsplit)call die('circleB boundary works with dimsplit=T only!')
if(gencoord)call die('circleB boundary is intended for Cartesian grids!'

oktest=index(teststr,'specialbound')>=1
if(oktest)write(*,*)'Circular boundary, idimsplit',idimsplit

if(.not.initialized)then
   initialized=.true.
   w0(1:nw)=w(ixMmin^D,1:nw)
   rcircle=half*(x(ixmax^D,1)-x(ixmin^D+1,1))
   write(*,*)'Radius of ragged circle:',rcircle
   r(ixG^S)=sqrt(^D&x(ixG^S,^D)**2+)
   do idim=1,ndim
      jx^L=ix^L+kr(^D,idim);
      edgeminB(ix^S,idim)= r(jx^S)<rcircle.and.r(ix^S)>=rcircle
      edgemaxB(ix^S,idim)= r(ix^S)<rcircle.and.r(jx^S)>=rcircle
   enddo
   outsideB(ixG^S)= r(ixG^S)>rcircle
endif

! The direction of the NEXT sweep will be idim based on idimsplit and it
if(idimsplit==0)then
  ! The initial getboundary call is before any advect1, idim is based on it
  if(it==2*(it/2))then
     idim=1
  else
     idim=ndim
  endif
else
  ! The getboundary calls after advect1, if last one, prepare for next sweep
  !!! We ignored ensure_bound here!!!
  if(istep<nstep)then
     idim=idimsplit
  else if(it==2*(it/2))then
     idim=min(idimsplit+1,ndim)
  else
     idim=max(idimsplit-1,1)
  endif
endif

hx^L=ix^L-kr(^D,idim);
jx^L=ix^L+kr(^D,idim);
kx^L=jx^L+kr(^D,idim);

! Outside circle use w0
where(outsideB(ixG^S)) w(ixG^S,iw)=w0(iw)

! Calculate two ghost cells from the inside using symmetry or anti-symmetry
! First do right/top semi circle then left/bottom semi circle

if(iw /= m0_+idim)then
   ! Symmetry
   where(edgemaxB(ix^S,idim))
      w(jx^S,iw)=w(ix^S,iw)
      w(kx^S,iw)=w(hx^S,iw)
   end where
   where(edgeminB(ix^S,idim))
      w(ix^S,iw)=w(jx^S,iw)
      w(hx^S,iw)=w(kx^S,iw)
   end where
else
   ! Anti-symmetry
   where(edgemaxB(ix^S,idim))
      w(jx^S,iw)=-w(ix^S,iw)
      w(kx^S,iw)=-w(hx^S,iw)
   end where
   where(edgeminB(ix^S,idim))
      w(ix^S,iw)=-w(jx^S,iw)
      w(hx^S,iw)=-w(kx^S,iw)
   end where
endif

return
end

!=============================================================================
subroutine specialbound_bilinear(qt,ix^L,iw,iB,w)

! Reflection on a circle, done for all variables at the same time
! Uses bilinear interpolation for reflected ghost cells

include 'vacdef.f'

integer:: ix^L,iw,iB,idist
double precision:: qt,w(ixG^T,nw)

double precision:: r(ixG^T),rCircle,rBmin,rBmax,rB
double precision:: x^D,xR^D,dx^D,coeff^D(0:1),w0(nw),wB(nw),m_rad
integer:: ix^D,ixR^D,dixR^D
logical:: oktest1,initialized
data initialized/.false./

save r,dx^D,w0
common/circle/rCircle,rBmin,rBmax
!-----------------------------------------------------------------------------

if(iw>1)return
if(gencoord)call die('CircleB boundary is intended for Cartesian grids!'

oktest=index(teststr,'specialbound')>=1
if(oktest)write(*,*)'Circular boundary with bilinear interpolations'

if(.not.initialized)then
   initialized=.true.

   r(ixG^S)=sqrt(^D&x(ixG^S,^D)**2+)
   dx^D=dx(ixMmin^DD,^D);

   rBmin=max(^D&x(ixGmax^DD,^D)-dx^D*1.0001,zero)
   rBmax=rBmin+max(2.0001*dx^D,zero)
   rCircle=rBmin-1.0001*sqrt(dx^D**2+)
   write(*,*)'Radius of bilinear circle:',rcircle

   w0(1:nw)=w(ixMmin^D,1:nw)

   if(oktest)write(*,*)'xmax^D,dx^D:',^D&x(ixGmax^DD,^D),dx^D
   if(oktest)write(*,*)'rBmin,rBmax:',rBmin,rBmax
endif

{do ix^D=ixGmin^D,ixGmax^D\}
   rB=r(ix^D)
   oktest1=oktest.and.(ixtest^D==ix^D|.and.)

   if(oktest1)write(*,*)'rB:',rB

   if(rB>rBmax)then
      w(ix^D,1:nw)=w0(1:nw)
      if(oktest1)write(*,*)'Outside w:',w(ixtest^D,iwtest)
   else if(rB>rBmin) then
      ! Reflect ghost cell center into inside, find closeby center
      x^D=x(ix^DD,^D);
      xR^D=(2*rCircle/rB-1)*x^D;
      ixR^D=(xR^D-x(ixGmin^DD,^D))/dx^D+ixGmin^D;

      if(oktest1)write(*,*)'x,xR,ixR:',x^D,xR^D,ixR^D

      ! Calculate bilinear interpolation coefficients
      coeff^D(1)=(xR^D-x(ixR^DD,^D))/dx^D;
      coeff^D(0)=1-coeff^D(1);

      if(oktest1)write(*,*)'coeff0:',coeff^D(0)
      if(oktest1)write(*,*)'coeff1:',coeff^D(1)

      ! Interpolate all variables into wB
      wB(1:nw)=zero
      {do dixR^D=0,1\}
         wB(1:nw)=wB(1:nw)+(coeff^D(dixR^D)*)*w(ixR^D+dixR^D,1:nw)
         if(oktest1)write(*,*)'dixR,w:',dixR^D,w(ixR^D+dixR^D,iwtest)
         if(it==itmin)then
           if(r(ixR^D+dixR^D)>one)then
              write(*,*)'Problem: ix,r:',ix^D,rB
              write(*,*)'Problem: ixR,dixR,r:',ixR^D,dixR^D,r(ixR^D+dixR^D)
              call die('Error in circleB')
           endif
         endif
      {enddo^D&\}

      if(oktest1)write(*,*)'w interp:',wB(iwtest)

      ! Reflect momentum: m = m - 2*(x.m)/(x.x)*x
      m_rad=two*(x^D*wB(m^D_)+)/(x^D**2+)
      ^D&wB(m^D_)=wB(m^D_)-m_rad*x^D;
      ! Put wB into ghost cell
      w(ix^D,1:nw)=wB(1:nw)
      if(oktest1)write(*,*)'Inside w:',w(ixtest^D,iwtest)
   endif
{enddo^D&\}

return
end

!=============================================================================
subroutine savefilelog_special(qunit,w,ix^L)

! Save averages for cells within the circle only
! Use typefilelog='special' in &filelist.

include 'vacdef.f'

integer:: qunit,ix^L
double precision:: w(ixG^T,nw)

logical:: fileopen
integer:: iw
double precision:: dvol(ixG^T),vol,wmean(nw),rCircle,rBmin,rBmax

save dvol,vol
common/circle/rCircle,rBmin,rBmax
!-----------------------------------------------------------------------------

inquire(qunit,opened=fileopen)
if(.not.fileopen)then
    open(qunit,file=filename(filelog_),status='unknown')
    write(qunit,'(a)')fileheadout
    write(qunit,'(a15,a64)')'it   t   dt   ',wnames

    ! Do not integrate for cells outside of the circle
    dvol(ixG^S)=dvolume(ixG^S)
    where((^D&x(ixG^S,^D)**2+)>rCircle**2) dvol(ixG^S)=zero
    vol=sum(dvol(ixG^S))
endif

do iw=1,nw 
   wmean(iw)=sum(dvol(ix^S)*w(ix^S,iw))/vol
end do

write(qunit,'(i7,100(1pe13.5))')it,t,dt,wmean
call flushunit(qunit)

return
end
!=============================================================================
subroutine readfileini_special(w)
double precision:: w(*)
call die('CircleB:Special readfileini is not defined')
end
!=============================================================================
subroutine savefileout_special(qunit,w,ix^L)
integer:: qunit,ix^L
double precision:: w(*)
call die('CircleB:Special savefileout is not defined')
end
!=============================================================================
! end module vacusr - onewave
!#############################################################################
