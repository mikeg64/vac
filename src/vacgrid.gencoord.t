! SUBROUTINES used by VACGRID in generalized coordinates
!=============================================================================
{INCLUDE:vacgrid.rotate3.t ^IFTHREED}
!=============================================================================
subroutine gridsetup2

! General grid. Calculate x at the boundaries.
! Determine often needed combinations of x, such as 
! dx, dvolume, surfaceC, normalC
!
! qx           - x with an extended index range for calculation of dx

include 'vacdef.f'

integer:: ix^L,hx^L,jx^L,h1x^L,h2x^L,h3x^L,j1x^L,j2x^L,j3x^L
integer:: ix,ix^D,dix^D,idim,jdim,kdim,idir,jdir,kdir
double precision:: qx(IXG^LL^LADD1:,ndim),xC(IXG^LL^LADD1:,ndim)
{^IFTHREED 
double precision:: norm(IXG^T,ndim)
common /gen3D/ norm
}
!-----------------------------------------------------------------------------

oktest=index(teststr,'gridsetup')>=1
if(oktest)write(*,*)'GridSetup2 x:',&
   ^D&x(ixtest^DD,^D)

{^IFMPI
! Fill in ghost cells for x from MPI neighbors
call mpibound(ndim,x)
}
qx(ixG^LL^LADD1:,1:ndim)=zero
qx(ixG^S,1:ndim) = x(ixG^S,1:ndim)
do idim=1,ndim
   ! First do the upper layers then the lower layers
   call extrapolate2(idim,.true. ,qx)
   call extrapolate2(idim,.false.,qx)
end do

if(oktest)write(*,*)'qx:',&
   ^D&qx(ixtest^DD,^D)

! Approximate dx(idim) by half the distance of x(jx) and x(hx) with shifts
! in direction idim. This is the cell-width in the local idim-th direction
do idim=1,ndim
   hx^L=ixG^L-kr(idim,^D); jx^L=ixG^L+kr(idim,^D);
   dx(ixG^S,idim)=half*sqrt(^D&(qx(jx^S,^D)-qx(hx^S,^D))**2+)
enddo

if(oktest)write(*,*)'dx:',^D&dx(ixtest^DD,^D)

! Calculate position of cell corners, average the 2**ndim closest cell centers
ixmax^D=ixGmax^D; ixmin^D=ixGmin^D-1;

xC(ix^S,1:ndim)=zero
{do dix^D=0,1\}
   jx^L=ix^L+dix^D;
   xC(ix^S,1:ndim)=xC(ix^S,1:ndim)+qx(jx^S,1:ndim)
{^D&enddo\}
xC(ix^S,1:ndim)=xC(ix^S,1:ndim)/2**ndim

if(oktest)write(*,*)'Corner xC:',^D&xC(ixtest^DD,^D)

! For EMMA_D we have the option to save the cell corners
if(index(teststr,'savecorners')>=1)then
   write(*,*)'Writing grid corners into data/GRID!!!'
   open(unitini-1,file='data/GRID',status='unknown')
   ixmax^D=ixMmax^D; ixmin^D=ixMmin^D-1;
   {do ix^D= ixmin^D,ixmax^D \}
      write(unitini-1,*)xC(ix^D,1:ndim)
   {enddo^D& \}
   close(unitini-1)
endif

! Calculate volume and normal vectors of cells 
{^IFONED call die('No generalized coordinates in 1D')}
{^IFTWOD
    ! This is the exact volume of the cell defined by four xC corners
    ix^L=ixG^L; hx^L=ix^L-1; h^Dx^L=ix^L-kr(^DD,^D);
    dvolume(ix^S)=((xC(ix^S,1)-xC(hx^S,1))*(xC(h1x^S,2)-xC(h2x^S,2))-&
                   (xC(ix^S,2)-xC(hx^S,2))*(xC(h1x^S,1)-xC(h2x^S,1)))/2

    ! Calculate the normal vectors for each side and each component
    ixmax^D=ixGmax^D-kr(^D,1); h2x^L=ix^L-kr(^D,2);
    normalC(ix^S,1,1)=   xC(ix^S,2)-xC(h2x^S,2)
    normalC(ix^S,1,2)=  -xC(ix^S,1)+xC(h2x^S,1)
    ixmax^D=ixGmax^D-kr(^D,2); h1x^L=ix^L-kr(^D,1);
    normalC(ix^S,2,1)=  -xC(ix^S,2)+xC(h1x^S,2)
    normalC(ix^S,2,2)=   xC(ix^S,1)-xC(h1x^S,1)
}
{^IFTHREED
    if(polargrid)then
       ! volume=area(r,z)*r*dphi, area(r,z)= abs(diag(hh,ii) x diag(hi,ih))/2
       ix^L=ixG^L; hx^L=ix^L-kr(^D,r_)-kr(^D,zz_);
       h1x^L=ix^L-kr(^D,r_); h2x^L=ix^L-kr(^D,zz_);
       dvolume(ix^S)=&
         abs((xC(ix^S,r_) -xC(hx^S,r_)) *(xC(h1x^S,zz_)-xC(h2x^S,zz_))-&
             (xC(ix^S,zz_)-xC(hx^S,zz_))*(xC(h1x^S,r_) -xC(h2x^S,r_)))/2&
         *dx(ix^S,pphi_)
    else
       ! This is an approximate volume calculated from the determinant of the 
       ! Jacobian J=d3x.(d1x X d2x) where dIx is the difference taken in 
       ! direction I
       ix^L=ixG^L; h^Dx^L=ix^L-kr(^DD,^D); j^Dx^L=ix^L+kr(^DD,^D);
       dvolume(ix^S)=zero
       do idim=1,ndim
          jdim=idim+1-ndim*(idim/ndim);kdim=jdim+1-ndim*(jdim/ndim)
          dvolume(ix^S)=dvolume(ix^S)+(qx(j3x^S,kdim)-qx(h3x^S,kdim))*(&
          (qx(j1x^S,idim)-qx(h1x^S,idim))*(qx(j2x^S,jdim)-qx(h2x^S,jdim))-&
          (qx(j1x^S,jdim)-qx(h1x^S,jdim))*(qx(j2x^S,idim)-qx(h2x^S,idim)))
       enddo
       dvolume(ix^S)=dvolume(ix^S)/8
    endif

    do idim=1,ndim
       ! Calculate the other two directions
       jdim=idim+1-ndim*(idim/ndim);kdim=jdim+1-ndim*(jdim/ndim)
       ! The normal vector is face centered in the idim direction
       ixmax^D=ixGmax^D-kr(^D,idim); 
       ! For the vectorproduct of the diagonals of the face
       ! we need differences in the jdim and kdim directions
       hx^L=ix^L-1+kr(^D,idim); h1x^L=ix^L-kr(^D,jdim); h2x^L=ix^L-kr(^D,kdim);
       ! The normal vector components for the side
       do idir=1,ndim
          jdir=idir+1-ndim*(idir/ndim);kdir=jdir+1-ndim*(jdir/ndim)
          normalC(ix^S,idim,idir)=&
          ((xC(ix^S,jdir)-xC(hx^S,jdir))*(xC(h1x^S,kdir)-xC(h2x^S,kdir))-&
           (xC(ix^S,kdir)-xC(hx^S,kdir))*(xC(h1x^S,jdir)-xC(h2x^S,jdir)))/2
       end do
    end do
}

if(minval(dvolume(ixM^S))<zero) &
     call die('left handed coordinate system!!!')
dvolume(ixG^S)=abs(dvolume(ixG^S))

{^IFMPI
! Correct volumes in 2nd ghost cells from neighboring processors
! This quantity is used in TVD scheme
call mpibound(1,dvolume)
}

if(oktest)write(*,*)'dvolume:',dvolume(ixtest^D)

! The area of the sides are the lengths of the normal vector
! The normal vectors are normalized to 1
do idim=1,ndim
   ixmax^D=ixGmax^D-kr(^D,idim);
   surfaceC(ix^S,idim)=sqrt(^D&normalC(ix^S,idim,^D)**2+)
   where(surfaceC(ix^S,idim)<smalldouble)
     ^D&normalC(ix^S,idim,^D)=one/sqrt(ndim+zero);
   elsewhere
     ^D&normalC(ix^S,idim,^D)=normalC(ix^S,idim,^D)/surfaceC(ix^S,idim);
   endwhere

   if(oktest)write(*,*)'idim Surface:',idim,surfaceC(ixtest^D,idim)
   if(oktest)write(*,*)'Normal:',^D&normalC(ixtest^DD,idim,^D)

   ! In 3D we will often need 1/sqrt(1-normal(3)**2) or 1/sqrt(1-normal(1)**2)
   {^IFTHREED
      where(abs(normalC(ix^S,idim,3))<half)
         norm(ix^S,idim)=1/sqrt(1-normalC(ix^S,idim,3)**2)
      elsewhere
         norm(ix^S,idim)=1/sqrt(1-normalC(ix^S,idim,1)**2)
      endwhere
   }
end do

if(typeaxial=='cylinder')then
   ! Multiply cell surfaces by the averaged radial coordinate xC(..,r_)
   ! except for the surface in the phi_ direction in 3D polar grid
   ! The averaging is done for the 2 vortices in the R-Z plane

   h1x^L=ixG^L-kr(^D,r_); h2x^L=ixG^L-kr(^D,zz_);

   surfaceC(ixG^S,r_) =surfaceC(ixG^S,r_) *half*abs(xC(ixG^S,r_)+xC(h2x^S,r_))
   surfaceC(ixG^S,zz_)=surfaceC(ixG^S,zz_)*half*abs(xC(ixG^S,r_)+xC(h1x^S,r_))

   if(oktest)write(*,*)'Multiplied surfaces:',&
      ^D&surfaceC(ixtest^DD,^D)

   ! Multiply cell volumes by the radial coordinate qx(..,1)
   dvolume(ixG^S)= dvolume(ixG^S)*abs(qx(ixG^S,r_))

   if(oktest)write(*,*)'Multiplied volume:',&
      dvolume(ixtest^D)

   ! For polar grid dx_phi=r*dphi. qphi is used to avoid compiler warnings
   if(polargrid)then
      dx(ixG^S,pphi_)=dx(ixG^S,pphi_)*abs(qx(ixG^S,r_))
      if(oktest)write(*,*)'Multiplied dx_phi:',&
         dx(ixtest^D,pphi_)
   endif

endif

x(ixG^S,1:ndim)=qx(ixG^S,1:ndim)

{^IFTHREED if(index(teststr,'correctgrid')>0)call correctgrid}

if(index(teststr,'savegrid')>0)call savegrid

volume=sum(dvolume(ixM^S))
{^IFMPI call mpiallreduce(volume,MPI_SUM)}

if(oktest)write(*,*)'volume:',volume

if(oktest)write(*,*)'Finish GridSetup2'

return 
end

!=============================================================================
subroutine savegrid

! Save surfaces and volumes into grid.dat for correctgrid

include 'vacdef.f'

integer:: qunit
!-----------------------------------------------------------------------------

print *,'Writing grid.dat!!!'

qunit=unitini-1
open(qunit,file='grid.dat',form='unformatted',status='unknown')
write(qunit) ixGmax^D
write(qunit)normalC(ixG^S,1:ndim,1:ndim)
write(qunit)surfaceC(ixG^S,1:ndim)
write(qunit)dvolume(ixG^S)
close(qunit)

return
end

{^IFTHREED
!=============================================================================
subroutine correctgrid

! Correct dvolume, surface to be consistent with 2D poloidal axisymmetry
!!! Assuming symmetry axis to be the Y axis

include 'vacdef.f'

double precision:: norm2D(ixGlo1:ixGhi1,ixGlo2:ixGhi2,2,2)
double precision:: surf2D(ixGlo1:ixGhi1,ixGlo2:ixGhi2,2)
double precision:: dvol2D(ixGlo1:ixGhi1,ixGlo2:ixGhi2),qphi,dphi,coef
integer:: qunit,qn1,qn2,iphi,idim
double precision:: norm(IXG^T,ndim)
common /gen3D/ norm
!-----------------------------------------------------------------------------

print *,'Correcting surfaces and volumes!'

! Calculate dphi from scalar product of 333 and 334 vectors
dphi=abs((x(3,3,3,1)*x(3,3,4,1)+x(3,3,3,3)*x(3,3,4,3)))&
     /sqrt(x(3,3,3,1)**2+x(3,3,3,3)**2)/sqrt(x(3,3,4,1)**2+x(3,3,4,3)**2)

dphi=acos(dphi)

print *,'dphi:',dphi,8*atan(one)/nx3

print *,'Reading grid.dat!!!'

qunit=unitini-1
open(qunit,file='grid.dat',form='unformatted',status='old')
read(qunit)qn1,qn2
if(qn1/=ixGmax1.or.qn2/=ixGmax2)call die('incorrect 2D grid size')
read(qunit)norm2D(ixG^LIM1:,ixG^LIM2:,1:2,1:2)
read(qunit)surf2D(ixG^LIM1:,ixG^LIM2:,1:2)
read(qunit)dvol2D(ixG^LIM1:,ixG^LIM2:)
close(qunit)

! Coefficietn between the axisymmetric 2D and full 3D volumes and surfaces
coef=two*sin(half*dphi)

print *,'Old normals:',normalC(ixtest^D,1:3,1:3)
print *,'Old surface:',surfaceC(ixtest^D,1:3)
print *,'Old volume :',dvolume(ixtest^D)

do iphi=1,ixGmax3
   qphi=atan2(x(3,3,iphi,3),x(3,3,iphi,1))

   ! Y components (parallel to axis)
   normalC(ixG^LIM1:,ixG^LIM2:,iphi,1,2)=norm2D(ixG^LIM1:,ixG^LIM2:,1,2)
   normalC(ixG^LIM1:,ixG^LIM2:,iphi,2,2)=norm2D(ixG^LIM1:,ixG^LIM2:,2,2)

   ! X and Z components
   normalC(ixG^LIM1:,ixG^LIM2:,iphi,1,1)=norm2D(ixG^LIM1:,ixG^LIM2:,1,1)&
       *cos(qphi)
   normalC(ixG^LIM1:,ixG^LIM2:,iphi,1,3)=norm2D(ixG^LIM1:,ixG^LIM2:,1,1)&
       *sin(qphi)
   normalC(ixG^LIM1:,ixG^LIM2:,iphi,2,1)=norm2D(ixG^LIM1:,ixG^LIM2:,2,1)&
       *cos(qphi)
   normalC(ixG^LIM1:,ixG^LIM2:,iphi,2,3)=norm2D(ixG^LIM1:,ixG^LIM2:,2,1)&
       *sin(qphi)

   surfaceC(ixG^LIM1:,ixG^LIM2:,iphi,1)=surf2D(ixG^LIM1:,ixG^LIM2:,1)*coef
   surfaceC(ixG^LIM1:,ixG^LIM2:,iphi,2)=surf2D(ixG^LIM1:,ixG^LIM2:,2)*coef
   dvolume(ixG^LIM1:,ixG^LIM2:,iphi)   =dvol2D(ixG^LIM1:,ixG^LIM2:)*coef
   surfaceC(ixG^LIM1:,ixG^LIM2:,iphi,3)=dvol2D(ixG^LIM1:,ixG^LIM2:)&
     /sqrt(x(ixG^LIM1:,ixG^LIM2:,iphi,1)**2+x(ixG^LIM1:,ixG^LIM2:,iphi,3)**2)
enddo

do idim=1,3
   where(abs(normalC(ixG^S,idim,3))<half)
      norm(ixG^S,idim)=1/sqrt(1-normalC(ixG^S,idim,3)**2)
   elsewhere
      norm(ixG^S,idim)=1/sqrt(1-normalC(ixG^S,idim,1)**2)
   endwhere
enddo

print *,'New normals:',normalC(ixtest^D,1:3,1:3)
print *,'New surface:',surfaceC(ixtest^D,1:3)
print *,'New volume :',dvolume(ixtest^D)

return
end
}

!=============================================================================
subroutine extrapolate2(idim,upper,qx)

! General grid. Calculate qx in the boundary layer in direction idim,
! side determined by upper, by circularly extrapolating x from
! the touching edge (xe) and the two inner edges (xf and xg).
! Circular extrapolation means that a circle is drawn trough xe,xf,xg and
! the boundary points are placed on this circle with spacing xe-xf.

include 'vacdef.f'

double precision:: qx(IXG^LL^LADD1:,ndim)
integer:: idim
logical:: upper
integer:: ix^L,ix^D,ixe,ixf,ixg,dix,idir
integer:: idm,ix^LIM(ndim),ixG^LIM(ndim),ixM^LIM(ndim)

! Vectors, squared lengths and dot products
double precision:: a(ndim),b(ndim),c(ndim),d(ndim),a2,b2,c2,d2,ab,ad,bd
! Linear combination coefficients
double precision:: alpha,beta
!-----------------------------------------------------------------------------

oktest=index(teststr,'extrapolate')>=1
if(oktest)write(*,*)'Extrapolate2 idim,upper:',idim,upper

^D&ixG^LIM(^D)=ixG^DL;
^D&ixM^LIM(^D)=ixM^DL;

select case(idim)
{case(^D)
   ! Setup limits for boundary with ixmin closest to mesh, ixmax further away
   ! For directions already done use the extra layers ixGmin-1..ixGmax+1
   if(upper)then
      do idm=1,idim
         ixmin(idm)=ixGmin(idm)-1; ixmax(idm)=ixGmax(idm)+1
      end do
      do idm=idim+1,ndim
         ixmin(idm)=ixMmin(idm); ixmax(idm)=ixMmax(idm)
         if(fullgridini .or.mpilowerB(idm)^IFMPI )ixmin(idm)=ixGmin(idm)
         if(fullgridini .or.mpiupperB(idm)^IFMPI )ixmax(idm)=ixGmax(idm)
      enddo
      ixmin(idim)=ixMmax(idim)+1
      if(fullgridini .or.mpiupperB(idim)^IFMPI )ixmin(idim)=ixGmax(idim)+1
      dix=1
   else
      do idm=1,idim
         ixmin(idm)=ixGmax(idm)+1; ixmax(idm)=ixGmin(idm)-1;
      end do
      do idm=idim+1,ndim
         ixmin(idm)=ixMmax(idm); ixmax(idm)=ixMmin(idm)
         if(fullgridini .or.mpilowerB(idm)^IFMPI )ixmax(idm)=ixGmin(idm)
         if(fullgridini .or.mpiupperB(idm)^IFMPI )ixmin(idm)=ixGmax(idm)
      end do
      ixmin(idim)=ixMmin(idim)-1
      if(fullgridini .or.mpilowerB(idim)^IFMPI )ixmin(idim)=ixGmin(idim)-1
      dix=-1
   endif
   ixmax^DD=ixmax(^DD);ixmin^DD=ixmin(^DD);

   if(oktest)write(*,*)'ixmin,ixmax:',ixmin^DD,ixmax^DD

   ! Indices for edge, inner and inner by two layers
   ixe=ixmin^D-dix; ixf=ixe-dix; ixg=ixf-dix
   {do ix^DD=ixmin^DD,ixmax^DD,dix;}
      ! a=xe->xf,b=xe->xg,d=xg->xf
      a2=zero; b2=zero; c2=zero; d2=zero; ab=zero; ad=zero; bd=zero;
      do idir=1,ndim
         a(idir)=qx(ixe^D%ix^DD,idir)-qx(ixf^D%ix^DD,idir)
         b(idir)=qx(ixe^D%ix^DD,idir)-qx(ixg^D%ix^DD,idir)
         d(idir)=a(idir)-b(idir)
         a2=a2+a(idir)**2
         b2=b2+b(idir)**2
         d2=d2+d(idir)**2
         ab=ab+a(idir)*b(idir)
         ad=ad+a(idir)*d(idir)
         bd=bd+b(idir)*d(idir)
         ! Distance to boundary point (xe->x)**2:=(xf->xprev)**2
         c(idir)=qx(ixf^D%ix^DD,idir)-qx(ix^D-dix^D%ix^DD,idir)
         c2=c2+c(idir)**2
      enddo
      ! Linear combination coefficients
      alpha=c2*(ad-sqrt(a2*b2*d2/c2-a2*b2+ab*ab))/a2/d2
      beta=(c2-alpha*a2)/b2
      do idir=1,ndim
         ! c=x->xe
         c(idir)=alpha*a(idir)+beta*b(idir)
         qx(ix^DD,idir)=qx(ixe^D%ix^DD,idir)-c(idir)
      enddo
   ^DD&enddo; 
\}
end select

if(oktest)write(*,*)'Gridsetup2 final qx:',^D&qx(ixtest^DD,^D)

return
end

!=============================================================================
subroutine rotatew(ix^L,idim,w)

! Rotate all vector variables of w into the normalC(.,idim,.) coord. system

include 'vacdef.f'

integer:: ix^L,idim
double precision:: w(ixG^T,1:nw)
!-----------------------------------------------------------------------------

{^IFTWOD   call rotatew2(ix^L,idim,w)}
{^IFTHREED 
if(polargrid.and.angmomfix)then
   ! Do 2D rotations in the R-Z plane
   if(idim/=phi_)call rotatew2(ix^L,idim,w)
else
   call rotatew3(ix^L,idim,w)
endif
}

return
end

!=============================================================================
subroutine rotatew2(ix^L,idim,w)

! Rotate all vector variables of w into the normalC(:,idim,:)  coord. system.
!
! This is the 2D version. The 2 base vectors are 
! ( nx, ny) ,for idim==1 and ( ny,-nx) for idim==2.  
! (-ny, nx)                  ( nx, ny) 

include 'vacdef.f'

integer:: ix^L,idim
double precision:: w(ixG^T,1:nw)
integer:: ivector,iw,jw,jdim
!-----------------------------------------------------------------------------

oktest=index(teststr,'rotatew')>=1
if(oktest)write(*,*)'Rotatew2 idim,wold:',idim,&
   w(ixtest^D,iwtest)

jdim=3-idim
do ivector=1,nvector
   iw=iw_vector(ivector)+idim; jw=iw_vector(ivector)+jdim
   if(oktest)write(*,*)'ivector,iw,jw:',ivector,iw,jw

   tmp(ix^S)=w(ix^S,iw)
   tmp2(ix^S)=w(ix^S,jw)
   w(ix^S,iw)=normalC(ix^S,idim,idim)*tmp(ix^S)&
             +normalC(ix^S,idim,jdim)*tmp2(ix^S)
   w(ix^S,jw)=normalC(ix^S,idim,idim)*tmp2(ix^S)&
             -normalC(ix^S,idim,jdim)*tmp(ix^S)

   if(oktest.and.(iw==iwtest.or.jw==iwtest))&
       write(*,*)'tmp,tmp2,wnewi,wnewj:',&
       tmp(ixtest^D),tmp2(ixtest^D),w(ixtest^D,iw),w(ixtest^D,jw)
enddo

if(oktest)write(*,*)'Rotatew2 iwtest,wnew:',iwtest,&
   w(ixtest^D,iwtest)

return
end

!=============================================================================
subroutine rotateback(ix^L,idim,idir,idirf,rotC)

! Calculate the RotC(idir,idirf) rotation matrix element for the idim side
! which is the idir component of the idirf-th base vector.

! Can not use tmp,tmp2!!!

include 'vacdef.f'

double precision:: rotC(ixG^T)
integer:: ix^L,idim,idir,idirf
!-----------------------------------------------------------------------------

{^IFTWOD   call rotateback2(ix^L,idim,idir,idirf,rotC)}
{^IFTHREED 
if(polargrid.and.angmomfix)then
   if(idim/=phi_.and.idir/=phi_)then
      ! Do 2D rotations in the R-Z plane
      call rotateback2(ix^L,idim,idir,idirf,rotC)
   else
      ! Identity matrix for the interface normal to phi
      rotC(ix^S)=kr(idir,idirf)
   endif
else
   call rotateback3(ix^L,idim,idir,idirf,rotC)
endif
}

return
end

!=============================================================================
subroutine rotateback2(ix^L,idim,idir,idirf,rotC)

! Calculate the RotC(idir,idirf) rotation matrix element for the idim side
! which is the idir component of the idirf-th base vector.
!
! This is the 2D version. The 2 base vectors are 
! b1=(nx, ny), b2=(-ny, nx) for idim==1
! b1=(ny,-nx), b2=( nx, ny) for idim==2

! Can not use tmp,tmp2!!!

include 'vacdef.f'

double precision:: rotC(ixG^T)
integer:: ix^L,idim,idir,idirf,jdir
!-----------------------------------------------------------------------------

jdir=3-idir

! When idirf==idim the base vector is normalC.
! When idirf/=idim the base vector is orthogonal to normalC, swapped components
if(idirf==idim)then
   rotC(ix^S)=  normalC(ix^S,idim,idir)
elseif(idir/=idim)then
   rotC(ix^S)=  normalC(ix^S,idim,jdir)
else
   rotC(ix^S)= -normalC(ix^S,idim,jdir)
endif

return
end

!=============================================================================
subroutine addflux_rotate(qdt,ix^L,iw,idim,fRC,ixR^L,fLC,ixL^L,wnew)

! Used in generalized coordinates when the fRC,fLC fluxes of a vector 
! variable are in the local coordinate system of the idim-th side, 
! while the variables in wnew are in the normal coordinate system. 
! For idir=1,ndim the contributions of fR and fL to variables in wnew are
!
! wnew(vect_idir)=wnew(vect_idir)-
!   (ProjRC(idir,idirf)*fRC-ProjLC(idir,idirf)*fLC)/dvolume
!
! where idirf is the direction of the flux, determined by the iw parameter.
!
! For qdt>=zero, the fluxes are physical and the ProjC contains qdt & surfaceC,
! for qdt< zero, ie. for TVD type dissip. flux, it does not & sign is opposite

include 'vacdef.f'

double precision:: qdt,fRC(ixG^T),fLC(ixG^T),wnew(ixG^T,nw)

integer:: ix^L,ixR^L,ixL^L,iw,idim

integer:: hx^L,ixC^L,iwv,idir,idirf
!-----------------------------------------------------------------------------

oktest= index(teststr,'addflux')>=1 .and. iw==iwtest
if(oktest)write(*,*)'AddFluxRotate wold:',wnew(ixtest^D,iwtest)

!SHIFT
!S ixRmin1=ixmin1;
!SHIFT MORE
!S ixLmin1=ixmin1-1
!SHIFT MORE
hx^L=ix^L-kr(idim,^D);
ixCmax^D=ixmax^D; ixCmin^D=hxmin^D;

iwv=vectoriw(iw)
idirf=iw-iwv

if(oktest)write(*,*)'idim,iw,iwv,idirf,ixL:',idim,iw,iwv,idirf

! Set flux to zero if required by boundary condition
! This is done in the local coordinate system, 
! this is how the boundary type is given
if(nofluxB(iw,idim)) call setnoflux(iw,idim,ix^L,fRC,ixR^L,fLC,ixL^L)

do idir=1,ndim

   call rotateback(ixC^L,idim,idir,idirf,tmp)

   if(oktest)write(*,*)'idir,projC(ix,hx):',idir,&
      tmp(ixtest^D),tmp(ixtest^D-kr(idim,^D))

   !SHIFT BEGIN
   if(qdt>=zero)then
      tmp(ixC^S)=tmp(ixC^S)*qdt*surfaceC(ixC^S,idim)
      wnew(ix^S,iwv+idir)=wnew(ix^S,iwv+idir) -&
         (tmp(ix^S)*fRC(ixR^S)-tmp(hx^S)*fLC(ixL^S))/dvolume(ix^S)
   else
      wnew(ix^S,iwv+idir)=wnew(ix^S,iwv+idir) +&
         (tmp(ix^S)*fRC(ixR^S)-tmp(hx^S)*fLC(ixL^S))/dvolume(ix^S)
   endif
   !SHIFT END
enddo

if((typeconstrain.eq.'fluxCT'.or.typeconstrain.eq.'fluxCD')&
   .and.iw>b0_.and.iw<=b0_+ndim.and.iw/=b0_+idim.and.istep==nstep)&
   call storeflux(qdt,fRC,ixLmin^D,ixRmax^D,idim,iw)

if(oktest)write(*,*)'qdt,fRC,fLC:',qdt,&
   fRC(ixRmin^D-ixmin^D+ixtest^D),fLC(ixLmin^D-ixmin^D+ixtest^D)
if(oktest)write(*,*)'AddFluxRotate wnew:',wnew(ixtest^D,iwtest)

return
end

