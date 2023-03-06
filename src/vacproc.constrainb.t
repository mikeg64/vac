!##############################################################################
! module vacproc.constrainb - Constrained transport for mhd(iso) in 2 or 3D

!=============================================================================
subroutine constrainb(w)

include 'vacdef.f'

double precision:: w(ixG^T,nw),w1(ixG^T,nw)
double precision:: coef{^IFTHREED, qtmp(ixG^T)}

integer:: ix^L,ixC^L,jxC^L,hxC^L,dix^D,j1x^L,j2x^L,h1x^L,h2x^L
integer:: idim,idir,iws(niw_)
{^IFTHREED 
integer:: j3x^L,h3x^L
double precision:: efield(ixG^T,ndim)
}

double precision:: qvol(ixG^T)
save qvol
logical:: initialized
data initialized/.false./
!-----------------------------------------------------------------------------

oktest=index(teststr,'constrainb')>=1
if(oktest)write(*,*)'ConstrainB wold, w:',&
    wold(ixtest^D,iwtest),w(ixtest^D,iwtest)

if(abs(eqpar(eta_))>smalldouble)&
   call die('ConstrainB currently does not work with resistivity!')

! Remove magnetic part of total energy if thermal pressure should be fixed
if(typephys=='mhd'.and.index(teststr,'keeppth')>0)&
    w(ixM^S,ee_)=w(ixM^S,ee_)-(^D&w(ixM^S,b^D_)**2+)/2

! Index limits for CD and CT fluxes
ix^L=ixM^L^LADD1;
ixCmin^D=ixmin^D; ixCmax^D=ixMmax^D;

! Calculate -dt/2 times the electric field components
! For CD schemes the 1/2 is for the /2dx, for CT schemes for the b->B average
select case(typeconstrain)
   case('fieldCT')
      ! Average v and B in time and then for cell corners/edges (2D/3D)
      do idim=1,ndim
         tmp(ix^S)=w(ix^S,m0_+idim)/w(ix^S,rho_)&
                   +wold(ix^S,m0_+idim)/wold(ix^S,rho_)
         tmp2(ix^S)=w(ix^S,b0_+idim)+wold(ix^S,b0_+idim)
         ! Average 4 neighbouring cells for cell corners/edges (2D/3D)
         w1(ix^S,v0_+idim)=zero
         w1(ix^S,b0_+idim)=zero
         {do dix^D=0,1\}
            jxC^L=ixC^L+dix^D;
            w1(ixC^S,v0_+idim)=w1(ixC^S,v0_+idim)+tmp(jxC^S)
            w1(ixC^S,b0_+idim)=w1(ixC^S,b0_+idim)+tmp2(jxC^S)
         {enddo^D&\}
       enddo
      ! Both v and B should be divided by 8 and another /2 for b->B average
      coef=dt/8**2/2
      fstore(ixC^S,ndim)=coef*&
         (w1(ixC^S,v1_)*w1(ixC^S,b2_)-w1(ixC^S,v2_)*w1(ixC^S,b1_))
      {^IFTHREED call die('fieldCT algorithm is not implemented for 3D yet')}
   case('EfieldCT')
      ! Average v x B in time and then for cell corners
      ! so it should be divided by 2*2**ndim and another /2 for b->B average
      coef=dt/8/2
      !(VxB)_z=m1*b2-m2*b1/rho
      idim=ndim
      fstore(ix^S,ndim)=coef*(&
         (w(ix^S,m1_)*w(ix^S,b2_)-w(ix^S,m2_)*w(ix^S,b1_))/w(ix^S,rho_)&
        +(wold(ix^S,m1_)*wold(ix^S,b2_)-wold(ix^S,m2_)*wold(ix^S,b1_))&
         /wold(ix^S,rho_))
      {^IFTHREED call die('EfieldCT algorithm is not implemented for 3D')}
      tmp(ix^S) =zero
      ! Average neighbouring cells for cell corners
      {do dix^D=0,1\}
         jxC^L=ixC^L+dix^D;
         tmp(ixC^S)=tmp(ixC^S)+fstore(jxC^S,idim)
      {enddo^D&\}
      fstore(ix^S,idim)=tmp(ixC^S)
   case('fieldCD1')
      ! Calculate v x B
      ! there is /2 for /2dx
      coef=dt/2
      fstore(ix^S,ndim)=coef*&
        (w(ix^S,m1_)*w(ix^S,b2_)-w(ix^S,m2_)*w(ix^S,b1_))/w(ix^S,rho_)
      {^IFTHREED
      !(VxB)_x=(m2*b3-m3*b2)/rho
      fstore(ix^S,1)=coef*&
         (w(ix^S,v2_)*w(ix^S,b3_)-w(ix^S,v3_)*w(ix^S,b2_))/w(ix^S,rho_)
      !(VxB)_y=m3*b1-m1*b3/rho
      fstore(ix^S,2)=coef*&
         (w(ix^S,v3_)*w(ix^S,b1_)-w(ix^S,v1_)*w(ix^S,b3_))/w(ix^S,rho_)
      }
   case('fieldCD')
      ! Average B and v in time
      do idim=1,ndim
         w1(ix^S,v0_+idim)=w(ix^S,m0_+idim)/w(ix^S,rho_)&
                          +wold(ix^S,m0_+idim)/wold(ix^S,rho_)
         w1(ix^S,b0_+idim)=w(ix^S,b0_+idim)+wold(ix^S,b0_+idim)
      enddo
      ! Both v and B should be divided by 2, and also there is 2 for /2dx
      coef=dt/8
      !(VxB)_z=v1*b2-v2*b1
      fstore(ix^S,ndim)=coef*&
         (w1(ix^S,v1_)*w1(ix^S,b2_)-w1(ix^S,v2_)*w1(ix^S,b1_))
      {^IFTHREED
      !(VxB)_x=v2*b3-v3*b2
      fstore(ix^S,1)=coef*&
         (w1(ix^S,v2_)*w1(ix^S,b3_)-w1(ix^S,v3_)*w1(ix^S,b2_))
      !(VxB)_y=v3*b1-v1*b3
      fstore(ix^S,2)=coef*&
         (w1(ix^S,v3_)*w1(ix^S,b1_)-w1(ix^S,v1_)*w1(ix^S,b3_))
      }
   case('EfieldCD1')
      ! Temporally first order variant
      ! There is /2 for /2dx
      coef=dt/2
      ! (VxB)_z=(m1*b2-m2*b1)/rho
      fstore(ix^S,ndim)=coef*&
            (w(ix^S,m1_)*w(ix^S,b2_)-w(ix^S,m2_)*w(ix^S,b1_))/w(ix^S,rho_)
      {^IFTHREED
      !(VxB)_x=(m2*b3-m3*b2)/rho
      fstore(ix^S,1)=coef*&
            (w(ix^S,m2_)*w(ix^S,b3_)-w(ix^S,m3_)*w(ix^S,b2_))/w(ix^S,rho_)

      !(VxB)_y=m3*b1-m1*b3/rho
      fstore(ix^S,2)=coef*&
            (w(ix^S,m3_)*w(ix^S,b1_)-w(ix^S,m1_)*w(ix^S,b3_))/w(ix^S,rho_)
      }
   case('EfieldWCD1')
      {^IFTHREED call die('EfieldWCD1 is not implemented for 3D')}
      ! Temporally first order weighted variant
      ! There is /2 for /2dx
      coef=dt/2
      ! Calculate the negative electric field v X B
      tmp(ixG^S)=coef*&
            (w(ixG^S,m1_)*w(ixG^S,b2_)-w(ixG^S,m2_)*w(ixG^S,b1_))/w(ixG^S,rho_)
      ! Do a spatial averaging with 1/2 weight for the center cell
      fstore(ix^S,ndim)=half*tmp(ix^S)
      ! The four neighbouring cells have 1/8 weight
      coef=one/8
      do idim=1,ndim
         j1x^L=ix^L+kr(idim,^D); h1x^L=ix^L-kr(idim,^D);
         fstore(ix^S,ndim)=fstore(ix^S,ndim)+coef*(tmp(j1x^S)+tmp(h1x^S))
      enddo
   case('EfieldCD')
      ! v x B is averaged in time and there is /2 for /2dx
      coef=dt/4
      ! (VxB)_z=m1*b2-m2*b1/rho
      fstore(ix^S,ndim)=coef*&
           ((w(ix^S,m1_)*w(ix^S,b2_)-w(ix^S,m2_)*w(ix^S,b1_))/w(ix^S,rho_)&
           +(wold(ix^S,m1_)*wold(ix^S,b2_)-wold(ix^S,m2_)*wold(ix^S,b1_))&
           /wold(ix^S,rho_))

      {^IFTHREED
      !(VxB)_x=(m2*b3-m3*b2)/rho
      fstore(ix^S,1)=coef*&
           ((   w(ix^S,m2_)*w(ix^S,b3_)-w(ix^S,m3_)*w(ix^S,b2_))/w(ix^S,rho_)&
           +(wold(ix^S,m2_)*wold(ix^S,b3_)-wold(ix^S,m3_)*wold(ix^S,b2_))&
           /wold(ix^S,rho_))

      !(VxB)_y=m3*b1-m1*b3/rho
      fstore(ix^S,2)=coef*&
           ((w(ix^S,m3_)*w(ix^S,b1_)-w(ix^S,m1_)*w(ix^S,b3_))/w(ix^S,rho_)&
           +(wold(ix^S,m3_)*wold(ix^S,b1_)-wold(ix^S,m1_)*wold(ix^S,b3_))&
           /wold(ix^S,rho_))
      }
   case('EfieldWCD')
      {^IFTHREED call die('EfieldWCD is not implemented for 3D')}
      ! v x B is averaged in time and there is /2 for /2dx
      coef=dt/4
      ! Weighted variant of field-CD scheme
      tmp(ixG^S)=coef*&
          ((w(ixG^S,m1_)*w(ixG^S,b2_)-w(ixG^S,m2_)*w(ixG^S,b1_))/w(ixG^S,rho_)&
          +(wold(ixG^S,m1_)*wold(ixG^S,b2_)-wold(ixG^S,m2_)*wold(ixG^S,b1_))&
          /wold(ixG^S,rho_))
      ! Do a spatial averaging with 0.5 weight for the center cell
      fstore(ix^S,ndim)=half*tmp(ix^S)
      ! Other four cells have 1/8 weight
      coef=one/8
      do idim=1,ndim
         j1x^L=ix^L+kr(idim,^D); h1x^L=ix^L-kr(idim,^D);
         fstore(ix^S,ndim)=fstore(ix^S,ndim)+coef*(tmp(j1x^S)+tmp(h1x^S))
      enddo
   case('fluxCT')
      call boundfluxb(fstore)
      ! Omega_z=(fx+fx-fy-fy)/4 and another /2 for b->B average
      j1x^L=ixC^L+kr(^D,1); j2x^L=ixC^L+kr(^D,2);
      tmp(ixC^S)=(fstore(ixC^S,1)+fstore(j2x^S,1)&
                 -fstore(ixC^S,2)-fstore(j1x^S,2))/8
      {^IFTHREED
      j3x^L=ixC^L+kr(^D,3);
      ! Omega_y=(fz+fz-fx-fx)/4 and another /2 for b->B average
      tmp2(ixC^S)=(fstore(ixC^S,3)+fstore(j1x^S,3)&
                  -fstore(ixC^S,1)-fstore(j3x^S,1))/8
      ! Omega_x=(fy+fy-fz-fz)/4 and another /2 for b->B average
      qtmp(ixC^S)=(fstore(ixC^S,2)+fstore(j3x^S,2)&
                  -fstore(ixC^S,3)-fstore(j2x^S,3))/8
      fstore(ixC^S,2)=tmp2(ixC^S)
      fstore(ixC^S,1)=qtmp(ixC^S)
      }
      fstore(ixC^S,ndim)=tmp(ixC^S)
   case('ryuCT','trfluxCT')
      {^IFTHREED call die('trfluxCT scheme is not implemented for 3D')}
      j1x^L=ixC^L+kr(^D,1); j2x^L=ixC^L+kr(^D,2);
      ! Add transport fluxes at time level n to the numerical fluxes in fstore
      fstore(ixC^S,1)=fstore(ixC^S,1)+half*dt*&
        (wold(ixC^S,b2_)*wold(ixC^S,m1_)/wold(ixC^S,rho_)+&
         wold(j1x^S,b2_)*wold(j1x^S,m1_)/wold(j1x^S,rho_))
      fstore(ixC^S,2)=fstore(ixC^S,2)+half*dt*&
        (wold(ixC^S,b1_)*wold(ixC^S,m2_)/wold(ixC^S,rho_)+&
         wold(j2x^S,b1_)*wold(j2x^S,m2_)/wold(j2x^S,rho_))

      call boundfluxb(fstore)
      ! Omega_3=(fx+fx-fy-fy)/2 and another /2 for b->B average
      tmp(ixC^S)=(fstore(ixC^S,1)+fstore(j2x^S,1)&
                 -fstore(ixC^S,2)-fstore(j1x^S,2))/4
      fstore(ixC^S,ndim)=tmp(ixC^S)

      {^IFTWOD
        if(oktest)write(*,*)'omega, fx ii,ij, fy ii,ji:',&
          tmp(ixtest1,ixtest2),&
          fstore(ixtest1,ixtest2,1),fstore(ixtest1,ixtest2+1,1),&
          fstore(ixtest1,ixtest2,2),fstore(ixtest1+1,ixtest2,2)
      }
   case('fluxCD')
      call boundfluxb(fstore)
      ! Omega_z=(fx+fx-fy-fy)/4 and another /2 for the /2dx
      h1x^L=ix^L-kr(^D,1); h2x^L=ix^L-kr(^D,2);
      tmp(ix^S)=(fstore(ix^S,1)+fstore(h1x^S,1)&
                -fstore(ix^S,2)-fstore(h2x^S,2))/8
      {^IFTHREED
      h3x^L=ix^L-kr(^D,3);
      ! Omega_y=(fz+fz-fx-fx)/4 and another /2 for the /2dx
      tmp2(ix^S)=(fstore(ix^S,3)+fstore(h3x^S,3)&
                 -fstore(ix^S,1)-fstore(h1x^S,1))/8
      ! Omega_x=(fy+fy-fz-fz)/4 and another /2 for the /2dx
      qtmp(ix^S)=(fstore(ix^S,2)+fstore(h2x^S,2)&
                 -fstore(ix^S,3)-fstore(h3x^S,3))/8
      fstore(ix^S,2)=tmp2(ix^S)
      fstore(ix^S,1)=qtmp(ix^S)
      }
      fstore(ix^S,ndim)=tmp(ix^S)
   case default
      call die('Unknown value for typeconstrain='//typeconstrain)
end select

if(typeaxial=='cylinder')then
   !Multiply B components and fstore=-E*dt by radial distance
   fstore(ix^S,ndim)=fstore(ix^S,ndim)*x(ix^S,r_)
   ^D&wold(ixM^S,b^D_)=wold(ixM^S,b^D_)*x(ix^S,r_);
endif

! Difference the -dt/2*electric field in fstore

! Index limits for differencing
j1x^L=ixM^L+kr(^D,1); j2x^L=ixM^L+kr(^D,2);
h1x^L=ixM^L-kr(^D,1); h2x^L=ixM^L-kr(^D,2);
{^IFTHREED j3x^L=ixM^L+kr(^D,3); h3x^L=ixM^L-kr(^D,1); }

{^IFGEN
if(.not.initialized .and.gencoord)then
   initialized=.true.
   ! Calculate Jacobian volume V eq(38) in Toth "The divB=0 constraint"
   ! Take half of it since the electric field is also divided by 2
   if(ndim==2)then
      qvol(ixM^S)=((x(j1x^S,1)-x(h1x^S,1))*(x(j2x^S,2)-x(h2x^S,2))-&
                   (x(j2x^S,1)-x(h2x^S,1))*(x(j1x^S,2)-x(h1x^S,2)))/2
      ! In cylindrical symmetry multiply volume by radial distance
      if(typeaxial=='cylinder')qvol(ixM^S)=qvol(ixM^S)*x(ixM^S,r_);
   else
      ! The 3D volume is almost fine, but there is a factor 4 missing
      qvol(ixM^S)=4*dvolume(ixM^S)
   endif
endif
}

if(index(typeconstrain,'CD')>0)then
   ! Calculate new B based on centered differences of v x B

   {^IFTWOD
   if(gencoord)then
      !dBx/dt eq(37) in Toth "The divB=0 constraint", fstore=-E_z*dt
      w(ixM^S,b1_)=wold(ixM^S,b1_)+&
         ((x(j1x^S,1)-x(h1x^S,1))*(fstore(j2x^S,2)-fstore(h2x^S,2))-&
          (x(j2x^S,1)-x(h2x^S,1))*(fstore(j1x^S,2)-fstore(h1x^S,2)))&
         /qvol(ixM^S)
      w(ixM^S,b2_)=wold(ixM^S,b2_)+&
         ((x(j1x^S,2)-x(h1x^S,2))*(fstore(j2x^S,2)-fstore(h2x^S,2))-&
          (x(j2x^S,2)-x(h2x^S,2))*(fstore(j1x^S,2)-fstore(h1x^S,2)))&
         /qvol(ixM^S)
   else
     !dBx/dt=d(VxB_z)/dy and dBy/dt=-d(VxB_z)/dx
     w(ixM^S,b1_)=wold(ixM^S,b1_)+(fstore(j2x^S,2)-fstore(h2x^S,2))/dx(ixM^S,2)
     w(ixM^S,b2_)=wold(ixM^S,b2_)-(fstore(j1x^S,2)-fstore(h1x^S,2))/dx(ixM^S,1)
   endif
   }

   {^IFTHREED
   if(gencoord)then
      ! Calculate (-dt/2*) curvilinear electric field components 
      ! Ecurl=Jinvtransp.E where Jinvtransp(i,j)=d_i x_j
      do idim=1,ndim
         ! Index limits for electric field component idim. One extra cell in 
         ! directions orthogonal to idim. ixL has extra cell in all directions
         ix^L=ixM^L^LADD1;
         ixC^L=ix^L^LSUBkr(idim,^D);
         jxC^L=ixC^L+kr(idim,^D);
         hxC^L=ixC^L-kr(idim,^D);
         efield(ixC^S,idim)={^D&(x(jxC^S,^D)-x(hxC^S,^D))*fstore(ixC^S,^D)+}
      enddo
      ! Calculate (dt/2*) dB_xi/dt, dB_eta/dt, dB_zeta/dt from efield=-dt/2*E
      fstore(ixM^S,1)=efield(j2x^S,3)-efield(h2x^S,3) &
                     -efield(j3x^S,2)+efield(h3x^S,2)
      fstore(ixM^S,2)=efield(j3x^S,1)-efield(h3x^S,1) &
                     -efield(j1x^S,3)+efield(h1x^S,3)
      fstore(ixM^S,3)=efield(j1x^S,2)-efield(h1x^S,2) &
                     -efield(j2x^S,1)+efield(h2x^S,1)
      ! Update the B field B_n+1=B_n+dt*det(J)Jinv.dBcurle/dt
      ! and detJ=1/(8*dvolume) and qvol=4*dvolume, and Jinv(i,j)=d_j x_i
      do idim=1,ndim
         w(ixM^S,b0_+idim)=wold(ixM^S,b0_+idim)+&
            (^D&(x(j^Dx^S,idim)-x(h^Dx^S,idim))*fstore(ixM^S,^D)+)/qvol(ixM^S)
      enddo
   else
      ! Cartesian grid
      !dBx/dt=d(VxB_z)/dy-d(VxB_y)/dz
      w(ixM^S,b1_)=wold(ixM^S,b1_)&
         +(fstore(j2x^S,3)-fstore(h2x^S,3))/dx(ixM^S,2)&
         -(fstore(j3x^S,2)-fstore(h3x^S,2))/dx(ixM^S,3)

      !dBy/dt=d(VxB_x)/dz-d(VxB_z)/dx
      w(ixM^S,b2_)=wold(ixM^S,b2_)&
         +(fstore(j3x^S,1)-fstore(h3x^S,1))/dx(ixM^S,3)&
         -(fstore(j1x^S,3)-fstore(h1x^S,3))/dx(ixM^S,1)

      !dBz/dt=d(VxB_y)/dx-d(VxB_x)/dy
      w(ixM^S,b3_)=wold(ixM^S,b3_)&
         +(fstore(j1x^S,2)-fstore(h1x^S,2))/dx(ixM^S,1)&
         -(fstore(j2x^S,1)-fstore(h2x^S,1))/dx(ixM^S,2)
   endif
   }
else
   ! Constrained Transport style differencing
   {^IFTWOD
   !average CT estimate of the electric field in direction x (b->B)
   tmp(ixM^LIM1:,ixC^LIM2:)=&
      fstore(ixM^LIM1:,ixC^LIM2:,2)+fstore(h1x^LIM1:,ixC^LIM2:,2)

   !average CT estimate of the electric field in direction y (b->B)
   tmp2(ixC^LIM1:,ixM^LIM2:)=&
      fstore(ixC^LIM1:,ixM^LIM2:,2)+fstore(ixC^LIM1:,h2x^LIM2:,2)

   if(gencoord)then
      !dBx/dt eq(37) in Toth "The divB=0 constraint", fstore=-E_z*dt
      w(ixM^S,b1_)=wold(ixM^S,b1_)+&
         ((x(j1x^S,1)-x(h1x^S,1))*(tmp(ixM^S) -tmp(h2x^S))&
         -(x(j2x^S,1)-x(h2x^S,1))*(tmp2(ixM^S)-tmp2(h1x^S)))&
         /qvol(ixM^S)
      w(ixM^S,b2_)=wold(ixM^S,b2_)+&
         ((x(j1x^S,2)-x(h1x^S,2))*(tmp(ixM^S) -tmp(h2x^S))&
         -(x(j2x^S,2)-x(h2x^S,2))*(tmp2(ixM^S)-tmp2(h1x^S)))&
         /qvol(ixM^S)
   else
      ! dBx/dt=+d(VxB_z)/dy 
      w(ixM^S,b1_)=wold(ixM^S,b1_)+(tmp(ixM^S)-tmp(h2x^S))/dx(ixM^S,2)
      ! dBy/dt=-d(VxB_z)/dx
      w(ixM^S,b2_)=wold(ixM^S,b2_)-(tmp2(ixM^S)-tmp2(h1x^S))/dx(ixM^S,1)

      if(oktest)write(*,*)'omega hh,ih,hi,ii:',&
         fstore(ixtest1-1,ixtest2-1,2),fstore(ixtest1,ixtest2-1,2),&
         fstore(ixtest1-1,ixtest2,2)  ,fstore(ixtest1,ixtest2,2)
      if(oktest)write(*,*)'flux ii,hi,ih:',tmp(ixtest1,ixtest2),&
         tmp(ixtest1-1,ixtest2),tmp(ixtest1,ixtest2-1)
   endif
   }

   {^IFTHREED
   !average CT fluxes in direction x then dBx/dt=d(VxB_z)/dy-d(VxB_y)/dz
   tmp(ixM^LIM1:,ixC^LIM2:,ixC^LIM3:)=&
      fstore(ixM^LIM1:,ixC^LIM2:,ixC^LIM3:,3)+&
      fstore(h1x^LIM1:,ixC^LIM2:,ixC^LIM3:,3)
   tmp2(ixM^LIM1:,ixC^LIM2:,ixC^LIM3:)=&
      fstore(ixM^LIM1:,ixC^LIM2:,ixC^LIM3:,2)+&
      fstore(h1x^LIM1:,ixC^LIM2:,ixC^LIM3:,2)
   w(ixM^S,b1_)=wold(ixM^S,b1_)&
     +( tmp(ixM^S)- tmp(h2x^S))/dx(ixM^S,2)&
     -(tmp2(ixM^S)-tmp2(h3x^S))/dx(ixM^S,3)

   !average CT fluxes in direction y then dBy/dt=d(VxB_x)/dz-d(VxB_z)/dx
   tmp(ixC^LIM1:,ixM^LIM2:,ixC^LIM3:)=&
      fstore(ixC^LIM1:,ixM^LIM2:,ixC^LIM3:,1)+&
      fstore(ixC^LIM1:,h2x^LIM2:,ixC^LIM3:,1)
   tmp2(ixC^LIM1:,ixM^LIM2:,ixC^LIM3:)=&
      fstore(ixC^LIM1:,ixM^LIM2:,ixC^LIM3:,3)+&
      fstore(ixC^LIM1:,h2x^LIM2:,ixC^LIM3:,3)
   w(ixM^S,b2_)=wold(ixM^S,b2_)&
     +( tmp(ixM^S)- tmp(h3x^S))/dx(ixM^S,3)&
     -(tmp2(ixM^S)-tmp2(h1x^S))/dx(ixM^S,1)

   !average CT fluxes in direction z then dBz/dt=d(VxB_y)/dx-d(VxB_x)/dy
   tmp(ixC^LIM1:,ixC^LIM2:,ixM^LIM3:)=&
      fstore(ixC^LIM1:,ixC^LIM2:,ixM^LIM3:,2)+&
      fstore(ixC^LIM1:,ixC^LIM2:,h3x^LIM3:,2)
   tmp2(ixC^LIM1:,ixC^LIM2:,ixM^LIM3:)=&
      fstore(ixC^LIM1:,ixC^LIM2:,ixM^LIM3:,1)+&
      fstore(ixC^LIM1:,ixC^LIM2:,h3x^LIM3:,1)
   w(ixM^S,b3_)=wold(ixM^S,b3_)&
     +( tmp(ixM^S)- tmp(h1x^S))/dx(ixM^S,1)&
     -(tmp2(ixM^S)-tmp2(h2x^S))/dx(ixM^S,2)
   }
endif

if(typeaxial=='cylinder')then
   !Divide back B components by radial distance
   ^D&wold(ixM^S,b^D_)=wold(ixM^S,b^D_)/x(ix^S,r_);
   ^D&w(ixM^S,b^D_)=w(ixM^S,b^D_)/x(ix^S,r_);
endif

if(typephys=='mhd'.and.index(teststr,'keeppth')>0)then
   ! Add magnetic part of total energy if thermal pressure should be fixed
   w(ixM^S,ee_)=w(ixM^S,ee_)+(^D&w(ixM^S,b^D_)**2+)/2
   ! Apply boundary conditions for energy and the divergence free field
   call getboundary(t+dt,e_,b^ND_,1,ndim,w)
else
   ! Apply boundary conditions to the divergence free field
   call getboundary(t+dt,b1_,b^ND_,1,ndim,w)
endif

{^IFRES
if(.not.sourcesplit.and.abs(eqpar(eta_))>smalldouble)then
   ! Time center old and new w, both have zero divergence
   w1(ixG^S,1:nw)=half*(wold(ixG^S,1:nw)+w(ixG^S,1:nw))
   ! Update magnetic field components 1..ndim only
   ^D&iws(^D)=b^D_; iws(niw_)=ndim
   ! Add resistive sources with centered differences, should keep div B=0
   call addsource2(dt,ixG^L,ixM^L,iws,t+dt/2,w1,t+dt,w)
endif
}

if(oktest)write(*,*)'New w:',w(ixtest^D,iwtest)

fstore(ixG^S,1:ndim)=zero

return
end

!=============================================================================
subroutine boundfluxb(f)

! Apply boundary conditions on the flux of the B field
!
! The interface centered fluxes are known for the mesh cells,
! but for the flux interpolated CT/CD algorithms, fluxes are required
! for all the grid cells.

include 'vacdef.f'

double precision:: f(ixG^T,ndim)

integer:: idim,idir,iB,ix^L,ixpair^L,ixe,ix,ixC^L
!-----------------------------------------------------------------------------

oktest=index(teststr,'boundfluxb')>=1
if(oktest)write(*,*)'BoundFluxB f:',f(ixtest^D,idimtest)

! idir is the direction of the flux component to be bounded
do idir=1,ndim

   ! Copy the correct part of f into tmp
   ixCmin^D=ixMmin^D-kr(idir,^D); ixCmax^D=ixMmax^D;
   tmp(ixC^S)=f(ixC^S,idir)

   {^IFMPI
   ! Fill in ghost cells from MPI neighbors
   call mpibound(1,tmp)
   }

   if(oktest)write(*,*)'idir,ixC:',idir,ixC^L

   do iB=1,nB
      ! Normal direction and size of boundary region
      idim=idimB(iB)
      ix^L=ixB^LIM(^D,iB);

      ! Boundary conditions for the flux are based 
      ! on the boundary conditions for B normal
      select case(typeB(b0_+idim,iB))
      case('periodic')
         ixpair^L=ixB^LIM(^D,ipairB(iB));
         select case(idim)
         {case(^D)
         if(upperB(iB))then
            ixpair^LIM^D=ixpair^LIM^D+dixB^LIM^D;
         else
            ixpair^LIM^D=ixpair^LIM^D-dixB^LIM^D;
         endif
         \}
         end select
         tmp(ix^S)=tmp(ixpair^S)
      case('cont','fixed')
         select case(idim)
         {case(^D)
            if(upperB(iB))then
               ixe=ixmin^D-1
            else if(idim==idir) then
               ixe=ixmax^D
            else
               ixe=ixmax^D+1
            endif
            !HPF$ INDEPENDENT
            do ix= ix^DL
               tmp(ix^D%ix^S)=tmp(ixe^D%ix^S)
            end do
            \}
         end select
      case('asymm')
         ! asymmetric orthogonal component means perfectly conducting boundary
         ! consequently v and B are parallel at the boundary, and Efield=0
         ! Flux for the ghost cells should be antisymmetric
         {^IFTHREED call die( &
            'Error in BoundFlux: asymm boundary type not implemented for 3D'}
         select case(idim)
         {case(^D)
            if(upperB(iB))then
               ixe=ixmin^D-1
            else if(idim==idir) then
               ixe=ixmax^D
            else
               ixe=ixmax^D+1
            endif
            if(idim==idir)then
               !HPF$ INDEPENDENT
               do ix= ix^DL
                  tmp(ix^D%ix^S)=-tmp(2*ixe-ix^D%ix^S)
               enddo
               tmp(ixe^D%ix^S)=zero
            else
               !HPF$ INDEPENDENT
               do ix= ix^DL
                  tmp(ix^D%ix^S)=-tmp(2*ixe+1-ix^D%ix^S)
               enddo
            endif
         \}
         end select
      {^IFMPI
      case('mpi','mpiperiod')
         ! This boundary is handled by MPI\}
      case default
         call die('Error in BoundFlux: Unimplemented boundary condition='//&
            typeB(b0_+idim,iB))
      end select
   enddo

   ! Copy the correct part of f into tmp again, since part of the
   ! boundary regions (defined for cell centered values) overlap it.

   tmp(ixC^S)=f(ixC^S,idir);

   ! Copy back the bounded tmp into f for the full grid
   f(ixG^S,idir)=tmp(ixG^S);

   if(oktest)write(*,*)'new f:',f(ixtest^D,idimtest)
enddo

return
end

!=============================================================================
! end module vacproc.constrainb
!##############################################################################
