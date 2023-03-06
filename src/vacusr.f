!##############################################################################
! module vacusr - EXAMPLE ! setvac -d=22 -g=104,104 -p=hdadiab -u=example 

! INCLUDE:vacnul.specialini.t
! INCLUDE:vacnul.specialbound.t
! INCLUDE:vacnul.specialsource.t
! INCLUDE:vacnul.specialio.t

!=============================================================================
subroutine specialini(ixmin1,ixmin2,ixmax1,ixmax2,w)

! This subroutine initializes w to have a high density blob in it

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmax1,ixmax2
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

double precision:: rhoin,xcent1,xcent2,radius
!-----------------------------------------------------------------------------

write(*,*)'DENSITY BLOB'
write(*,*)'Density inside blob:'
read(*,*) rhoin
write(*,*)'Radius and center coordinates:'
read(*,*) radius,xcent1,xcent2

where(radius**2 > (x(ixmin1:ixmax1,ixmin2:ixmax2,1)-xcent1)**2&
   +(x(ixmin1:ixmax1,ixmin2:ixmax2,2)-xcent2)**2 )
   w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)=rhoin
endwhere

return
end

!=============================================================================
subroutine specialbound(qt,ixmin1,ixmin2,ixmax1,ixmax2,iw,iB,w)

! qt is the time, the iB-th boundary region is in the ix^L region,
! and the subroutine should provide the values for w(ix^S,iw).

! Normal velocity and density are driven at the left boundary by a 
! rho=rho0+drho*sin(omega*t); v=dv*sin(omega*t) perturbation.

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmax1,ixmax2,iw,iB
double precision:: qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),drho

double precision, parameter:: rho0=1.D0, dv=0.01D0, omega=20.D0
!-----------------------------------------------------------------------------

drho=dv/sqrt(eqpar(adiab_)*eqpar(gamma_)/rho0)

select case(iB)
case(1)
   select case(iw)
   case(rho_)
      w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)= rho0 + drho*sin(qt*omega)
   case(m1_)
      w(ixmin1:ixmax1,ixmin2:ixmax2,m1_) = rho0 *   dv*sin(qt*omega)
   case default
      call die('Special boundary is not defined for iw=m2_')
   end select
case default
   call die('Special boundary is defined for iB=1 region only')
end select

return
end

!=============================================================================
subroutine specialsource(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iws,qtC,wCT,qt,w)

! Coriolis forces: SOURCE(m1)=f*m2; SOURCE(m2)=-f*m1 where f=eqpar(coriolis_)
! This is a second order point-implicit trapezoidal scheme meant for split
! source terms (sourcesplit=T) when wCT==w at input.

include 'vacdef.f'

integer:: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
   iws(niw_)
double precision:: qdt,qtC,qt,wCT(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),&
   w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
double precision:: fdt,fdthalf2

integer:: iw,iiw
!-----------------------------------------------------------------------------

if(eqpar(coriolis_)==zero)return

if(.not.sourcesplit)call die&
   ('Coriolis SpecialSource is designed for sourcesplit=T')

if(it==itmin)write(*,*)'Coriolis forces are added'

fdt     = eqpar(coriolis_)*qdt
fdthalf2=(fdt*half)**2

do iiw=1,iws(niw_); iw=iws(iiw)
   select case(iw)
   case(m1_)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)=((1-fdthalf2)&
          *wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)+fdt*wCT(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,m2_))/(1+fdthalf2)
   case(m2_)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)=((1-fdthalf2)&
          *wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)-fdt*wCT(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,m1_))/(1+fdthalf2)
   end select
end do

return
end

!=============================================================================
subroutine getdt_special(w,ixmin1,ixmin2,ixmax1,ixmax2)

! If the Coriolis force is made very strong it may require time step limiting,
! but this is not implemented here.

include 'vacdef.f'
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
integer:: ixmin1,ixmin2,ixmax1,ixmax2
!-----------------------------------------------------------------------------
return
end

!=============================================================================
subroutine readfileini_special(w)
include 'vacdef.f'
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
call die('Special readfileini is not defined')
end
!=============================================================================
subroutine savefileout_special(qunit,w,ixmin1,ixmin2,ixmax1,ixmax2)
include 'vacdef.f'
integer:: qunit,ixmin1,ixmin2,ixmax1,ixmax2
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
call die('Special savefileout is not defined')
end
!=============================================================================
subroutine savefilelog_special(qunit,w,ixmin1,ixmin2,ixmax1,ixmax2)

! The average potential enstrophy (1/V)*Integral dV*(curl v)**2/rho and
! the average ABSOLUTE transverse momentum (1/V)*Integral dV*|m2| 
! are calculated and saved.

include 'vacdef.f'

integer:: qunit,ixmin1,ixmin2,ixmax1,ixmax2,jx1min1,jx1min2,jx1max1,jx1max2,&
   jx2min1,jx2min2,jx2max1,jx2max2,hx1min1,hx1min2,hx1max1,hx1max2,hx2min1,&
   hx2min2,hx2max1,hx2max2
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),denstrophy(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2),enstrophy,meanabsm2
logical:: fileopen
!-----------------------------------------------------------------------------

if(gencoord.or.ndim/=2)call die('Enstrophy for 2D Cartesian grid only')


inquire(qunit,opened=fileopen)
if(.not.fileopen)then
    open(qunit,file=filename(filelog_),status='unknown')
    write(qunit,'(a)')fileheadout
    write(qunit,'(a)')'it   t    dt     enstrophy      abs_m2'
endif


! Shift indices and calculate curl v = dv2/dx1-dv1/dx2, velocity=momentum/rho
jx1min1=ixmin1+kr(1,1);jx1min2=ixmin2+kr(2,1);jx1max1=ixmax1+kr(1,1)
jx1max2=ixmax2+kr(2,1);jx2min1=ixmin1+kr(1,2);jx2min2=ixmin2+kr(2,2)
jx2max1=ixmax1+kr(1,2);jx2max2=ixmax2+kr(2,2);
hx1min1=ixmin1-kr(1,1);hx1min2=ixmin2-kr(2,1);hx1max1=ixmax1-kr(1,1)
hx1max2=ixmax2-kr(2,1);hx2min1=ixmin1-kr(1,2);hx2min2=ixmin2-kr(2,2)
hx2max1=ixmax1-kr(1,2);hx2max2=ixmax2-kr(2,2);
denstrophy(ixmin1:ixmax1,ixmin2:ixmax2)=half*((w(jx1min1:jx1max1,&
   jx1min2:jx1max2,m2_)/w(jx1min1:jx1max1,jx1min2:jx1max2,rho_)&
   -w(hx1min1:hx1max1,hx1min2:hx1max2,m2_)/w(hx1min1:hx1max1,hx1min2:hx1max2,&
   rho_))/dx(ixmin1:ixmax1,ixmin2:ixmax2,1)-(w(jx2min1:jx2max1,&
   jx2min2:jx2max2,m1_)/w(jx2min1:jx2max1,jx2min2:jx2max2,rho_)&
   -w(hx2min1:hx2max1,hx2min2:hx2max2,m1_)/w(hx2min1:hx2max1,hx2min2:hx2max2,&
   rho_))/dx(ixmin1:ixmax1,ixmin2:ixmax2,2))

! Calculate local potential enstrophy density
denstrophy(ixmin1:ixmax1,ixmin2:ixmax2)=denstrophy(ixmin1:ixmax1,&
   ixmin2:ixmax2)**2/w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)

! Average potential enstrophy over the volume
enstrophy=sum(dvolume(ixmin1:ixmax1,ixmin2:ixmax2)*denstrophy(ixmin1:ixmax1,&
   ixmin2:ixmax2))/volume


! Calculate average absolute value of m2
meanabsm2=sum(dvolume(ixmin1:ixmax1,ixmin2:ixmax2)*abs(w(ixmin1:ixmax1,&
   ixmin2:ixmax2,m2_)))/volume


 write(qunit,'(i7,4(1pe13.5))')it,t,dt,enstrophy,meanabsm2

 call flushunit(qunit)

return
end
!=============================================================================
! end module vacusr - EXAMPLE
!##############################################################################
