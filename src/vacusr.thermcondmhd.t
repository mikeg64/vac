! -*- mode: f90 -*-
!
! PURPOSE: ADDS THE THERMAL CONDUCTION SOURCE TO THE ENERGY EQUATION IN MHD
! 
! S=DIV(KAPPA_i,j . GRAD_j T)
!
! where KAPPA_i,j is a matrix proportional to the B_i B_j diadic matrix.
! Therefore heat conduction accross magnetic field lines is not allowed. 
! The coefficent is determined by the eqpar(kappa_) parameter, and it is
! proportional to temperarture to the 2.5th exponent:
!
! KAPPA_i,j = (B_i B_j / B**2) * eqpar(kappa_) * (p/rho)**2.5
! 
! AUTHOR : Sander Belien (sbelien@esa.nascom.nasa.gov or belien@rijnh.nl)
! Modifications: Gabor Toth
!
! To include these subroutines into the VACUSR module at the top write 
!
!INCLUDE: 'vacusr.thermcondmhd.t'
!
! then call addsource_tcond from the specialsource subroutine like
!
!if(eqpar(kappa_)>smalldouble)&
!   call addsource_tcond(qdt,ixI^L,ixO^L,qtC,wCT,qt,w)
!
! If there are no other special parameters vacusrpar.t might look like this:
!
!INTEGER,PARAMETER :: kappa_  =neqpar+1   ! heat conduction coefficient
!INTEGER,PARAMETER :: nspecialpar = 1
!CHARACTER*5,PARAMETER:: specialparname='kappa'
!
! To ensure numerical stability call the getdt_tcond subroutine from 
! getdt_special:
!
!if(eqpar(kappa_)>smalldouble)&
!   call getdt_tcond(w,ix^L)
!
!##############################################################################
subroutine addsource_tcond(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)

include 'vacdef.f'

integer          :: ixI^L,ixO^L,iws(niw_)
double precision :: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)

double precision :: B2(ixG^T),BgradT(ixG^T)
integer          :: ix^L,idim
!------------------------------------------------------------------------------

oktest=index(teststr,'addsource')>=1
if(oktest)write(*,*)'addsource_tcond_mhd with kappa=',eqpar(kappa_)

if(typephys/='mhd')&
   call die('MHD thermal conduction can only be used for MHD equations!')

! Calculating thermal conduction sources 
! involves second derivatives, two extra layers
call ensurebound(2,ixI^L,ixO^L,qtC,wCT)

ix^L=ixO^L^LADD1;

! Determine temperature on ixI^L
call getpthermal(.true.,wCT,ixI^L,tmp)
tmp(ixI^S) = tmp(ixI^S) / wCT(ixI^S,rho_)

! Determine B**2
B2(ix^S)=( ^C&wCT(ix^S,b^C_)**2+ )

BgradT(ix^S)=zero
do idim=1,ndim
   ! idirth component of gradient of temperature at cell center 
   call gradient(.true.,tmp,ix^L,idim,tmp2)

   ! B dot gradT at cell center
   BgradT(ix^S) = BgradT(ix^S) +wCT(ix^S,b0_+idim)*tmp2(ix^S)
enddo

! Calculate dt*kappa*(B.gradT)/B**2 with kappa=eqpar(kappa)*T**2.5
where(B2(ix^S)>smalldouble)
   BgradT(ix^S)  = BgradT(ix^S)*qdt*eqpar(kappa_)*tmp(ix^S)**2.5/B2(ix^S)
elsewhere
   ! No heat conduction where B=0, which is a poor approximation, but B**2
   ! should not be zero anywhere
   BgradT(ix^S)  = zero
endwhere

do idim=1,ndim 
   ! dt * kappa_i,j . grad T= dt * B_i kappa B.grad T / B.B
   tmp2(ix^S) =wCT(ix^S,b0_+idim)*BgradT(ix^S)

   ! Contribution to div (kappa grad T) from direction idim
   call gradient(.false.,tmp2,ixO^L,idim,tmp)

   ! update w ... 
   w(ixO^S,ee_)=w(ixO^S,ee_)+tmp(ixO^S)
enddo

return
end

!=============================================================================
subroutine getdt_tcond(w,ix^L)

! Check diffusion time limit dt < dtdiffpar * dx_i**2 / ((gamma-1)*kappa_i/rho)
! where                      kappa_i=eqpar(kappa_)*T**2.5*B_i**2/B**2
! and                        T=p/rho

! Baased on C. Hirsch volume 2, p.631, eq.23.2.17 and the remarks below:

include 'vacdef.f'

double precision :: w(ixG^T,nw)
integer          :: ix^L

double precision :: dtdiff_tcond
integer          :: idim
!-----------------------------------------------------------------------------

oktest=index(teststr,'getdt')>=1

! Determine pressure and then the direction independent part of 
! (gamma-1)*kappa_i/rho: tmp=(gamma-1)*eqpar(kappa_)*(p/rho)**2.5/rho
call getpthermal(.true.,w,ix^L,tmp)
tmp(ix^S) = (eqpar(gamma_)-1) * eqpar(kappa_) * &
            (tmp(ix^S) / w(ix^S,rho_))**2.5 / w(ix^S,rho_) 

! Divide by tmp2 = B**2 where the magnetic field is non-vanishing
tmp2(ix^S)=( ^C&w(ix^S,b^C_)**2+ )
where(tmp2(ix^S)>smalldouble)tmp(ix^S) = tmp(ix^S)/tmp2(ix^S)

do idim=1,ndim
   ! (gamma-1)*kappa_i/rho = tmp*B_i**2
   tmp2(ix^S)=tmp(ix^S)*w(ix^S,b0_+idim)**2

   ! dt< dtdiffpar * dx_idim**2/((gamma-1)*kappa_idim/rho)
   dtdiff_tcond=dtdiffpar/maxval(tmp2(ix^S)/dx(ix^S,idim)**2)
   {^IFMPI call mpiallreduce(dtdiff_tcond,MPI_MIN)}

   ! limit the time step
   dt=min(dt,dtdiff_tcond)
   if(oktest)write(*,*)'Thermal cond. limit on dt for direction:',&
       idim,dtdiff_tcond
enddo

return
end
!=============================================================================

