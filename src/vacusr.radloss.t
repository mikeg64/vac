!##############################################################################
subroutine addsource_rloss(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)
!
! PURPOSE: ADDS THE RADIATIVE LOSSES -dt* rho^2 * Q(T) TO ENERGY
!
!  where 
!
! T=T0*p/rho
!
! de/dt=Q0*rho**2*khi(T)*T**alpha(T)
!
!  The constants Q0 and T0 depend on the average particle mass and the density,
!  energy density, and temperature units.
!  The functions khi and alpha depend on the temperature.
!
! AUTHOR : Sander Belien (sbelien@esa.nascom.nasa.gov or belien@rijnh.nl)
! MODIFIED: Gabor Toth and Rony Keppens
!
! To include these subroutines into the VACUSR module at the top write 
!
!INCLUDE: 'vacusr.radloss.t'
!
! then call addsource_rloss from the specialsource subroutine like
!
!if(eqpar(qhow_)>smalldouble)&
!   call addsource_rloss(qdt,ixI^L,ixO^L,wCT,w)
!
! If there are no other special parameters vacusrpar.t might look like this:
!
!INTEGER,PARAMETER :: qcoef_ =neqpar+1   ! radiative loss coefficient Q0
!INTEGER,PARAMETER :: tunit_ =neqpar+2   ! temperature unit T0
!INTEGER,PARAMETER :: qhow_  =neqpar+3   ! parameter for the method =1 or 2
!INTEGER,PARAMETER :: nspecialpar = 3    !   use 0 for no thermal conduction
!
! To ensure numerical stability call the getdt_rloss subroutine from 
! getdt_special.

include 'vacdef.f'

integer          :: ixI^L,ixO^L,iws(niw_)
double precision :: qdt,qtC,qt,wCT(ixG^T,nw),w(ixG^T,nw)
!------------------------------------------------------------------------------

oktest=index(teststr,'radloss')>=1
if(oktest)write(*,*)'AddSource_RadLoss, energy:',w(ixtest^D,ee_)

if(e_<1)call die('Radiative losses require energy variable indexed by e_>0!')

! de/dt = -n_e**2 * Q(T)
call getQ(wCT,ixO^L,tmp2)
w(ixO^S,e_) = w(ixO^S,e_) - qdt * wCT(ixO^S,rho_)**2 * tmp2(ixO^S)

if(oktest)write(*,*)'Q, energy:',tmp2(ixtest^D),w(ixtest^D,ee_)

return 
end

!==============================================================================
subroutine getQ(w,ix^L,Q)

!!! This subroutine cannot use tmp2 variable!!!

!
! Calculates the optically thin loss function Q according to eqpar(Qhow_).
! Qhow_=1 or Qhow_=2
! (1) See Priest, "Solar MHD", (Reidel, Dordrecht, 1982), p. 88
! (2) See Klimchuk & Gary (Ap. J., 448:925-937,1995
!

include 'vacdef.f'

integer:: ix^L
double precision :: w(ixG^T,nw),Q(ixG^T)
!----------------------------------------------------------------------------
oktest=index(teststr,'getQ')>=1
if(oktest)write(*,*)'GetQ: Q0, T0, Qhow:',&
   eqpar(qcoef_),eqpar(tunit_),eqpar(qhow_)

! Calculate temperature in units given by eqpar(tunit_)
! --> back to dimensionfull temperature in Kelvin
call getpthermal(.true.,w,ix^L,tmp)
tmp(ix^S) = eqpar(tunit_) * tmp(ix^S) / w(ix^S,rho_) 

! Calculate radiative loss based on the method defined by eqpar(Qhow_)
! --> note: Q-values with factor 10**-40 taken out for compilation
!     pusposes
select case(nint(eqpar(Qhow_)))
case(1)
     where (tmp(ix^S)<= 10**3.89063)
        Q(ix^S) = 1.D-40*10**(-2.9) * tmp(ix^S)**(11.7)
     endwhere
     where (tmp(ix^S)>= 10**3.89063.and.tmp(ix^S)<= 10**4.30195)
        Q(ix^S) = 10**(-21.307) * tmp(ix^S)**(6.15)
     endwhere
     where (tmp(ix^S)>= 10**4.30195.and.tmp(ix^S)<= 10**4.575)
        Q(ix^S) = 10**(5.15)
     endwhere
     where (tmp(ix^S)>= 10**4.575.and.tmp(ix^S)<= 10**4.9)
        Q(ix^S) = 10**(-4.0) * tmp(ix^S)**(2.0)
     endwhere
     where (tmp(ix^S)>= 10**4.9.and.tmp(ix^S)<= 10**5.4)
        Q(ix^S) = 10**(5.8)
     endwhere
     where (tmp(ix^S)>= 10**5.4.and.tmp(ix^S)<= 10**5.77)
        Q(ix^S) = 10**(16.6) * tmp(ix^S)**(-2.0)
     endwhere
     where (tmp(ix^S)>= 10**5.77.and.tmp(ix^S)<= 10**6.315)
        Q(ix^S) = 10**(5.06)
     endwhere
     where (tmp(ix^S)>= 10**6.315.and.tmp(ix^S)<= 10**7.60457)
        Q(ix^S) = 10**(9.27) * tmp(ix^S)**(-0.66667)
     endwhere
     where (tmp(ix^S)>= 10**7.60457)
        Q(ix^S) = 10**(0.4) * tmp(ix^S)**(0.5)
     endwhere
case(2)
     where (tmp(ix^S)<= 10**6.18)
        Q(ix^S) = 10**(5.28)
     endwhere
     where (tmp(ix^S)>= 10**6.18.and.tmp(ix^S)<= 10**6.55)
        Q(ix^S) = 10**(14.55) * tmp(ix^S)**(-1.5)
     endwhere
     where (tmp(ix^S)>= 10**6.55)
        Q(ix^S) = 10**(2.54) * tmp(ix^S)**(0.33333)
     endwhere
case default
   write(*,*)'Error in GetQ: Unknown value for eqpar(qhow_):',eqpar(qhow_)
   call die('Error in vacusr.radloss')
end select

! Multiply the result with the eqpar(qcoef_) coefficient
Q(ix^S)=Q(ix^S)*eqpar(qcoef_)*1.D-40

return
end 

!==============================================================================
subroutine getdt_rloss(w,ix^L)
!
! Check diffusion time limit for dt
!

include 'vacdef.f'

double precision :: w(ixG^T,nw),dtdiff_rloss
integer          :: ix^L
!----------------------------------------------------------------------------

oktest=index(teststr,'getdt')>=1

! Calculate Q
call getQ(w,ix^L,tmp2)

! dt< e/(Q*rho**2)
dtdiff_rloss = dtdiffpar/maxval(w(ix^S,rho_)**2*tmp2(ix^S)/w(ix^S,e_))
{^IFMPI call mpiallreduce(dtdiff_rloss,MPI_MIN)}

if (oktest)write(*,*)'Thin radiative losses dt:',dtdiff_rloss

dt = min(dt,dtdiff_rloss)

return
end
!=============================================================================
