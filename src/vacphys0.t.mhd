!##############################################################################
! module vacphys0 - mhd

!=============================================================================
subroutine conserve(ix^L,w)

! Transform primitive variables into conservative ones

include 'vacdef.f'

integer:: ix^L
double precision:: w(ixG^T,nw)
!-----------------------------------------------------------------------------

oktest=index(teststr,'conserv')>=1
if(oktest)write(*,*)'Conserve w:',w(ixtest^D,iwtest)

! Calculate total energy from pressure, kinetic and magnetic energy
w(ix^S,e_)=w(ix^S,p_)/(eqpar(gamma_)-1)+&
   half*(w(ix^S,rho_)*(^C&w(ix^S,v^C_)**2+)+(^C&w(ix^S,b^C_)**2+))

! Convert velocity to momentum
^C&w(ix^S,m^C_)=w(ix^S,rho_)*w(ix^S,v^C_);

if(oktest)write(*,*)'New w:',w(ixtest^D,iwtest)

return
end

!=============================================================================
subroutine primitive(ix^L,w)

! Transform conservative variables into primitive ones

include 'vacdef.f'

integer:: ix^L
double precision:: w(ixG^T,nw)
!-----------------------------------------------------------------------------

oktest=index(teststr,'primitive')>=1
if(oktest)write(*,*)'Primitive w:',w(ixtest^D,iwtest)

! Calculate pressure
call getpthermal(.true.,w,ix^L,tmp)
w(ix^S,p_)=tmp(ix^S)

! Convert momentum to velocity
^C&w(ix^S,v^C_)=w(ix^S,m^C_)/w(ix^S,rho_);

if(oktest)write(*,*)'New w:',w(ixtest^D,iwtest)

return
end

!=============================================================================
subroutine getv(w,ix^L,idim,v)

! Calculate v_idim=m_idim/rho within ix

include 'vacdef.f'

integer:: ix^L,idim
double precision:: w(ixG^T,nw),v(ixG^T)
!-----------------------------------------------------------------------------

oktest=index(teststr,'getv')>=1
if(oktest)write(*,*)'GetV w:',w(ixtest^D,iwtest)

v(ix^S)=w(ix^S,m0_+idim)/w(ix^S,rho_)

if(oktest)write(*,*)'GetV v:',v(ixtest^D)

return 
end


!=============================================================================
subroutine getcmax(new_cmax,w,ix^L,idim,cmax)

! Calculate cmax_idim=cfast_i+abs(v_idim) within ix^L
! where cfast_i=sqrt(0.5*(cf**2+sqrt(cf**4-4*cs**2*b_i**2/rho)))
! and cf**2=b**2/rho+cs**2/rho is the square of the speed of the fast wave 
! perpendicular to the magnetic field, and cs is the sound speed.

include 'vacdef.f'

logical:: new_cmax
integer:: ix^L,idim
double precision:: w(ixG^T,nw),cmax(ixG^T)
double precision:: csound2(ixG^T),cfast2(ixG^T)
save csound2,cfast2
!-----------------------------------------------------------------------------

oktest=index(teststr,'getcmax')>=1

!Direction independent part of getcmax:
if(new_cmax)then
   new_cmax=.false.
   call getcsound2(w,ix^L,csound2)
   if(oktest)write(*,*)'csound2:',csound2(ixtest^D)
   cfast2(ix^S)=( ^C&w(ix^S,b^C_)**2+ )/w(ix^S,rho_)+csound2(ix^S)
end if
if(oktest)write(*,*)'cfast2:',cfast2(ixtest^D)
!if(oktest.and.index(teststr,'discriminant')>=1)write(*,*)'Discriminant:',&
!   cfast2(ix^S)**2-4*csound2(ix^S)*(w(ix^S,b0_+idim)**2)/w(ix^S,rho_)

cmax(ix^S)=sqrt(half*(cfast2(ix^S)+ &
   sqrt(cfast2(ix^S)**2-4*csound2(ix^S)* &
   (w(ix^S,b0_+idim)**2)/w(ix^S,rho_)))) &
   +abs(w(ix^S,m0_+idim)/w(ix^S,rho_))

if(oktest) write(*,*)'cmax:',cmax(ixtest^D)

return 
end

!=============================================================================
subroutine getcsound2prim(w,ix^L,csound2)

! Calculate the square of the thermal sound speed csound2 within ix^L
! from the primitive variables in w.
! csound2=gamma*p/rho

include 'vacdef.f'

double precision:: w(ixG^T,nw),csound2(ixG^T)
integer:: ix^L
!-----------------------------------------------------------------------------

if(eqpar(gamma_)<=zero)&
   call die('Correct Getcsound2prim for NONIDEAL gas in vacphys.t.mhd')

csound2(ix^S)=eqpar(gamma_)*w(ix^S,p_)/w(ix^S,rho_)

return 
end

!=============================================================================
subroutine getcsound2(w,ix^L,csound2)

! Calculate the square of the thermal sound speed csound2 within ix^L.
! csound2=gamma*p/rho

include 'vacdef.f'

double precision:: w(ixG^T,nw),csound2(ixG^T)
integer:: ix^L
!-----------------------------------------------------------------------------

if(eqpar(gamma_)<=zero)&
   call die('Correct Getcsound2 for NONIDEAL gas in vacphys.t.mhd')

oktest=index(teststr,'getcsound2')>=1
if(oktest) write(*,*)'Getcsound2'

call getpthermal(.true.,w,ix^L,csound2)
if(oktest) write(*,*)'p(ixtest)=',csound2(ixtest^D)
csound2(ix^S)=eqpar(gamma_)*csound2(ix^S)/w(ix^S,rho_)

return 
end

!=============================================================================
subroutine getpthermal(clipping,w,ix^L,p)

!!! This subroutine should not use tmp,tmp2

! Calculate thermal pressure within ix^L
! p=(g-1)*(e-0.5*(m**2/r+b**2))
! where g is the adiabatic index gamma, e is the total energy density,
! m is the momentum density and b is the magnetic field
! If clipping is .true. clip off negative pressures

include 'vacdef.f'

double precision:: w(ixG^T,nw),p(ixG^T)
integer:: ix^L
logical:: clipping
!-----------------------------------------------------------------------------

if(eqpar(gamma_)<=zero)&
   call die('Correct GetPthermal for NONIDEAL gas in vacphys.t.mhd')

oktest=index(teststr,'getpthermal')>=1

if(oktest) write(*,*)'GetPthermal'

! First calculate kinetic energy*2=m**2/rho

p(ix^S)=( ^C&w(ix^S,m^C_)**2+ )/w(ix^S,rho_)
if(oktest) write(*,*)'p(ixtest)=',p(ixtest^D)

! Add magnetic energy*2=b**2
p(ix^S)=p(ix^S)+ ^C&w(ix^S,b^C_)**2+
if(oktest) write(*,*)'p(ixtest)=',p(ixtest^D)

! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
p(ix^S)=(eqpar(gamma_)-1)*(w(ix^S,e_)-half*p(ix^S))
if(oktest) write(*,*)'p(ixtest)=',p(ixtest^D)

if(smallp<zero)then
   smallp=minval(p(ix^S))*smallpcoeff
   {^IFMPI call mpiallreduce(smallp,MPI_MIN)}
   if(oktest) write(*,*)'smallp, smallpcoeff:',smallp,smallpcoeff
   if(oktest) write(*,*)'minval(p):',minval(p(ix^S))
   if(oktest) write(*,*)'ix-limits:',ix^L
   if(smallp<zero)then
      write(*,*)'Initial condition contains negative pressure!'
      smallp=maxval(p(ix^S))*smallpcoeff
      {^IFMPI call mpiallreduce(smallp,MPI_MAX)}
      write(*,*)'Using smallp=',smallp
   endif
endif

! Clip off negative pressure if clipping is set
if(clipping)p(ix^S)=max(p(ix^S),smallp)

!!!if(clipping)p(ix^S)=abs(p(ix^S)) Hans-s solution for negative pressure

return 
end

!=============================================================================
subroutine getptotal(w,ix^L,p)

! Calculate total pressure within ix^L including magnetic pressure
! p=(g-1)*e-0.5*(g-1)*m**2/rho+(1-0.5*g)*b**2

!!! This subroutine should not use tmp

include 'vacdef.f'

double precision::  w(ixG^T,nw),p(ixG^T),gamma
integer:: ix^L
!-----------------------------------------------------------------------------

if(eqpar(gamma_)<=zero)&
   call die('Correct GetPtotal for NONIDEAL gas in vacphys.t.mhd')

oktest=index(teststr,'getptotal')>=1
if(oktest) write(*,*)'GetPtotal'

gamma=eqpar(gamma_)
p(ix^S)=(one-half*gamma)*( ^C&w(ix^S,b^C_)**2+ )

! Add contribution of total energy=(g-1)*e
p(ix^S)=p(ix^S)+(gamma-one)*w(ix^S,e_)
if(oktest) write(*,*)'p(ixtest)=',p(ixtest^D)

! Subtract contribution of kinetic energy=0.5*(g-1)*m**2/rho
p(ix^S)=p(ix^S)-half*(gamma-one)*(^C&w(ix^S,m^C_)**2+)/w(ix^S,rho_)

if(oktest) write(*,*)'p(ixtest)=',p(ixtest^D)

return 
end

