!=============================================================================
subroutine setnozzle(r,rC)

! Set global geometric factors for a nozzle:
!
!   area(ixG^T)     cross section at cell centers.
!   areaC(ixG^T)    cross section at cell boundaries.
!   areadx(ixG^T)   cell volume (integrated cross section).
!   areaside(ixG^T) =(areaC_i+1/2 - areaC_i-1/2)/dr
!
! The input arrays 'r' and 'rC' are the locations of cell centers and
! cell boundaries, respectively.

! This is Yee's nozzle with a cross section 1.398+0.347*tanh(0.8*x-4) 

include 'vacdef.f'

double precision:: r(IXGLO1-1:IXGHI1+1),rC(IXGLO1-1:IXGHI1+1)
double precision:: qareaC(IXGLO1-1:IXGHI1+1)
!----------------------------------------------------------------------------

oktest=index(teststr,'setnozzle')>=1

! ADAPTOR does not know tanh function, so tanh(x)=1-2/(exp(2*x)+1)
area(ixG^LIM1:) =1.398+0.347*(1-2/(exp(1.6*r(ixG^LIM1:) -8)+1))
qareaC(ixGmin1-1:ixGmax1)=&
                 1.398+0.347*(1-2/(exp(1.6*rC(ixGmin1-1:ixGmax1)-8)+1))
areaC(ixG^LIM1:)=qareaC(ixG^LIM1:)
! We need to integrate areaC for exact cell volume, e.g. with Mathematica
areadx(ixG^LIM1:)=1.398*(rC(ixG^LIM1:)-rC(ixG^LIM1-1:))+0.347/0.8* &
   log(cosh(0.8*rC(ixG^LIM1:)-4.)/cosh(0.8*rC(ixG^LIM1-1:)-4.))
! We use qareaC to calculate areaside
areaside(ixG^LIM1:)=(qareaC(ixG^LIM1:)-qareaC(ixG^LIM1-1:))/areadx(ixG^LIM1:)

if(oktest)write(*,*)'SetNozzle area,areaC,areadx,areaside:',&
   area(ixtest1),areaC(ixtest1),areadx(ixtest1),areaside(ixtest1)

return 
end
!=============================================================================
