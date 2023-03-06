!##############################################################################
! module vacproc.projectb - Projection of B for mhd(iso) in 2 or 3D

!=============================================================================
subroutine projectb(w)

! Project B according to B'=B-grad phi, where Laplace phi=div B within ixM

include 'vacdef.f'

double precision:: w(ixG^T,nw)

double precision:: divb(ixG^T),phi(ixG^T),divbmax,qdivbmax,qdivbmin
integer:: ix^L,iB,idim,info,matvecmax
character*3:: typestop

save divbmax,matvecmax,typestop
data divbmax /0.D0/
!-----------------------------------------------------------------------------

oktest=index(teststr,'projectb')>=1
if(oktest)write(*,*)'ProjectB w:',w(ixtest^D,iwtest)

! Determine the div B error in ix=ixM+1
call getdivb(w,ixM^L,divb)

! Calculate and save divbmax, the maximum div B allowed
if(divbmax<smalldouble)then
   divbmax=maxval(abs(divb(ixM^S)))
   {call mpiallreduce(divbmax,MPI_MAX) ^IFMPI}
   if(verbose)write(*,*)'ProjectB it,divbmax:',it,divbmax
   ! Default values if procpar was not set
   if(procpar(divbcoeff_)==-one)procpar(divbcoeff_)=one/10
   if(procpar(divbconst_)==-one)procpar(divbconst_)=zero

   ! Determine relative or absolute limit for div B
   if(procpar(divbcoeff_)<zero.and.procpar(divbcoeff_)>-one)then
      ! 1 > -divbcoeff > 0 means relative stopping criterion
      divbmax=-procpar(divbcoeff_)
      typestop='rel'
      if(verbose)write(*,*)'Reduction factor for div B:',divbmax
   else
      divbmax=divbmax*max(zero,procpar(divbcoeff_))+&
                      max(zero,procpar(divbconst_))
      typestop='max'
      if(verbose)write(*,*)'Allowed maximum for div B:',divbmax
      if(divbmax<smalldouble) call die(&
         'Error in ProjectB: Too small value for divbmax')
   endif

   ! -divbcoeff>1 or -divbconst>1 gives maximum number of iterations
   matvecmax=1000
   if(procpar(divbconst_)<-one)matvecmax=nint(-procpar(divbconst_))
   if(procpar(divbcoeff_)<-one)matvecmax=nint(-procpar(divbcoeff_))

   if(verbose)write(*,*)'Maximum number of matvecs:',matvecmax

endif

if(oktest)then
   call getdivb(w,ixM^L,divb)
   qdivbmax=maxval(divb(ixM^S))
   qdivbmin=minval(divb(ixM^S))
   if(verbose)write(*,*)'Max and min of divb:',qdivbmax,qdivbmin
   ! if(oktest)write(*,*)'Indices for max:',maxloc(divb(ix^S))
endif

! Determine boundary condition for the Poisson solver based on 
! the bpundary condition for the normal component of the magnetic field
do iB=1,nB
   select case(typeB(b0_+idimB(iB),iB))
   case('periodic')
     typeBscalar(iB)='periodic'
   case('symm','symm0')
     typeBscalar(iB)='asymm'
   case('asymm')
     typeBscalar(iB)='symm'
   case('fixed','fixed1')
     typeBscalar(iB)='grad0'
   {^IFMPI
   case('mpi')
     typeBscalar(iB)='mpi'
   case('mpiperiod')
     typeBscalar(iB)='mpiperiod' \}
   case default
     typeBscalar(iB)='nul'
   end select
   if(it==itmin.and.oktest)write(*,*)'iB,idim,typeB,typeBscalar:',&
                        iB,idimB(iB),typeB(b0_+idimB(iB),iB),typeBscalar(iB)
enddo 


! Initial guess for phi is zero for the iterative solvers (nonzero='false')
phi(ixM^S)=zero

! Solve the Poisson equation
{^IFPOISSON
call poisson('project B ',divb,divbmax,typestop,matvecmax,info,.false.,phi)
!} call die('Error: Poisson solver is OFF! setvac -on=poisson; make vac')

if(oktest)write(*,*)'Poisson solver info:',info

! Do not do anything if the initial guess satisfied the stopping criterion
if(info==3)return

! Do not subtract grad(phi) from B if the iterations did not reduce the error
if(info<0)return

! Subtract tmp=grad(phi) from the first ndim components of the B field in ixM+1

! First get the ghost cell values for the solution
call boundscalar(phi);

! For the full MHD equations correct total energy according to the value of
! |procpar(3)|=1, 2, or 3
if(fourthorder)then
   ix^L=ixM^L;
else
   ix^L=ixM^L^LADD1;
endif
do idim=1,ndim
   if(fourthorder)then
      call gradient4(.true.,phi,ix^L,idim,tmp)
   else
      call gradient(.true.,phi,ix^L,idim,tmp)
   endif
   if(typephys=='mhd')then
      select case(nint(abs(procpar(divbbound_))))
      case(2)
         ! Correct total energy so that thermal pressure is kept constant
         w(ix^S,ee_)=w(ix^S,ee_)-w(ix^S,b0_+idim)*tmp(ix^S)+tmp(ix^S)**2/2
      case(3)
         ! Correct total energy so that total pressure is kept constant
         if(eqpar(gamma_)/=two.and.eqpar(gamma_)/=one)&
         w(ix^S,ee_)=w(ix^S,ee_)+(eqpar(gamma_)-two)/(eqpar(gamma_)-one)*&
                          (-w(ix^S,b0_+idim)*tmp(ix^S)+tmp(ix^S)**2/2)
      end select
   endif
   ! B'=B-grad(Phi)
   w(ix^S,b0_+idim)=w(ix^S,b0_+idim)-tmp(ix^S)
end do

! Recalculate boundaries for the first ndim components of B if 
! procpar(divbbound_) is positive
if(procpar(divbbound_)>zero) call getboundary(t,b1_,b^ND_,1,ndim,w)

if(oktest)then
   call getdivb(w,ixM^L,divb)
   qdivbmax=maxval(divb(ixM^S))
   qdivbmin=minval(divb(ixM^S))
   if(verbose) write(*,*)'New   extrema of divb:',qdivbmax,qdivbmin
endif

return
end

!=============================================================================
! end module vacproc.projectb
!##############################################################################
