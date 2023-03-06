!##############################################################################
! module vacpoisson 
! Poisson solver and iterative schemes 
! for a scalar quantity defined on the mesh
!=============================================================================
subroutine poisson(purpose,rhs,tolerance,typestop,matvecmax,info,nonzero,phi)

! Solve grad div phi = rhs for the purpose defined in that string
! The right hand side of the Poisson equation is in rhs
! The initial guess for the solution is given by phi and nonzero
!    nonzero=.true. : the initial guess is non-zero and it is in phi
!    nonzero=.false.: the initial guess is zero and phi=0
! The required accuracy of the solution is given by tolerance and typestop:
!    typestop='max': resid=max(abs(Laplace(phi)-rhs)) < tolerance
!    typestop='abs': resid=sum((Laplace(phi)-rhs)**2) < tolerance
!    typestop='rel': resid=sum((Laplace(phi)-rhs)**2) < tolerance*resid_init
! At most matvecmax matrix vector multiplications can be performed
! The success of the iteration is given by returned value of info:
!     abs(info)=  0 - solution found satisfying given tolerance.
!                 1 - iteration aborted due to division by very small value.
!                 2 - no convergence within maximum number of iterations.
!                 3 - initial guess satisfies the stopping criterion.
!    sign(info)=  + - residual decreased
!                 - - residual did not reduce

include 'vacdef.f'

double precision:: rhs(ixG^T),phi(ixG^T),tolerance
integer:: matvecmax,info
character*^LENTYPE :: purpose
character*3:: typestop
logical:: nonzero

double precision:: resid
integer:: ix^L,idim,iter
logical:: okprint

external matvec_poisson
!-----------------------------------------------------------------------------

oktest=index(teststr,'poisson')>=1
okprint=index(teststr,'printiter')>=1
if(oktest)write(*,*)'Poisson solver called for ',purpose

if(typepoisson=='default')then
   ! Select method for iterative solver depending on the grid type
   if(typeaxial=='slab'.and. .not.gencoord)then
      typepoisson='cg'
   else
      typepoisson='bicgstab'
   endif
endif
if(it==itmin.and.verbose) write(*,*) &
  'Using ',typepoisson,' for solving the Poisson equation for ',purpose

! Initialize parameters and phi for the iterative solvers
iter=matvecmax
resid=tolerance

! Solve the Poisson equation
select case(typepoisson)
case('cg')
   call cgscalar(okprint,rhs,ixM^L,nonzero,phi,matvec_poisson, &
      iter,resid,typestop,info)
case('bicgstab')
   call bicgstabscalar(okprint,rhs,ixM^L,nonzero,phi,matvec_poisson, &
      iter,resid,typestop,info)
case('minres')
   call minresscalar(okprint,rhs,ixM^L,nonzero,phi,matvec_poisson, &
      iter,resid,typestop,info)
case default
   call die('Error in Poisson: Unknown type of iterative method:'// &
      typepoisson)
end select

if(oktest)write(*,*)'Poisson info,nmatvec,resid',info,iter,resid

if(info/=0.and.info/=3)then
   nerror(poissonerr_)=nerror(poissonerr_)+1
   if(nerror(poissonerr_)==1)then
      write(*,*)'No convergence for Poisson eq. for ',purpose
      write(*,"(a,i2,a,i5,a,i2,a,i5,a,g10.3)")&
         'Error code=',poissonerr_,' it=',it,' info=',info,' iter=',iter,&
         ' resid=',resid
      select case(abs(info))
      case(1)
         write(*,*)'Breakdown due to division by a small value'
      case(2)
         write(*,*)'No convergence within maximum number of iterations'
      end select
      if(info>0)write(*,*)'The residual decreased'
      if(info<0)write(*,*)'The residual did not decrease'
   end if
end if

return
end

!=============================================================================
subroutine matvec_poisson(qx,qy)

! Calculate qy=laplace(qx) for the Poisson solvers

include 'vacdef.f'

double precision:: qx(ixG^T),qy(ixG^T)

double precision:: qdx(ixG^T),qddx(ixG^T)

integer:: ix,ix^D,ix^L,idim
!-----------------------------------------------------------------------------

oktest=index(teststr,'matvec')>=1
if(oktest)write(*,*)'Matvec_Poisson'

if(oktest)write(*,*)'before bound qx:',qx(ixtest^D)

call boundscalar(qx);

if(oktest)write(*,*)'after  bound qx:',qx(ixtest^D)

if(index(teststr,'uselap4')>0)then
   call laplace4(qx,ixM^L,qy)
   return
endif

if(fourthorder)then
   qy(ixM^S)=zero
   do idim=1,ndim
      call gradient4(.true. ,qx ,ixM^L ,idim,qdx )
      call boundgradient(idim,qdx)
      call gradient4(.false.,qdx,ixM^L,idim,qddx)
      qy(ixM^S)=qy(ixM^S)+qddx(ixM^S)
   enddo
   return
endif

! Calculate y=Laplace(phi)
qy(ixM^S)=zero

ix^L=ixM^L^LADD1;
do idim=1,ndim
   ! qdx=d qx/dx_idim (gradient)
   call gradient(.true. ,qx ,ix^L ,idim,qdx )
   ! qddx=d qdx/dx_idim (divergence of gradient)
   call gradient(.false.,qdx,ixM^L,idim,qddx)
   ! qy=Laplace(qx)=Sum_idim(qddx)
   qy(ixM^S)=qy(ixM^S)+qddx(ixM^S)

   if(oktest)write(*,*)'idim,qdx,qddx,qy:',idim,&
                        qdx(ixtest^D),qddx(ixtest^D),qy(ixtest^D)
enddo

return
end

!=============================================================================
subroutine boundscalar(phi)

! Calculate boundary for phi based on typeBscalar

include 'vacdef.f'

double precision:: phi(ixG^T)
integer:: ix,ixe,ixf,ix^L,ixpair^L,idim,iB
!-----------------------------------------------------------------------------

oktest=index(teststr,'boundscalar')>=1
if(oktest)write(*,*)'BoundScalar phi:',phi(ixtest^D)

{^IFMPI 
! Get boundaries from other PE-s
call mpibound(1,phi)
}

do iB=1,nB
   idim=idimB(iB)
   ix^L=ixB^LIM(^D,iB);

   if(oktest)write(*,*)'iB,idim,ixL,typeB',iB,idim,ix^L,typeBscalar(iB)

   select case(typeBscalar(iB))
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
      phi(ix^S)=phi(ixpair^S)
   case('cont')
      ! ghost cells = edge
      select case(idim)
      {case(^D)
         if(upperB(iB))then
           ixe=ixmin^D-1
         else
           ixe=ixmax^D+1
         endif
         !HPF$ INDEPENDENT
         do ix= ix^DL
            phi(ix^D%ix^S)=phi(ixe^D%ix^S)
         end do 
      \}
      end select
   case('cont1')
      ! ghost cells are extrapolated from edge
      select case(idim)
      {case(^D)
         if(upperB(iB))then
           ixe=ixmin^D-1; ixf=ixe-1
         else
           ixe=ixmax^D+1; ixf=ixe+1
         endif
         !HPF$ INDEPENDENT
         do ix= ix^DL
            phi(ix^D%ix^S)=(abs(ix-ixe)+1)*phi(ixe^D%ix^S)-&
                            abs(ix-ixe)   *phi(ixf^D%ix^S)
         end do 
      \}
      end select
   case('symm')
      select case(idim)
      {case(^D)
         if(upperB(iB))then
            ixe=2*ixmin^D-1
         else
            ixe=2*ixmax^D+1
         endif
         !HPF$ INDEPENDENT
         do ix= ix^DL
            phi(ix^D%ix^S)=+phi(ixe-ix^D%ix^S)
         end do
      \}
      end select
   case('asymm')
      select case(idim)
      {case(^D)
         if(upperB(iB))then
            ixe=2*ixmin^D-1
         else
            ixe=2*ixmax^D+1
         endif
         !HPF$ INDEPENDENT
         do ix= ix^DL
            phi(ix^D%ix^S)=-phi(ixe-ix^D%ix^S)
         end do
      \}
      end select
   case('grad0')
      ! phi should have 0 gradient in all directions at the first ghost cells,
      ! therefore phi is 0 at 1st, and copy of mesh edge at 2nd ghost cells
      ! In generalized coordinates and/or in axial symmetry, multiply by 
      ! the ratio of surfaces on the two sides of the 1st ghost cell.
      select case(idim)
      {case(^D)
         if(upperB(iB))then
            phi(ixmin^D^D%ix^S)=zero
            if(gencoord)then
               phi(ixmin^D+1^D%ix^S)=phi(ixmin^D-1^D%ix^S)*&
                  surfaceC(ixmin^D-1^D%ix^S,idim)/surfaceC(ixmin^D^D%ix^S,idim)
            else if(typeaxial/='slab'.and.idim==r_)then
               phi(ixmin^D+1^D%ix^S)=phi(ixmin^D-1^D%ix^S)*&
                  areaC(ixmin1-1)/areaC(ixmin1)
            else
               phi(ixmin^D+1^D%ix^S)=phi(ixmin^D-1^D%ix^S)
            endif
         else
            phi(ixmax^D^D%ix^S)=zero
            if(gencoord)then
               phi(ixmax^D-1^D%ix^S)=phi(ixmax^D+1^D%ix^S)*&
                  surfaceC(ixmax^D^D%ix^S,idim)/surfaceC(ixmax^D-1^D%ix^S,idim)
            else if(typeaxial/='slab'.and.idim==1)then
               phi(ixmax^D-1^D%ix^S)=phi(ixmax^D+1^D%ix^S)*&
                  areaC(ixmax1)/areaC(ixmax1-1)
            else
               phi(ixmax^D-1^D%ix^S)=phi(ixmax^D+1^D%ix^S)
            endif
         endif
      \}
      end select
   case('nul')
       phi(ix^S)=zero
   {^IFMPI
   case('mpi','mpiperiod')
       ! This boundary is handled by MPI\}
   case default
       write(*,*)'Error in BoundScalar, unknown boundary type:', &
                typeBscalar(iB),' iB=',iB
       call die('Correct parameter file')
   end select
end do

if(oktest)write(*,*)'final phi:',phi(ixtest^D)

return
end

!=============================================================================
subroutine boundgradient(idim,phi)

! Calculate boundary for gradient of phi based on typeBscalar 
! and the direction idir in which the gradient was taken
! This is only needed for fourth order scheme. Note that the symm and antisymm
! conditions are reversed for the gradient of phi. 

include 'vacdef.f'

double precision:: phi(ixG^T)
integer:: ix,ixe,ix^L,ixpair^L,idim,iB
!-----------------------------------------------------------------------------

oktest=index(teststr,'boundgrad')>=1
if(oktest)write(*,*)'BoundGradient grad phi:',phi(ixtest^D)

{^IFMPI 
! Get boundaries from other PE-s
call mpibound(1,phi)
}

do iB=1,nB
   if(idim==idimB(iB))then
      ix^L=ixB^LIM(^D,iB);

      if(oktest)write(*,*)'iB,idim,ixL,typeB',iB,idim,ix^L,typeBscalar(iB)

      select case(typeBscalar(iB))
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
         phi(ix^S)=phi(ixpair^S)
      case('symm')
         select case(idim)
         {case(^D)
         if(upperB(iB))then
            ixe=2*ixmin^D-1
         else
            ixe=2*ixmax^D+1
         endif
         !HPF$ INDEPENDENT
         do ix= ix^DL
            phi(ix^D%ix^S)=-phi(ixe-ix^D%ix^S)
         end do
         \}
         end select
      case('asymm')
         select case(idim)
         {case(^D)
         if(upperB(iB))then
            ixe=2*ixmin^D-1
         else
            ixe=2*ixmax^D+1
         endif
         !HPF$ INDEPENDENT
         do ix= ix^DL
            phi(ix^D%ix^S)=+phi(ixe-ix^D%ix^S)
         end do
         \}
         end select
      case('cont','nul')
         phi(ix^S)=zero
      {^IFMPI
      case('mpi','mpiperiod')
         ! This boundary is handled by MPI\}
      case default
         write(*,*)'Error in BoundGradient, unknown boundary type:', &
                typeBscalar(iB),' iB=',iB
         call die('Correct parameter file')
      end select
   endif
end do

if(oktest)write(*,*)'final grad phi:',phi(ixtest^D)

return
end

!=============================================================================
subroutine cgscalar(okprint,rhs,ix^L,nonzero,qx,matvec,iter,tol,typestop,info)

! The CG-algorithm is implemented as shown on page 12 of the thesis
! "Preconditioning for sparse matrices with applications."
! Auke van der Ploeg, University of Groningen, 1994.
! Rewritten to F90 by G. Toth based on the F77 subroutine in src/conjgrad.f 

! This subroutine determines the solution of A.QX=RHS, where
! the matrix-vector multiplication with A is performed by 
! the subroutine 'matvec'.

!!! If the matrix is not symmetric positive definite, CG is likely to fail.

!     Description of arguments:

!     okprint: (input) (boolean) 
!        Determines whether of not output is printed.
!     matvec: external routine for matrix-vector multiplication.
!     rhs: (input/output)
!        on input:  right-hand side vector.
!        on output: residual vector.
!     ixL:
!        Region of unknowns within the computational grid ixG
!     nonzero: (input) (boolean)
!        Tells CG if initial guess in qx is zero or not. 
!        If nonzero is .FALSE., one MATVEC call is saved.
!     qx: (input/output)
!        on input:  initial guess for the solution vector.
!        on output: solution vector.
!     matvec: (subroutine)
!        performes the action of the A matrix
!     iter: (input/output)
!       on input:  maximum number of iterations to be performed.
!       on output: actual  number of iterations done.
!     tol: (input/output) 
!       on input:  required (relative) 2-norm or maximum norm of residual
!       on output: achieved (relative) 2-norm or maximum norm of residual
!     typestop (input) (character*3)
!       Determine stopping criterion (||.|| denotes the 2-norm):
!       typestop='rel'    -- relative stopping crit.: ||res|| <= tol*||res0||
!       typestop='abs'    -- absolute stopping crit.: ||res|| <= tol
!       typestop='max'    -- maximum  stopping crit.: max(abs(res)) <= tol
!     info (output)
!         Gives reason for returning:
!     abs(info)=  0 - solution found satisfying given tolerance.
!                 1 - iteration aborted due to division by very small value.
!                 2 - no convergence within maximum number of iterations.
!                 3 - initial guess satisfies the stopping criterion.
!    sign(info)=  + - residual decreased
!                 - - residual did not reduce

include 'vacdef.f'

integer:: ix^L,iter,info
double precision:: qx(ixG^T),rhs(ixG^T),tol
logical:: okprint,nonzero
character*3:: typestop

integer:: i,itr,matv
double precision:: rho,rhonew,res,res0,bet,alf,assumedzero

external matvec
!----------------------------------------------------------------------------
if(okprint)write(*,*)'CGscalar tol,mxmv,ixL:',tol,iter,ix^L

assumedzero=1.D-16; itr=0; matv=0

if (typestop/='rel'.and.typestop/='abs'.and.typestop/='max') then
   write(*,*)'Error in CG:'
   call die('Parameter typestop='//typestop//' should be one of rel/abs/max')
end if

if(okprint) write(*,*)'n gives the number of CG-iterations.'

! Calculate the initial residual R:=RHS-A*X and its 2-norm.

if (nonzero) then
   call matvec(qx,tmp)
   matv = matv + 1
   rhs(ix^S)=rhs(ix^S)-tmp(ix^S)
endif

! rho=||rhs||
rho=sum(rhs(ix^S)**2)
{call mpiallreduce(rho,MPI_SUM) ^IFMPI}
rho=sqrt(rho)

res0=rho
if (typestop=='max')then
   res0=maxval(abs(rhs(ix^S)))
   {call mpiallreduce(res0,MPI_MAX) ^IFMPI}
endif
res=res0

assumedzero = assumedzero*res0

if (okprint) then
   IF (typestop=='max') THEN
      write(*,*)'n:',itr,' Maximum norm initial residual:',res0
   ELSE
      write(*,*)'n:',itr,' 2-norm intial residual:',res0
   END IF
end if

if (res0<smalldouble.or.(typestop/='rel'.and.res0<=tol)) then
   info = 3
else
   ! Initialize rho and tmp=Z
   rho=rho*rho
   tmp(ix^S)=rhs(ix^S)

   ! Do iteration
   do
      ! AZ=A.Z
      call matvec(tmp,tmp2) 
      matv=matv+1

      ! alf=A.AZ
      alf=sum(tmp(ix^S)*tmp2(ix^S))
      {call mpiallreduce(alf,MPI_SUM) ^IFMPI}
      if(okprint)write(*,*)'alf=',alf

      if (abs(alf)<=assumedzero**2) then
         info = 1
         exit
      end if
      alf=rho/alf
      if(okprint)write(*,*)'alf=',alf

      qx(ix^S) =qx(ix^S)  + alf*tmp(ix^S)
      rhs(ix^S)=rhs(ix^S) - alf*tmp2(ix^S)

      ! rhonew=||rhs||
      rhonew=sum(rhs(ix^S)**2)
      {call mpiallreduce(rhonew,MPI_SUM) ^IFMPI}
      rhonew=sqrt(rhonew)
      if(okprint)write(*,*)'rhonew=',rhonew

      select case(typestop)
      case('max')
         res=maxval(abs(rhs(ix^S)))
         {call mpiallreduce(res,MPI_MAX) ^IFMPI}
      case('rel')
         res=rhonew/res0
      case('abs')
         res=rhonew
      end select
      rhonew=rhonew*rhonew
      itr=itr+1
      IF (okprint) write(*,*)'n:',itr,' ',typestop,'. norm of residual:',res

      if (res<=tol) then
         info = 0
         exit
      end if
      if (itr>=iter) then
         info = 2
         exit
      end if
      if (rho<=assumedzero**2) then
         info = 1
         exit
      end if 

      bet=rhonew/rho
      if(okprint)write(*,*)'bet=',bet

      tmp(ix^S)=rhs(ix^S)+bet*tmp(ix^S)

      IF (okprint) write(*,*)'alf,bet,rho,rhonew:',alf,bet,rho,rhonew
      rho=rhonew
   enddo
endif

! return number of iterations and achieved residual
iter=itr
tol =res
if((typestop=='rel'.and.res>one).or.(typestop/='rel'.and.res>res0))&
   info=-info

! report results
if(okprint)then
   write(*,*)'Total Number of CG-iterations:',itr
   write(*,*)'Number of matrix-vector mult.:',matv
   select case(abs(info))
   case(0)
      write(*,*)'Successful iteration, norm of res. is:',tol
   case(1)
      write(*,*)'Iteration aborted due to division by a'
      write(*,*)'very small value.' 
   case(2)
      write(*,*)'Stopping crit. not fulfilled within given'
      write(*,*)'maximum number of iterations.'
   case(3)
      write(*,*)'Initial guess for the solution satisfies'
      write(*,*)'given stopping criterion.' 
   case default
      write(*,*)'Impossible value for info:',info
   end select
   if(info<0)write(*,*)'The residual did not reduce'
endif

return
! -------------------------- end of cgscalar -----------------------------
end

!=============================================================================
subroutine bicgstabscalar(okprint,rhs,ix^L,nonzero,qx,matvec,mxmv,tol,&
         typestop,info)

! Simple BiCGstab(\ell=1) iterative method
! Modified by G.Toth from the \ell<=2 version written
! by M.A.Botchev, Jan.'98
!
! This is the "vanilla" version of BiCGstab(\ell) as described
! in PhD thesis of D.R.Fokkema, Chapter 3.  It includes two enhancements 
! to BiCGstab(\ell) proposed by G.Sleijpen and H.van der Vorst in
! 1) G.Sleijpen and H.van der Vorst "Maintaining convergence 
!    properties of BiCGstab methods in finite precision arithmetic",
!    Numerical Algorithms, 10, 1995, pp.203-223
! 2) G.Sleijpen and H.van der Vorst "Reliable updated residuals in
!    hybrid BiCG methods", Computing, 56, 1996, 141-163
!
! {{ This code is based on:
! subroutine bistbl v1.0 1995
!
! Copyright (c) 1995 by D.R. Fokkema.
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.  }}
!
! ix^L    == (input) INTEGER sizes of the system to solve 
! qx      == (input/output) DOUBLE PRECISION array dimension n
!            initial guess on input, solution on output
! rhs     == (input) DOUBLE PRECISION array dimension n
!            right-hand side (rhs) vector
!            (output) changed to initial residual=rhs-A.qx
! matvec  == (input) EXTERNAL name of matrix vector subroutine
!            to deliver y:=A*x by CALL matvec(n,x,y)
! nonzero == (input) LOGICAL tells
!            BiCGstab(\ell) if initial guess in x is zero or not. 
!            If nonzero is .FALSE., one MATVEC call is saved.
! tol     == (input/output) DOUBLE PRECISION tolerance for all possible
!            stopping criterions (see the 'typestop' parameter)
!            On output, if info=0 or 1, tol is actually achieved
!            residual reduction or whatever (see the 'typestop' parameter)
! typestop== (input) CHARACTER*3 stopping criterion (||.|| denotes 
!            the 2-norm):
!            typestop='rel' -- relative stopping crit.: ||res|| <= tol*||res0||
!            typestop='abs' -- absolute stopping crit.: ||res|| <= tol
!            typestop='max' -- maximum  stopping crit.: max(abs(res)) <= tol
! NOTE(for typestop='rel' and 'abs'): To save computational work, the value of 
!            residual norm used to check the convergence inside the main 
!            iterative loop is computed from 
!            projections, i.e. it can be smaller than the true residual norm
!            (it may happen when e.g. the 'matrix-free' approach is used).
!            Thus, it is possible that the true residual does NOT satisfy
!            the stopping criterion ('rel' or 'abs').
!            The true residual norm (or residual reduction) is reported on 
!            output in parameter TOL -- this can be changed to save 1 MATVEC
!            (see comments at the end of the subroutine)
! mxmv   ==  (input/output) INTEGER.  On input: maximum number of matrix 
!            vector multiplications allowed to be done.  On output: 
!            actual number of matrix vector multiplications done
! info   ==  (output) INTEGER. 
!     abs(info)=  0 - solution found satisfying given tolerance.
!                 1 - iteration aborted due to division by very small value.
!                 2 - no convergence within maximum number of iterations.
!                 3 - initial guess satisfies the stopping criterion.
!    sign(info)=  + - residual decreased
!                 - - residual did not reduce
!
! Memory requirement:
! nwrk (see vacdef.t) >= 2*lmax+5
! wrk    ==  (workspace) DOUBLE PRECISION array dimension (ixG^T,nwrk)
!            

! ----------------------------------------------------------
include 'vacdef.f'
!
!     .. Work Aliases ..
!
integer, parameter :: lmax=1,l=1,r=1,u=r+(l+1),xp=u+(l+1)
integer, parameter :: z=1,zz=z+(l+1),y0=zz+(l+1),yl=y0+1,y=yl+1
!
!     .. Parameters ..
! 
double precision:: qx(ixG^T),rhs(ixG^T),tol
integer :: ix^L, mxmv, info
logical :: okprint,nonzero
character*3 :: typestop
!
!     .. Matrix ..
!
external matvec
!
!     .. Local ..
!
double precision :: rwork(lmax+1,2*lmax+5)

logical GoOn, rcmp, xpdt
integer i, j, k, nmv
double precision alpha, beta, omega, rho0, rho1, sigma
double precision varrho, hatgamma
double precision assumedzero, rnrm0, rnrm, rnrmMax0, rnrmMax
double precision mxnrmx, mxnrmr, kappa0, kappal
double precision qtmp

!----------------------------------------------------------------------------
if(okprint)write(*,*)'BiCGSTABscalar tol,mxmv:',tol,mxmv

if(nwrk<2*lmax+3) then
  write(*,*) 'Error in BiCGSTABscalar: in vacdef.t adjust nwrk>=',2*lmax+3
  call die('Edit src/vacdef.t or the parameter file')
endif

info = 0

if (tol<=zero) call die('Error in BiCGSTABscalar: tolerance < 0')
if (mxmv<=1)   call die('Error in BiCGSTABscalar: maxmatvec < 1')

!
!     --- Initialize first residual
!
assumedzero = 1.d-16
if (nonzero) then
   call matvecxV ( qx, wrk,r, matvec )
   wrk(ix^S,r) = rhs(ix^S) - wrk(ix^S,r)
   nmv = 1
else
   wrk(ix^S,r) = rhs(ix^S)
   nmv = 0
endif
!
!     --- Initialize iteration loop
!

rhs(ix^S) = wrk(ix^S,r)
wrk(ix^S,xp) = qx(ix^S)

if (nonzero) then
   qx(ix^S) = zero
endif

rnrm0 = sum( wrk(ix^S,r)**2 )
{call mpiallreduce(rnrm0,MPI_SUM) ^IFMPI}
rnrm0 = sqrt(rnrm0)
rnrm = rnrm0
if(okprint) print *,'initial rnrm:',rnrm

mxnrmx = rnrm0
mxnrmr = rnrm0
rcmp = .false.
xpdt = .false.

alpha = zero
omega = one
sigma = one
rho0 =  one
!
!     --- Iterate
!
select case(typestop)
case('rel')
   GoOn = rnrm>tol*rnrm0 .and. nmv<mxmv
   assumedzero = assumedzero*rnrm0
   rnrmMax = 0
   rnrmMax0 = 0
case('abs')
   GoOn = rnrm>tol       .and. nmv<mxmv
   assumedzero = assumedzero*rnrm0
   rnrmMax = 0
   rnrmMax0 = 0
case('max')
   rnrmMax0 = maxval( abs( wrk(ix^S,r) ) )
   {call mpiallreduce(rnrmMax0,MPI_MAX) ^IFMPI}
   rnrmMax  = rnrmMax0
   if(okprint) print *,'initial rnrmMax:',rnrmMax
   GoOn = rnrmMax>tol    .and. nmv<mxmv
   assumedzero = assumedzero*rnrmMax
case default
   call die('Error in BiCGSTABScalar: unknown typestop value:'//typestop)
end select

if (.not.GoOn) then
   if(okprint) print *,'BiCGSTABScalar: nothing to do. info = ',info
   mxmv = nmv
   info = 3
   return
endif

do while (GoOn)
!
!     =====================
!     --- The BiCG part ---
!     =====================
!
   rho0 = -omega*rho0
   do k=1,l
      rho1 = sum ( rhs(ix^S)*wrk(ix^S,r+k-1) )
      {call mpiallreduce(rho1,MPI_SUM) ^IFMPI}

      if (abs(rho0)<assumedzero**2) then
         info = 1
         return
      endif
      beta = alpha*(rho1/rho0)
      rho0 = rho1
      do j=0,k-1
         wrk(ix^S,u+j) = wrk(ix^S,r+j) - beta*wrk(ix^S,u+j)
      enddo
      
      call matvecVV (wrk,u+k-1,u+k, matvec)
      nmv = nmv+1
      
      sigma = sum ( rhs(ix^S)*wrk(ix^S,u+k) )
      {call mpiallreduce(sigma,MPI_SUM) ^IFMPI}

      if (abs(sigma)<assumedzero**2) then
         info = 1
         return
      endif

      alpha = rho1/sigma
      qx(ix^S) = alpha*wrk(ix^S,u) + qx(ix^S)
      
      do j=0,k-1
         wrk(ix^S,r+j) = -alpha*wrk(ix^S,u+j+1)    &
              + wrk(ix^S,r+j)
      enddo

      call matvecVV(wrk,r+k-1,r+k, matvec)
      nmv = nmv+1
      
      rnrm = sum ( wrk(ix^S,r)**2 )
      {call mpiallreduce(rnrm,MPI_SUM) ^IFMPI}
      rnrm = sqrt( rnrm )
      
      mxnrmx = max (mxnrmx, rnrm)
      mxnrmr = max (mxnrmr, rnrm)
   enddo

   !
   !  ==================================
   !  --- The convex polynomial part ---
   !  ================================== 
   !
   !    --- Z = R'R
   !
   do i=1,l+1 
      do j=i-1,l
         qtmp = sum( wrk(ix^S,r+j)*wrk(ix^S,r+i-1) )
         {call mpiallreduce(qtmp,MPI_SUM) ^IFMPI}
         rwork(j+1,z+i-1) = qtmp
         rwork(z+i-1,j+1) = rwork(j+1,z+i-1) 
      enddo
   enddo

   !
   !   --- tilde r0 and tilde rl (small vectors)
   !
   rwork(1:l+1,zz:zz+l)   = rwork(1:l+1,z:z+l) 
   rwork(1,y0) = -one
   rwork(2,y0) = rwork(2,z)
   rwork(2,y0) = rwork(2,y0) / rwork(2,zz+1)
   rwork(l+1,y0) = zero

   rwork(1,yl) = zero
   rwork(2,yl) = rwork(2,z+l)
   rwork(2,yl) = rwork(2,yl) / rwork(2,zz+1)
   rwork(l+1,yl) = -one
   !
   !   --- Convex combination
   !
   rwork(1:l+1,y) = zero
   do j=1,l+1
      rwork(1:l+1,y) = rwork(1:l+1,y) + &
           rwork(j,yl)*rwork(1:l+1,z+j-1)
   enddo
   kappal = sqrt( sum( rwork(1:l+1,yl)*rwork(1:l+1,y) ) )
   rwork(1:l+1,y) = zero
   do j=1,l+1
      rwork(1:l+1,y) = rwork(1:l+1,y) + &
           rwork(j,y0)*rwork(1:l+1,z+j-1)
   enddo
   kappa0 = sqrt( sum( rwork(1:l+1,y0)*rwork(1:l+1,y)  ) )

   varrho = sum( rwork(1:l+1,yl)*rwork(1:l+1,y) )  
   varrho = varrho / (kappa0*kappal)
   
   hatgamma = sign(1d0,varrho)*max(abs(varrho),7d-1) * (kappa0/kappal)

   rwork(1:l+1,y0) = -hatgamma*rwork(1:l+1,yl) + rwork(1:l+1,y0)

   !
   !    --- Update
   !
   omega = rwork(l+1,y0)

   do j=1,l
      wrk(ix^S,u) = wrk(ix^S,u) - rwork(1+j,y0)*wrk(ix^S,u+j)
      qx(ix^S)    = qx(ix^S)    + rwork(1+j,y0)*wrk(ix^S,r+j-1)
      wrk(ix^S,r) = wrk(ix^S,r) - rwork(1+j,y0)*wrk(ix^S,r+j)
   enddo

   rwork(1:l+1,y) = zero
   do j=1,l+1
      rwork(1:l+1,y) = rwork(1:l+1,y) + rwork(j,y0)*rwork(1:l+1,z+j-1)
   enddo

   rnrm = sqrt( sum( rwork(1:l+1,y0)*rwork(1:l+1,y) ) )

   select case(typestop)
   case('rel')
      GoOn = rnrm>tol*rnrm0 .and. nmv<mxmv
      if(okprint) print *, nmv,' matvecs, ', &
           ' ||rn||/||r0|| =',rnrm/rnrm0
   case('abs')
      GoOn = rnrm>tol       .and. nmv<mxmv
      if(okprint) print *, nmv,' matvecs, ||rn|| =',rnrm
   case('max')
      rnrmMax = maxval( abs( wrk(ix^S,r) ) )
      {call mpiallreduce(rnrmMax,MPI_MAX) ^IFMPI}
      GoOn = rnrmMax>tol    .and. nmv<mxmv
      if(okprint) print *, nmv,' matvecs, max(rn) =',rnrmMax
   end select

enddo
!
!     =========================
!     --- End of iterations ---
!     =========================
!

qx(ix^S) = wrk(ix^S,xp) + qx(ix^S)

! --------------------- One matvec can be saved by commenting out this:
!
!     --- Check stopping criterion
!
!!call matvecxV (qx, wrk,r, matvec)
!!nmv = nmv + 1
!!wrk(ix^S,r) = rhs(ix^S) - wrk(ix^S,r)   
!!rnrm = sqrt( sum( wrk(ix^S,r)**2 ) )
! --------------------- One matvec can be saved by commenting out this

select case(typestop)
case('rel')
   if (rnrm>tol*rnrm0) info = 2
   tol = rnrm/rnrm0
case('abs')
   if (rnrm>tol) info = 2
   tol = rnrm
case('max')
   if (rnrmMax>tol) info = 2
   tol = rnrmMax
end select

if((typestop/='max'.and.rnrm>rnrm0).or.&
   (typestop=='max'.and.rnrmMax>rnrmMax0)) info=-info

mxmv = nmv

return
end  
!=======================================================================
subroutine minresscalar(okprint,rhs,ix^L,nonzero,qx,matvec,mxmv,tol,  &
                   typestop,info)
! By M.A.Botchev, Apr.'98 
! Matlab implementation of MINRES by G.Sleijpen was used as a base 
!
! parameters are standard for VAC
! info = 2 -- no convergence achieved and residual decreased
! info =-2 -- no convergence achieved and residual did not decrease
! info = 3 -- nothing to do

include 'vacdef.f'

character*3:: typestop
integer:: ix^L,info,mxmv
double precision:: rhs(ixG^T),qx(ixG^T),tol
logical:: okprint,nonzero
external matvec
!-----------------------------------------------
integer, parameter:: v=1,vold=2,vtld=3,qw=4,wtld=5,wtld2=6,r=7
integer:: nmv
double precision:: beta,btld,alpha,atld,el0,el1,el2,ci,es, &
                   rnrm,rnrm0,rnrmMax,rnrmMax0
logical:: GoOn
!-----------------------------------------------------------------
if(okprint)write(*,*)'MINRESscalar tol,mxmv,ixL:',tol,mxmv,ix^L

if (typestop/='rel'.and.typestop/='abs'.and.typestop/='max') then
   write(*,*)'Error in MINRES:'
   write(*,*)'Parameter typestop=',typestop,' should be one of rel/abs/max.'
   call die('Error in vacpoisson.t')
end if

if(nwrk<7) write(*,*) 'Error in MINRESscalar: in vacdef.t adjust nwrk>=7'

info=0

! initial residual:
if (nonzero) then
   call matvecxV ( qx, wrk,r, matvec )
   wrk(ix^S,r) = rhs(ix^S) - wrk(ix^S,r)
   nmv = 1
else
   wrk(ix^S,r) = rhs(ix^S)
   nmv = 0
endif
rnrm0 = sum( wrk(ix^S,r)**2 )
{call mpiallreduce(rnrm0,MPI_SUM) ^IFMPI}
rnrm0 = sqrt(rnrm0)
rnrm = rnrm0
rnrmMax = zero
rnrmMax0= zero

select case(typestop)
case('rel')
   GoOn = rnrm>tol*rnrm0 .and. nmv<mxmv
case('abs')
   GoOn = rnrm>tol       .and. nmv<mxmv
case('max')
   rnrmMax = maxval( abs( wrk(ix^S,r) ) )
   {call mpiallreduce(rnrmMax,MPI_MAX) ^IFMPI}
   rnrmMax0=rnrmMax
   GoOn = rnrmMax>tol    .and. nmv<mxmv
end select

if (.not.GoOn) then
   info = 3
   if(okprint) print *,'VACMINRES90: nothing to do'
   return
endif

beta = zero
btld = zero
ci = - one
es = zero
wrk(ix^S,v)     = wrk(ix^S,r)/rnrm0
wrk(ix^S,wtld2) = wrk(ix^S,v)
wrk(ix^S,vold)  = zero
wrk(ix^S,qw)     = zero
     
! main loop:
do while (GoOn)

   call matvecVV (wrk,v,vtld, matvec)
   nmv = nmv + 1
   wrk(ix^S,vtld) = wrk(ix^S,vtld) - beta*wrk(ix^S,vold)

   alpha = sum ( wrk(ix^S,v)*wrk(ix^S,vtld) )
   {call mpiallreduce(alpha,MPI_SUM) ^IFMPI}
   wrk(ix^S,vtld) = wrk(ix^S,vtld) - alpha*wrk(ix^S,v)

   beta = sum( wrk(ix^S,vtld)**2 )
   {call mpiallreduce(beta,MPI_SUM) ^IFMPI}
   beta = sqrt(beta)
   wrk(ix^S,vold) = wrk(ix^S,v)
   wrk(ix^S,v) = wrk(ix^S,vtld)/beta

   el1 = es*alpha - ci*btld
   el2 = es*beta
   atld = -es*btld-ci*alpha
   btld = ci*beta
   el0 = sqrt(atld**2 + beta**2)
   ci = atld/el0
   es = beta/el0

   wrk(ix^S,wtld)  = wrk(ix^S,wtld2) - el1*wrk(ix^S,qw)
   wrk(ix^S,wtld2) = wrk(ix^S,v)     - el2*wrk(ix^S,qw)
   wrk(ix^S,qw) = wrk(ix^S,wtld)/el0

   qx(ix^S) = qx(ix^S) + (rnrm*ci)*wrk(ix^S,qw)

   rnrm = rnrm*es

   select case(typestop)
   case('rel')
      GoOn = rnrm>tol*rnrm0 .and. nmv<mxmv
      if(okprint) print *, nmv,' matvecs, ', &
           ' ||rn||/||r0|| =',rnrm/rnrm0
   case('abs')
      GoOn = rnrm>tol       .and. nmv<mxmv
      if(okprint) print *, nmv,' matvecs, ||rn|| =',rnrm
   case('max')
      if(mod(nmv,3).eq.0) then
         call matvecxV ( qx, wrk,r, matvec )
         wrk(ix^S,r) = rhs(ix^S) - wrk(ix^S,r)
         rnrmMax = maxval( abs( wrk(ix^S,r) ) )
         {call mpiallreduce(rnrmMax,MPI_MAX) ^IFMPI}
         GoOn = rnrmMax>tol    .and. nmv<mxmv
         if(okprint) print *, nmv,' matvecs, max(rn) =',rnrmMax
      endif
   end select

enddo

select case(typestop)
case('rel')
   if (rnrm>tol*rnrm0) info = 2
   tol = rnrm/rnrm0
case('abs')
   if (rnrm>tol) info = 2
   tol = rnrm
case('max')
   if (rnrmMax>tol) info = 2
   tol = rnrmMax
end select

if((typestop/='max'.and.rnrm>rnrm0).or.&
   (typestop=='max'.and.rnrmMax>rnrmMax0)) info=-info

mxmv = nmv

return
end  
!=============================================================================
subroutine matvecxV (q, workVv,k, matvec)
! Performs v_k := A . q for BICGSTAB
include 'vacdef.f'
integer, parameter :: lmax=1
integer:: k
double precision :: q(ixG^T),workVv(ixG^T,2*lmax+5)
external matvec
!----------------------------------------------------------------------
call matvec(q,tmp)
workVv(ixM^S,k)=tmp(ixM^S)
return
end

!=======================================================================
subroutine matvecVV (workVv,k1,k2, matvec)
! Performs v_k2 := A . v_k1
include 'vacdef.f'
integer, parameter :: lmax=1
integer:: k1,k2
double precision :: workVv(ixG^T,2*lmax+5)
external matvec
!-----------------------------------------------------------------------
tmp2(ixM^S)=workVv(ixM^S,k1)
call matvecxV(tmp2, workVv,k2, matvec)
return
end
!=======================================================================
! end module vacpoisson
!##############################################################################
