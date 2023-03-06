!##############################################################################
! module vaciter

subroutine solveiter(qnimpl,rhs,dw,nonzero,method,matvec,&
   iter,resid,typestop,info,qwork,qnwork)

! Solves A.dw=rhs with the method given, the matvec subroutine defines A.
!
! resid is the dimensionless relative error
! resid

include 'vacdef.f'

integer:: qnimpl,iter,info,qnwork
double precision:: rhs(qnimpl),dw(qnimpl),resid,qwork(qnwork)
logical:: nonzero
character*^LENTYPE :: method
character*3:: typestop
external matvec

integer:: nworkmin
logical:: okprint,initialized
data initialized/.false./
!-----------------------------------------------------------------------------

oktest=index(teststr,'solveiter')>=1
okprint=index(teststr,'printiter')>=1
if(oktest)write(*,*)'SolveIter method,typestop:',method,' ',typestop

if(.not.initialized)then
   initialized=.true.
   ! Check for the size of the needed work array
   select case(method)
   case('vac_cg')
      nworkmin=2
   case('vac_bicg')
      nworkmin=7
   case('vac_bicgl')
      nworkmin=5+2*implrestart2
   case('vac_gmres','vac_mrpc')
      nworkmin=1+implrestart
   case('vac_gcr')
      nworkmin=2+2*implrestart
   case('vac_gmresr')
      nworkmin=2+2*implrestart+implrestart2
   end select

   nworkmin=nworkmin*qnimpl

   if(nworkmin>qnwork)then
       write(*,*)'Not enough work space left for iterative method ',method
       write(*,*)'Change parameters or increase nwork in src/vacdef.t and'
       write(*,*)'recompile VAC. Minimum value for nwork:',&
           nwork+nworkmin-qnwork
       call die(' ')
   endif

endif

if(oktest)write(*,*)'Before ',method,' iter,resid:',iter,resid

select case(method)
case('vac_cg')
   call cg77(okprint,matvec,rhs,dw,nonzero,qnimpl,resid,typestop,&
      iter,qwork,qwork(qnimpl+1),info)
case('vac_bicg')
   call bicgstabl(okprint,matvec,rhs,dw,nonzero,qnimpl,resid,typestop,&
      1,iter,qwork,info)
case('vac_bicgl')
   call bicgstabl(okprint,matvec,rhs,dw,nonzero,qnimpl,resid,typestop,&
      implrestart2,iter,qwork,info)
case('vac_gmresr')
   call vacgmresr(okprint,qnimpl,rhs,dw,qwork,implrestart,implrestart2,&
      iter,resid,matvec,typestop,info)
case('vac_gcr')
   call vacgmresr(okprint,qnimpl,rhs,dw,qwork,implrestart,0,iter,resid,&
          matvec,typestop,info)
case('vac_gmres')
   call vacgmres(okprint,qnimpl,rhs,dw,nonzero,qwork,implrestart,iter,&
          resid,matvec,typestop,info)
case('vac_mrpc')
   dtmrpc=dt
   call vacmrpc(okprint,qnimpl,rhs,dw,nonzero,qwork,implrestart,resid,&
          matvec,typestop,dtmrpc,implmrpcpar,typeimplinit,info)
   if(implmrpcpar>=0.and.typeimplinit/='explicit2')dt=dtmrpc
case default
      call die('Unknown type of iterative method:'//method)
end select

if(oktest)write(*,*)'iter,info,resid:',iter,info,resid

return
end

!=======================================================================

! end module vaciter
!##############################################################################


