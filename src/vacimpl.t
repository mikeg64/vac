!##############################################################################
! module vacimpl

INCLUDE:vaciter.t

!=============================================================================
subroutine advance_part(qdt,ixI^L,w)

! Advance w with F_imp only by qdt. Explicit time dependece is taken at
! t+implpar*dt, however.

include 'vacdef.f'

integer:: ixI^L
double precision:: qdt,w(ixG^T,nw),w1(ixG^T,nw)
logical:: sourceunsplit0
!-----------------------------------------------------------------------------

select case(typeimpl1)
   case('source')
      w1(ixG^S,1:nw)=w(ixG^S,1:nw);
      call addsource2(qdt,ixI^L,ixM^L,iw_impl,&
         t+implpar*dt,w1,t+implpar*dt,w)
   case('special')
      write(*,*)'special semi-implicit not implemented yet!!!'
      stop
   case default
      ! If the source is not to be treated implicitly pretend split sources
      if(.not.implsource)then
         sourceunsplit0=sourceunsplit
         sourceunsplit=.false.
      endif

      w1(ixG^S,1:nw)=w(ixG^S,1:nw);
      call advect1(typeimpl1,qdt,ixI^L,iw_impl,1,ndim,&
         t+implpar*dt,w1,t+implpar*dt,w)

      if(.not.implsource)sourceunsplit=sourceunsplit0
end select

return
end

!=============================================================================
subroutine advance_impl(ixI^L,w1,w)

! Within ixI^L solve the (semi-)implicit problem:
!
! w_n+1 = w_n + F_n - beta*F_n_imp + beta*F_n+1_imp
!
! where F=F(w)=dw/dt is the 'residual'. A Newton-Raphson iteration is used:
!
! w_(k=0)=w_n
!
! (I/dt-beta*dF_imp/dw).dw = (w_n-w_k)/dt + F_n - beta*F_n_imp + beta*F_k_imp
! w_k+1=w_k+dw
!
! where the linearized problem for k=0,1,2.. is solved by a direct method (1D)
! or an iterative method (2 and 3D).
!
! When the Newton-Raphson iteration is off, only the k=0 problem is solved,
! this is sufficient for a linear F_imp or for steady state problems.
!
! We use advance_expl() to calculate F
! and    advance_part() to calculate F_imp and dF_imp/dw.dw
!
! Boundaries are kept updated by advance_expl() and getboundary()

include 'vacdef.f'

integer:: ixI^L
double precision:: w(ixG^T,nw),w1(ixG^T,nw)

integer:: n,ix^D,iw,iiw,iter,info,iter_nr,matsize
double precision:: REStot(ixG^T,nw),resid,dwnrm,wnrm2,coef1,coef2

double precision:: rhs0(nhi),dw(nhi),rhs(nhi),dtold,qwold(ixG^T,nw)

save dtold,qwold

logical:: non0dw,converged,oktime
character*3 :: typestop

external matvec_impl

double precision, dimension(ixG^T,nw):: w_k,RESimp
double precision:: wnrm(nw),dtexpl,dtcoeff
common/impld/w_k,RESimp,wnrm,dtexpl,dtcoeff

double precision:: time1,timeiter,timeprec,timejac,cputime
data timeiter,timeprec,timejac /3*0.D0/
!-----------------------------------------------------------------------------

oktest=index(teststr,'implicit')>=1
if(oktest)write(*,*)'Advance_Impl'

oktime=index(teststr,'timeimpl')>=1

nwimpl=iw_impl(niw_)
nimpl=nwimpl*{nx^D*}

if(typeimplmat=='free')then
   ! Reserve memory for weps1 and weps2 in the WORK array
   matsize=2*nhi
else
   ! Reserve memory for the full implicit matrix MAT
   matsize=nwimpl*nimpl*(2*ndim+1)
endif

if(it==itmin)then
   if(residmin>zero.and.implconserv.and.it==itmin)write(*,*)&
      'Warning: implconserv=T may be inaccurate for steady state!'
endif

wnrm(1:nw)=-one
! Norm of current w_n, only needed for implicitly treated variables
do iiw=1,nwimpl; iw=iw_impl(iiw)
   wnrm(iw)=sqrt(sum(w(ixM^S,iw)**2)/(nimpl/nwimpl))
   if(wnrm(iw)<smalldouble) wnrm(iw)=one
end do
if(oktest)write(*,*)'wnrm:',wnrm

! Scale down dt for explicit advances to satisfy the CFL condition
if(dtpar>zero)then
   courantpar=half
   call getdt_courant(w,ixM^L)
   call getdt(w,ixM^L)
   call getdt_special(w,ixM^L)
   dtexpl=dt
   dt=dtpar
   dtcoeff=dtexpl/dt
else
   dtcoeff=half/courantpar
   dtexpl =dt*dtcoeff
endif

if(implpred)then
   dt    =half*dt
   dtcoeff=two*dtcoeff
endif

! For 3 level BDF2 scheme set beta=implpar for it>1.
! For it==1 the original implpar is used.
if(impl3level.and.it>itmin)implpar=(dt+dtold)/(2*dt+dtold)

if(oktest)write(*,*)'dtcoeff,dtexpl,dt:',dtcoeff,dtexpl,dt

if(nproc(1)/=0)then
   if(it-itmin==((it-itmin)/nproc(1))*nproc(1)) call process(1,1,ndim,w)
endif 

! w_(k=0) = w_n (passed through common/impld/)
w_k(ixI^S,1:nw)=w(ixI^S,1:nw)

! advance the explicitly treated variables in w_k
if(iw_semi(niw_)>0)call advance_expl(typefull1,ixI^L,iw_semi,w1,w_k)

! initialization for NR
! RESimp= dtexpl * F_n_imp (for the whole grid)
RESimp(ixI^S,1:nw)=w_k(ixI^S,1:nw)
call advance_part(dtexpl,ixI^L,RESimp)
RESimp(ixM^S,1:nw)=RESimp(ixM^S,1:nw)-w_k(ixM^S,1:nw)

if(typeimpl1=='source'.and.residmin<zero)then
   ! Time accurate for semi-implicit (eq. 13 and 14 in Paper implvac)
   ! REStot = dtexpl * F_n_expl(w_n+beta*dt_n*F_n) + dtexpl * F_n_imp

   ! Calculate predictor step w_n+beta*dt_n*F_n
   w1(ixI^S,1:nw)=w(ixI^S,1:nw)
   call advect1(typefull1,dt*implpar,ixI^L,iw_full,1,ndim,&
      t,w,t+dt*implpar,w1)
   call getboundary(t+implpar*dt,1,nw,1,ndim,w1)

   ! Calculate corrector step excluding source terms
   sourceunsplit=.false.
   REStot(ixI^S,1:nw)=w(ixI^S,1:nw)
   call advect1(typefull1,dt,ixI^L,iw_impl,1,ndim,&
      t+dt*implpar,w1,t+dt,REStot)
   sourceunsplit=.true.
   ! Add first order source terms dt_expl*F_impl
   REStot(ixM^S,1:nw)=dtcoeff*(REStot(ixM^S,1:nw)-w(ixM^S,1:nw))&
                       +RESimp(ixM^S,1:nw)

   ! Calculate (1-implpar)*F_n used on RHS of BDF2 scheme (eq 14)
   coef1=(1-implpar)/(implpar*dt)
   w1(ixM^S,1:nw)=coef1*(w1(ixM^S,1:nw)-w(ixM^S,1:nw))

else if(typefull1/=typeimpl1)then
   ! REStot = dtexpl * (F_n_exp+F_n_imp) 
   REStot(ixI^S,1:nw)=w(ixI^S,1:nw)
   call advect1(typefull1,dtexpl,ixI^L,iw_impl,1,ndim,&
      t+dt*implpar,w,t+dt*implpar,REStot)
   REStot(ixM^S,1:nw)=REStot(ixM^S,1:nw)-w(ixM^S,1:nw)
else
   ! REStot = dtexpl * (F_n_exp+F_n_imp) = RESimp
   REStot(ixM^S,1:nw)=RESimp(ixM^S,1:nw)
endif

! Calculate rhs used for iter_nr=1
n=0
if(impl3level.and.it>itmin)then
   if(typeimpl1=='source'.and.residmin<zero)then
      ! Collect RHS terms from Eq 14 in Paper implvac
      coef2=(1-implpar)/dtold
      {do ix^DB=ixMmin^DB,ixMmax^DB\}
      do iiw=1,nwimpl; iw=iw_impl(iiw)
         n=n+1
         ! RHS = (beta*F_n+alpha*(w_n-w_n-1)/dt_n-1)/wnrm for the first
         ! Newton-Raphson iteration. The BDF2 scheme implies
         ! beta+alpha=1 and beta=(dt_n+dt_n-1)/(2*dt_n+dt_n-1)
         rhs(n)=(REStot(ix^D,iw)/dtexpl+(w(ix^D,iw)-qwold(ix^D,iw))*coef2&
                -w1(ix^D,iw))/wnrm(iw)
      end do
      {enddo^D&\}
   else
      ! Collect RHS terms from Eq 8 in Paper implvac
      coef1=implpar/dtexpl
      coef2=(1-implpar)/dtold
      {do ix^DB=ixMmin^DB,ixMmax^DB\}
      do iiw=1,nwimpl; iw=iw_impl(iiw)
         n=n+1
         ! RHS = (beta*F_n+alpha*(w_n-w_n-1)/dt_n-1)/wnrm for the first
         ! Newton-Raphson iteration. The BDF2 scheme implies
         ! beta+alpha=1 and beta=(dt_n+dt_n-1)/(2*dt_n+dt_n-1)
         rhs(n)=(REStot(ix^D,iw)*coef1+(w(ix^D,iw)-qwold(ix^D,iw))*coef2)&
                /wnrm(iw)
      end do
      {enddo^D&\}
   endif
else
   {do ix^DB=ixMmin^DB,ixMmax^DB\}
      do iiw=1,nwimpl; iw=iw_impl(iiw)
         n=n+1
         ! RHS = F_n/wnrm for the first iteration
         rhs(n)=REStot(ix^D,iw)/dtexpl/wnrm(iw)
      end do
   {enddo^D&\}
endif

! Save previous timestep for 3 level scheme
if(impl3level)then
   dtold=dt
   do iiw=1,nwimpl; iw=iw_impl(iiw)
      qwold(ixM^S,iw)=w(ixM^S,iw)
   enddo
endif

if(implnewton.or.implconserv)then
   ! Calculate RHS0 used for RHS when iter_nr>1
   n=0
   {do ix^DB=ixMmin^DB,ixMmax^DB\}
      do iiw=1,nwimpl; iw=iw_impl(iiw)
         n=n+1
         !RHS0= F_n - beta*F_n_imp + w_n/dt = RHS-beta*F_n_impl+w_n/dt
         rhs0(n)=rhs(n)+&
                 (-implpar*RESimp(ix^D,iw)/dtexpl+w(ix^D,iw)/dt)/wnrm(iw)
      end do
   {enddo^D&\}
endif

if(oktest)write(*,*)'Sum REStot,RESimp,rhs:',&
   sum(RESimp(ixM^S,1:nw)),sum(REStot(ixM^S,1:nw)),sum(rhs(1:nimpl)**2)

! Initial guess for dw = w_n+1 - w_n
non0dw=.true.
select case(typeimplinit)
case('explicit2')
   ! For time accurate MR-PC the 2nd order "predictor" step.
   ! The method is defined by "typeadvance" in the par file.
   REStot(ixI^S,1:nw)=w(ixI^S,1:nw)
   call advance_expl(typefull1,ixG^L,iw_full,w1,REStot)
   n=0
   {do ix^DB=ixMmin^DB,ixMmax^DB\}
      do iiw=1,nwimpl; iw=iw_impl(iiw)
         n=n+1
         dw(n)=(REStot(ix^D,iw)-w(ix^D,iw))/wnrm(iw)
      end do
   {enddo^D&\}
case('explicit')
   ! w_n+1-w_n = dt * F_n
   dw(1:nimpl)=rhs(1:nimpl)*dt
case('scaled')
   ! Like explicit, but amplitude reduced
   ! w_n+1-w_n = dtexpl * F_n
   dw(1:nimpl)=rhs(1:nimpl)*dtexpl
case('old')
   !!! Leave dw as it is from previous timestep. Save dw if this works.
   !HPF_ if(.false.)write(*,*)'Avoiding SP2 HPF compiler bug'
case('nul')
   ! w_n+1-w_n = 0
   dw(1:nimpl)=zero
   non0dw=.false.
case default
   stop 'Unknown type for typeimplinit'
end select

if(index(teststr,'savedwrhs')>=1.and.it==itmin)then
   write(*,*)'Saving dw and rhs into binary file fort.99'
   write(99)dw(1:nimpl),rhs(1:nimpl)
   close(99)
endif

! Newton-Raphson iteration and iterative linear solver
dwnrm=bigdouble
iter_nr=0
do
   iter_nr=iter_nr+1;
   if(oktest)write(*,*)'iter_nr:',iter_nr
   if(iter_nr>impliternr)then
      write(*,*)'Newton-Raphson failed to converge iter_nr=',iter_nr
      if(residmin<smalldouble)stop
      exit
   endif
   nnewton=nnewton+1

   if(typeimplmat/='free'.and.(iter_nr==1.or.implnewmat))then
      if(oktime)time1=cputime()
      call getjacobian(matsize,work)
      if(oktime)then
         time1=cputime()-time1; timejac=timejac+time1;
         write(*,*)'GetJ.time:',time1,timejac
      endif
   endif

   if (iter_nr>1) then
      ! Caculate new RHS
      n=0
      {do ix^DB=ixMmin^DB,ixMmax^DB\}
         do iiw=1,nwimpl; iw=iw_impl(iiw)
            n=n+1
            ! RHS = (F_n - beta*F_n_imp + w_n/dt + beta*F_k_imp - w_k/dt)/wnrm
            ! use: RHS0 and RESimp = dtexpl * F_k_imp
            rhs(n)= rhs0(n)&
              +(implpar*RESimp(ix^D,iw)/dtexpl - w_k(ix^D,iw)/dt)/wnrm(iw)
         enddo
      {enddo^D&\}

      ! Initial guess for dw is always zero in later NR iterations
      dw(1:nimpl)=zero
      non0dw=.false.
   endif

   if(oktest)write(*,*)'norm of rhs:',sqrt(sum(rhs(1:nimpl)**2)/nimpl)

   if(implnewton)then
      ! link the inner iterative solver with the outer NR
      ! require inner loop more accurate than outer loop
      resid=one/10
      typestop='rel'
   else if(typeimplmat=='prec')then
      ! No normalization is needed for preconditioned solvers (???!!!)
      resid=implerror
      typestop='rel'
   else
      ! Normalize resid by dt(expl) for non-preconditioned solvers
      ! Distinguish between steady state and time accurate cases
      if(residmin<=zero)then
         resid=implerror/dt
      else
         resid=implerror/dtexpl
      endif
      typestop='abs'
   endif

   iter=impliter

   if(oktest)write(*,*)'Before ',typeimpliter,' iter,resid:',iter,resid

   if(typeimpliter=='tridiag')then
      if(oktime)time1=cputime()
      call tridiag(dt,dw,rhs,oktest)
      if(oktime)then
         time1=cputime()-time1; timeiter=timeiter+time1
         write(*,*)'Trid.time:',time1,timeiter
      endif
      info=0
      resid=zero
      iter=1
   else
      if(typeimplmat=='prec')then
         ! Preconditioning
         if(oktime)time1=cputime()
         call prepost(.true.,dw,rhs,work,matsize,iter,resid,info,oktest)
         if(oktime)then
            time1=cputime()-time1; timeprec=timeprec+time1
            write(*,*)'Prec.time:',time1,timeprec
         endif
      endif
      if(oktime)time1=cputime()
      call solveiter(nimpl,rhs,dw,non0dw,typeimpliter,matvec_impl,&
         iter,resid,typestop,info,work(matsize+1),nwork-matsize)
      if(oktime)then
         time1=cputime()-time1; timeiter=timeiter+time1
         write(*,*)'Iter.time:',time1,timeiter
      endif
      ! Postprocessing
      if(typeimplmat=='prec')&
         call prepost(.false.,dw,rhs,work,matsize,iter,resid,info,oktest)
   endif

   ! MR-PC may change the value of dt
   if(typeimpliter=='vac_mrpc')dtcoeff=dtexpl/dt

   if(oktest)write(*,*)'After iter,info,resid:',iter,info,resid

   niter=niter+iter
   if(oktest.and.info/=0)write(*,*) &
       'Advance_Impl warning: no convergence, info:',info

   ! Update w: w_k+1 = w_k + spar*dw  with spar=1 or spar<1 from
   ! backtracking (for steady state only) based on reducing the residual 
   ! F(w_k+1)<=F(w_k). Also calculates RESimp=F_imp_k+1*dtexpl and converged.
   call backtrack(dwnrm, dw, converged)
   if(index(teststr,'newtonsave')>=1)call savefile(fileout_,w_k)
   if(converged)exit
enddo

if(implpred)then
   ! w_n+1 = w_n + dt*F(w_k)
   dt=two*dt
   call advect1(typefull1,dt,ixG^L,iw_full,1,ndim,t+dt/2,w_k,t,w)
   call getboundary(t+dt,1,nw,1,ndim,w)
else if(implconserv)then
   ! To be conservative we need to recalculate w from conservative fluxes
   ! w = w_n + dt*F_n - beta*dt*F_n_imp + beta*dt*F_n+1_imp 
   !   = rhs0*dt*wnrm + beta*dt*RESimp/dtexpl
   n=0
   {do ix^DB=ixMmin^DB,ixMmax^DB\}
      do iiw=1,nwimpl; iw=iw_impl(iiw)
         n=n+1
         w(ix^D,iw)=rhs0(n)*dt*wnrm(iw)+implpar/dtcoeff*RESimp(ix^D,iw)
      enddo
   {enddo^D&\}
   call getboundary(t+dt,1,nw,1,ndim,w)
   if(index(teststr,'newtonsave')>=1)call savefile(fileout_,w)
else
   w(ixI^S,1:nw)=w_k(ixI^S,1:nw)
endif

if(oktest)write(*,*)'nmatvec=',nmatvec
if(implnewton.and.oktest)write(*,*)'Final iter_nr=',iter_nr

return
end

!=============================================================================
subroutine backtrack(dwnrm, dw, converged)

! Update w: w_k+1 = w_k + spar*dw  with spar from backtracking
! such that F(w_k+1) <= F(w_k) if possible

include 'vacdef.f'

double precision:: dwnrm,dw(nhi)
logical:: converged,compactres0

integer:: iw, iiw, n, ix^D, ires
double precision:: REStot(ixG^T,nw),spar,resold,wnrm2

double precision, dimension(ixG^T,nw):: w_k,RESimp
double precision:: wnrm(nw),dtexpl,dtcoeff
common/impld/w_k,RESimp,wnrm,dtexpl,dtcoeff
!-----------------------------------------------------------------------------

oktest=index(teststr,'backtrack')>=1
if(oktest)write(*,*)'Backtrack'

if(implnewton)then
   ! Calculate progress in NR scheme to set linear solver accuracy
   ! dwnrm=||w_k+1 - w_k||/||w_n||
   dwnrm=sqrt(sum(dw(1:nimpl)**2)/nimpl)
   converged = dwnrm<implerror
   if(oktest)write(*,*)'dwnrm:',dwnrm
else
   converged = .true.
endif

! Initial guess for spar is limited to avoid non-physical w for pseudo-time
if(residmin>zero.and.impldwlimit<bigdouble)then
  spar=min(one,impldwlimit/(maxval(abs(dw(1:nimpl)))+smalldouble))
else
  spar=one
endif

ires=0
resold=residual
do
   ires=ires+1
   if(oktest)write(*,*)'ires, spar:',ires,spar

   n=0
   {do ix^DB=ixMmin^DB,ixMmax^DB\}
      do iiw=1,nwimpl; iw=iw_impl(iiw)
         n=n+1
         w_k(ix^D,iw)=w_k(ix^D,iw)+spar*dw(n)*wnrm(iw)
      enddo
   {enddo^D&\}
   call getboundary(t+dt,1,nw,1,ndim,w_k)

   ! Exit if Newton-Raphson converged or no Newton-Raphson is done
   ! but still calculate RESimp if it is needed for implconserv
   if(converged.and..not.implconserv) exit

   !calculate implicit residual RESimp = dtexpl*F_k+1_imp

   if(converged.and.implconserv)then
      compactres0=compactres
      compactres=.false.
   endif

   RESimp(ixG^S,1:nw)=w_k(ixG^S,1:nw)
   call advance_part(dtexpl,ixG^L,RESimp)
   RESimp(ixM^S,1:nw)=RESimp(ixM^S,1:nw)-w_k(ixM^S,1:nw)

   if(converged.and.implconserv)compactres=compactres0

   if(converged) exit

   ! Do not backtrack if not a steady state calculation
   if(residmin<smalldouble)exit

   ! calculate residual ||dt*F(w_k+1)||
   if(typefull1/=typeimpl1)then
      REStot(ixG^S,1:nw)=w_k(ixG^S,1:nw)
      call advect1(typefull1,dtexpl,ixG^L,iw_impl,1,ndim,&
         t+dt,w_k,t+dt,REStot)
      REStot(ixM^S,1:nw)=REStot(ixM^S,1:nw)-w_k(ixM^S,1:nw)
   else
      REStot(ixM^S,1:nw)=RESimp(ixM^S,1:nw)
   endif
   residual=zero
   do iw=1,nw
      wnrm2=sum(w_k(ixM^S,iw)**2)
      if(wnrm2<smalldouble)wnrm2=one
      residual = residual + sum(REStot(ixM^S,iw)**2)/wnrm2
   enddo
   residual=sqrt(residual/nw)
   if(oktest)write(*,*)'resold,residual:',resold,residual

   ! Exit if backtracked towards steady-state or giving up
   if(residual<=resold .or. ires>3)exit
   spar=spar*half
end do

return
end

!=============================================================================
subroutine matvec_impl(dw,qy,ipar)

! Calculate qy= L.dw = (I/dt - beta*dF/dw).dw
! Do it matrix-free or with the matrix or preconditioned
!     depending on "typeimplmat"
! If matvec is called from the VAC library, the ipar array has one element, 
!     and it contains the number of elements in dw and qy.
! If matvec is called from the PIM library, then all the PIM integer 
!     parameters are there.

include 'vacdef.f'

double precision:: dw(nhi),qy(nhi)
integer:: ipar(*)
!-----------------------------------------------------------------------------

nmatvec=nmatvec+1

select case(typeimplmat)
case('free')
   call matvec_free(dw,qy,work,work(nhi+1))
case('with')
   call matvec_with(dw,qy,work)
case('prec')
   call matvec_prec(dw,qy,work)
case default
   stop 'Unknown value for typeimplmat'
end select

return
end

!=============================================================================
subroutine matvec_free(dw,qy,weps1,weps2)

! Calculate qy=L.dw for the iterative solver, matrix-free 
! where L= I/dt - beta*dF/dw   (dt=dt_implicit)
!
! One sided derivative:
!----------------------
! weps = w + eps*dw
!
! dF/dw.dw = (F(w+eps*dw)-F(w))/eps = [(weps'-weps) - (w'-w)]/eps/dtexpl
!
!                                   = (weps'-w')/eps/dtexpl - dw/dtexpl
!
! L.dw = dw/dt - beta*dF/dw.dw = (1/dt+beta/dtexpl)*dw - 
!                                beta*(weps'-w')/eps/dtexpl
!
! where w=w_k, w'=w+F_imp, beta=implpar, 
!       eps=(impldiffpar)^(1/2)*(w_k.dw)/(dw.dw)
! ----------------------------------------------------------------------------
! Centered derivative:
! --------------------
! weps1 = w + eps*dw, weps2 = w - eps*dw
!
! dF/dw.dw = (F(weps1)-F(weps2))/2/eps = 
!   [(weps1'-weps1) - (weps2'-weps2)]/2/eps/dtexpl  =
!   (weps1'-weps2')/2/eps/dtexpl - dw/dtexpl
!
! L.dw = dw/dt-beta*dF/dw.dw = (1/dt+beta/dtexpl)*dw - 
!                              beta*(weps1'-weps2')/2/eps/dtexpl
!
! where w=w_k, beta=implpar, eps=(impldiffpar)^(1/3)*(w_k.dw)/(dw.dw)

include 'vacdef.f'

double precision:: dw(nhi),qy(nhi),weps1(ixG^T,nw),weps2(ixG^T,nw)

integer:: n,ix^D,idim,iw,iiw,ntest(ndim+1)
double precision:: spar,dwnrm

double precision:: w_k(ixG^T,nw),RESimp(ixG^T,nw),wnrm(nw),dtexpl,dtcoeff
common/impld/w_k,RESimp,wnrm,dtexpl,dtcoeff
!-----------------------------------------------------------------------------

oktest=index(teststr,'matvec')>=1
if(oktest)write(*,*)'MatVec_Free'

dwnrm=sqrt(sum(dw(1:nimpl)**2)/nimpl)
if(dwnrm<smalldouble)dwnrm=one
if(implcentered) then
   spar=impldiffpar**(one/3)/dwnrm
else
   spar=sqrt(impldiffpar)/dwnrm
endif
if(oktest)write(*,*)'spar, dwnrm =',spar,dwnrm

weps1(ixM^S,1:nw)=w_k(ixM^S,1:nw)
if(implcentered)weps2(ixM^S,1:nw)=w_k(ixM^S,1:nw)

n=0
{do ix^DB=ixMmin^DB,ixMmax^DB\}
   do iiw=1,nwimpl; iw=iw_impl(iiw)
      n=n+1
      weps1(ix^D,iw)=weps1(ix^D,iw)+spar*dw(n)*wnrm(iw)
      if(implcentered)weps2(ix^D,iw)=weps2(ix^D,iw)-spar*dw(n)*wnrm(iw)
   enddo
{enddo^D&\}

if(oktest)write(*,*)'impldiffpar,spar:',impldiffpar,spar

if(index(teststr,'matvecsave')>=1)call savefile(fileout_,weps1)

call getboundary(t+dt*implpar,1,nw,1,ndim,weps1)
call advance_part(dtexpl,ixG^L,weps1)

if(implcentered)then
  call getboundary(t+dt*implpar,1,nw,1,ndim,weps2)
  call advance_part(dtexpl,ixG^L,weps2)
else
  weps2(ixM^S,1:nw)=w_k(ixM^S,1:nw)+RESimp(ixM^S,1:nw)
endif

! Calculate qy = L.dw

if(implcentered)then
    spar=implpar*half/spar/dtexpl
else
    spar=implpar/spar/dtexpl
endif

n=0
{do ix^DB=ixMmin^DB,ixMmax^DB\}
   do iiw=1,nwimpl; iw=iw_impl(iiw)
      n=n+1
      qy(n)=(dtcoeff+implpar)*dw(n)/dtexpl-&
         spar*(weps1(ix^D,iw)-weps2(ix^D,iw))/wnrm(iw)
   enddo
{enddo^D&\}

return
end
!=============================================================================
subroutine matvec_with(dw,qy,MAT)

! Calculate qy=MAT.dw for the iterative solver, non-matrix-free
! where MAT= I/dt - beta*dF/dw   (dt=dt_implicit)
!
! and F from dw/dt = F(w)
!
! explicit calculation of the Jacobian dF/dw based on
! 3-point (1D); 5-point (2D); 7-point (3D) stencils

include 'vacdef.f'

integer, parameter:: nstencil=2*ndim+1

double precision:: dw(nhi),qy(nhi),MAT(nwimpl,nwimpl,ixM^S,nstencil)

logical:: inside
integer:: n,ndw,ndwsave(nstencil),ix^D,jstencil,jshift,iw,jw,iiw
double precision:: spar, dwnrm

double precision:: w_k(ixG^T,nw),RESimp(ixG^T,nw),wnrm(nw),dtexpl,dtcoeff
common/impld/w_k,RESimp,wnrm,dtexpl,dtcoeff

integer:: jshift^D(nstencil)
logical:: periodic(ndim)
common/shifts/ jshift^D
common/period/periodic
!-----------------------------------------------------------------------------

oktest=index(teststr,'matvec')>=1
if(oktest)write(*,*)'MatVec_With'

! Calculate qy = MAT.dw = (I/dt - beta dF/dw).dw 
n=0
{do ix^DB=ixMmin^DB,ixMmax^DB\}
   do iw=1,nwimpl
      n=n+1
      qy(n)=zero

      do jstencil=1,nstencil
         if(iw==1)then
            ndw=n-iw+ jshift^D(jstencil)+
            inside=.true.
            {
            jshift=jshift^D(jstencil)
            if(jshift/=0.and.inside)then
               if((jshift>0.and.ix^D==ixMmax^D).or.&
                  (jshift<0.and.ix^D==ixMmin^D))then
                  if(periodic(^D))then 
                     ndw=ndw-nx^D*jshift
                  else
                     inside=.false.
                  endif
               endif
            endif
            \}
            if(.not.inside)ndw=-1
            ndwsave(jstencil)=ndw
         else
            ndw=ndwsave(jstencil)
         endif
         
         if(ndw>=0)then
            do jw=1,nwimpl
! write(*,*)'iw,jw,ix,jx,jstencil,n,m:',iw,jw,ix^D,jx^D,jstencil,n,ndw+jw
               qy(n)=qy(n)+MAT(iw,jw,ix^D,jstencil)*dw(ndw+jw)
            enddo
         endif
      enddo
   enddo
{enddo^D&\}

return
end

!=============================================================================
subroutine matvec_prec(dw,qy,MAT2)

! Calculate qy=MAT2.dw for the iterative solver, preconditioned
! where MAT2 is the preconditioned I/dt - beta*dF/dw matrix
!
! and F from dw/dt = F(w)
!
! explicit calculation of the Jacobian dF/dw based on
! 5-point (2D) or 7-point (3D) stencils and block-preconditioner

include 'vacdef.f'

double precision:: dw(nhi),qy(nhi),scratch(nhi)
double precision:: MAT2(nwimpl*nwimpl*{nx^D*},2*ndim+1)
!-----------------------------------------------------------------------------

! Calculate qy=L.dw

{^IFONED stop 'No preconditioning in 1D'}
{^IFTWOD
call pentasoleis(MAT2(1,1),MAT2(1,2),MAT2(1,3),MAT2(1,4),MAT2(1,5),&
    qy,dw,scratch,nx1*nx2,nwimpl,nx1)}

{^IFTHREED
call heptasoleis(MAT2(1,1),MAT2(1,2),MAT2(1,3),MAT2(1,4),MAT2(1,5),&
    MAT2(1,6),MAT2(1,7),qy,dw,scratch,nx1*nx2*nx3,nwimpl,nx1,nx1*nx2)}

return^NOONED

end

!=============================================================================
subroutine getjacobian(matsize,MAT)

include 'vacdef.f'

integer, parameter:: nstencil=2*ndim+1

! The MAT array is actually the first MATSIZE element of the WORK array
integer:: matsize
double precision:: MAT(nwimpl,nwimpl,ixM^S,nstencil)

! Index names for temporary matrices stored in the WORK array
integer:: flux_,weps_,w1_,w2_,cmax_,powellsrc_,nworkmin

! Local variables
integer:: n,ix^D,jx^D,jshift,iw,jw,iiw,jjw,jstencil,idim,iB
integer:: nrow,ipoint,ncol
logical:: compactres0

! Input common variables
double precision:: w_k(ixG^T,nw),RESimp(ixG^T,nw),wnrm(nw),dtexpl,dtcoeff
common/impld/w_k,RESimp,wnrm,dtexpl,dtcoeff

! Output common variables
integer:: jshift^D(nstencil)
logical:: periodic(ndim)
common/shifts/ jshift^D
common/period/periodic

logical:: initialized
data initialized/.false./
!-----------------------------------------------------------------------------

! Initialize the PERIODIC array and the JSHIFT array. The latter will be
!!! replaced by an indirect index array !!! 
if(.not.initialized)then
   initialized=.true.
   ! Check for periodic boundaries !!! shifted will not work ???
   periodic(1:ndim)=.false.
   do iB=1,nB
      if(typeB(1,iB)=='periodic') then
         idim=idimB(iB)
         periodic(idim)=.true.
         if((.not.impljacfast).and.(nx(idim)/nstencil)*nstencil/=nx(idim))then
            write(*,*)'Number of grid points in direction ',idimB(iB)
            write(*,*)'should be a multiple of ',nstencil
            write(*,*)'for impljacfast=F and periodic boundaries!'
            stop
         endif 
      endif 
   enddo
   
   jshift^D(1:nstencil)=0;
   jshift1(1)=0; jshift1(2)=-1; jshift1(3)=+1
   {^NOONED      jshift2(4)=-nx1; jshift2(5)=+nx1}
   {^IFTHREED    jshift3(6)=-nx1*nx2; jshift3(7)=+nx1*nx2}
   jshift^D(1:nstencil)=jshift^D(1:nstencil)*nwimpl;
endif

if(impljacfast)then
   weps_=matsize+1; w1_=weps_+nhi; w2_=w1_+nhi; flux_=w2_+nhi;
   cmax_=flux_+ndim*nwimpl*{ixGmax^D*}
   if(typeimpl1.eq.'tvdlf1')then
      powellsrc_=cmax_+ndim*{ixGmax^D*}
   else
      powellsrc_=cmax_
   endif
   if(divbfix.and.ndim>1.and.(typephys.eq.'mhd'.or.typephys.eq.'mhdiso'))then
      nworkmin=powellsrc_+nwimpl*{ixGmax^D*}
   else
      nworkmin=powellsrc_
   endif
else
   weps_=matsize+1; w1_=weps_+nhi; nworkmin=w1_+nhi;
endif
if(nwork<nworkmin)then
   write(*,*)'Not enough work space left for getjacobian(fast=',impljacfast,')'
   write(*,*)'Change parameters or increase nwork in src/vacdef.t and'
   write(*,*)'recompile VAC. Minimum value for nwork:',nworkmin
   stop
end if
if(impljacfast)then
   call getjacobian_new(sqrt(impldiffpar),periodic,MAT,&
      work(weps_),work(w1_),work(w2_),work(flux_),work(cmax_),work(powellsrc_))
else
   compactres0=compactres
   compactres=.true.
   call getjacobian_old(sqrt(impldiffpar),MAT,work(weps_),work(w1_))
   compactres=compactres0
endif

if(index(teststr,'savejacobian')>=1.and.it==itmin)then
   if(impljacfast)then
      open(99,file='jacnew',status='UNKNOWN')
      write(*,*)'Saving Jacobian into jacnew'
   else
      open(99,file='jacold',status='UNKNOWN')
      write(*,*)'Saving Jacobian into jacold'
   endif
   {do ix^DB=ixMmin^DB,ixMmax^DB\}
      do iw=1,nwimpl
         do jstencil=1,nstencil
            do jw=1,nwimpl
               write(99,'(g15.7,10i5)')MAT(iw,jw,ix^D,jstencil),&
                                       iw,jw,ix^D,jstencil
            enddo
         enddo
      enddo
   {enddo^D&\}
   close(99)
endif

!!! Testing the influence of index order. 
!!! TST(ixM^S,jstencil,iiw,jjw) is faster than MAT(iiw,jjw,ixM^S,jstencil) !!!
!   do iiw=1,nwimpl; iw=iw_impl(iiw); do jjw=1,nwimpl; jw=iw_impl(jjw)
!   do jstencil=1,nstencil
!      TST(ixM^S,jstencil,iiw,jjw)=MAT(iiw,jjw,ixM^S,jstencil)
!   end do; end do; end do
!   do iiw=1,nwimpl; iw=iw_impl(iiw); do jjw=1,nwimpl; jw=iw_impl(jjw)
!   do jstencil=1,nstencil
!      MAT(iiw,jjw,ixM^S,jstencil)=TST(ixM^S,jstencil,iiw,jjw)
!   end do; end do; end do
!endif

if(index(teststr,'savemat')>=1.and.it==itmin)then
   write(*,*)'Saving Matrix into binary file fort.98'
   do jstencil=1,nstencil
      write(98){^D&(|}((MAT(iw,jw,ix^D,jstencil),iw=1,nwimpl),jw=1,nwimpl),&
         {ix^D=ixMmin^D,ixMmax^D)}
   enddo
   close(98)
endif

if(index(teststr,'vsmjac')>=1.and.it==itmin)then
   open(20,file='coJ.mat' ,status='UNKNOWN')
   open(21,file='jcoJ.mat',status='UNKNOWN')
   open(22,file='diJ.mat' ,status='UNKNOWN')
   nrow=0
   ipoint=1
   {do ix^DB=ixMmin^DB,ixMmax^DB\}
      do iw=1,nwimpl
         nrow=nrow+1
         write(22,'(i8)') ipoint
         do jstencil=1,nstencil
            ncol=nrow-iw+(jshift^D(jstencil)+)
            {
            jshift=jshift^D(jstencil)
            if(jshift/=0)then
               jx^D=ix^D+sign(1,jshift)
               if(jx^D>ixMmax^D.or.jx^D<ixMmin^D)ncol=ncol-nx^D*jshift
            endif
            \}

            do jw=1,nwimpl
               write(20,'(g12.5)') MAT(iw,jw,ix^D,jstencil) 
               write(21,'(i8)') ncol+jw
               ipoint=ipoint+1
            enddo
         enddo
      enddo
   {enddo^D&\}
   write(22,'(i8)') ipoint
   close(20)
   close(21)
   close(22)
endif 

return
end
!=============================================================================
subroutine getjacobian_new(spar,periodic,MAT,weps,w1,w2,flux,cmaxC,powellsrc)

! calculate Jacobian matrix for TVDLF1 or EULER

include 'vacdef.f'

integer, parameter:: nstencil=2*ndim+1

logical:: periodic(ndim)
double precision:: spar,MAT(nwimpl,nwimpl,ixM^S,nstencil)
double precision, dimension(ixG^T,nw):: weps,w1,w2
double precision:: flux(ixG^S,nwimpl,ndim),cmaxC(ixG^S,ndim)
double precision:: powellsrc(ixG^S,nwimpl)

double precision, dimension(ixG^T)   :: f0,feps,v,dfdw,divb
double precision:: coeff,qareaC(ixGlo1:ixGhi1),qareadx(ixGlo1:ixGhi1)
logical:: transport,new_cmax,tvdflux,divbsrc,divbfix0
integer:: ixI^L,jxM^L,hxM^L,hxC^L,ixC^L,jxC^L,ix^D,jstencil
integer:: ix,iw,iiw,jw,jjw,idim,idir,iB

! Index names for Powell sources
integer, parameter:: qrho_=1, qm0_=1, qm^C_=qm0_+^C
integer:: qe_,qb0_,qb^C_

! Input common arrays, w_k, wnrm, dtexpl are used
double precision:: w_k(ixG^T,nw),RESimp(ixG^T,nw),wnrm(nw),dtexpl,dtcoeff
common/impld/w_k,RESimp,wnrm,dtexpl,dtcoeff
!-----------------------------------------------------------------------------

oktest=index(teststr,'getjacobian')>=1
if(oktest)write(*,*)'GetJacobian_new'

tvdflux= typeimpl1=='tvdlf1'
divbsrc= divbfix .and. ndim>1 .and. (typephys=='mhd'.or.typephys=='mhdiso')

ixI^L=ixM^L^LADD1;
ixCmax^D=ixMmax^D;

do idim=1,ndim
   ! Calculate reference flux for unperturbed w_k in each idim and iw
   call getv(w_k,ixI^L,idim,v)
   do iiw=1,nwimpl; iw=iw_impl(iiw)
      call getflux(w_k,ixI^L,iw,idim,f0,transport)
      if(transport)then
         flux(ixI^S,iiw,idim)=f0(ixI^S)+v(ixI^S)*w_k(ixI^S,iw)
      else
         flux(ixI^S,iiw,idim)=f0(ixI^S)
      endif
      if(oktest)write(*,*)'idim,iw,f0:',idim,iw,flux(ixtest^D,iiw,idim)
   enddo

   ! Calculate cmax for each direction in advance
   if(tvdflux)then
      ixCmin^D=ixMmin^D-kr(idim,^D); jxC^L=ixC^L+kr(^D,idim);

      ! Average w for interface into w1, calculate cmax and store it in v
      w1(ixC^S,1:nw)=(w_k(ixC^S,1:nw)+w_k(jxC^S,1:nw))/2

      ! In generalized coordinates we need cmax orthogonal to the surface
      {^IFGEN if(gencoord)call rotatew(ixC^L,idim,w1)}
      new_cmax=.true.
      call getcmax(new_cmax,w1,ixC^L,idim,v)
      coeff=-half*implpar
      if(gencoord)then
         ! In generalized coordinates cmax occurs as -implpar*cmax*surface/2
         cmaxC(ixC^S,idim)=coeff*v(ixC^S)*surfaceC(ixC^S,idim)
      else if(typeaxial=='slab'.or.idim/=r_)then
         ! For slab symmetric Cartesian grids only -implpar*cmax/2 occurs
         cmaxC(ixC^S,idim)=coeff*v(ixC^S)
      else
         ! For axial symmetry on Cartesian grid -implpar*cmax*area/2 occurs
         forall(ix=ixCmin1:ixCmax1)cmaxC(ix,ixC^SE,idim)=&
                           coeff*v(ix,ixC^SE)*areaC(ix)
      endif
   endif
enddo

if(divbsrc)then
   ! Determine index names based on physics
   if(typephys=='mhdiso')then
      qe_=0;     qb0_=1+^NC; qb^C_=qb0_+^C;
   else
      qe_=2+^NC; qb0_=qe_;   qb^C_=qb0_+^C;
   endif
   if(oktest)write(*,*)'Divbsrc rho,m^C,e,b^C:',qrho_,qm^C_,qe_,qb^C_
   ! Switch off divbfix for addsource
   divbfix0=divbfix
   divbfix=.false.
   ! Calculate div B for middle cell contribution to Powell's source terms
   divb(ixM^S)=zero
   do idim=1,ndim
      tmp(ixI^S)=w_k(ixI^S,qb0_+idim)
      call gradient(.false.,tmp,ixM^L,idim,tmp2)
      divb(ixM^S)=divb(ixM^S)+tmp2(ixM^S)
   enddo
   ! Calculate the coefficients that multiply div B
   do iiw=1,nwimpl; iw=iw_impl(iiw)
     {if(iw==qm^C_)powellsrc(ixI^S,iiw)=w_k(ixI^S,qb^C_) \}
     {if(iw==qb^C_)powellsrc(ixI^S,iiw)=w_k(ixI^S,qm^C_)/w_k(ixI^S,qrho_) \}
     if(iw==qe_)powellsrc(ixI^S,iiw)=(^C&w_k(ixI^S,qb^C_)*w_k(ixI^S,qm^C_)+)&
                                     /w_k(ixI^S,qrho_)
   enddo
endif

! For center cells initialize MAT here since that may not be set at all
coeff=one/dt
do iiw=1,nwimpl
   do jjw=1,nwimpl
      ! The diagonal elements are initialized to 1/dt others to 0
      if(iiw==jjw)then
         MAT(iiw,iiw,ixM^S,1)=coeff
      else
         MAT(iiw,jjw,ixM^S,1)=zero
      endif
   enddo
enddo

if(typeaxial/='slab'.or.implsource)then
   ! Put into w1 w_k+dtexpl*S(w_k)
   w1(ixM^S,1:nw)=w_k(ixM^S,1:nw)
   if(typeaxial/='slab')call addgeometry(dtexpl,ixM^L,iw_impl,w_k,w1)
   if(implsource)call addsource2(dtexpl,ixG^L,ixM^L,iw_impl,&
       t+implpar*dt,w_k,t+implpar*dt,w1)
endif

! The w to be perturbed and the index for the perturbed variable
weps(ixI^S,1:nw)=w_k(ixI^S,1:nw)
jw=0
do jjw=1,nwimpl; 
   ! Remove perturbation from previous jw if there was a previous one
   if(nwimpl>1.and.jw>0)weps(ixM^S,jw)=w_k(ixM^S,jw)
   jw=iw_impl(jjw)
   ! Perturb new jw variable
   coeff=spar*wnrm(jw)
   weps(ixM^S,jw)=w_k(ixM^S,jw) + coeff
   ! Apply boundary condition for the perturbation
   call getboundary(t+dt*implpar,1,nw,1,ndim,weps)

   do idim=1,ndim
      ! Interface indices (ixC) and shifted indices for average (jxC)
      ixCmin^D=ixMmin^D-kr(idim,^D); jxC^L=ixC^L+kr(^D,idim);

      ! Calculate velocity for perturbed transport flux in advance
      call getv(weps,ixI^L,idim,v)

      ! Calculate dfdw for each iw variable
      do iiw=1,nwimpl; iw=iw_impl(iiw)

         call getflux(weps,ixI^L,iw,idim,feps,transport)
         ! dfdw is multiplied by -implpar/2*wnrm(jw)/wnrm(iw) in all formulae
         coeff=-implpar*half/spar/wnrm(iw)
         if(transport)then
            dfdw(ixI^S)=coeff*(v(ixI^S)*weps(ixI^S,iw)+&
                               feps(ixI^S)-flux(ixI^S,iiw,idim))
         else
            dfdw(ixI^S)=coeff*(feps(ixI^S)-flux(ixI^S,iiw,idim))
         endif

         if(oktest)write(*,*)'jw,feps,dfdw:',jw,feps(ixtest^D),dfdw(ixtest^D)

         if(gencoord)then
            ! Add contribution of dfdw and cmaxC*dwdw for all interfaces
            do idir=1,ndim
               ! Interface indices (ixC) and shifted indices for average (jxC)
               ixCmin^D=ixMmin^D-kr(idir,^D); jxC^L=ixC^L+kr(^D,idir);
               hxM^L=ixM^L-kr(^D,idir); jxM^L=ixM^L+kr(^D,idir); 

               if(divbsrc.and.iw/=qrho_)then
                  call adddivbsrcdw(powellsrc,w_k,weps,&
                       -implpar*half*wnrm(jw)/wnrm(iw),&
                       -implpar*half/spar/wnrm(iw),hxM^L,jxM^L,&
                       iiw,iw,jw,qb0_+idim,idir,dfdw,tmp,tmp2)
                  tmp(ixC^S) =tmp(ixC^S)*&
                     normalC(ixC^S,idir,idim)*surfaceC(ixC^S,idir)
                  tmp2(jxC^S)=tmp2(jxC^S)*&
                     normalC(ixC^S,idir,idim)*surfaceC(ixC^S,idir)
               else
                  ! Multiply dfdw with interface normals for both sides
                  tmp(ixC^S) =dfdw(ixC^S)*&
                     normalC(ixC^S,idir,idim)*surfaceC(ixC^S,idir)
                  tmp2(jxC^S)=dfdw(jxC^S)*&
                     normalC(ixC^S,idir,idim)*surfaceC(ixC^S,idir)
               endif
               if(tvdflux.and.idim==idir)then
                  call addcmaxdwdw(cmaxC,w_k,weps,one/spar/wnrm(iw),&
                     hxM^L,jxM^L,iw,jw,idim,tmp,tmp2)

                  ! TVD flux due to cell center is added
                  if(iw==jw)MAT(iiw,jjw,ixM^S,1)=MAT(iiw,jjw,ixM^S,1)&
                        -(cmaxC(ixM^S,idim)+cmaxC(hxM^S,idim))/dvolume(ixM^S)
               endif
               ! Fluxes due to left and right neighbors
               if(idim==1)then
                  MAT(iiw,jjw,ixM^S,2*idir)  =+tmp(hxM^S)/dvolume(ixM^S)
                  MAT(iiw,jjw,ixM^S,2*idir+1)=-tmp2(jxM^S)/dvolume(ixM^S)
               else
                  MAT(iiw,jjw,ixM^S,2*idir)=MAT(iiw,jjw,ixM^S,2*idir)&
                    +tmp(hxM^S)/dvolume(ixM^S)
                  MAT(iiw,jjw,ixM^S,2*idir+1)=MAT(iiw,jjw,ixM^S,2*idir+1)&
                    -tmp2(jxM^S)/dvolume(ixM^S)
               endif
            enddo ! idir
            ! For axial symmetry the contribution of the middle cell is
            ! -implpar*wnrm(jw)/wnrm(iw)*half*dF_r/dw * (1/R)
            if(typeaxial/='slab'.and.idim==r_)then
               if(divbsrc.and.jw==qb1_.and.iw/=qrho_)then
                  coeff=-implpar*half*wnrm(jw)/wnrm(iw)
                  MAT(iiw,jjw,ixM^S,1)=MAT(iiw,jjw,ixM^S,1)&
                     -(dfdw(ixM^S)+coeff*powellsrc(ixM^S,iiw))/x(ixM^S,r_)
               else
                  MAT(iiw,jjw,ixM^S,1)=MAT(iiw,jjw,ixM^S,1)&
                     -dfdw(ixM^S)/x(ixM^S,r_)
               endif
            endif
         else if(typeaxial=='slab'.or.idim/=r_)then
            ! Put contribution of dfdw and cmax*dwdw through 
            ! the left and right interfaces in direction idim
            hxM^L=ixM^L-kr(^D,idim); jxM^L=ixM^L+kr(^D,idim);

            if(divbsrc.and.iw/=qrho_)then
                call adddivbsrcdw(powellsrc,w_k,weps,&
                    -implpar*half*wnrm(jw)/wnrm(iw),&
                    -implpar*half/spar/wnrm(iw),hxM^L,jxM^L,&
                     iiw,iw,jw,qb0_+idim,idim,dfdw,tmp,tmp2)
            else
               tmp(ixC^S)= dfdw(ixC^S)
               tmp2(jxC^S)=dfdw(jxC^S)
            endif

            if(tvdflux)then
               call addcmaxdwdw(cmaxC,w_k,weps,one/spar/wnrm(iw),&
                  hxM^L,jxM^L,iw,jw,idim,tmp,tmp2)

               ! Contribution of TVD flus due to middle cell
               if(iw==jw)MAT(iiw,jjw,ixM^S,1)=MAT(iiw,jjw,ixM^S,1)-&
                      (cmaxC(ixM^S,idim)+cmaxC(hxM^S,idim))/dx(ixM^S,idim)
            endif
            MAT(iiw,jjw,ixM^S,2*idim)=   tmp(hxM^S) /dx(ixM^S,idim)
            MAT(iiw,jjw,ixM^S,2*idim+1)=-tmp2(jxM^S)/dx(ixM^S,idim)
         else
            hxM^L=ixM^L-kr(^D,r_); jxM^L=ixM^L+kr(^D,r_);

            ! Correct coefficients for angmomfix
            if(angmomfix.and.iw==mphi_)then
               if(tvdflux)forall(ix=ixCmin1:ixCmax1)&
                  cmaxC(ix,ixC^SE,idim)=cmaxC(ix,ixC^SE,idim)*areaC(ix)
               qareaC(ixC^LIM1:)=areaC(ixC^LIM1:)
               areaC(ixC^LIM1:)=qareaC(ixC^LIM1:)**2
               qareadx(ixM^LIM1:)=areadx(ixM^LIM1:)
               areadx(ixM^LIM1:)=qareadx(ixM^LIM1:)*area(ixM^LIM1:)
            endif

            ! In axial symmetry on Cartesian grid multiply df/dw by area
            ! but first add Powell source contributions to df/dw if any
            if(divbsrc.and.iw/=qrho_)then
               call adddivbsrcdw(powellsrc,w_k,weps,&
                    -implpar*half*wnrm(jw)/wnrm(iw),&
                    -implpar*half/spar/wnrm(iw),hxM^L,jxM^L,&
                    iiw,iw,jw,qb1_,r_,dfdw,tmp,tmp2)
               forall(ix=ixCmin1:ixCmax1) tmp(ix,ixM^SE)= &
                  tmp(ix,ixM^SE)*areaC(ix)
               forall(ix=jxCmin1:jxCmax1) tmp2(ix,ixM^SE)= &
                  tmp2(ix,ixM^SE)*areaC(ix-1)
            else
               forall(ix=ixCmin1:ixCmax1) tmp(ix,ixM^SE)= &
                  dfdw(ix,ixM^SE)*areaC(ix)
               forall(ix=jxCmin1:jxCmax1) tmp2(ix,ixM^SE)= &
                  dfdw(ix,ixM^SE)*areaC(ix-1)
            endif

            if(tvdflux)call addcmaxdwdw(cmaxC,w_k,weps,one/spar/wnrm(iw),&
               hxM^L,jxM^L,iw,jw,r_,tmp,tmp2)

            ! Put contributions due to left and right neighbors and middle
            forall(ix=ixMmin1:ixMmax1) MAT(iiw,jjw,ix,ixM^SE,2)= &
               +tmp(ix-1,ixM^SE)/areadx(ix)
            forall(ix=ixMmin1:ixMmax1) MAT(iiw,jjw,ix,ixM^SE,3)= &
               -tmp2(ix+1,ixM^SE)/areadx(ix)
            if(divbsrc.and.jw==qb1_.and.iw/=qrho_)then
               ! Powell's source is not conservative, we use 
               ! the identity div_R(const)/Volume = const/R
               coeff=-implpar*half*wnrm(jw)/wnrm(iw)
               MAT(iiw,jjw,ixM^S,1)=MAT(iiw,jjw,ixM^S,1)&
                  -(dfdw(ixM^S)+coeff*powellsrc(ixM^S,iiw))/x(ixM^S,r_)

               ! Here TVD flux is added separately
               if(tvdflux.and.iw==jw)forall(ix=ixMmin1:ixMmax1) &
               MAT(iiw,jjw,ix,ixM^SE,1)=MAT(iiw,jjw,ix,ixM^SE,1)&
                  -(cmaxC(ix,ixM^SE,idim)+cmaxC(ix-1,ixM^SE,idim))/areadx(ix)
            else
               ! This is the general (ang.mom.) conservative way (with tvdflux)
               forall(ix=ixMmin1:ixMmax1) MAT(iiw,jjw,ix,ixM^SE,1)= &
               MAT(iiw,jjw,ix,ixM^SE,1)+&
                  (-tmp(ix,ixM^SE)+tmp2(ix,ixM^SE))/areadx(ix)
            endif

            ! Correct coefficients back after angmomfix
            if(angmomfix.and.iw==mphi_)then
               areaC(ixC^LIM1:)=qareaC(ixC^LIM1:)
               areadx(ixM^LIM1:)=qareadx(ixM^LIM1:)
               if(tvdflux)forall(ix=ixCmin1:ixCmax1)&
                 cmaxC(ix,ixC^SE,idim)=cmaxC(ix,ixC^SE,idim)/areaC(ix)
            endif
         endif

      enddo ! iw
   enddo ! idim
   if(oktest)then
      write(*,*)'After fluxes jw=',jw,' stencil,   MAT'
      do jstencil=1,nstencil
         write(*,*)jstencil,':',MAT(1:nwimpl,1:nwimpl,ixtest^D,jstencil)
      enddo
   endif

   !Derivatives of local and geometrical source terms 

   if(implsource.or.typeaxial/='slab')then
      if(oktest)write(*,*)'Adding dS/dw'

      ! Put into w2 w_k+dtexpl*S(w_k+eps)
      w2(ixM^S,1:nw)  =w_k(ixM^S,1:nw)
      if(typeaxial/='slab') &
         call addgeometry(dtexpl,ixM^L,iw_impl,weps,w2)
      if(implsource)call addsource2(dtexpl,ixG^L,ixM^L,iw_impl,&
         t+implpar*dt,weps,t+implpar*dt,w2)
      do iiw=1,nwimpl; iw=iw_impl(iiw)
         coeff=-implpar/spar/wnrm(iw)/dtexpl
         MAT(iiw,jjw,ixM^S,1)=  MAT(iiw,jjw,ixM^S,1)+&
            coeff*(w2(ixM^S,iw)-w1(ixM^S,iw))
      enddo
   endif
enddo

if(oktest)write(*,*)'After fluxes and sources:  MAT(...,1):', &
   MAT(1:nwimpl,1:nwimpl,ixtest^D,1)

! Contribution of middle to Powell's source terms
if(divbsrc)then
   do iiw=1,nwimpl; iw=iw_impl(iiw)
   do jjw=1,nwimpl; jw=iw_impl(jjw)
      coeff=-implpar*wnrm(jw)/wnrm(iw)

      ! S(mom)= -B*divB
      {if(iw==qm^C_.and.jw==qb^C_)MAT(iiw,jjw,ixM^S,1)=MAT(iiw,jjw,ixM^S,1)&
         -coeff*divb(ixM^S) \}
      ! S(B)= -mom/rho*divB
      {if(iw==qb^C_)then
         if(jw==qrho_)MAT(iiw,jjw,ixM^S,1)=MAT(iiw,jjw,ixM^S,1)&
            +coeff*divb(ixM^S)*w_k(ixM^S,qm^C_)/w_k(ixM^S,qrho_)**2
         if(jw==qm^C_)MAT(iiw,jjw,ixM^S,1)=MAT(iiw,jjw,ixM^S,1)&
            -coeff*divb(ixM^S)/w_k(ixM^S,qrho_) 
      endif \}

      ! S(e)= -mom.B/rho*divB
      if(iw==qe_)then
         if(jw==qrho_)MAT(iiw,jjw,ixM^S,1)=MAT(iiw,jjw,ixM^S,1)&
            +coeff*divb(ixM^S)*(^C&w_k(ixM^S,qm^C_)*w_k(ixM^S,qb^C_)+)&
             /w_k(ixM^S,qrho_)**2
         {if(jw==qm^C_)MAT(iiw,jjw,ixM^S,1)=MAT(iiw,jjw,ixM^S,1)&
            -coeff*divb(ixM^S)*w_k(ixM^S,qb^C_)/w_k(ixM^S,qrho_) \}
         {if(jw==qb^C_)MAT(iiw,jjw,ixM^S,1)=MAT(iiw,jjw,ixM^S,1)&
            -coeff*divb(ixM^S)*w_k(ixM^S,qm^C_)/w_k(ixM^S,qrho_) \}
      endif
   enddo
   enddo

   if(oktest)write(*,*)'After divb sources:  MAT(...,1):', &
      MAT(1:nwimpl,1:nwimpl,ixtest^D,1)
endif

! Add contributions from ghost cells to non-periodic boundaries
{
  if(.not.periodic(^D))then
      ixI^L=ixM^L; ixImax^D=ixMmin^D;
      MAT(1:nwimpl,1:nwimpl,ixI^S,1)=MAT(1:nwimpl,1:nwimpl,ixI^S,1)&
                                    +MAT(1:nwimpl,1:nwimpl,ixI^S,2*^D)
      MAT(1:nwimpl,1:nwimpl,ixI^S,2*^D)=zero
      ixI^L=ixM^L; ixImin^D=ixMmax^D;
      MAT(1:nwimpl,1:nwimpl,ixI^S,1)=MAT(1:nwimpl,1:nwimpl,ixI^S,1)&
                                    +MAT(1:nwimpl,1:nwimpl,ixI^S,2*^D+1)
      MAT(1:nwimpl,1:nwimpl,ixI^S,2*^D+1)=zero
   endif 
\}

if(oktest)write(*,*)'After boundary correction: MAT:',&
   MAT(1:nwimpl,1:nwimpl,ixtest^D,1)

! Restore divbfix
if(divbsrc)divbfix=divbfix0

return
end

!=============================================================================
subroutine adddivbsrcdw(powellsrc,w0,weps,coef1,coef2,hxM^L,jxM^L,&
                        iiw,iw,jw,qb_,idim,dfdw,dfdwL,dfdwR)

! Add dfdw and contributions from div B in Powell's source terms amd put it
! into dfdwL,dfdwR.
! Div B behaves the same way as div f but it is mutiplied by extra sources 
! which are associated with the central cell, thus shifted relative to the cell
! that contributs to div B via its qb_ component of B.

include 'vacdef.f'

double precision:: powellsrc(ixG^S,nwimpl)
double precision:: w0(ixG^T,nw),weps(ixG^T,nw),coef1,coef2
double precision:: dfdw(ixG^T),dfdwL(ixG^T),dfdwR(ixG^T)
integer:: hxM^L,jxM^L,iiw,iw,jw,qb_,idim
!-----------------------------------------------------------------------------

! In the mesh dB/dw is a simple Kronecker delta. The source terms are 
! multiplied by the general coef1=-implpar*half*wnrm(jw)/wnrm(iw) coefficient.

if(jw==qb_)then
   dfdwL(ixM^S)=dfdw(ixM^S)+coef1*powellsrc(jxM^S,iiw)
   dfdwR(ixM^S)=dfdw(ixM^S)+coef1*powellsrc(hxM^S,iiw) 
else
   dfdwL(ixM^S)=dfdw(ixM^S)
   dfdwR(ixM^S)=dfdw(ixM^S)
endif

! For the ghost cells calculate dB/dw explicitly from the perturbed weps and 
! the original w0. The perturbation was wnrm(jw)*spar, therefore we use
! the coefficient coef2=-implpar*half/spar/wnrm(iw).
! dfdwL and dfdwR have elements at the lower and upper edges, respectively.
select case(idim)
{case(^D)
    dfdwL(hxMmin^D^D%hxM^S)=dfdw(hxMmin^D^D%hxM^S)+&
       coef2*powellsrc(ixMmin^D^D%ixM^S,iiw)*&
       (weps(hxMmin^D^D%hxM^S,qb_)-w0(hxMmin^D^D%hxM^S,qb_))
    dfdwR(jxMmax^D^D%jxM^S)=dfdw(jxMmax^D^D%jxM^S)+&
       coef2*powellsrc(ixMmax^D^D%ixM^S,iiw)*&
       (weps(jxMmax^D^D%jxM^S,qb_)-w0(jxMmax^D^D%jxM^S,qb_))
    \}
end select

return
end

!=============================================================================
subroutine addcmaxdwdw(cmaxC,w0,weps,coeff,hxM^L,jxM^L,iw,jw,idim,dfdwL,dfdwR)

! Add contribution from cmax*dwdw to the dfdwL and dfdwR arrays

include 'vacdef.f'

double precision:: dfdwL(ixG^T),dfdwR(ixG^T),cmaxC(ixG^S,ndim)
double precision:: w0(ixG^T,nw),weps(ixG^T,nw),coeff
integer:: hxM^L,jxM^L,iw,jw,idim

!-----------------------------------------------------------------------------

if(iw==jw)then
   dfdwL(ixM^S)=dfdwL(ixM^S)+cmaxC(ixM^S,idim)
   dfdwR(ixM^S)=dfdwR(ixM^S)-cmaxC(hxM^S,idim)
endif
! For the ghost cells calculate dwdw explicitly from the perturbed weps and the
! original w0. The perturbation was wnrm(jw)*spar. Therefore we use
! coeff=one/spar/wnrm(iw), while cmaxc contains -implpar/2 etc.
! dfdwL and dfdwR have elements at the lower and upper edges, respectively.
select case(idim)
{case(^D)
    dfdwL(hxMmin^D^D%hxM^S) =dfdwL(hxMmin^D^D%hxM^S) +&
        cmaxC(hxMmin^D^D%hxM^S,idim)*coeff*&
       (weps(hxMmin^D^D%hxM^S,iw)-w0(hxMmin^D^D%hxM^S,iw))
    dfdwR(jxMmax^D^D%jxM^S)=dfdwR(jxMmax^D^D%jxM^S)-&
        cmaxC(ixMmax^D^D%ixM^S,idim)*coeff*&
       (weps(jxMmax^D^D%jxM^S,iw)-w0(jxMmax^D^D%jxM^S,iw))
    \}
end select

return
end

!=============================================================================
subroutine getjacobian_old(spar,MAT,weps,weps1)

! calculates Jacobian matrix

include 'vacdef.f'

integer, parameter:: nstencil=2*ndim+1

double precision:: spar,MAT(nwimpl,nwimpl,ixM^S,nstencil)
double precision, dimension(ixG^T,nw):: weps,weps1

integer:: n,ix^D,jx^D,iw,jw,iiw,jjw,istencil,jstencil,idiag

double precision:: coeff

double precision:: w_k(ixG^T,nw),RESimp(ixG^T,nw),wnrm(nw),dtexpl,dtcoeff
common/impld/w_k,RESimp,wnrm,dtexpl,dtcoeff

integer:: initialized,ixmask(ixG^T),icycle(nstencil,nstencil)
save initialized,ixmask,icycle
!-----------------------------------------------------------------------------

oktest=index(teststr,'getjacobian')>=1
if(oktest)write(*,*)'GetJacobian_old'

if(initialized/=12345)then
   initialized=12345
   ! Set masks, shifts and cycles for implicit solver
   call gridmasks(ixmask,icycle)     
endif

! buildup of jacobian per jw-row: dF_imp/dw(_jw)
do jjw=1,nwimpl; jw=iw_impl(jjw)
   ! calculation of residual when jw-variable changes 
   ! within stencil, using boxed changes
   ! buildup of jw-row in Jacobian dF_imp/dw(_jw)

   do jstencil=1,nstencil
      weps(ixG^S,1:nw) = w_k(ixG^S,1:nw)
      where(ixmask(ixM^S)==jstencil) &
         weps(ixM^S,jw) = w_k(ixM^S,jw) + spar*wnrm(jw)
      weps1(ixM^S,1:nw) = weps(ixM^S,1:nw)
      call getboundary(t+dt*implpar,1,nw,1,ndim,weps1)
      call advance_part(dtexpl,ixG^L,weps1)

      if(oktest)write(*,*)'weps,weps1,RESimp,jjw,jstencil:',&
         weps(ixtest^D,iwtest),weps1(ixtest^D,iwtest),&
         RESimp(ixtest^D,iwtest),jjw,jstencil

      ! store dF(_ix,iw)/dw(_jx,jw)=(F'-F)/spar in weps
      do iiw=1,nwimpl; iw=iw_impl(iiw)
         coeff= -implpar/spar/wnrm(iw)/dtexpl
         weps(ixM^S,iw)=coeff*(weps1(ixM^S,iw)-weps(ixM^S,iw)-RESimp(ixM^S,iw))
      enddo
      ! Update all elements of Jacobian for which information is available
      do istencil=1,nstencil
         idiag= icycle(istencil,jstencil)
         do iiw=1,nwimpl; iw=iw_impl(iiw)
            where(ixmask(ixM^S)==istencil)
               ! dF(_ix,iw)/dw(_jx,jw) = 
               !   (F_ix,iw(w+epsilon_jx,jw) - F_ix,iw(w))/epsilon
               MAT(iiw,jjw,ixM^S,idiag)=weps(ixM^S,iw)
            endwhere
         enddo
      enddo
   enddo
enddo

! Add 1/dt to the diagonal matrix elements
coeff=one/dt
do iiw=1,nwimpl
   MAT(iiw,iiw,ixM^S,1)=MAT(iiw,iiw,ixM^S,1)+coeff
enddo

return
end
!=============================================================================
subroutine gridmasks(ixmask,icycle)

! calculates grid masks and sets stencil cycles needed when calculating
! Jacobian matrix in implicit solver with getjacobian_old

include 'vacdef.f'

integer, parameter:: nstencil=2*ndim+1
integer:: ixmask(ixG^T),icycle(nstencil,nstencil)
integer:: ix^D,istencil,jstencil
!-----------------------------------------------------------------------------

oktest=index(teststr,'gridmask')>=1
if(oktest)write(*,*)'Gridmask'

ixmask(ixG^S)=0 

! 1D grid mask for three point stencil
!      _____________
!      | 1 | 2 | 3 |
!      |___|___|___|
! cycles for Jacobian diagonals MAIN (i) - LOWER (i-1) - UPPER (i+1)

{^IFONED

do ix1=ixMmin1,ixMmax1
   ixmask(ix1)=mod(ix1-ixMmin1,nstencil)+1
enddo

icycle(1,1)=1;icycle(1,2)=3;icycle(1,3)=2

if(oktest)then
 write(*,*)'nstencil =',nstencil
 write(*,*)'mask on grid =',ixmask(ixG^S)
endif

}

! 2D grid mask for five point stencil
!      _____________________
!      | 1 | 2 | 3 | 4 | 5 |
!      |___|___|___|___|___|
!      | 3 | 4 | 5 | 1 | 2 |
!      |___|___|___|___|___|
!      | 5 | 1 | 2 | 3 | 4 |
!      |___|___|___|___|___|
!      | 2 | 3 | 4 | 5 | 1 |
!      |___|___|___|___|___|
!      | 4 | 5 | 1 | 2 | 3 |
!      |___|___|___|___|___|
! cycles for Jacobian diagonals in order
! MAIN (i,j) - LOWER   (i-1,j) - UPPER  (i+1,j) 
!            - LOWMOST (i,j-1) - UPMOST (i,j+1)

{^IFTWOD

ix2=ixMmin2
do ix1=ixMmin1,ixMmax1
   ixmask(ix1,ix2)=mod(ix1-ixMmin1,nstencil)+1
enddo
do ix2=ixMmin2+1,ixMmax2
   ixmask(ixMmin1:ixMmax1,ix2)=mod(ixmask(ixMmin1:ixMmax1,ix2-1)+1,nstencil)+1
enddo

icycle(1,1)=1;icycle(1,2)=3;icycle(1,3)=5;icycle(1,4)=4;icycle(1,5)=2

if(oktest)then
 write(*,*)'nstencil =',nstencil
 write(*,*)'mask on grid =',ixmask(ixG^S)
endif

}

! 3D grid mask for seven point stencil
!      _____________________________
! x :  | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
!      |___|___|___|___|___|___|___|
!      _____________________________
! y :  | 1 | 3 | 5 | 7 | 2 | 4 | 6 |
!      |___|___|___|___|___|___|___|
!      _____________________________
! z :  | 1 | 4 | 7 | 3 | 6 | 2 | 5 |
!      |___|___|___|___|___|___|___|
! cycles for Jacobian diagonals in order
! MAIN (i,j,k) - LOWER   (i-1,j,k) - UPPER  (i+1,j,k)
!              - 2LOW    (i,j-1,k) - 2UP    (i,j+1,k)
!              - LOWMOST (i,j,k-1) - UPMOST (i,j,k+1)

{^IFTHREED

ix3=ixMmin3
ix2=ixMmin2
do ix1=ixMmin1,ixMmax1
   ixmask(ix1,ix2,ix3)=mod(ix1-ixMmin1,nstencil)+1
enddo
do ix2=ixMmin2+1,ixMmax2
   ixmask(ixMmin1:ixMmax1,ix2,ix3)= &
      mod(ixmask(ixMmin1:ixMmax1,ix2-1,ix3)+1,nstencil)+1
enddo
do ix3=ixMmin3+1,ixMmax3
   ixmask(ixMmin1:ixMmax1,ixMmin2:ixMmax2,ix3)= &
      mod(ixmask(ixMmin1:ixMmax1,ixMmin2:ixMmax2,ix3-1)+2,nstencil)+1
enddo

icycle(1,1)=1;icycle(1,2)=3;icycle(1,3)=5;
icycle(1,4)=7;icycle(1,5)=6;icycle(1,6)=4;icycle(1,7)=2

if(oktest)then
 write(*,*)'nstencil =',nstencil
 write(*,*)'mask on grid =',ixmask(ixG^S)
endif

}

! Fill in the rest of ICYCLE with cyclic permutations of the first row
do istencil=1,nstencil-1
   icycle(istencil+1,1)=icycle(istencil,nstencil)
   do jstencil=1,nstencil-1
      icycle(istencil+1,jstencil+1)=icycle(istencil,jstencil)
   enddo
enddo

return
end

!=======================================================================
subroutine tridiag(qdt,dw,rhs,oktest)

include 'vacdef.f'

double precision:: qdt,dw(nhi),rhs(nhi)
!-----------------------------------------------------------------------------

! The 3 diagonals of the matrix are stored in the common work array
! The size of each diagonal is nwimpl*nwimpl*nx1=nwimpl*nimpl
{^IFONED 
call tridiag1(oktest,dw,rhs,work,work(nwimpl*nimpl+1),work(2*nwimpl*nimpl+1))
}
{^IFTWOD
!!!call tridiag2
stop 'ADI in 2D is not implemented yet!!!'
}
{^IFTHREED
! ------ 3D approximate factorisation
stop 'ADI in 3D is not implemented yet!!!'
}

return
end

!=======================================================================
subroutine tridiag1(oktest,dw,rhs,qd,qe,qf)

! Use Auke van der Ploeg's block tridiagonal solver to solve MAT*dw=rhs in 1D.
! We set Auke's STEP parameter to NBLOCK, which is best for serial computers.

include 'vacdef.f'

double precision:: dw(nhi),rhs(nhi)
double precision, dimension(nwimpl,nwimpl,nx1):: qd,qe,qf

double precision:: scratch(nw,nw),scratchv(nw)
integer:: pivot(nw,ixGhi1),pivotv(nw),nworkmin

logical:: periodic(ndim)
common/period/periodic
!-----------------------------------------------------------------------------

if(oktest)write(*,*)'TRIDIAG1'

dw(1:nimpl)=rhs(1:nimpl)

if(periodic(1))then
   nworkmin=4*nwimpl*nimpl
else
   nworkmin=3*nwimpl*nimpl
endif
if(nwork<nworkmin)then
   write(*,*)'Mot enough work space left for tridiag(cyclic=',periodic(1),')'
   write(*,*)'Change parameters or increase nwork in src/vacdef.t and'
   write(*,*)'recompile VAC. Minimum value for nwork:',nworkmin
   stop
endif
if(periodic(1))then
   call precyclic(qe,qd,qf,pivot,nx1,nwimpl)
   call solcyclic(qe,qd,qf,dw,pivot,work(3*nwimpl*nimpl+1),&
      scratch,pivotv,nx1,nwimpl,scratchv)
else
   call preprocLO(qe,qd,qf,pivot,nx1,nwimpl)
   call solutionLO(qe,qd,qf,dw,nx1,nwimpl,pivot,scratchv)
endif

return
end

!=======================================================================
!subroutine tridiag2(qdt,dw,rhs,oktest,MAT,qd,qe,qf)

! In 2D use approximate factorization of MAT:
!
! MAT = (I/dt-JAC) = (I/dt-JAC_x)*(I/dt-JAC_y) + O(dt**2)
!
! where JAC = JAC_x + JAC_y + JAC_z, and JAC_x, JAC_y, and JAC_z are 
! (possibly cyclic) block tridiagonal matrices. The diagonal part of
! JAC is distributed evenly between JAC_x, JAC_y, and JAC_z.

!include 'vacdef.f'

!integer, parameter:: nstencil=2*ndim+1

!------

!double precision:: qdt,dw(nhi),rhs(nhi),MAT(nwimpl,nwimpl,ixM^S,nstencil)
!double precision, dimension(nwimpl,nwimpl,nx1):: qe,qd,qf

!integer:: ix1,ix2,nmin,iw

!double precision:: qx(nw*ixGhi1),scratch(nw,nw),scratchv(nw)
!integer:: pivot(nw,ixGhi1),pivotv(nw)

!logical:: periodic(ndim)
!common/period/periodic
!-----------------------------------------------------------------------------

!if(oktest)write(*,*)'TRIDIAG2'

!dw(1:nimpl)=rhs(1:nimpl)

!if(nx2>ixGhi1)stop 'nx2>ixGhi1 is not allowed for ADI scheme!!!'

! X-sweep row by row, solve (I/dt-JAC_x)*dw_x=rhs
!do ix2=1,nx2
   ! Take one half of the diagonal part but keep the full I/dt
!   qd(1:nwimpl,1:nwimpl,1:nx1)=MAT(1:nwimpl,1:nwimpl,ixMmin1:ixMmax1,ix2,1)/2
!   do iw=1,nwimpl
!      qd(iw,iw,1:nx1)=qd(iw,iw,1:nx1)+half/qdt
!   enddo
!   ! Take the lower and upper diagonals for the X direction
!   qe(1:nwimpl,1:nwimpl,1:nx1)=MAT(1:nwimpl,1:nwimpl,ixMmin1:ixMmax1,ix2,2)
!   qf(1:nwimpl,1:nwimpl,1:nx1)=MAT(1:nwimpl,1:nwimpl,ixMmin1:ixMmax1,ix2,3)

!   call preprocLO(qe,qd,qf,pivot,nx1,nwimpl,nx1,scratch)
!   nmin=(ix2-1)*nx1*nwimpl+1 ! The first element in dw for the ix2-th row
!   call solutionLO(qe,qd,qf,dw(nmin),nx1,nwimpl,pivot,nx1,scratchv)
!enddo

! Y-sweep column by column, solve (I/dt-JAC_y)*dw=dw_x
!do ix1=1,nx1
!   ! Take one half of the diagonal part but keep the full I/dt
!   qd(1:nwimpl,1:nwimpl,1:nx2)=MAT(1:nwimpl,1:nwimpl,ix1,ixMmin2:ixMmax2,1)/2
!   do iw=1,nwimpl
!      qd(iw,iw,1:nx2)=qd(iw,iw,1:nx2)+half/qdt
!   enddo
!   ! Take the lower and upper diagonals for the Y direction
!   qe(1:nwimpl,1:nwimpl,1:nx2)=MAT(1:nwimpl,1:nwimpl,ix1,ixMmin2:ixMmax2,4)
!   qf(1:nwimpl,1:nwimpl,1:nx2)=MAT(1:nwimpl,1:nwimpl,ix1,ixMmin2:ixMmax2,5)

!   call preprocLO(qe,qd,qf,pivot,nx2,nwimpl,nx2,scratch)
!   do ix2=0,nx2-1
!      nmin=nwimpl*(ix1-1+ix2*nx1)
!      qx(1+nwimpl*ix2:nwimpl+nwimpl*ix2)=dw(1+nmin:nwimpl+nmin)
!   enddo
!   call solutionLO(qe,qd,qf,qx,nx2,nwimpl,pivot,nx2,scratchv)
!   do ix2=0,nx2-1
!      nmin=nwimpl*(ix1-1+ix2*nx1)
!      dw(1+nmin:nwimpl+nmin)=qx(1+nwimpl*ix2:nwimpl+nwimpl*ix2)
!   enddo
!enddo

!return
!end

!=============================================================================
subroutine prepost(preproc,dw,rhs,MAT2,matsize,iter,resid,info,oktest)

! Preconditioning the matrix or postprocessing the solution depending on the 
! PREPROC logical parameter.

include 'vacdef.f'

logical:: preproc
double precision:: dw(nhi),rhs(nhi),MAT2(nwimpl*nimpl,2*ndim+1),resid
integer:: matsize,iter,info

integer:: pivot(nw,ixGhi^D*),extradiag_,extrawork_,nworkmin

!-----------------------------------------------------------------------------

{^IFONED   stop 'No preconditioning in 1D, use typeimpliter=tridiag.'}

if(preproc)then
    extradiag_=matsize+1; extrawork_=extradiag_+nwimpl*nimpl
    nworkmin=extrawork_+nimpl
    if(nwork<nworkmin)then
       write(*,*)'Not enough workspace left for pentadiagonal preconditioner.'
       write(*,*)'Change parameters or increase nwork in src/vacdef.t and'
       write(*,*)'recompile VAC. Minimum value for nwork:',nworkmin
       stop
    endif
    {^IFTWOD
    call prepenta(MAT2(1,1),MAT2(1,2),MAT2(1,3),MAT2(1,4),MAT2(1,5),&
       work(extradiag_),&
       pivot,nx1*nx2,nwimpl,nx1,implrelax,dw,rhs,work(extrawork_))}
    {^IFTHREED
    call prehepta(MAT2(1,1),MAT2(1,2),MAT2(1,3),MAT2(1,4),MAT2(1,5),&
       MAT2(1,6),MAT2(1,7),work(extradiag_),&
       pivot,nx1*nx2*nx3,nwimpl,nx1,nx1*nx2,implrelax,dw,rhs,work(extrawork_))}
else
    {^IFTWOD
    call postpenta(MAT2(1,3),MAT2(1,4),MAT2(1,5),&
       nx1*nx2,nwimpl,nx1,dw)}
    {^IFTHREED
    call posthepta(MAT2(1,3),MAT2(1,5),MAT2(1,6),MAT2(1,7),&
       nx1*nx2*nx3,nwimpl,nx1,nx1*nx2,dw)}
endif

return^NOONED
end

!=======================================================================

! end module vacimpl
!##############################################################################
