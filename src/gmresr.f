C=============================================================================
C This library contains a GMRESR and a plain GMRES algorithm, furthermore
C the MR-PC algorithm with built in time-step control. All was provided by
C M. Botchev, Mathematics Dept., Utrecht University, December 1996, July 1997.
C Some minor changes by G. Toth (Feb. 1997, July 1997).
C
C=============================================================================
C subroutines: vacgmresr,gmresr,gmres0          ! GMRESR
C              vacgmres,gmres                   ! GMRES
C              vacmrpc,ritzval,solvels,handv    ! MRPC
C
C GMRES and GMRESR use BLAS calls: DCOPY, DDOT, DNRM2, DAXPY, DSCAL, DTSRV
C MRPC uses lots of extra BLAS and LAPACK calls, directly or indirectly.
C All the needed routines can be found in src/blas.f and src/lapack.f.
C=============================================================================

C=============================================================================
C********************************* GMRESR ************************************
C=============================================================================
      subroutine vacgmresr(oktest,n,rhs,x,work,truncate,restartinner,
     &     iter,resid,matvec,stc,info)

C Solves L*x=rhs where L is defined by the matvec subroutine
C-----------------------------------------------------------------------------
C oktest is TRUE if intermediate residuals should be printed
C
C N     is the number of unknowns
C
C RHS   is the right hand side
C
C X     is the initial guess on input, and the solution on output
C
C WORK  is a work array for extra vectors needed by GMRESR
C
C TRUNCATE is the number of vectors in the truncated series
C
C RESTARTINNER is the dimension of the inner Krylov subspace
C       RESTARTINNER = 0 results in the GCR method
C       which seems to be promising for the matrix-free approach
C
C ITER  ITER*TRUNCATE is the maximum number of outer iterations on input and 
C       ITER is the actual number of iterations on output inculding inner
C       iterations as well.
C
C RESID is the required accuracy on input, and the achieved accuracy on output
C
C STC   determines stopping criterion (||.|| denotes the 2-norm):
C       stc='rel'    -- relative stopping crit.: ||res|| < RESID*||res0||
C       stc='abs'    -- absolute stopping crit.: ||res|| < RESID
C
C INFO  is 0 for success and a different number for each kinds of errors:
C
C INFO  is 1 if the required accuracy has not been achieved for given ITER 
C
C NOTE ON MEMORY: the number of vectors required is
C
C     NVECTOR = 2*TRUNCATE + RESTARTINNER + 2  
C-----------------------------------------------------------------------------
      logical oktest
      integer n,truncate,restartinner,iter,info
      double precision x(n),rhs(n),work(*),resid
      character*3 stc

C internal variables
      double precision eps

      external matvec
C-----------------------------------------------------------------------------
 
C Set internal parameters for GMRESR 

C     Save the input value of resid into eps
      eps = resid

      call gmresr(oktest,n,truncate,restartinner,rhs,x,work,eps,stc,
     &            iter,resid,matvec,info)

      return
      end

C=============================================================================
      subroutine gmresr(oktest,n,j,mgmres,b,x,work,eps,stc,
     &                  maxits,resid,matvec,iflag)
C*********************************************************
C GMRESR algorithm to solve linear system Ax = b
C
C  m.botchev, utrecht university, december 1996
C
C Details to the algorithm may be found in
C  H.A. van der Vorst, C. Vuik, "GMRESR: a Family of Nested GMRES
C  Methods", Num. Lin. Alg. Appl., vol. 1(4), 369--386 (1994)
C
C parameter list:
C oktest is TRUE if intermediate residuals should be printed
C n      == INTEGER size of problem
C j      == INTEGER truncation parameter (work with j last vectors)
C mgmres == INTEGER dimension of the envoked GMRES
C           if mgmres.eq.0 then we get GCR algorithm, the simplest
C           version of GMRESR 
C b      == DOUBLE PRECISION righthand side vector
C x      == DOUBLE PRECISION initial guess on input,
C           (approximate) solution on output
C work   == DOUBLE PRECISION work space 1 of size n x (2*j+mgmres+2)
C eps    == DOUBLE PRECISION tolerance of stopping criterion. 
C           process is stopped as soon as 
C           (1) residual norm has been dumped by factor eps, 
C           i.e.  ||res|| / ||res0|| <= eps   OR
C           (2) maximum number of iterations maxit has been performed
C stc    == CHARACTER*3
C           Determine stopping criterion (||.|| denotes the 2-norm):
C           stc='rel'    -- relative stopping crit.: ||res|| < eps*||res0||
C           stc='abs'    -- absolute stopping crit.: ||res|| < eps
C maxits == INTEGER max. no. outer_iterative_steps/truncation_length on input
C           on output it is the actual number of total iterative steps   
C resid  == DOUBLE PRECISION residual measure (depends on stopping criterion)
C           achieved on output 
C iflag  == INTEGER on output 0 - solution found within tolerance
C                             1 - no convergence within maxits
C ----------------------------------------------------------
C subroutines used
C   matvec   == matrix-vector product y <- A*x
C  blas subroutines:
C   dscal
C   daxpy
C  blas functions:
C   ddot
C   dnrm2 
C**********************************************************

      external matvec

C     list of variables: arrays in alphabetical order,
C     then other variables in alphabetical order

      logical oktest
      character*3 stc

      integer i,iflag,its,nits,itsinn,j,k,maxits,mgmres,n

      double precision b(n),x(n),work(n,0 : 2*j+mgmres+2 -1),
     &     alpha,alphai,cknorm,ckres,ddot,dnrm2,eps,epsinn,
     &     res0,resnor,resid

C distribute work space work(n,*) among some virtual variables;
C namely, we think of columns of work as being occupied by 
C c(n,0:j-1), u(n,0:j-1), resid(n), workgmr(n,mgmres+1)
C therefore we define "shifts"

      integer c, u, workres, workgmre
C ----------------------------------------------------------------------------

      if((stc.NE.'rel').and.(stc.NE.'abs'))then
         PRINT *,'Error in VACGMRESR:'
         PRINT *,'PARAMETER STC=',stc,' SHOULD BE rel OR abs.'
         STOP
      endif

C     c occupies j columns 0...j-1:
      c = 0 
C     u occupies j columns j...2*j-1:
      u = j
C     resid occupies 1 column No. 2*j:
      workres = 2*j    
C     workgmre occupies mgmres+1 columns 2*j+1...2*j+mgmres+1:
      workgmre = 2*j+1 
C     so that we can access, say, to the (k)th column of the virtual
C     array c(n,0:j-1) as work(1,c+k),
C     virtual residual vector resid(n) is work(1,workres) and so on ...

C ***Furthermore, we build sequences c_k and u_k, k = 0,...,m-1
C but we store only the last j vectors in the following way:
C Assume j=3, then
C --------------------------------------------------------------
C  k    |  number of column of work | columns of work which are vectors
C       |  in which we store c_k    |  we actually store and use
C  0    |           0               |   c_0             u_0            ...
C  1    |           1               |   c_0  c_1        u_0  u_1       ...
C  2    |           2               |   c_0  c_1  c_2   u_0  u_1  u_2  ...
C  3    |           0               |   c_3  c_1  c_2   u_3  u_1  u_2  ...
C  4    |           1               |   c_3  c_4  c_2   u_3  u_4  u_2  ... 
C  5    |           2               |   c_3  c_4  c_5   u_3  u_4  u_5  ...
C  6    |           0               |   c_6  c_4  c_5   u_6  u_4  u_5  ...
C ...   |           ...             |      ...               ...
C This mapping is performed by function mod(k,j)

C     Reset iteration counter
      nits= 0
      its = 0
      
C     Calculate (initial) residual norm
      call matvec(x,work(1,workres),n)
      alpha = -1
      call daxpy(n,alpha,b,1,work(1,workres),1)
      call dscal(n,alpha,work(1,workres),1)

C     Calculate residual norm and quit if it is zero
      res0 = dnrm2(n,work(1,workres),1)
      resnor = res0
      resid = 0

      if ( res0 .eq. 0.0d0 ) then
         iflag = 0
         maxits = 0
         return  
      end if

      if (stc.eq.'abs') then
         resid=resnor
      else
         resid=resnor/res0
      endif

      if ( resid .le. eps ) then
         iflag = 0
         maxits = 0
         return
      end if 
 
C     Main iterative loop ============================
      k = -1
      do while (.true.)

         if(oktest)write(*,199)its,resid
 199     format('   its =', i4, ' resid =', d20.6)

C        Loop to increment dimension of projection subspace
         k=k+1

C        Number of step (not restart) to be done
         its = its + 1
C        write(*,'(A,i3)') '+++++++++++++++++ its ',its 

C        - - - - - - - - - - - - - - - - - - - - - - - - -
C        This part should deliver 
C        u(1,k) <-- invA * resid
C        where u(1,k) is the k-th column of array u(1:n,0:m) and
C        invA is some reasonable approximation to the inverse of A
C
C        If mgmres=0 then no inner iterations are performed 
C        to get invA, so that invA is just the identity matrix. 
C        In this case algorithm GMRESR is nothing but GCR
C
C        Otherwise for inner iterations we perform ONE restart of GMRES
C        ATTN: for this implementation it is crucial to perform only
C        ONE restart of GMRES
         if (mgmres.eq.0) then
C           u(1,k) := resid  
            call dcopy(n,work(1,workres),1,work(1,u+mod(k,j)),1)
            call matvec(work(1,u+mod(k,j)),work(1,c+mod(k,j)),n)
            nits=nits+1
         else
C           Solve linear system A*u(1,k)=resid by GMRES
C           The stopping criterion for inner iterations is 
C           always absolute but it is automatically adjusted
C           not to be stricter than the stopping criterion for the 
C           outer iterations.  For example, if stop.criterion for
C           the outer iterations is relative than absolute stop.
C           criterion for the inner iterations is (eps*res0)
C           Accuracy for inner iteration:

            if(stc.eq.'abs')then
               epsinn = eps
            else
               epsinn = eps*res0
            endif

C           After envoking gmres0 epsinn and itsinn contain actual achieved
C           accuracy and number of performed iterations respectively

            itsinn=mgmres

            call gmres0(oktest,n,mgmres,
     &           work(1,workres),work(1,u+mod(k,j)),
     &           work(1,c+mod(k,j)),work(1,workgmre),
     &           epsinn,itsinn,matvec)

            nits=nits+itsinn
         end if           
C - - - - - - - - - - - - - - - - - - - - - - - - 
      
C        Inner loop to orthogonalize 
C        c(1,k) with respect to c(1,k-j),...,c(1,k-1)
C        and to update correspondingly 
C        u(1,k) with respect to u(1,k-j),...,u(1,k-1)
C        parameter j is used only here
         do i = max0(0,k-j),k-1
            alphai = ddot(n,work(1,c+mod(i,j)),1,work(1,c+mod(k,j)),1)
            call daxpy(n,-alphai,work(1,c+mod(i,j)),1,
     &           work(1,c+mod(k,j)),1)
            call daxpy(n,-alphai,work(1,u+mod(i,j)),1,
     &           work(1,u+mod(k,j)),1)
         end do

C        Normalize c(1,k) and "normalize" u(1,k)
         cknorm = dnrm2(n,work(1,c+mod(k,j)),1)
         cknorm = 1 / cknorm
         call dscal(n,cknorm,work(1,c+mod(k,j)),1)
         call dscal(n,cknorm,work(1,u+mod(k,j)),1)

C        Update current solution and residual
         ckres = ddot(n,work(1,c+mod(k,j)),1,work(1,workres),1)
         call daxpy(n, ckres,work(1,u+mod(k,j)),1,x,          1)
         call daxpy(n,-ckres,work(1,c+mod(k,j)),1,work(1,workres),1)

C        call show(n,10,x,'GMRESR       ')  

C        Calculate residual norm, check convergence
         resnor = dnrm2(n,work(1,workres),1)

         if (stc.eq.'abs') then
            resid=resnor
         else
            resid=resnor/res0
         endif

         if ( resid .le. eps ) then
            iflag = 0
            maxits = nits
            return
         end if
         if (its .ge. maxits*j) then
            iflag = 1
            maxits = nits
            return
         end if

C        print 11, '            ||res|| = ',resnor 
C 11     format(A,d)
C 13     format(i4,A,d)
      
      end do
C End of inifinite iterative loop =================
C End of GMRESR subroutine      
      end 

C=============================================================================
       subroutine gmres0(oktest,n,im,rhs,uu,cc,work0,eps,maxits,matvec)

C This is the modified GMRES routine gmres0 adapted for GMRESR by 
C Mike Botchev, Utrecht University, Dec. 1996
C For detail on how to make GMRES (for GMRESR) cheaper see 
C the above-mentioned paper on GMRESR 
c*************************************************************
C This code was initially written by Youcef Saad
C then revised by Henk A. van der Vorst  
C and Mike Botchev (oct. 1996)
C ************************************************************ 
c gmres algorithm . simple version .  (may 23, 1985)
c parameter list:
c oktest == TRUE for printing intermediate results
c n      == size of problem
c im     == size of krylov subspace:  should not exceed 50 in this
c          version (can be reset in code. looking at comment below)
c rhs    == right hand side
c uu     == initial guess for vector u (see above-mentioned paper on GMRESR)
c           on input, approximate solution on output
c cc     == initial guess for vector c (see above-mentioned paper on GMRESR)
c           on input, approximate solution on output
c work0  == work space of size n x (im+1)
c eps    == tolerance for stopping criterion. process is stopped
c           as soon as ( ||.|| is the euclidean norm):
c           || current residual || <= eps  
c maxits == maximum number of iterations allowed
c           on OUTPUT: actual number of iterations performed
c ----------------------------------------------------------------
c subroutines 
c matvec      == matrix vector multiplication y <- A*x
c
c BLAS:
c dcopy       == y <-- x routine
c ddot        == dot product function
c dnrm2       == euclidean norm function
c daxpy       == y <-- y+ax routine
c dscal       == x <-- ax routine
c dtsrv       == to solve linear system with a triangular matrix
c*************************************************************
c-------------------------------------------------------------
c arnoldi size should not exceed 10 in this version..
c to reset modify maxdim. BUT:             ----------------
c maxdim was set to 10 because of two reasons:
c (1) it is assumed in this implementation that it is cheaper to
c make maxdim vector updates than to make 1 matrix-vector
c multiplication;
c (2) for large maxdim we may lose the orthogonality property
c on which this cheap implementation is based.
c Please keep it in mind changing maxdim
c-------------------------------------------------------------
      integer maxdim,maxd1,md1max
      parameter (maxdim=10, maxd1=maxdim+1, md1max=maxdim*maxd1)
      external matvec

      logical oktest
      integer jjj,jj1
      integer i,i1,im,its,j,k,k1,maxits,n
      double precision cc(n),coeff,coef1,dabs,ddot,dnrm2,dsqrt,eps,epsmac,
     &                 gam,rhs(n),ro,uu(n),work0(n,im+1),t     

      double precision hh(maxd1,maxdim),hh1(maxd1,maxdim),c(maxdim),
     &                 s(maxdim),rs(maxd1),rs1(maxd1)

      data (( hh(jj1,jjj), jj1=1,maxd1), jjj=1,maxdim) / md1max*0.0 / ,
     &      epsmac / 1.d-16 / 
C-----------------------------------------------------------------------------

      if (im .gt. maxdim) then
         im = maxdim
         write (*,'(A,i2)') 'GMRES0: dimension has been reduced to ',im
         write (*,'(A)') ' => reset MAXDIM if you want it to be more'
         write (*,'(A)') ' BUT read comments near MAXDIM before'
      end if

      its = 0

C     ----------------------------
C     Outer loop starts here.. 
C     BUT here (for GMRESR) only 1 outer loop is allowed
C     Compute initial residual vector 
C     ----------------------------
 10   continue
C        do not calculate initial residual first restart because 
C        initial guess is always zero. 
C        make initial guess zero:
         coeff = 0.0
         call dscal(n,coeff,uu,1)
C        make initial residual right-hand side:
         call dcopy(n,rhs,1,work0,1)

	 ro = dnrm2 (n, work0, 1)
	 if ((ro .eq. 0.0d0).or.(ro .le. eps)) then
            call matvec(uu, cc, n)
            eps = ro
            maxits = its 
            return
         end if

         coeff = 1 / ro
         call dscal(n,coeff,work0,1)

	 if (oktest) write(*, 199) its, ro

c        initialize 1-st term  of rhs of hessenberg system..
	 rs(1) = ro
	 i = 0

 4       continue
            i=i+1
            its = its + 1
            i1 = i + 1
            call  matvec(work0(1,i), work0(1,i1), n)
c           -----------------------------------------
c           modified gram - schmidt...
c           -----------------------------------------
            do j=1, i
               t = ddot(n, work0(1,j),1,work0(1,i1),1)
               hh(j,i) = t
               call daxpy(n, -t, work0(1,j), 1, work0(1,i1), 1)
            end do
            t = dnrm2(n, work0(1,i1), 1)
            hh(i1,i) = t
            if (t .ne. 0.0d0)then
               t = 1 / t
               call dscal(n, t, work0(1,i1), 1)
C              save new column of hh in hh1 to reproduce vector cc later on
               call dcopy(maxd1,hh(1,i),1,hh1(1,i),1)
            endif
c           done with modified gram schmidt and arnoldi step..

c           now  update factorization of hh
            if (i .ne. 1) then
c              perform previous transformations  on i-th column of h
               do k=2,i
                  k1 = k-1
                  t = hh(k1,i)
                  hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
                  hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
               end do
            endif
            gam = dsqrt(hh(i,i)**2 + hh(i1,i)**2)
            if (gam .eq. 0.0d0) gam = epsmac

c           determine next plane rotation
            c(i) = hh(i,i)/gam
            s(i) = hh(i1,i)/gam
            rs(i1) = -s(i)*rs(i)
            rs(i) =  c(i)*rs(i)

c           determine residual norm and test for convergence-
            hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
            ro = dabs(rs(i1))
            if (oktest) write(*, 199) its, ro
	 if ((i .lt. im) .and. (ro .gt. eps))  goto 4
c
c        now compute solution. first solve upper triangular system.
c
C        rs := hh(1:i,1:i) ^-1 * rs

         call dtrsv('U','N','N',i,hh,maxd1,rs,1)   
c        done with back substitution..

c        now form linear combination to get vector uu
	 do j=1, i
            t = rs(j)
	    call daxpy(n, t, work0(1,j), 1, uu,1)
         end do
C        DO NOT restart outer loop EVEN when necessary (that is crucial
C        for this implementation of GMRESR):  NEVER goto 10 !  
C     if (ro .gt. eps .and. its .lt. maxits) goto 10

C     Finally, reproduce vector cc as cc = A*uu = work0*hh1*rs:
C     rs := hh1(1:i1,1:i) * rs
      coeff = 1
      coef1 = 0
      call dgemv('N',i1,i,coeff,hh1,maxd1,rs,1,coef1,rs1,1)

C     now multiply Krylov basis vectors work0 by rs:
C     cc := work0*rs
      call dscal(n,coef1,cc,1)
      do j=1, i1
         t = rs1(j)
	 call daxpy(n, t, work0(1,j), 1, cc,1)
      end do        

 199  format('itsinn =', i4, ' res. norm =', d20.6)

      maxits=its
      eps=ro 
      return
c------------------------------- end of gmres0 ----------------------
      end

C=============================================================================
C********************************* GMRES *************************************
C=============================================================================
      subroutine vacgmres(oktest,n,rhs,x,nonzero,work,restart,iter,
     &                    resid,matvec,stc,info)

C Solves L*x=rhs where L is defined by the matvec subroutine
C-----------------------------------------------------------------------------
C oktest == TRUE for printing intermediate results
C
C N     is the number of unknowns
C
C RHS   is the right hand side
C
C X     is the initial guess on input, and the solution on output
C
C WORK  is a work array for extra vectors needed by GMRES
C
C RESTART is the dimension of the Krylov subspace
C
C ITER  ITER*RESTART is the maximum number of iterations on input and 
C       ITER is the actual number of iterations on output
C
C RESID is the required accuracy on input, and the achieved accuracy on output
C
C
C STC   is the stopping criterion (||.|| denotes the 2-norm):
C           stc='rel'    -- relative stopping crit.: ||res|| < RESID*||res0||
C           stc='abs'    -- absolute stopping crit.: ||res|| < RESID
C
C INFO  is 0 for success and a different number for each kinds of errors:
C
C INFO  is 1 if required accuracy has not been achieved
C
C NOTE ON MEMORY: the number of vectors required is
C
C     NVECTOR = RESTART + 1  
C-----------------------------------------------------------------------------

      logical nonzero,oktest
      character*3 stc
      integer n,restart,iter,info
      double precision x(n),rhs(n),work(*),resid

C internal variables
      double precision eps
      external matvec
C-----------------------------------------------------------------------------

C     Save input value of resid into eps
      eps = resid
      info = 0

      call gmres(oktest,n,restart,rhs,x,nonzero,work,resid,stc,iter,
     &           matvec)

      if (resid .gt. eps) info = 1

      return
      end

C==================================================================
      subroutine gmres(oktest,n,im,rhs,sol,nonzerox,work,eps,stc,
     &                 maxits,matvec)
C ************************************************************
C This code was initially written by Youcef Saad
C then revised by Henk A. van der Vorst  
C and Mike Botchev (Oct. 1996)
C Some minor changes by G. Toth (Feb. 1997)
C ************************************************************ 
c gmres algorithm . simple version .  (may 23, 1985)
c parameter list:
c n      == size of problem
c im     == size of krylov subspace:  should not exceed 50 in this
c          version (can be reset in code. looking at comment below)
c rhs    == right hand side
c sol    == initial guess on input, approximate solution on output
c nonzerox == tells gmres if initial guess in sol is zero or not. 
c           if nonzerox is .FALSE., one MATVEC call is saved.
c work   == work space of size n x (im+1)
c eps    == on input:  tolerance for stopping criterion
c           on output: actual achieved (relative) norm of the residual
C stc    == stopping criterion (||.|| denotes the 2-norm):
C           stc='rel'    -- relative stopping crit.: ||res|| < eps*||res0||
C           stc='abs'    -- absolute stopping crit.: ||res|| < eps
c maxits == maximum number of iterations allowed
c oktest == TRUE for printing intermediate results
c ----------------------------------------------------------------
c subroutines used =
c matvec      == matrix vector multiplication y <- A*x
c
c BLAS:
c ddot        == dot product function
c dnrm2       == euclidean norm function
c daxpy       == y <-- y+ax routine
c dscal       == x <-- ax routine
c dtsrv       == to solve linear system with a triangular matrix
c*************************************************************
c-------------------------------------------------------------
c arnoldi size should not exceed 50 in this version..
c to reset modify maxdim             ----------------
c-------------------------------------------------------------
      integer maxdim, maxd1
      parameter (maxdim=50, maxd1=maxdim+1)
      external matvec

      logical nonzerox,oktest
      character*3 stc
      integer i,i1,im,its,j,k,k1,maxits,n,n1
      double precision coeff,dabs,ddot,dnrm2,dsqrt,eps,eps1,epsmac,
     &                 gam,rhs(n),ro,sol(n),work(n,im+1),t

      double precision hh(maxd1,maxdim),c(maxdim),s(maxdim),rs(maxd1)

      data epsmac/1.d-16/

      integer iii
C ===========================================================================

      if((stc.NE.'rel').and.(stc.NE.'abs'))then
         PRINT *,'Error in VACGMRES:'
         PRINT *,'PARAMETER STC=',stc,' SHOULD BE rel OR abs.'
         STOP
      endif

      if (im .gt. maxdim) then
         im = maxdim
         write (*,'(A,i2)') 'GMRES: dimension has been reduced to ',im
         write (*,'(A)') ' => reset MAXDIM if you want it to be more'
      end if

      n1 = n + 1
      its = 0
c     -------------------------------------------------------------
c     outer loop starts here..
c     -------------- compute initial residual vector --------------
 10   continue
         if (nonzerox.or.its.gt.0) then
	    call matvec(sol, work, n)
            coeff = -1
            call daxpy(n,coeff,rhs,1,work,1)
            call dscal(n,coeff,work,1)
         else
c           do not perform MATVEC since initial residual = rhs       
            call dcopy(n,rhs,1,work,1)
         end if
c-------------------------------------------------------------
	 ro = dnrm2 (n, work, 1)
	 if (ro .eq. 0.0d0) then
            eps = ro
            maxits = its 
            return
         end if

	 if (its .eq. 0) then      
C           set eps1 for stopping criterion
            if (stc.eq.'abs') then
               eps1=eps
               if (ro .le. eps1) then
                  eps = ro
                  maxits = its
                  return
               end if
            else
               eps1=eps*ro
            end if
         end if

         coeff = 1 / ro
         call dscal(n,coeff,work,1)

	 if (oktest) write(*, 199) its, ro

c        initialize 1-st term  of rhs of hessenberg system..
	 rs(1) = ro
	 i = 0
 4       continue
            i=i+1
            its = its + 1
            i1 = i + 1
            call matvec(work(1,i), work(1,i1), n)

c           -----------------------------------------
c           modified gram - schmidt...
c           -----------------------------------------
            do j=1, i
               t = ddot(n, work(1,j),1,work(1,i1),1)
               hh(j,i) = t
               call daxpy(n, -t, work(1,j), 1, work(1,i1), 1)
            end do
            t = dnrm2(n, work(1,i1), 1)
            hh(i1,i) = t
            if (t .ne. 0.0d0) then
               t = 1 / t
               call dscal(n, t, work(1,i1), 1)
            endif
c           done with modified gram schmidt and arnoldi step..

c           now  update factorization of hh
            if (i .ne. 1) then

c              perform previous transformations  on i-th column of h
               do k=2,i
                  k1 = k-1
                  t = hh(k1,i)
                  hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
                  hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
               end do
            endif
            gam = dsqrt(hh(i,i)**2 + hh(i1,i)**2)
            if (gam .eq. 0.0d0) gam = epsmac

c           determine next plane rotation
            c(i) = hh(i,i)/gam
            s(i) = hh(i1,i)/gam
            rs(i1) = -s(i)*rs(i)
            rs(i) =  c(i)*rs(i)

c           determine residual norm and test for convergence
            hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
            ro = dabs(rs(i1))
            if (oktest) write(*, 199) its, ro
 199        format('   its =', i4, ' res. norm =', d20.6)
	 if (i .lt. im .and. (ro .gt. eps1))  goto 4

c        now compute solution. first solve upper triangular system.  
C        rs := hh(1:i,1:i) ^-1 * rs

         call dtrsv('U','N','N',i,hh,maxd1,rs,1)
c        done with back substitution

c        now form linear combination to get solution
	 do j=1, i
            t = rs(j)
	    call daxpy(n, t, work(1,j), 1, sol,1)
         end do

c        restart outer loop  when necessary

      if (ro .gt. eps1 .and. its .lt. maxits*im) goto 10

      maxits=its
      eps=eps/eps1*ro

      return
c------------------------------- end of gmres ----------------------
      end 

C=============================================================================
C**************************** MR-PC ******************************************
C=============================================================================

      subroutine vacmrpc(oktest,n,rhs,x,nonzero,workV,iter0,
     &           resnorm,matvec,stc,tau,implmrpc,typepred,info)
C ******************************************************  
C M.Botchev, Math. Dept., Utrecht University, July 1997
C ******************************************************
C Approximately solves L*x=rhs where L is defined by the matvec subroutine
C-----------------------------------------------------------------------------
C OKTEST is TRUE for printing intermediate results
C
C N     is the number of unknowns
C
C RHS   is the right hand side
C
C X     is the initial guess (by predictor) on input, the solution on output
C
C NONZERO is TRUE if the RHS array contains any nonzero elements
C
C WORKV  is a work array for extra vectors needed by VACMRPC
C
C ITER0  is the dimension of the Krylov subspace, i.e.,
C        ITER0 (nonrestarted) steps of GMRES is performed 
C
C RESNORM is the required accuracy on input, and the achieved 
C       accuracy on output
C
C STC   analogue of the stopping criterion. Output parameter INFO
C       (see below) is set to 0 if (||*|| denotes the 2-norm):
C       ||res|| < RESNORM*||res0|| when stc='rel' -- relative stop.criterion
C       ||res|| < RESNORM          when stc='abs' -- absolute stop.criterion
C
C TAU   on input, is the current time step
C       on output, time step chosen by the MR-PC on output, if IMPLMRPC>=0
C                  -1.0D0 if IMPLMRPC=-1       
C
C IMPLMRPC is integer parameter for step-size control:
C       implmrpc<0 -- do not change time step tau (may be useful e.g.
C                     if vacmrpc is used like GMRES linear solver)
C       implmrpc=0 -- gentle   control, without extra-corrector
C       implmrpc=1 -- gentle   control, with extra-corrector
C       implmrpc=2 -- mild     control, with extra-corr.
C       implmrpc=3 -- standard control, without extra-corrector
C       implmrpc=4 -- standard control, with extra-corrector
C       implmrpc=5 -- strict   control, with extra-corrector
C
C The standard setting is 4. Switching off the extra-corrector 
C allows to decrease computational work, but it is prone to instability 
C (when the step size is to be decreased and actual change of the
C step size is postponed till the next time step)  
C
C TYPEPRED is a character*10 parameter for the explicit predictor:
C    nul       --  zeroth order predictor step, time step will be 
C                  changed immediately, rhs=0, solution rescaled.
C    explicit  --  first order predictor step, time step will be 
C                  changed immediately, rhs and solution rescaled.
C    explicit2 --  second order predictor step, time step will be
C                  only be changed in next time step
C
C INFO  is 0 for success and a different number for each kinds of errors:
C
C INFO  is 1 if required accuracy has not been achieved
C INFO  is 2 if there is nothing to do because initial guess
C            is accurate enough    
C
C NOTE ON MEMORY: the number of vectors required is
C
C     NVECTOR = ITER + 1  
C
C Uses BLAS and LAPACK calls
C-----------------------------------------------------------------------------
      integer itMax
      parameter (itMax=130)
      integer gentle,gentlecorr,mildcorr
      integer standard,standardcorr,strictcorr
      parameter (gentle=0,gentlecorr=1,mildcorr=2)
      parameter (standard=3,standardcorr=4,strictcorr=5)
 
      logical nonzero,oktest,TauBad
      character*3 stc
      character*10 typepred
      integer i,iii,info,iter,iter0,iterAct,n,implmrpc
C Debug:
      integer j 
      double precision al,be,bLS(itMax+1),coeff,dmax1,dnrm2,dsqrt,eps,
     & resnorm,resnorm0,tau,tau0,tau1,tauNew,HRVabs,HRVmi0,HRVmi1,
     & HRVmin,Hhh(itMax+1,itMax),HRVa,HRVb,HRVr(itMax),HRVi(itMax),
     & Hwork(itMax+2,itMax+2),rhs(n),workV(n,iter0+1),x(n),xLS(itMax+1)

      external matvec
C-----------------------------------------------------------------------------
c      oktest = .true.
C-----------------------------------------------------------------------------

C Initialize Krylov dimension
      iter = iter0
C     (NB: iter0 is introduced because otherwise, when iter is changed,
C      the corresponding parameter in VAC implrestart also) 
  
C Save required residual value resnorm into eps
      eps = resnorm
      info = 0

      if (iter.gt.itMax) then
         iter = itMax
         write (*,'(A,i2)') 'VACMRPC: iter has been reduced to ',itMax
         write (*,'(A)') ' => reset itMax if you want it to be more'
      end if

C Calculate initial residual and store it in workV(1,1)
      if (nonzero) then
C (nonzero) means not 'nul', so that the residual has to be computed
	    call matvec(x, workV(1,1), n)
            coeff = -1
            call daxpy(n,coeff,rhs,1,workV(1,1),1)
            call dscal(n,coeff,workV(1,1),1)
      else
C  otherwise residual = rhs 
            call dcopy(n,rhs,1,workV(1,1),1)
      endif 

      resnorm0 =  dnrm2( n, workV(1,1),1 )
      if ((stc.eq.'abs').and.(resnorm0.le.eps)) then
         resnorm = resnorm0
         if(oktest) write(*,'(A)') 'VACMRPC: nothing to do'
         info = 2
         if(implmrpc.lt.0) tau = -1.D0 
         return
      endif

C Form Krylov subspace
      call   handv(n,itMax+1,iter,iterAct,workV(1,1),Hhh,workV, matvec)
C     sub... handv(n,LDH,dimGMR,dimGMRact,res0,    Hh, workV,matvec)

C Prescribed Krylov dimension iter may be decreased to iterAct 
C if the next vector of Krylov subspace is zero
      if (oktest.and.(iter.gt.iterAct)) then
         write(*,'(A,i2)')
     &         'VACMRPC: Krylov dimension is decreased to ',iterAct
      endif

      if (iterAct.eq.0) then
C     -- initial guess appears to be an exact solution
         resnorm = resnorm0
         if(oktest) write(*,'(A)') 'VACMRPC: nothing to do'
         info = 2 
         if(implmrpc.lt.0) tau = -1.D0 
         return
      endif

      if (implmrpc.ge.0) then
         if (implmrpc.eq.mildcorr) then
C           "between gentle and standard" step size control
            HRVa = 3.7
            HRVb = 4.5
         else if ((implmrpc.eq.gentle).or.(implmrpc.eq.gentlecorr)) then
C           gentle step size control
            if(oktest) write(*,'(A)') 'gentle step size control'
            HRVa = 9.5  !7.0 
            HRVb = 12.5 !12.0
         else if(implmrpc.eq.strictcorr) then
C           strict step size control
            if(oktest) write(*,'(A)') 'strict step size control' 
            HRVa = 1.8 
            HRVb = 2.5 
         else
C           standard step size control
            if(oktest) write(*,'(A)') 'standard step size control' 
            HRVa = 2.95
            HRVb = 3.15
         endif    
      endif

C Find harmonic Ritz values:
      call ritzval(itMax,iterAct,Hhh,HRVr,HRVi,'HRV',Hwork)
C     sub..ritzval(kk,kkAct,Hh,HRVr,HRVi,HRVorRV,Hwork)

C Scale harmonic Ritz values and find closest to origin in
C the left-hand complex half-plane
      HRVmi0 = HRVmin(tau,HRVr,HRVi,iterAct)
      if(oktest) then
         write(*,'(f10.7,A,f13.4)') tau,' >>>>>>> ',HRVmi0
         do i=1,iterAct
            write(*,'(2f12.4)') HRVr(i),HRVi(i)
         enddo
      endif

      tauNew = tau
C Do not make step size control if implmrpc<0, consider step size to be good
      if (implmrpc.ge.0) then
         TauBad = ((HRVmi0.lt.HRVa).or.(HRVmi0.gt.HRVb))
      else
         TauBad = .false.
      endif   

      if(TauBad) then
         tau0 = tau
         if (HRVmi0.lt.HRVa) then
C           Time step is too small
            tau1 = 1.8*tau0
         else
C           Time step is too large
            tau1 = 0.5*tau0
         endif
C        Update the projection matrix w.r.t. tau1
         do i=1,iterAct
            Hhh(i,i)=Hhh(i,i) - 1./tau0 + 1./tau1
         enddo
C        Find harmonic Ritz values:
         call ritzval(itMax,iterAct,Hhh,HRVr,HRVi,'HRV',Hwork)
C        Scale harmonic Ritz values and find closest to origin in
C        the left-hand complex half-plane
         HRVmi1 = HRVmin(tau1,HRVr,HRVi,iterAct)
         if(oktest) then
           write(*,'(f10.7,A,f13.4)') tau1,' >>>>>>> ',HRVmi1
           do i=1,iterAct
              write(*,'(2f12.4)') HRVr(i),HRVi(i)
           enddo
         endif

C        Check whether tau1 is a good time step
         TauBad = ((HRVmi1.lt.HRVa).or.(HRVmi1.gt.HRVb))
         tauNew = tau1
      endif
      
      iii = 1
      do while (TauBad)
C        The secant method to find good step size
         tauNew = tau1 - (HRVmi1-0.5*(HRVa+HRVb))*
     &                               (tau1-tau0)/(HRVmi1-HRVmi0)
         tauNew = dmax1(tauNew,1.d-16)
         tau0 = tau1
         HRVmi0 = HRVmi1
         tau1 = tauNew
C        Update the projection matrix w.r.t. tau1
         do i=1,iterAct
            Hhh(i,i)=Hhh(i,i) - 1./tau0 + 1./tau1
         enddo
C        Find harmonic Ritz values:
         call ritzval(itMax,iterAct,Hhh,HRVr,HRVi,'HRV',Hwork)
C        Scale harmonic Ritz values and find closest to origin in
C        the left-hand complex half-plane
         HRVmi1 = HRVmin(tau1,HRVr,HRVi,iterAct)

C        Check whether tau1 is a good time step
         TauBad = ((HRVmi1.lt.HRVa).or.(HRVmi1.gt.HRVb))
         if(oktest) then 
            write(*,'(f13.4,A,f13.4)') tau1,' >>>>>>> ',HRVmi1
            do i=1,iterAct
               write(*,'(2f12.4)') HRVr(i),HRVi(i)
            enddo
         endif
         iii = iii + 1
         if ((iii.ge.9).and.TauBad) then
C           Quit the secant method loop
            if (HRVmi1.gt.HRVb) then
               tauNew = 0.5*tau
            else
               tauNew = tau
            endif
            TauBad = .false.
         endif   
      enddo

      if(oktest) write(*,'(A,f12.4,A,f12.4)') ' dt =',tau,' dtNew =',
     &           tauNew

      if (tauNew.ne.tau) then
         if (typepred.ne.'explicit2') then
C           Step size can (and will) be changed immediately
            if(oktest) write(*,'(A,f12.4,A,f12.4)') 'VACMRPC: ',HRVmi1,
     &                 ' Ok: dt is changed by factor >>>>>',(tauNew/tau)
            if (nonzero) then
C              adjust rhs of the LS-problem for the new time step: 
               if(nonzero) resnorm0 = (tauNew/tau)*resnorm0
             
C              rescale Predictor guess:
               coeff = tauNew/tau
               call dscal(n,coeff,x,1)
            endif 
            tau = tauNew
         else
C           Step size will be changed next time-step

C           Recreate the projection matrix for the former time step
            do i=1,iterAct
               Hhh(i,i)=Hhh(i,i) - 1./tauNew + 1./tau
            enddo

C           Do not allow growth factor of step size be more than 1+sqrt(2)
            if(tauNew.gt.(2.41*tau)) tauNew = 2.41*tau
            if(oktest) write(*,'(A,f12.4,A,f12.4)') 'VACMRPC: ',HRVmi1,
     &            ' Ok: dt will be changed by factor >>>>>',(tauNew/tau)
         endif
      endif

C     Solve the LS-problem
      xLS(1)=resnorm0
      do i=2,iterAct+1
         xLS(i)=0.0
      enddo
      call solvels(itMax,iterAct,Hhh,xLS,resnorm,Hwork)
C     subroutine solvels(kk,kkAct,Hh,bb,Resid,work)      
  
C Check if the required residual reduction was not achieved
      if (stc.eq.'abs') then
         if (resnorm .gt. eps) info = 1
      elseif(stc.eq.'rel') then
         if ((resnorm/resnorm0).gt.eps) info = 1
      endif 

C Update solution vector x := x + V*xLS
      al = 1.0
      be = 1.0
      call dgemv('N',n,iterAct,al,workV,n,xLS,1,be,x,1)

      if ((implmrpc.eq.gentlecorr.or.implmrpc.eq.mildcorr.or.
     &     implmrpc.eq.standardcorr.or.implmrpc.eq.strictcorr)
     &   .and.(tauNew.lt.(0.8*tau))) then
C        time step is "unstable"--------------------- make extra-corrector

            call matvec(x, workV(1,1), n)
            coeff = -1
            call daxpy(n,coeff,rhs,1,workV(1,1),1)
            call dscal(n,coeff,workV(1,1),1)

         resnorm0 = dnrm2( n, workV(1,1),1 )
         if(oktest) write(*,*)'Resnorm before extra-corrector = ',
     &                                                          resnorm0
  
C        number of MINRES steps for extra-corrector can be set here:
         iter = iter0 - 2
c         if(tauNew.lt.(0.7*tau)) iter = iter0

C        Form new Krylov subspace
         call handv(n,itMax+1,iter,iterAct,workV(1,1),Hhh,workV,matvec)

C        Solve the LS-problem
         xLS(1)=resnorm0
         do i=2,iterAct+1
            xLS(i)=0.0
         enddo
         call solvels(itMax,iterAct,Hhh,xLS,resnorm,Hwork)
C  subroutine solvels(kk,kkAct,Hh,bb,Resid,work)      
         if(oktest) write(*,*)'Resnorm after  extra-corrector = ',       
     &                                                          resnorm

C Update solution vector x := x + V*xLS
         al = 1.0
         be = 1.0
         call dgemv('N',n,iterAct,al,workV,n,xLS,1,be,x,1) 
C        -------------------------------------------- extra-corrector END
      endif

      tau = tauNew

      if(implmrpc.lt.0) tau = -1.D0 
      return
      end

C --------------------------------------------------------------
      double precision function HRVmin(tau,HRVr,HRVi,iterAct)
C     Scales harmonic Ritz values and finds modulus of 
C     the one closest to 0; we ignore positive values
      integer i,iterAct
      double precision dsqrt,HRVabs,HRVr(iterAct),HRVi(iterAct),tau

      HRVmin=1.d16
      do i=1,iterAct
         HRVr(i)=1-tau*HRVr(i)
         HRVi(i)=tau*HRVi(i)
         HRVabs = dsqrt(HRVr(i)*HRVr(i)+HRVi(i)*HRVi(i))
         if((HRVabs.lt.HRVmin).and.(HRVr(i).lt.0.0)) HRVmin=HRVabs
      end do
      return
      end


C --------------------------------------------------------
      subroutine solvels(kk,kkAct,Hh,bb,Resid,work)
C ******************************************************  
C M.Botchev, Math. Dept., Utrecht University, May 1997
C ******************************************************
C solvels solves the least squares (LS) problem Hh*xx=bb
C where Hh is (kkAct+1)-by-kkAct matrix of full rank (kkAct)
C Parameter list:
C kk      == (input) INTEGER. kk+1 is leading dimension of array Hh
C kkAct   == (input) INTEGER size of the problem, i.e.,
C            number of columns in array Hh, kkAct <= kk
C Hh      == (input) DOUBLE PRECISION array dimension (kk+1,kk)
C bb      == (input/output) DOUBLE PRECISION array dimension (kk+1)
C            On entry, right hand side vector of kk+1 elements
C            On exit, the first kk elements contains the LS-solution
C Resid   == (output) DOUBLE PRECISION l2-norm of residual
C            ||bb-Hh*xx|| where xx is the LS-solution
C work    == (workspace) DOUBLE PRECISION array dimension (kk+2,kk+2)
C ------------------------------------------------------
C Subroutines used:
C LAPACK subroutines
C      dgeqrf == QR factorization of general matrix
C      dormqr == multiplication by orthogonal matrix 
C      dtrtrs == solution of the linear system with triangular matrix
C ------------------------------------------------------

      integer kk,kkAct,kwork,Info,i,j
      double precision dabs,Resid,bb(kk+1),Hh(kk+1,kk),work(kk+2,kk+2)

C Place Hh the first kk+1 rows and kk columns of work
      do i=1,kkAct+1
         do j=1,kkAct
            work(i,j)=Hh(i,j)
         end do
      end do

C QR factorization of Hh, the upper triangular matrix R is in work
      kwork=kkAct+1
      call dgeqrf(kkAct+1,kkAct,work,kk+2,work(1,kk+1),work(1,kk+2),
     &            kkAct,Info)
      if(Info.ne.0)then
           write(*,'(A,i2)') 'solvels - Error in DGEQRF: info = ',Info
           return
      endif

C Multiplying bb=Q^T*bb we reduce the LS-problem to R*xx=bb
      call dormqr('L','T',kkAct+1,1,kkAct,work,kk+2,work(1,kk+1),bb,
     &            kk+1,work(1,kk+2),kkAct,Info)
      if(Info.ne.0)then
           write(*,'(A,i2)') 'solvels - Error in DORMQR: info = ',Info
           return
      endif

C The last component of bb is residual norm of the LS problem
          Resid = dabs(bb(kkAct+1))

C Solve system with triangular matrix R*xx=bb
      call dtrtrs('U','N','N',kkAct,1,work,kk+2,bb,kk+1,Info)        
      if(Info.ne.0)then
           write(*,'(A,i2)') 'solvels - Error in DTRTRS: info = ',Info
           return
      end if

      return
      end

C --------------------------------------------------------
      subroutine ritzval(kk,kkAct,Hh,HRVr,HRVi,HRVorRV,Hwork)
C ******************************************************  
C M.Botchev, Math. Dept., Utrecht University, May 1997
C ******************************************************
C For a given upper-Hessenberg projection matrix Hh,
C ritzval calculates Ritz values or harmonic Ritz values
C of the original matrix
C This algorithm is introduced in:
C M.B. van Gijzen, "A Polynomial Preconditioner for the 
C GMRES algorithm", J. of Computational and Applied Math.,
C 1995, vol 59, pp. 91-107.
C Parameter list:
C kk      == (input) INTEGER; kk+1 is leading dimension of Hh
C            kk can not be more than 100(=kkMax - enlarge if necessary) 
C kkAct   == (input) INTEGER ACTual dimension kkAct <= kk
C            i.e., the projection matrix is Hh(kkAct+1,kkAct)
C Hh      == (input) DOUBLE PRECISION array dimension (kk+1,kk)
C HRVr    == (output) DOUBLE PRECISION array dimension (kkAct) 
C            the computed real parts of (harmonic) Ritz values
C HRVi    == (output) DOUBLE PRECISION array dimension (kkAct) 
C            the computed imaginary parts of (harmonic) Ritz values
C HRVorRV == (input) CHARACTER*3
C            set HRVorRV = 'HRV' to get harmonic Ritz values
C            otherwise Ritz values will be calculated
C Hwork   == (workspace) DOUBLE PRECISION array dimension 
C            (kk+2,kk+2)
C ------------------------------------------------------
C Subroutines used:
C BLAS subroutines
C      dcopy  == makes a copy of vector y <-- x
C      dger   == rank-one update of a matrix
C LAPACK subroutines
C      dgbsv  == solves banded system of linear equations
C      dhseqr == computes eigenvalues of upper Hessenberg matrix 
C -------------------------------------------------------
      integer kkMax
      parameter (kkMax=100)
  
      character*3 HRVorRV
      integer kk,kkAct,i,j,Info,kl,ku,Ipiv(kkMax)
      double precision al,Hh(kk+1,kk),
     &                HRVr(kkAct),HRVi(kkAct),Hwork(kk+2,kk+2)

      if (kkAct.gt.kkMax) then
         write(*,'(A)') 
     &   'ERROR in ritzval: input parameter kkAct cannot me more than'
         write(*,'(i4,A)') kkMax,'=kkMax.  Enlarge constant kkMax'
         return
      endif

C -----------------------------------------------if (HRVorRV.eq.'HRV') then
      if (HRVorRV.eq.'HRV') then
C Rewrite Hh^T into the first kkAct columns of Hwork in the 
C LAPACK storage format for band matrices
C Set column kkAct+1 of Hwork to the kkAct-th basis vector 
C e_kkAct=(0,...,0,1)^T 
         kl = 1
         ku = kkAct-1
         do i=1,kkAct
            do j=1,kkAct
            Hwork(kl+ku+1+i-j,j) = Hh(j,i)
            end do
            Hwork(i,kkAct+1)=0.0
            HRVr(i) = 0.0
            HRVi(i) = 0.0
         end do
         Hwork(kkAct,kkAct+1)=1.0

C Solve linear system Hh(1:kkAct,1:kkAct) f_kkAct = e_kkAct  -->
C the solution is in vector Hwork(1,kkAct+2)
         call dgbsv(kkAct,kl,ku,1,Hwork,kk+2,Ipiv,Hwork(1,kkAct+1),
     &           kk+2,Info)
         call dcopy(kkAct,Hwork(1,kkAct+1),1,Hwork(1,kkAct+2),1)
         if(Info.ne.0)then
           write(*,'(A,i2)') 'Error in DGBSV: info = ',Info
           return
         endif
      endif
C --------------------------------------------------------- endif

C Clean up first kkAct columns of Hwork and place in its
C first kkAct rows and first kkAct columns those of Hh 
      do j=1,kkAct
         do i=1,kkAct
         Hwork(i,j) = Hh(i,j)
         end do
         Hwork(kkAct+1,j) = 0.0
         Hwork(kkAct+2,j) = 0.0
      end do

C -------------------------------------- if (HRVorRV.eq.'HRV') then
      if (HRVorRV.eq.'HRV') then  
C Set kkAct+1 column of Hwork to the kkAct-th basis vector e_kkAct=(0,...,0,1)^T
         do i=1,kkAct
            Hwork(i,kkAct+1) = 0.0
         end do
         Hwork(kkAct,kkAct+1) = 1.0

C Make rank-one update of Hwork(1:kkAct,1:kkAct) 
C Hwork = Hwork + (Hh(kkAct+1,kkAct)^2)*f_kkAct*e_kkAct^T
C then Hwork contains inverse of the generalized inverse - 
C eigenvalues of this matrix are harmonic Ritz values
         al = Hh(kkAct+1,kkAct)*Hh(kkAct+1,kkAct)
         call dger(kkAct,kkAct,al,Hwork(1,kkAct+2),1,Hwork(1,kkAct+1),
     &             1,Hwork,kk+2)
      endif
C --------------------------------------------------------- endif

C Calculate eigenvalues of square upper Hessenberg matrix
C Hwork(1:kkAct,1:kkAct) - they are (harmonic) Ritz values
C Workspace variable al is not referenced by the DHSEQR 
C Column kkAct+1 of Hwork is used like a workspace
      call dhseqr('E','N',kkAct,1,kkAct,Hwork,kk+2,HRVr,HRVi,
     &            al,1,Hwork(1,kkAct+1),i,Info)

      if(Info.ne.0)then
        write(*,'(A,i2)') 'Error in DHSEQR: info = ',Info
        return
      end if

      return
      end

C --------------------------------------------------------
      subroutine handv(n,LDH,dimGMR,dimGMRact,res0,Hh,workV,matvec)
C ******************************************************  
C M.Botchev, Math. Dept., Utrecht University, March 1997
C ******************************************************
C handv delivers matrices H and V as defined by the Arnoldi
C algorithm in GMRES framework
C Parameter list:
C n         == (input) INTEGER leading dimension of workV,
C              i.e., size of the problem
C LDH       == (input) INTEGER leading dimension of Hh
C dimGMR    == (input) desired dimension of Krylov subspace
C dimGMRact == (output) actual dimension of Krylov subspace
C              is less than dimGMR if the next vector of 
C              Krylov basis ( workV(1,i1) ) occurs to be zero
C res0      == (input) DOUBLE PRECISION array dimension (n),
C              initial vector for Krylov subspace (typically,
C              (scaled) initial residual)
C Hh        == (output) DOUBLE PRECISION array dimension (LDH,LDH-1)
C              upper Hessenberg projection matrix 
C workV     == (workspace/output) DOUBLE PRECISION array 
C              dimension (n,dimGMR+1) 
C              On exit, workV(n,dimGMRact+1) is matrix V
C matvec    == (input) EXTERNAL name of the MatVec subroutine   
C ------------------------------------------------------
C subroutines used =
C ope(n,x,y,matvec) == matrix vector multiplication y <- A*x
C                      envoking subroutine matvec
C BLAS:
C ddot        == dot product function
C dnrm2       == euclidean norm function
C daxpy       == y <-- y+ax routine
C dcopy       == y <-- x  routine 
C dscal       == x <-- ax routine
C ------------------------------------------------------
      external matvec

      integer dimGMR,dimGMRact,i,i1,j,LDH,n
      double precision ddot,dnrm2,epsmax,Hh(LDH,*),res0(n),t,
     &                 workV(n,dimGMR+1)
      
      data epsmax / 1.d-16 /

      do i=1,dimGMR+1
         do j=1,dimGMR
         Hh(i,j) = 0.0
         end do
      end do

      t =  dnrm2( n, res0,1 )
      if (t.gt.0) then
         t = 1.0 / t    
      else
         dimGMRact = 0
         return
      end if
      call dcopy( n, res0,1, workV(1,1),1 ) 
      call dscal( n, t, workV(1,1),1 )

C      call show(n,1,workV(1,1),' V1:')                 ! debug     
 
      dimGMRact = dimGMR 
      do i=1,dimGMR
         i1 = i + 1
C         call ope( n, workV(1,i), workV(1,i1), matvec )
         call matvec(workV(1,i), workV(1,i1), n)
C        modified Gram - Schmidt... ------------------ BEGIN 
	 do j=1,i
	    t = ddot( n, workV(1,j),1, workV(1,i1),1 )
	    Hh(j,i) = t
	    call daxpy( n, -t, workV(1,j),1, workV(1,i1),1 )
         end do
	 t = dnrm2(n, workV(1,i1),1 )
	 Hh(i1,i) = t
	 if ( (1.0+epsmax).ge.(1.0+t) ) then
            t = 0.0
            dimGMRact = i
C            write(*,'(A,i3)') 'Krylov dimension has been decreased to ',
C     &            i
         else 
            t = 1.0 / t
         end if
         call dscal(n, t, workV(1,i1), 1)
C        modified Gram - Schmidt... ------------------ END
         if ( t.eq.0.0 ) then
            return
         end if
      end do

      return

      end

C=================================================================
