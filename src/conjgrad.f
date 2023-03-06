C This library contains CG and and BiCGSTAB(l) written by Auke van der Ploeg.

C*****************************************************************************

      SUBROUTINE CG77(outp,matvec,R,X,nonzerox,N,tol,stc,iter,Z,AZ,INFO)
C
C       February 12, 1997, written by Auke van der Ploeg
C       February 16, 1997, Gabor Toth added NONZEROX parameter and 
C                                     improved diagnostics
C
C     This subroutine determines the solution of A(X)=R, where
C     the matrix-vector multiplication with A is performed by
C     the subroutine 'matvec'.
C !!!
C !!! If the matrix is not symmetric positive definite, CG 
C !!! is likely to fail.
C !!!
C
C     Description of arguments:
C
C     outp: (input) (boolean) 
C        Determines whether of not output is printed.
C     matvec: external routine for matrix-vector multiplication.
C     R: (input/output)
C        on input:  right-hand side vector.
C        on output: residual vector.
C     X: (input/output)
C        on input:  initial quess for the solution vector.
C        on output: solution vector.
C     nonzerox: (input) (boolean)
C        Tells CG if initial guess in X is zero or not. 
C        If nonzerox is .FALSE., one MATVEC call is saved.
C     N: (input)
C        Number of unknowns.
C     stc (input) (character*3)
C       Determine stopping criterion (||.|| denotes the 2-norm):
C       stc='rel'    -- relative stopping crit.: ||res|| < tol*||res0||
C       stc='abs'    -- absolute stopping crit.: ||res||<tol
C       stc='max'    -- maximum  stopping crit.: max(abs(res))<tol
C     tol: (input/output) 
C       on input:  required (relative) 2-norm or maximum norm of residual
C       on output: achieved (relative) 2-norm or maximum norm of residual
C     iter: (input/output)
C       on input:  maximum number of iterations to be performed.
C       on output: actual  number of iterations done.
C     Z,AZ: work arrays.
C     INFO (output)
C         Gives reason for returning:
C                 0 - solution found satisfying given tolerance.
C                 1 - iteration aborted due to division by very small value.
C                 2 - no convergence within maximum number of iterations.
C                 3 - Initial guess satisfies the stopping criterion.
C                 4 - Initial residual smaller than parameter 'assumedzero'.
C
C     The CG-algorithm is implemented as shown on page 12 of the thesis
C      "Preconditioning for sparse matrices with applications."
C      Auke van der Ploeg
C      University of Groningen, 1994.
C
      IMPLICIT NONE
      CHARACTER*3 stc
      LOGICAL nonzerox,outp
      Integer N,i,itr,matv,iter,INFO
      REAL*8 X(N),R(N),Z(N),AZ(N),rho,rhonew,res,res0,bet,dnrm2,
     *       alf,assumedzero,ddot,tol
C
      external matvec
C ----------------------------------------------------------------------------
      assumedzero=1.D-16
      itr=0
      matv=0
      IF ((stc.NE.'rel').and.(stc.NE.'abs').and.(stc.NE.'max')) THEN
         PRINT *,'Error in CG:'
         PRINT *,'Parameter STC=',stc,' should be one of rel/abs/max.'
         STOP
      END IF
      IF (outp) PRINT *,'n gives the number of CG-iterations.'
C
C       Calculate the initial residual R:=RHS-A*X and its 2-norm.
C
      IF (nonzerox) THEN
         CALL matvec(X,Z,N) 
         matv = matv + 1
         DO i=1,N
            R(i)=R(i)-Z(i)
         ENDDO
      ENDIF
      rho=dnrm2(N,R,1)
      IF (stc.EQ.'max') THEN
         res0=0.D0
         DO i=1,N
            IF (ABS(R(i)).GT.res0) res0=ABS(R(i))
         ENDDO
      ELSE
         res0=rho
      END IF
      assumedzero = assumedzero*res0
      IF (outp) THEN
         IF (stc.EQ.'max') THEN
            PRINT *,'n:',itr,' Maximum norm initial residual:',res0
         ELSE
            PRINT *,'n:',itr,' 2-norm intial residual:',res0
         END IF
      END IF
      IF ((stc.NE.'rel').AND.(res0.LT.tol)) THEN
         INFO = 3
         GO TO 3000
      END IF
      IF (res0.LT.assumedzero) THEN
         INFO = 4
         GO TO 3000
      END IF
      rho=rho*rho
C     Initialize Z:=R-A*X
      call dcopy(N,R,1,Z,1)
C      DO i=1,N
C         Z(i)=R(i)
C      ENDDO

C     Beginning of the REPEAT-UNTIL loop:
 200  CONTINUE
         CALL matvec(Z,AZ,N) 
         matv=matv+1
         alf=ddot(N,Z,1,AZ,1)
         IF (ABS(alf).GT.assumedzero**2) THEN
            alf=rho/alf
         ELSE
            INFO = 1
            GO TO 3000
         END IF
         CALL daxpy(N,alf,Z,1,X,1)
         CALL daxpy(N,-alf,AZ,1,R,1)
         rhonew=dnrm2(N,R,1)
         IF (stc.EQ.'max') THEN
            res=0.D0
            DO i=1,N
               IF (ABS(R(i)).GT.res) res=ABS(R(i))
            ENDDO
         ELSE IF (stc.EQ.'rel') THEN
            res=rhonew/res0
         ELSE
            res=rhonew
         END IF
         rhonew=rhonew*rhonew
         itr=itr+1
         IF (outp) PRINT *,'n:',itr,' ',stc,'. norm of residual:',res
         IF (res.LE.tol) THEN
            INFO = 0
            GO TO 3000
         END IF
         IF (itr.GE.iter) THEN
            INFO = 2
            GOTO 3000
         END IF

         IF (rho.GT.assumedzero**2) THEN
            bet=rhonew/rho
         ELSE
            INFO = 1
            GO TO 3000
         END IF 
         DO i=1,N
            Z(i)=R(i)+bet*Z(i)
         ENDDO
         IF (outp) PRINT *,'alf,bet,rho,rhonew:',alf,bet,rho,rhonew
         rho=rhonew
C     End of REPEAT-UNTIL loop:
      GO TO 200

 3000 CONTINUE
      iter=itr
      tol =res
C
C      Print the results:
C
      IF(outp)then
         IF (INFO.LE.2) THEN
            PRINT *, 'Total Number of CG-iterations:',itr
            PRINT *, 'Number of matrix-vector mult.:',matv
         ENDIF
         IF (INFO.EQ.0) THEN
            PRINT *,'Successful iteration, norm of res. is:',tol
         ELSE IF (INFO.EQ.1) THEN
            PRINT *,'Iteration aborted due to division by a'
            PRINT *,'very small value.' 
         ELSE IF (INFO.EQ.2) THEN
            PRINT *,'Stopping crit. not fulfilled within given'
            PRINT *,'maximum number of iterations.'
         ELSE IF (INFO.EQ.3) THEN
            PRINT *,'Initial guess for the solution satisfies'
            PRINT *,'given stopping criterion.' 
         ELSE IF (INFO.EQ.4) THEN
            PRINT *,'Initial residual is very small.'
         END IF
      ENDIF
C
C       {End of subroutine CG.}
      RETURN
      END

C*****************************************************************************
      SUBROUTINE BICGSTABl(oktest,matvec,R,X,nonzerox,N,tol,stc,
     &     l,iter,work,info)

C     This interface routine chops work array into smaller arrays

      IMPLICIT NONE
      INTEGER N,l,iter,info
      CHARACTER*3 stc
      LOGICAL oktest,nonzerox
      REAL*8 X(N),tol,R(N),work(*)

      external matvec
C----------------------------------------------------------------------------
      CALL BICGSTABl1(oktest,matvec,R,X,nonzerox,N,tol,stc,l,iter,
     &   work,work(N+1),work(2*N+1),work(3*N+1),work((3+l+1)*N+1),info)

      return
      end

C=============================================================================
      SUBROUTINE BICGSTABl1(outp,matvec,R,X,nonzerox,N,tol,stc,l,iter,
     &                      Rs,U,Y,Ud,Rd,INFO)
C
C
C       Original version, 1994, Auke van der Ploeg
C       OUTP feature added on February 11, 1997, Gabor Toth
C       Revised, February 13, 1997, Auke van der Ploeg
C       NONZEROX added and diagnostics improved, February 19, 1977, Gabor Toth
C
C     This subroutine determines the solution of A(X)=R, where
C     the matrix-vector multiplication with A is performed by
C     the subroutine 'matvec'.
C
C     Description of arguments:
C
C     outp: (boolean) determines whether of not output is 
C            printed.
C     matvec: external routine for matrix-vector multiplication.
C     R: (input/output)
C        on input:  right-hand side vector.
C        on output: residual vector.
C     X: (input/output)
C        on input:  initial quess for the solution vector.
C        on output: solution vector.
C     nonzerox: (input) (boolean)
C        Tells CG if initial guess in X is zero or not. 
C        If nonzerox is .FALSE., one MATVEC call is saved.
C     N: (input)
C        number of unknowns.
C     stc (input) (character*3)
C       Determine stopping criterion (||.|| denotes the 2-norm):
C       stc='rel'    -- relative stopping crit.: ||res|| < tol*||res0||
C       stc='abs'    -- absolute stopping crit.: ||res||<tol
C       stc='max'    -- maximum  stopping crit.: max(abs(res))<tol
C     tol: (input/output) 
C       on input:  required (relative) 2-norm or maximum norm of residual.
C       on output: achieved (relative) 2-norm or maximum norm of residual.
C     l: (input) 
C       number of 'inner-'iterations.
C     iter: (input/output)
C       on input:  maximum number of iterations to be performed.
C       on output: actual  number of iterations done.
C     Rs, U, Y, Ud, Rd
C       work arrays.
C     INFO (output)
C         Gives reason for stopping:
C                 0 - solution found satisfying given tolerance.
C                 1 - iteration aborted due to division by very small value.
C                 2 - no convergence within maximum number of iterations.
C                 3 - Initial guess satisfies the stopping criterion.
C                 4 - Initial residual smaller than parameter 'assumedzero'.
C
C     The Bicgstab(l)-algorithm is implemented as shown on page
C     10 of the REPORT 772
C      "BICGSTAB(l) for linear equations involving unsymmetric matrices
C       with complex spectrum." of Gerard Sleijpen and Diederik Fokkema.
C
      IMPLICIT NONE
      CHARACTER*3 stc
      Integer N,i,j,p,itr,k,matv,l,iter,lg,l1,INFO
      LOGICAL outp,nonzerox,more
      PARAMETER (lg=32)
      REAL*8 X(N),rho0,rho1,alf,w,bet,gam,tol,
     *       sig,gam1,tau,gam0,gam2,U,Rs,Ud,Rd,Y,
     *       R(N),res,res0,ff1,ff2,ff,ddot,dnrm2,assumedzero
      DIMENSION sig(lg),gam1(lg),tau(lg,lg),gam0(lg),gam2(lg),
     +          U(N),Rs(N),Ud(N,0:l),Rd(N,0:l),Y(N)

      external matvec
C ---------------------------------------------------------------------------
      assumedzero=1.D-16
      itr=0
      matv=0
      IF ((stc.NE.'rel').and.(stc.NE.'abs').and.(stc.NE.'max')) THEN
         PRINT *,'Error in BiCGSTAB(l):'
         PRINT *,'Parameter STC=',stc,' should be one of rel/abs/max.'
         STOP
      ELSE IF (l.LT.1.OR.l.GT.lg) THEN
         PRINT *,'Error in BiCGSTAB(l):'
         PRINT *,'Parameter L=',l,' should be L>0 and L<LG=',lg
         STOP
      END IF
C
C       Calculate the initial residual R=RHS-A*X and its 2-norm.
C
      IF (nonzerox) THEN
         CALL matvec(X,Rs,N) 
         matv = matv + 1
         DO i=1,N
            R(i)=R(i)-Rs(i)
         ENDDO
      ENDIF
      IF (stc.EQ.'max') THEN
         res0=0.D0
         DO i=1,N
            IF (ABS(R(i)).GT.res0) res0=ABS(R(i))
         ENDDO
      ELSE
         res0=dnrm2(N,R,1)
      END IF
      assumedzero = assumedzero*res0
      IF (outp) THEN
         IF (stc.EQ.'max') THEN
            PRINT *,'n:',itr,' Maximum norm initial residual:',res0
         ELSE
            PRINT *,'n:',itr,' 2-norm intial residual:',res0
         END IF
      END IF
      IF ((stc.NE.'rel').AND.(res0.LT.tol)) THEN
         INFO = 3
         GO TO 3000
      END IF
      IF (res0.LT.assumedzero) THEN
         INFO = 4
         GO TO 3000
      END IF

C     Initialize Rs:=RHS-A*X and U:=0
      DO i=1,N
         Rs(i)=R(i)
         U(i)=0.D0
      ENDDO

      rho0=1.D0
      alf=0.D0
      w=1.D0
      k=-l
      IF (outp) PRINT *,'n gives the number of CGSTAB(l)-iterations.'

C     Beginning of the REPEAT-UNTIL loop:
 200  CONTINUE
        k=k+l
        CALL dcopy(N, U, 1, Ud(1,0), 1)
        CALL dcopy(N, R, 1, Rd(1,0), 1)
        rho0=-w*rho0

C       Bi-CG PART:
        j=0
        more=.TRUE.

 300    CONTINUE
          rho1=ddot(N,Rd(1,j),1,Rs,1)
          IF (ABS(rho0).GT.assumedzero) THEN
            bet=alf*rho1/rho0
          ELSE
            more=.FALSE.
            INFO = 1
            GO TO 1000
          END IF 
          rho0=rho1
          DO i=0,j
            DO p=1,N
              Ud(p,i)=Rd(p,i)-bet*Ud(p,i)
            ENDDO
          ENDDO
          CALL dcopy(N, Ud(1,j), 1, Y, 1)
C
C           U:=A*Y
C
          CALL matvec(Y,U,n) 
          matv=matv+1
          CALL dcopy(N, U, 1, Ud(1,j+1), 1)
          gam=ddot(N,Ud(1,j+1),1,Rs,1)
          IF (ABS(gam).GT.assumedzero) THEN
            alf=rho0/gam
          ELSE
            more=.FALSE.
            INFO = 1
            GO TO 1000
          END IF
          CALL daxpy(N*(j+1),-alf,Ud(1,1),1,Rd(1,0),1)
          CALL dcopy(N, Rd(1,j), 1, Y, 1)
C
C           U:=A*Y
C
          CALL matvec(Y,U,n) 
          matv=matv+1
          CALL dcopy(N, U, 1, Rd(1,j+1), 1)
          CALL daxpy(N,alf,Ud(1,0),1,X,1)
          j=j+1
 1000   IF ((j.LE.l-1).and.more) GO TO 300
        IF (j.EQ.0) GO TO 3000
        l1=j
C
C          END Bi-CG PART,    begin MR PART:
C
        DO j=1,l1
          DO i=1,j-1
            ff=ddot(N,Rd(1,j),1,Rd(1,i),1)
            IF (ABS(sig(i)).GT.(assumedzero*tol)**2) THEN
              ff=ff/sig(i)
            ELSE
              more=.FALSE.
              INFO = 1
              ff=0.D0
            END IF
            tau(i,j)=ff
            CALL daxpy(N,-ff,Rd(1,i),1,Rd(1,j),1)
          ENDDO
          ff1=dnrm2(N,Rd(1,j),1)
          ff1=ff1*ff1
          ff2=ddot(N,Rd(1,j),1,Rd(1,0),1)
          sig(j)=ff1
          IF (ff1.GT.(assumedzero*tol)**2) THEN
            gam1(j)=ff2/ff1
          ELSE
            more=.FALSE.
            INFO = 1
            gam1(j)=0.D0
          END IF
        ENDDO
        gam0(l1)=gam1(l1)
        w=gam0(l1)
        DO j=l1-1,1,-1
          ff=gam1(j)
          DO i=j+1,l1
            ff=ff-tau(j,i)*gam0(i)
          ENDDO
          gam0(j)=ff
        ENDDO
        DO j=1,l1-1
          ff=gam0(j+1)
          DO i=j+1,l1-1
            ff=ff+tau(j,i)*gam0(i+1)
          ENDDO
          gam2(j)=ff
        ENDDO
        DO p=1,N
          X(p)=X(p)+gam0(1)*Rd(p,0)
          R(p)=Rd(p,0)-gam1(l1)*Rd(p,l1)
          U(p)=Ud(p,0)-gam0(l1)*Ud(p,l1)
        ENDDO
        DO j=1,l1-1
          CALL daxpy(N,-gam0(j),Ud(1,j),1,U,1)
          CALL daxpy(N,gam2(j),Rd(1,j),1,X,1)
          CALL daxpy(N,-gam1(j),Rd(1,j),1,R,1)
        ENDDO
        IF (stc.EQ.'max') THEN
          res=0.D0
          DO j=1,N
            IF (ABS(R(j)).GT.res) res=ABS(R(j))
          ENDDO
        ELSE IF (stc.EQ.'abs') THEN
          res=dnrm2(N,R,1)
        ELSE
          res=dnrm2(N,R,1)/res0
        END IF
        itr=itr+1
        IF (outp) PRINT *,'n:',itr,' ',stc,'. norm of residual:',res
        IF (res.LE.tol) THEN
           INFO = 0
           GO TO 3000
        END IF
        IF (itr.GE.iter) THEN
           INFO = 2
           GOTO 3000
        END IF

C     End of the REPEAT-UNTIL loop
      IF (more) GO TO 200

 3000 CONTINUE
      iter=itr
      tol=res
C
C      Print the results:
C
      IF(outp)then
         IF (INFO.LE.2) THEN
            PRINT *, 'Number of BiCGSTAB(l)-iterations:',itr
            PRINT *, 'Number of matrix-vector mult.:',matv
         ENDIF
         IF (INFO.EQ.0) THEN
            PRINT *,'Successful iteration, norm of res. is:',tol
         ELSE IF (INFO.EQ.1) THEN
            PRINT *,'Iteration aborted due to division by a'
            PRINT *,'very small value.' 
         ELSE IF (INFO.EQ.2) THEN
            PRINT *,'Stopping crit. not fulfilled within given'
            PRINT *,'maximum number of iterations.'
         ELSE IF (INFO.EQ.3) THEN
            PRINT *,'Initial guess for the solution satisfies'
            PRINT *,'given stopping criterion.' 
         ELSE IF (INFO.EQ.4) THEN
            PRINT *,'Initial residual is very small.'
         END IF
      ENDIF
C
C       {End of subroutine BiCGSTAB(l).}
      RETURN
      END
C*****************************************************************************
