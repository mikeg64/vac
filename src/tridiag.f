C=============================================================================
C This library contains the subroutines for the block tridiagonal solver
C and the pentadiagonal solver, both written by Auke van der Ploeg.
C The cyclic block tridiagonal solver was implemented by 
C Gabor Toth and Auke van der Ploeg.
C
C subroutines: preprocLO, solutionLO                          ! TRIDIAG
C              precyclic, solcyclic                           ! CYCLIC TRIDIAG
C              prepenta, pentapreceis, pentasoleis, postpenta ! PENTADIAG
C              prehepta, heptapreceis, heptasoleis, posthepta ! HEPTADIAG
C              solveb                                         ! modified LAPACK
C
C=============================================================================
C********************** TRIDIAGONAL ******************************************
C=============================================================================
      subroutine preprocLO(e,d,f,pivot,nblock,N)
C
C     This routine performs the preprocessing necessary for solving
C     a block tridiagonal system. It constructs
C     a complete block LU-decomposition of the tri-diagonal matrix. 
C     The coefficient matrix is overwritten by the factors L and U.
C 
      IMPLICIT NONE
      INTEGER N,nblock,pivot(N,nblock)
      REAL*8 e(N,N,nblock),d(N,N,nblock),f(N,N,nblock)
C
C=======================================================================
C         Description of arguments:
C=======================================================================
C
C     e, d, f, nblock, N:
C               The matrix is assumed to be block tridiagonal,
C               The size of the blocks is N*N, and the total number
C               of unknowns is N*nblock.
C               The blocks are stored as follows:
C
C               e(j-1)  d(j-1)  f(j-1)
C                       e(j)    d(j)    f(j)
C                               e(j+1)  d(j+1)  f(j+1)
C
C              for j=2..nblock-1, where e(1) and f(nblock) are 
C              zero blocks. It is assumed that the
C              blocks are not very sparse, so the sparsity pattern
C              within the separate blocks is not exploited.
C              For example, the (i,k)-element of the j-th block on the 
C              main diagonal is stored in d(i,k,j).
C              On output, the matrix L+U is stored as follows:
C              The blocks on the main diagonal of U are unity blocks,
C              the blocks on the main diagonal of L are stored in d,
C              e contains the subdiagonal blocks of L,
C              f contains the superdiagonal blocks of U.
C     pivot:   integer array which contains the sequence generated
C              by partial pivoting. used in subroutine 'DGetrf'.
C
C=======================================================================
C         Local variables:
C=======================================================================
C
      INTEGER j,INFO
      REAL*8 one
      PARAMETER (one=1.D0)
C
C     j       :integers necessary for implementing loops.
C     INFO:    (output) integer for LAPACK subroutines.
C     one:
C              variable necessary for BLAS-routine dgemm. 
C
C=======================================================================
C         External functions
C=======================================================================
C
C     dgemm,   BLAS level three Matrix-Matrix Product.
C              See 'man dgemm' for description.
C
C     DGetrf,  LAPACK routine, computes an LU factorization of a general 
C              M-by-N matrix A using partial pivoting with 
C              row interchanges. The factorization has the form
C                  A = P * L * U
C              where P is a permutation matrix, L is lower triangular 
C              with unit diagonal elements (lower trapezoidal if m > n),
C              and U is upper triangular (upper trapezoidal if m < n).
C              This is the right-looking Level 3 BLAS version of the
C              algorithm.
C
C     DGetrs,  LAPACK routine, solves a system of linear equations
C              A * X = B,  A**T * X = B,  or  A**H * X = B
C              with a general N-by-N matrix A using the LU factorization
C              computed by DGETRF.
C
C         Compute the LU-decomposition of d(1):
C
      CALL DGETRF( N, N, d(1,1,1), N, pivot(1,1), INFO )
C
C        Calculate inv(1)*f(1), where inv(1) is the inverse of
C        the block on the main diagonal.
C
      CALL DGETRS('n',N,N,d(1,1,1),N,pivot(1,1),f(1,1,1),N,INFO)
      DO j=2,nblock-1
C
C         Calculate d(j)-e(j)*f(j-1) and put it in d(j):
C  
        CALL dgemm('n','n',N,N,N,-one,e(1,1,j),N,f(1,1,j-1),N,
     &              one,d(1,1,j),N)
C
C         Compute the LU-decomposition of d(j):
C
        CALL DGETRF( N, N, d(1,1,j), N, pivot(1,j), INFO )
C
C        Calculate inv(j)*f(j), where inv(j) is the inverse of
C        the block on the main diagonal.
C
        CALL DGETRS('n',N,N,d(1,1,j),N,pivot(1,j),f(1,1,j),N,INFO)
      ENDDO
      j=nblock
      CALL dgemm('n','n',N,N,N,-one,e(1,1,j),N,f(1,1,j-1),N,
     &            one,d(1,1,j),N)
      CALL DGETRF( N, N, d(1,1,j), N, pivot(1,j), INFO )
      RETURN
      END
C=============================================================================
      subroutine solutionLO(e,d,f,x,nblock,N,pivot,scratch)
C
C     This routine solves a block tri-diagonal system with a block
C     LU-decomposition.
C     It is assumed that the preprocessing has already been performed
C     in subroutine 'preprocLO.' 
C
      IMPLICIT NONE
      INTEGER N,nblock,pivot(N,nblock)
      REAL*8 e(N,N,nblock),d(N,N,nblock),
     *           f(N,N,nblock),x(N,nblock),scratch(N)
C
C=======================================================================
C         Description of arguments:
C=======================================================================
C
C     e, d, f, nblock, N:
C              The matrix L+U is assumed to be block tridiagonal,
C              The size of the blocks is N*N, and the total number
C              of unknowns is N*nblock.
C              It is assumed that the
C              blocks are not very sparse, so the sparsity pattern
C              within the separate blocks is not exploited.
C              For example, the (i,k)-element of the j-th block on the 
C              main diagonal is stored in d(i,k,j).
C              The blocks on the main diagonal of U are unity blocks,
C              the blocks on the main diagonal of L are stored in d,
C              e contains the subdiagonal blocks of L,
C              f contains the superdiagonal blocks of U.
C     x:       On input, the right-hand side is stored in x.
C              This vector is overwritten by the solution. 
C              The vector x is partitioned according to
C              the blockstructure of the matrices. For example, the 
C              (j-1)*N+i-th component is stored in x(i,j), for
C              i=0,1,...N-1 and j=1,2,...,nblock. 
C
C=======================================================================
C         Local variables:
C=======================================================================
C
      INTEGER i,k,j, Bg
      REAL*8 one, ff
      PARAMETER (one=1.D0, Bg = 7)
C
C     Bg   :   For small values of the blocksize (<Bg) 
C              BLAS-routines for matrix-vector multiplication
C              are not called. The user is advised to experiment 
C              with several values for Bg for optimal performance.
C     i,k,j:   integers necessary for implementing loops.
C     one  :   variable necessary for BLAS-routine DGEMV. 
C     ff   :   work variable.
C
C=======================================================================
C         External function
C=======================================================================
C
C     DGEMV,   BLAS level two Matrix-Vector Product.
C              See 'man DGEMV' for description.
C
C     solveb,  solves a system of linear equations
C              A * X = B,  A**T * X = B,  or  A**H * X = B
C              with a general N-by-N matrix A using the LU factorization
C              computed by DGETRF.
C
C        First solve Ly=x, put the solution again in x:
C
      IF (N.LT.Bg) GO TO 1000
      CALL solveb( N, d(1,1,1), pivot(1,1), x(1,1))
      DO j=2,nblock
        CALL DGEMV('n',N,N,-one,e(1,1,j),N,x(1,j-1),1,one,x(1,j),1)
        CALL solveb( N, d(1,1,j), pivot(1,j), x(1,j))
      ENDDO
C
C        Now solve Uy=x, put the solution again in x:
C
      DO j=nblock-1,1,-1
        CALL DGEMV('n',N,N,-one,f(1,1,j),N,x(1,j+1),1,one,x(1,j),1)
      ENDDO
      RETURN
C
C================================================================
C
C      (N.LT.Bg)
C
C================================================================
C
 1000 CONTINUE
      CALL solveb( N, d(1,1,1), pivot(1,1), x(1,1))
      DO j=2,nblock
        DO i=1,N
          ff=x(i,j)
          DO k=1,N
            ff=ff-e(i,k,j)*x(k,j-1)
          ENDDO
          scratch(i)=ff
        ENDDO
        DO i=1,N
          ff=scratch(pivot(i,j))
          scratch(pivot(i,j))=scratch(i)
          DO k=1,i-1
            ff=ff-d(i,k,j)*scratch(k)
          ENDDO
          scratch(i)=ff
        ENDDO
        DO i=N,1,-1
          ff=scratch(i)
          DO k=i+1,N
            ff=ff-d(i,k,j)*scratch(k)
          ENDDO
          scratch(i)=ff/d(i,i,j)
          x(i,j)=scratch(i)
        ENDDO
      ENDDO
C
C        Now solve Uy=x, put the solution again in x:
C
      DO j=nblock-1,1,-1
        DO i=1,N
          ff=x(i,j)
          DO k=1,N
            ff=ff-f(i,k,j)*x(k,j+1)
          ENDDO
          x(i,j)=ff
        ENDDO
      ENDDO

      RETURN
      END
C=============================================================================
C********************** CYCLIC TRIDIAGONAL ***********************************
C=============================================================================
      subroutine precyclic(e,d,f,pivot,nblock,N)
C
C     This routine modifies the block-cyclic matrix to tridiagonal 
C     and does the LU-decomposition for it.
C 
      IMPLICIT NONE
      INTEGER N,nblock,pivot(N,nblock)
      REAL*8 e(N,N,nblock),d(N,N,nblock),f(N,N,nblock)
C
C=======================================================================
C         Description of arguments:
C=======================================================================
C
C     e, d, f, nblock, N:
C               The matrix is assumed to be block tridiagonal,
C               The size of the blocks is N*N, and the total number
C               of unknowns is N*nblock.
C               The blocks are stored as follows:
C
C               d(1)    f(1)                      e(1)
C               e(j-1)  d(j-1)  f(j-1)
C                       e(j)    d(j)    f(j)
C                               e(j+1)  d(j+1)    f(j+1)
C              f(nblock)                e(nblock) d(nblock)
C
C              for j=2..nblock-1, where e(1) and f(nblock) are 
C              blocks coming from periodic boundary conditions.
C     pivot:   integer array which contains the sequence generated
C              by partial pivoting. used in subroutine 'DGetrf'.
C
C=======================================================================
C         Local variables:
C=======================================================================
C
      INTEGER i,k

C----------------------------------------------------------------------

      do i=1,N
         do k=1,N
            d(i,k,1)=d(i,k,1)-e(i,k,1)
            d(i,k,nblock)=d(i,k,nblock)-f(i,k,nblock)
         enddo
      enddo

      call preprocLO(e,d,f,pivot,nblock,N)

      return
      end

C======================================================================
      subroutine solcyclic(e,d,f,x,pivot,z,h,pivoth,
     &     nblock,N,scratch)
C
C     This routine solves a block cyclic tri-diagonal system with a block
C     LU-decomposition.
C     It is assumed that the preprocessing has already been performed
C     in subroutine 'precyclic.' 
C
      IMPLICIT NONE
      INTEGER N,nblock,pivot(N,nblock),pivoth(N)
      REAL*8 e(N,N,nblock),d(N,N,nblock),f(N,N,nblock),
     *       x(N,nblock),z(N,nblock,N),h(N,N),scratch(N)
C
C=======================================================================
C         Description of arguments:
C=======================================================================
C
C     e, d, f, nblock, N:
C              The matrix L+U is assumed to be block tridiagonal,
C              with cyclic blocks e(1) and f(nblock) coming from
C              periodic boundary conditions. 
C              The size of the blocks is N*N, and the total number
C              of unknowns is N*nblock.
C              It is assumed that the
C              blocks are not very sparse, so the sparsity pattern
C              within the separate blocks is not exploited.
C              For example, the (i,k)-element of the j-th block on the 
C              main diagonal is stored in d(i,k,j).
C              The blocks on the main diagonal of U are unity blocks,
C              the blocks on the main diagonal of L are stored in d,
C              e contains the subdiagonal blocks of L,
C              f contains the superdiagonal blocks of U.
C     x:       On input, the right-hand side is stored in x.
C              This vector is overwritten by the solution. 
C              The vector x is partitioned according to
C              the blockstructure of the matrices. For example, the 
C              (j-1)*N+i-th component is stored in x(i,j), for
C              i=0,1,...N-1 and j=1,2,...,nblock. 
C
C=======================================================================
C         Local variables:
C=======================================================================
C
      INTEGER j, i, k, INFO
      REAL*8 one
      PARAMETER (one=1.D0)
C
C     j,i,k:   integers necessary for implementing loops.
C     INFO :   used for error indication by DGETRS.
C     one  :   variables necessary for BLAS-routine DGEMV. 
C
C=======================================================================
C         External function
C=======================================================================
C
C     DGEMV,   BLAS level two Matrix-Vector Product.
C              See 'man DGEMV' for description.
C
C     solveb,  solves a system of linear equations
C              A * X = B,  A**T * X = B,  or  A**H * X = B
C              with a general N-by-N matrix A using the LU factorization
C              computed by DGETRF.
C
C        First solve Ly=x, put the solution again in x:
C

C-----------------------------------------------------------------------


C     Identity part of H
      do i=1,N
         do k=1,N
             h(i,k)=0
         enddo
         h(i,i)=1
      enddo

C     Put 0-s into z
      do k=1,N
         do j=2,nblock-1
            do i=1,N
               z(i,j,k) = 0
            end do
         end do
      end do

C     Cycle over columns in U and form RHS-s
      do i=1,N
         do k=1,N
            z(k,1,i)      = e(k,i,1)
            z(k,nblock,i) = f(k,i,nblock)
         end do

C        Calculate Z(*,*,i)
         call solutionLO(e,d,f,z(1,1,i),nblock,N,pivot,scratch)

C        Form H=I+V^T*Z

         do k=1,N
            h(k,i)=h(k,i)+z(k,1,i)+z(k,nblock,i)
         enddo
      enddo

C     decompose H
      CALL DGETRF( N, N, h, N, pivoth, INFO )
      
C     Solve A'x=rhs

      call solutionLO(e,d,f,x,nblock,N,pivot,scratch)

C     Calculate V^T.y
      do i=1,N
         scratch(i)=x(i,1)+x(i,nblock)
      end do

C     Calculate H.(V^T.y)
      call solveb(N,h,pivoth,scratch)

C     Calculate x --> x-Z.(H.(V^T.y)))
      CALL DGEMV('n',N*nblock,N,-one,z,N*nblock,scratch,1,one,x,1)

      return
      end
C=============================================================================
C********************** PENTADIAGONAL ****************************************
C=============================================================================
      subroutine prepenta(dm,e,f,e1,f1,d,pivot,nblock,N,M,alf,x,b,work)
C
C       Auke van der Ploeg, Januari 8, 1997.
C
C     This routine performs the preprocessing necessary for constructing
C     a preconditioner for a block penta-diagonal system. It constructs
C     an incomplete block LU-decomposition of the penta-diagonal matrix,
C     in such a way that L+U has the same blockstructure as A.
C     Further, it generates a new right-hand side and starting vector
C     belonging to the preconditioned system
C
C         L^{-1}AU^{-1}x=L^{-1}b.
C 
C !!! This routine also scales the original matrix A, right-hand side
C     vector b and the starting vector x.
C
C     It can also deal with blocks in the
C     upper-right and lower-left corner that come from the treatment
C     of periodic boundary conditions.
C
C     Only the blocks on the main diagonal change, which enables the
C     efficient Eisenstat implementation.
*
*          See, for example, page 19 of the Phd thesis:
*          'Preconditioning for sparse matrices with applications.'
*          1994, University of Groningen.
*
C
C===================================================================
C         Gustafsson modification
C===================================================================
C
C     It is possible to encorporate the so-called Gustafsson
C     modification for the blocks on the main diagonal.
C     In this appoach, a splitting A=LU+R, is constructed in such a way
C     that the block row sums of R are zero. For systems of
C     linear equations coming from problems with an elliptic character
C     this can strongly improve the convergence behaviour of 
C     iterative methods. See page 22 of Phd thesis 'Preconditioning 
C     for sparse ..... ' for an illustration of this phenomenon.
C 
      IMPLICIT NONE
      INTEGER N,M,nblock,pivot(N,nblock)
      REAL*8 alf
      REAL*8 e(N,N,nblock), d(N,N,nblock), f(N,N,nblock),
     *       e1(N,N,nblock), f1(N,N,nblock), dm(N,N,nblock), 
     *       x(N*nblock), b(N*nblock), work(N*nblock)
C
C=======================================================================
C         Description of arguments:
C=======================================================================
C
C     alf       The parameter for Gustafsson modification. alf>=0 means 
C               no modification, 0> alf >= -1 is the valid parameter range
C     e, d, f, e1, f1, nblock, N, M:
C               The matrix L+U is assumed to be block pentadiagonal
C               with a blockstructure corrosponding with a five-point
C               stencil. The distance from the outermost blocks
C               (except those coming from periodic boundary conditions)
C               to the diagonal blocks is M.
C               The size of the subblocks is N*N.
C               The total number of diangonal blocks is nblock, hence  
C               the total number of unknowns is N*nblock.
C               In addition, the matrix contains blocks coming from
C               the treatment of periodic boundary conditions, M blocks
C               in the Lower-Left corner, and M blocks in the 
C               Upper-Right corner. The distance from these blocks
C               to the blocks on the main diagonal is nblock-M.
C
C               The blocks of L+U are stored as follows:
C               d(j): j=1..nblock        scaling blocks. L and U are
C                                        scaled in such a way that 
C                                        their blocks on the main
C                                        diagonal are unity blocks.
C               e(j): j=2..nblock        sub diagonal blocks, e(1)=0.
C               f(j): j=1..nblock-1      super diagonal blocks,
C                                        f(nblock)=0.
C               e1(j): j=M+1..nblock     blocks in the lower-triangular
C                                        part with distance M from 
C                                        the main diagonal.
C               f1(j): j=1..nblock-M     blocks in the upper-triangular
C                                        part with distance M from 
C                                        the main diagonal.
C               e1(j): j=1..M            blocks coming from periodic
C                                        boundary conditions in the 
C                                        upper-right corner, distance
C                                        nblock-M from blocks on the
C                                        main diagonal.
C               f1(j): j=1..M            blocks coming from periodic
C                                        boundary conditions in the 
C                                        lower-left corner, distance
C                                        nblock-M from blocks on the
C                                        main diagonal.
C
C              It is assumed that the
C              blocks are not very sparse, so the sparsity pattern
C              within the separate blocks is not exploited.
C              For example, the (i,k)-element of the j-th block on the 
C              main diagonal is stored in d(i,k,j).
C     dm:      on entrance: main diagonal DA of the coefficient
C              matrix A.
C              on exit: DA-2I, scaled in such a way that the main
C              diagonal of L is one.
C     pivot:   integer array which contains the sequence generated
C              by partial pivoting. used in subroutine 'DGetrf'.
C
C=======================================================================
C         Local variables:
C=======================================================================
C
      INTEGER i
      REAL*8 one
      PARAMETER (one=1.D0)
C
C     i        integer necessary for implementing loops.
C     one:     variable necessary for BLAS-routine dgemv. 
C
C=======================================================================
C         External functions
C=======================================================================
C
C     dgemv,   BLAS level three Matrix-Vector Product.
C              See 'man dgemv' for description.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      EXECUTABLE STATEMENTS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       Construct the incomplete block LU-decomposition:
C
      CALL pentapreceis(d,e,f,e1,f1,dm,pivot,nblock,N,M,alf)
C
C     b:=D^{-1}b
C
      DO i=1,nblock
        CALL solveb(N, d(1,1,i), pivot(1,i), b(1+(i-1)*N))
      ENDDO
C
C     b:=L^{-1}b
C
      DO i=2,nblock
        IF (i.GT.M) THEN
          CALL DGEMV('n',N,N,-one,e1(1,1,i),N,
     &           b(1+(i-M-1)*N),1,one,b(1+(i-1)*N),1)
        END IF
        CALL DGEMV('n',N,N,-one,e(1,1,i),N,
     &           b(1+(i-2)*N),1,one,b(1+(i-1)*N),1)
C
C        Deal with blocks in the lower-left corner coming from
C        periodic boundary conditions:
C
        IF (i+M.GT.nblock) 
     $     CALL DGEMV('n',N,N,-one,f1(1,1,i),N,
     $       b(1+(i-nblock+M-1)*N),1,one,b(1+(i-1)*N),1)
      ENDDO
C
C      x:=Ux
C
      CALL DCOPY(N*nblock, x, 1, work, 1)
      DO i=1,nblock-1
        CALL DGEMV('n',N,N,one,f(1,1,i),N,
     $             work(1+i*N),1,one,x(1+(i-1)*N),1)
        IF (i+M.LE.nblock) 
     $   CALL DGEMV('n',N,N,one,f1(1,1,i),N,
     $        work(1+(i+M-1)*N),1,one,x(1+(i-1)*N),1)
C
C        Deal with blocks in the upper-right corner coming from
C        periodic boundary conditions:
C
        IF (i.LE.M) 
     $   CALL DGEMV('n',N,N,one,e1(1,1,i),N,
     $    work(1+(i+nblock-M-1)*N),1,one,x(1+(i-1)*N),1)
      ENDDO

      RETURN
      END
C=============================================================================
      subroutine pentapreceis(d,e,f,e1,f1,dm,pivot,nblock,N,M,alf)
C
C     Auke van der Ploeg, December 18 1996
C
C     This routine performs the preprocessing necessary for constructing
C     a preconditioner for a block penta-diagonal system. It constructs
C     an incomplete block LU-decomposition of the penta-diagonal matrix,
C     in such a way that L+U has the same blockstructure as A.
C 
C !!! This routine also scales the original matrix A.
C
C     It can also deal with blocks in the
C     upper-right and lower-left corner that come from the treatment
C     of periodic boundary conditions.
C
C     Only the blocks on the main diagonal change, which enables the
C     efficient Eisenstat implementation.
*
*          See, for example, page 19 of the Phd thesis:
*          'Preconditioning for sparse matrices with applications.'
*          1994, University of Groningen.
*
C
C===================================================================
C         Gustafsson modification
C===================================================================
C
C     It is possible to encorporate the so-called Gustafsson
C     modification for the blocks on the main diagonal.
C     In this appoach, a splitting A=LU+R, is constructed in such a way
C     that the block row sums of R are zero. For systems of
C     linear equations coming from problems with an elliptic character
C     this can strongly improve the convergence behaviour of 
C     iterative methods. See page 22 of Phd thesis 'Preconditioning 
C     for sparse ..... ' for an illustration of this phenomenon.
C 
      IMPLICIT NONE
      INTEGER N,M,nblock,pivot(N,nblock)
      REAL*8 alf
      REAL*8 e(N,N,nblock), d(N,N,nblock), f(N,N,nblock),
     *           e1(N,N,nblock), f1(N,N,nblock), dm(N,N,nblock)
C
C=======================================================================
C         Description of arguments:
C=======================================================================
C
C     alf       Parameter for the Gustafsson modification. alf>=0 
C               means no modification, 0>alf>=-1 is the valid parameter range.
C     e, d, f, e1, f1, nblock, N, M:
C               The matrix L+U is assumed to be block pentadiagonal
C               with a blockstructure corrosponding with a five-point
C               stencil. The distance from the outermost blocks
C               (except those coming from periodic boundary conditions)
C               to the diagonal blocks is M.
C               The size of the subblocks is N*N.
C               The total number of diangonal blocks is nblock, hence  
C               the total number of unknowns is Nt=N*nblock.
C               In addition, the matrix contains blocks coming from
C               the treatment of periodic boundary conditions, M blocks
C               in the Lower-Left corner, and M blocks in the 
C               Upper-Right corner. The distance from these blocks
C               to the blocks on the main diagonal is nblock-M.
C
C               The blocks of L+U are stored as follows:
C               d(j): j=1..nblock        scaling blocks. L and U are
C                                        scaled in such a way that 
C                                        their blocks on the main
C                                        diagonal are unity blocks.
C               e(j): j=2..nblock        sub diagonal blocks, e(1)=0.
C               f(j): j=1..nblock-1      super diagonal blocks,
C                                        f(nblock)=0.
C               e1(j): j=M+1..nblock     blocks in the lower-triangular
C                                        part with distance M from 
C                                        the main diagonal.
C               f1(j): j=1..nblock-M     blocks in the upper-triangular
C                                        part with distance M from 
C                                        the main diagonal.
C               e1(j): j=1..M            blocks coming from periodic
C                                        boundary conditions in the 
C                                        upper-right corner, distance
C                                        nblock-M from blocks on the
C                                        main diagonal.
C               f1(j): j=1..M            blocks coming from periodic
C                                        boundary conditions in the 
C                                        lower-left corner, distance
C                                        nblock-M from blocks on the
C                                        main diagonal.
C
C              It is assumed that the
C              blocks are not very sparse, so the sparsity pattern
C              within the separate blocks is not exploited.
C              For example, the (i,k)-element of the j-th block on the 
C              main diagonal is stored in d(i,k,j).
C     dm:      on entrance: main diagonal DA of the coefficient
C              matrix A.
C              on exit: DA-2I, scaled in such a way that the main
C              diagonal of L is one.
C     pivot:   integer array which contains the sequence generated
C              by partial pivoting. used in subroutine 'DGetrf'.
C
C=======================================================================
C         Local variables:
C=======================================================================
C
      INTEGER i,j,INFO
      REAL*8 one, mone
      PARAMETER (one=1.D0, mone=-1.D0)
      LOGICAL gus
C
C     j        :integers necessary for implementing loops.
C     INFO:    (output) integer for LAPACK subroutines.
C     one, mone:
C              variables necessary for BLAS-routine dgemm. 
C     gus:     true if -1<=alf<0, i.e. Gustaffson modification is used
C=======================================================================
C         External functions
C=======================================================================
C
C     dgemm,   BLAS level three Matrix-Matrix Product.
C              See 'man dgemm' for description.
C
C     DGetrf,  LAPACK routine, computes an LU factorization of a general 
C              M-by-N matrix A using partial pivoting with 
C              row interchanges. The factorization has the form
C                  A = P * L * U
C              where P is a permutation matrix, L is lower triangular 
C              with unit diagonal elements (lower trapezoidal if m > n),
C              and U is upper triangular (upper trapezoidal if m < n).
C              This is the right-looking Level 3 BLAS version of the
C              algorithm.
C
C     DGetrs,  LAPACK routine, solves a system of linear equations
C              A * X = B,  A**T * X = B,  or  A**H * X = B
C              with a general N-by-N matrix A using the LU factorization
C              computed by DGETRF.
C
      CALL dcopy(N*N*nblock,dm(1,1,1),1,d(1,1,1),1) 
C
C       Relaxed form of Gustafsson modification:
C
      gus= alf.lt.0.D0 .and. alf.ge.-1.D0
C
C      CALL dscal( N*N*nblock, alf, d, 1)
C
C         Compute the LU-decomposition of d(1):
C
      CALL DGETRF( N, N, d(1,1,1), N, pivot(1,1), INFO )
      CALL DGETRS('n',N,N,d(1,1,1),N,pivot(1,1),dm(1,1,1),N,INFO)
      DO i=1,N
        dm(i,i,1)=dm(i,i,1)-2.D0
      ENDDO
      CALL DGETRS('n',N,N,d(1,1,1),N,pivot(1,1),f(1,1,1),N,INFO)
      CALL DGETRS('n',N,N,d(1,1,1),N,pivot(1,1),f1(1,1,1),N,INFO)
      CALL DGETRS('n',N,N,d(1,1,1),N,pivot(1,1),e1(1,1,1),N,INFO)
      DO j=2,nblock
        IF (j.GT.M) THEN
C
C         Interaction boundary conditions lower-left corner
C         (blocks stored in f1(nblock-M+1..nblock) and
C         upper-right corner (blocks stored in e1(1..M).
C
          IF (j+M.GT.nblock)
     $      CALL dgemm('n','n',N,N,N,mone,f1(1,1,j),N,
     $                e1(1,1,j-nblock+M),N,one,d(1,1,j),N)
          CALL dgemm('n','n',N,N,N,mone,e1(1,1,j),N,f1(1,1,j-M),N,
     $              one,d(1,1,j),N)
        END IF
        CALL dgemm('n','n',N,N,N,mone,e(1,1,j),N,f(1,1,j-1),N,
     $             one,d(1,1,j),N)
        IF (gus) THEN
          IF (j.GT.M) THEN
            IF (j+M.GT.nblock) THEN
               CALL dgemm('n','n',N,N,N,alf,f1(1,1,j),N,
     $                f(1,1,j-nblock+M),N,one,d(1,1,j),N)
               CALL dgemm('n','n',N,N,N,alf,f1(1,1,j),N,
     $                f1(1,1,j-nblock+M),N,one,d(1,1,j),N)
            END IF
            CALL dgemm('n','n',N,N,N,alf,e1(1,1,j),N,f(1,1,j-M),N,
     $                one,d(1,1,j),N)
            IF (j.LE.2*M)
     $        CALL dgemm('n','n',N,N,N,alf,e1(1,1,j),N,e1(1,1,j-M),
     $                N,one,d(1,1,j),N)
          END IF
          IF (j.LE.M+1)
     $     CALL dgemm('n','n',N,N,N,alf,e(1,1,j),N,e1(1,1,j-1),N,
     $               one,d(1,1,j),N)
          IF (j-1+M.LE.nblock)
     $     CALL dgemm('n','n',N,N,N,alf,e(1,1,j),N,f1(1,1,j-1),N,
     $               one,d(1,1,j),N)
        END IF
C
C         Compute the LU-decomposition of d(j):
C
        CALL DGETRF( N, N, d(1,1,j), N, pivot(1,j), INFO )
        CALL DGETRS('n',N,N,d(1,1,j),N,pivot(1,j),e1(1,1,j),N,INFO)
        CALL DGETRS('n',N,N,d(1,1,j),N,pivot(1,j),e(1,1,j),N,INFO)
        CALL DGETRS('n',N,N,d(1,1,j),N,pivot(1,j),dm(1,1,j),N,INFO)
        DO i=1,N
          dm(i,i,j)=dm(i,i,j)-2.D0
        ENDDO
        IF (j.LT.nblock)
     $   CALL DGETRS('n',N,N,d(1,1,j),N,pivot(1,j),f(1,1,j),N,INFO)
        CALL DGETRS('n',N,N,d(1,1,j),N,pivot(1,j),f1(1,1,j),N,INFO)
      ENDDO
      RETURN
      END

C=============================================================================
      subroutine pentasoleis(d,e,f,e1,f1,x,y,v,nblock,N,M)
C
C     Auke van der Ploeg, December 18 1996
C
C     This routine applies the preconditonering which  
C     must have been constructed in subroutine 'pentapreceis.' 
*          Contains the matrix-vector multiplication
*
*               x := L^{-1}AU^{-1}y
*  
*          implemented according to the efficient
*          Eisenstat implementation.
*
*          See, for example, page 19 of the Phd thesis:
*          'Preconditioning for sparse matrices with applications.'
*          1994, University of Groningen.
*
*          Vector y remains unchanged.
C
      IMPLICIT NONE
      INTEGER N,M,nblock
      REAL*8 e(N,N,nblock),d(N,N,nblock),e1(N,N,nblock),f1(N,N,nblock),
     *           f(N,N,nblock),x(N,nblock),y(N,nblock),v(N,nblock)
C
C=======================================================================
C         Description of arguments:
C=======================================================================
C
C     e, d, f, e1, f1, nblock, N, M:
C               The matrix L+U is assumed to be block pentadiagonal
C               with a blockstructure corrosponding with a five-point
C               stencil. The distance from the outermost blocks
C               (except those coming from periodic boundary conditions)
C               to the diagonal blocks is M.
C               The size of the subblocks is N*N.
C               The total number of diangonal blocks is nblock, hence  
C               the total number of unknowns is Nt=N*nblock.
C               In addition, the matrix contains blocks coming from
C               the treatment of periodic boundary conditions, M blocks
C               in the Lower-Left corner, and M blocks in the 
C               Upper-Right corner. The distance from these blocks
C               to the blocks on the main diagonal is nblock-M.
C
C               The blocks of L+U are stored as follows:
C               d(j): j=1..nblock        contains D-2I 
C                                        in which D contains the blocks
C                                        on the main diagonal of A.
C               e(j): j=2..nblock        sub diagonal blocks, e(1)=0.
C               f(j): j=1..nblock-1      super diagonal blocks,
C                                        f(nblock)=0.
C               e1(j): j=M+1..nblock     blocks in the lower-triangular
C                                        part with distance M from 
C                                        the main diagonal.
C               f1(j): j=1..nblock-M     blocks in the upper-triangular
C                                        part with distance M from 
C                                        the main diagonal.
C               e1(j): j=1..M            blocks coming from periodic
C                                        boundary conditions in the 
C                                        upper-right corner, distance
C                                        nblock-M from blocks on the
C                                        main diagonal.
C               f1(j): j=1..M            blocks coming from periodic
C                                        boundary conditions in the 
C                                        lower-left corner, distance
C                                        nblock-M from blocks on the
C                                        main diagonal.
C
C              It is assumed that the
C              blocks are not very sparse, so the sparsity pattern
C              within the separate blocks is not exploited.
C              For example, the (i,k)-element of the j-th block on the 
C              main diagonal is stored in d(i,k,j).
C     y:       right-hand side
C     x:       On output, the solution is stored in x.
C     v:       scratch variable
C              The vector x, y and v are partitioned according to
C              the blockstructure of the matrices. For example, the 
C              (j-1)*N+i-th component is stored in x(i,j), for
C              i=0,1,...N-1 and j=1,2,...,nblock. 
C
C=======================================================================
C         Local variables:
C=======================================================================
C
      INTEGER j, i, k, Ng
      REAL*8 one, vj1, vj2, vj3
      PARAMETER (one=1.D0, Ng=7)
C
C     Ng   :   for large values of N, each matrix-vector with a
C              subblock is implemented with a  BLAS-routine. For small 
C              values of N (N<Ng) this gives too much overhead, so it 
C              is better to write out the loops, which enables
C              more optimalizations by the compiler.
C     j,i  :   integers necessary for implementing loops.
C     INFO :   used for error indication by DGETRS.
C     one, :   variable necessary for BLAS-routine DGEMV. 
C
C=======================================================================
C         External function
C=======================================================================
C
C     DGEMV,   BLAS level two Matrix-Vector Product.
C              See 'man DGEMV' for description.
C
C=======================================================================

      IF (N.LT.Ng) GO TO 1000
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       Code for N >= Ng:
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C      v:=(I-U)^{-1}y
C
      CALL DCOPY( N*nblock, y(1,1), 1, v(1,1), 1)
      DO j=nblock-1,1,-1
        CALL DGEMV('n',N,N,-one,f(1,1,j),N,v(1,j+1),1,one,v(1,j),1)
        IF (j+M.LE.nblock) 
     $   CALL DGEMV('n',N,N,-one,f1(1,1,j),N,v(1,j+M),1,one,v(1,j),1)
C
C        Deal with blocks in the upper-right corner coming from
C        periodic boundary conditions:
C
        IF (j.LE.M) 
     $   CALL DGEMV('n',N,N,-one,e1(1,1,j),N,v(1,j+nblock-M),1,
     $               one,v(1,j),1)
      ENDDO
C
C      x:=y+(D-2I)v, D-2I is stored in d.
C
      CALL DCOPY( N*nblock, y(1,1), 1, x(1,1), 1)
      DO j=1,nblock
        CALL DGEMV('n',N,N,one,d(1,1,j),N,v(1,j),1,one,x(1,j),1)
      ENDDO
C
C        solve (I-L)y=x, put the solution again in x:
C
      DO j=2,nblock
        IF (j.GT.M)
     $     CALL DGEMV('n',N,N,-one,e1(1,1,j),N,x(1,j-M),1,one,x(1,j),1)
        CALL DGEMV('n',N,N,-one,e(1,1,j),N,x(1,j-1),1,one,x(1,j),1)
C
C        Deal with blocks in the lower-left corner coming from
C        periodic boundary conditions:
C
        IF (j+M.GT.nblock) 
     $     CALL DGEMV('n',N,N,-one,f1(1,1,j),N,x(1,j-nblock+M),1,
     $                one,x(1,j),1)
      ENDDO
      CALL daxpy( nblock*N, one, v(1,1), 1, x(1,1), 1 )
      RETURN
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       Code for 1 < N < Ng:
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      v:=(I-U)^{-1}y
C
 1000 IF (N.EQ.1) GO TO 2000
      CALL DCOPY( N*nblock, y(1,1), 1, v(1,1), 1)
      DO j=nblock-1,1,-1
        IF (j.LE.M) THEN
C
C        Deal with blocks in the upper-right corner coming from
C        periodic boundary conditions:
C
          DO k=1,N
            vj1 = v(k,j+M)
            vj2 = v(k,j+1)
            vj3 = v(k,j+Nblock-M)
            DO i=1,N
              v(i,j) = v(i,j) - f1(i,k,j)*vj1 -
     &          f(i,k,j)*vj2 - e1(i,k,j)*vj3
            ENDDO
          ENDDO
        ELSE
          IF (j+M.LE.nblock) THEN
            DO k=1,N
              vj1 = v(k,j+M)
              vj2 = v(k,j+1)
              DO i=1,N
                v(i,j) = v(i,j) - f1(i,k,j)*vj1 - f(i,k,j)*vj2
              ENDDO
            ENDDO
          ELSE
            DO k=1,N
              vj2 = v(k,j+1)
              DO i=1,N
                v(i,j) = v(i,j) - f(i,k,j)*vj2
              ENDDO
            ENDDO
          END IF
        END IF
      ENDDO
C
C      x := y + (D-2I)v, D-2I is stored in d.
C
      CALL DCOPY( N*nblock, y(1,1), 1, x(1,1), 1)
      DO j=1,nblock
        DO k=1,N
          vj1 = v(k,j)
          DO i=1,N
            x(i,j) = x(i,j) + d(i,k,j)*vj1
          ENDDO
        ENDDO
      ENDDO
C
C        solve (I-L)z = x, put the solution again in x:
C
      DO j=2,nblock
        IF (j+M.GT.nblock) THEN
C
C        Deal with blocks in the lower-left corner coming from
C        periodic boundary conditions:
C
          DO k=1,N
            vj1 = x(k,j-M)
            vj2 = x(k,j-1)
            vj3 = x(k,j-nblock+M)
            DO i=1,N
              x(i,j) = x(i,j) - e1(i,k,j)*vj1 - 
     &                e(i,k,j)*vj2 - f1(i,k,j)*vj3
            ENDDO
          ENDDO
        ELSE
          IF (j.GT.M) THEN
            DO k=1,N
              vj1 = x(k,j-M)
              vj2 = x(k,j-1)
              DO i=1,N
                x(i,j) = x(i,j) - e1(i,k,j)*vj1 - e(i,k,j)*vj2
              ENDDO
            ENDDO
          ELSE
            DO k=1,N
              vj2 = x(k,j-1)
              DO i=1,N
                x(i,j) = x(i,j) - e(i,k,j)*vj2
              ENDDO
            ENDDO
          END IF
        END IF
      ENDDO
      CALL daxpy( nblock*N, one, v(1,1), 1, x(1,1), 1 )
      RETURN
C
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       Code for N=1:
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      v:=(I-U)^{-1}y
C
 2000 v(1,nblock)=y(1,nblock)
      DO j=nblock-1,1,-1
        IF (j.LE.M) THEN
C
C        Deal with blocks in the upper-right corner coming from
C        periodic boundary conditions:
C
          v(1,j) = y(1,j) - f1(1,1,j)*v(1,j+M) -
     &        f(1,1,j)*v(1,j+1) - e1(1,1,j)*v(1,j+nblock-M)
        ELSE
          IF (j+M.LE.nblock) THEN
            v(1,j) = y(1,j) - f1(1,1,j)*v(1,j+M) - f(1,1,j)*v(1,j+1)
          ELSE
            v(1,j) = y(1,j) - f(1,1,j)*v(1,j+1) 
          END IF
        END IF
      ENDDO
C
C      x := y + (D-2I)v, D-2I is stored in d.
C
      DO j=1,nblock
        x(1,j) = y(1,j) + d(1,1,j)*v(1,j)
      ENDDO
C
C        solve (I-L)z = x, put the solution again in x:
C
      DO j=2,nblock
        IF (j+M.GT.nblock) THEN
C
C        Deal with blocks in the lower-left corner coming from
C        periodic boundary conditions:
C
          x(1,j) = x(1,j) - e1(1,1,j)*x(1,j-M) - 
     &            e(1,1,j)*x(1,j-1) - f1(1,1,j)*x(1,j-nblock+M)
        ELSE
          IF (j.GT.M) THEN
            x(1,j) = x(1,j) - e1(1,1,j)*x(1,j-M) - e(1,1,j)*x(1,j-1) 
          ELSE
            x(1,j) = x(1,j) - e(1,1,j)*x(1,j-1) 
          END IF
        END IF
      ENDDO
      CALL daxpy( nblock, one, v(1,1), 1, x(1,1), 1 )
      RETURN
      END

C=============================================================================
      SUBROUTINE postpenta(f,e1,f1,nblock,N,M,x)
C
C     Original by Auke van der Ploeg, December 18 1996
C     Modified by Gabor Toth,         January 14, 1997
C
C     Postprocessing pentadiagonal solver: x:=U^{-1}x
C
      IMPLICIT NONE
      INTEGER N,M,nblock
      REAL*8 f(N,N,nblock),e1(N,N,nblock),f1(N,N,nblock),x(N,nblock)
C
C=======================================================================
C         Description of arguments:
C=======================================================================
C
C     f, e1, f1, nblock, N, M:
C               The matrix U is the upper part of a block pentadiagonal
C               with a blockstructure corrosponding with a five-point
C               stencil. The distance from the outermost blocks
C               (except those coming from periodic boundary conditions)
C               to the diagonal blocks is M.
C               The size of the subblocks is N*N.
C               The total number of diangonal blocks is nblock, hence  
C               the total number of unknowns is Nt=N*nblock.
C               In addition, the matrix contains blocks coming from
C               the treatment of periodic boundary conditions, M blocks
C               in the Lower-Left corner, and M blocks in the 
C               Upper-Right corner. The distance from these blocks
C               to the blocks on the main diagonal is nblock-M.
C
C               The blocks of U are stored as follows:
C               f(j): j=1..nblock-1      super diagonal blocks,
C                                        f(nblock)=0.
C               e1(j): j=M+1..nblock     blocks in the lower-triangular
C                                        part with distance M from 
C                                        the main diagonal.
C               f1(j): j=1..nblock-M     blocks in the upper-triangular
C                                        part with distance M from 
C                                        the main diagonal.
C               e1(j): j=1..M            blocks coming from periodic
C                                        boundary conditions in the 
C                                        upper-right corner, distance
C                                        nblock-M from blocks on the
C                                        main diagonal.
C               f1(j): j=1..M            blocks coming from periodic
C                                        boundary conditions in the 
C                                        lower-left corner, distance
C                                        nblock-M from blocks on the
C                                        main diagonal.
C
C              It is assumed that the
C              blocks are not very sparse, so the sparsity pattern
C              within the separate blocks is not exploited.
C              For example, the (i,k)-element of the j-th block on the 
C              superdiagonal is stored in f(i,k,j).
C     x:       on input:  solution of preconditioned system,
C              on output: solution of the original system.
C
C=======================================================================
C         Local variables:
C=======================================================================
C
      INTEGER i
      REAL*8 one
      PARAMETER (one=1.D0)

C     i    :   integer necessary for implementing loops.
C     one  :   variable necessary for BLAS-routine DGEMV. 
C
C=======================================================================

      DO i=nblock-1,1,-1
         CALL DGEMV('n',N,N,-one,f(1,1,i),N,
     $        x(1,i+1),1,one,x(1,i),1)
         IF (i+M.LE.nblock)
     $   CALL DGEMV('n',N,N,-one,f1(1,1,i),N,
     $        x(1,i+M),1,one,x(1,i),1)
         IF (i.LE.M)
     $   CALL DGEMV('n',N,N,-one,e1(1,1,i),N,
     $        x(1,i+nblock-M),1,one,x(1,i),1)
      ENDDO

      RETURN
      END
C=============================================================================
C********************** HEPTADIAGONAL ****************************************
C=============================================================================

      subroutine prehepta(dm,e,f,e1,f1,e2,f2,d,pivot,nblock,N,
     &                         M1,M2,alf,x,b,work)
C
C       Auke van der Ploeg, Januari 8, 1997.
C
C     This routine performs the preprocessing necessary for constructing
C     a preconditioner for a block hepta-diagonal system. It constructs
C     an incomplete block LU-decomposition of the hepta-diagonal matrix,
C     in such a way that L+U has the same blockstructure as A.
C     Further, it generates a new right-hand side and starting vector
C     belonging to the preconditioned system
C
C         L^{-1}AU^{-1}x=L^{-1}b.
C 
C !!! This routine also scales the original matrix A, right-hand side
C     vector b and the starting vector x.
C
C     It can also deal with blocks in the
C     upper-right and lower-left corner that come from the treatment
C     of periodic boundary conditions.
C
C     Only the blocks on the main diagonal change, which enables the
C     efficient Eisenstat implementation.
*
*          See, for example, page 19 of the Phd thesis:
*          'Preconditioning for sparse matrices with applications.'
*          1994, University of Groningen.
*
C
C===================================================================
C         Gustafsson modification
C===================================================================
C
C     It is possible to encorporate the so-called Gustafsson
C     modification for the blocks on the main diagonal.
C     In this appoach, a splitting A=LU+R, is constructed in such a way
C     that the block row sums of R are zero. For systems of
C     linear equations coming from problems with an elliptic character
C     this can strongly improve the convergence behaviour of 
C     iterative methods. See page 22 of Phd thesis 'Preconditioning 
C     for sparse ..... ' for an illustration of this phenomenon.
C 
      IMPLICIT NONE
      INTEGER N,M1,M2,nblock,pivot(N,nblock)
      REAL*8 alf
      REAL*8 e(N,N,nblock), d(N,N,nblock), f(N,N,nblock),
     *       e1(N,N,nblock), f1(N,N,nblock), dm(N,N,nblock), 
     *       x(N*nblock), b(N*nblock), work(N*nblock),
     *       e2(N,N,nblock), f2(N,N,nblock)
C=======================================================================
C         Description of arguments:
C=======================================================================
C
C     dm:      on entrance: main diagonal DA of the coefficient
C              matrix A.
C              on exit: DA-2I, scaled in such a way that the main
C              diagonal of L is one.
C     d, e, f, e1, f1, e2, f2:
C               The matrix L+U is assumed to be block heptadiagonal.
C               In addition, the matrix contains blocks coming from
C               the treatment of periodic boundary conditions, M2 blocks
C               in the Lower-Left corner, and M2 blocks in the 
C               Upper-Right corner. The distance from these blocks
C               to the blocks on the main diagonal is nblock-M2.
C
C               The blocks are stored as follows:
C               d(j): j=1..nblock        scaling blocks. L and U are
C                                        scaled in such a way that 
C                                        their blocks on the main
C                                        diagonal are unity blocks.
C               e(j): j=2..nblock     sub diagonal blocks.
C               f(j): j=1..nblock-1   super diagonal blocks.
C               e1(j): j=M1+1..nblock     blocks in the lower-triangular
C                                        part with distance M1 from 
C                                        the main diagonal.
C               f1(j): j=1..nblock-M1     blocks in the upper-triangular
C                                        part with distance M1 from 
C                                        the main diagonal.
C               e2(j): j=M2+1..nblock     blocks in the lower-triangular
C                                        part with distance M2 from 
C                                        the main diagonal.
C               f2(j): j=1..nblock-M2     blocks in the upper-triangular
C                                        part with distance M2 from 
C                                        the main diagonal.
C               e2(j): j=1..M2            blocks coming from periodic
C                                        boundary conditions in the 
C                                        upper-right corner, distance
C                                        nblock-M2 from blocks on the
C                                        main diagonal.
C               f2(j): j=nblock-M2+1..nblock:
C                                        blocks coming from periodic
C                                        boundary conditions in the 
C                                        lower-left corner, distance
C                                        nblock-M2 from blocks on the
C                                        main diagonal.
C
C              It is assumed that the
C              blocks are not very sparse, so the sparsity pattern
C              within the separate blocks is not exploited.
C              For example, the (i,k)-element of the j-th block on the 
C              main diagonal is stored in mat(i,k,j,3).
C
C     pivot:   integer array which contains the sequence generated
C              by partial pivoting. used in subroutine 'DGetrf'.
C     nblock:  Number of diagonal blocks.
C     N:        the size of the blocks.
C     M1:       distance of blocks to the main diagonal blocks.
C     M2:       distance of outer-most blocks to main diagonal blocks.
C               (except those coming from periodic boundary conditions)
C               1 < M1 < M2.
C               The matrix has a blockstructure corresponding with
C               a seven-point stencil on a three-dimensional, 
C               rectangular grid. The blocks corresonding with the
C               direction in which grid points are numbered
C               first, are the sub- and superdiagonal blocks.
C               The blocks corresponding with the
C               direction in which grid points are numbered
C               secondly, have distance M1 from the main diagonal
C               blocks. Finally, The blocks corresponding with the
C               direction in which grid points are numbered
C               last, have distance M2 from the main diagonal blocks.
C     alf:      The parameter for Gustafsson modification. alf>=0 means 
C               no modification, 0> alf >= -1 is the valid parameter range.
C     x:      the solution is stored in x.
C     b:      right-hand side vector.
C     work:   work array
C
C=======================================================================
C         Local variables:
C=======================================================================
C
      INTEGER i
      REAL*8 one
      PARAMETER (one=1.D0)
C
C     i        integer necessary for implementing loops.
C     one:     variable necessary for BLAS-routine dgemv. 
C
C=======================================================================
C         External functions
C=======================================================================
C
C     dgemv,   BLAS level three Matrix-Vector Product.
C              See 'man dgemv' for description.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      EXECUTABLE STATEMENTS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       Construct the incomplete block LU-decomposition:
C
      CALL heptapreceis(d,e,f,e1,f1,e2,f2,dm,pivot,nblock,N,M1,M2,alf)
C
C     b:=D^{-1}b
C
      DO i=1,nblock
        CALL solveb(N, d(1,1,i), pivot(1,i), b(1+(i-1)*N))
      ENDDO
C
C     b:=L^{-1}b
C
      DO i=2,nblock
        IF (i.GT.M2) THEN
          CALL DGEMV('n',N,N,-one,e2(1,1,i),N,
     &           b(1+(i-M2-1)*N),1,one,b(1+(i-1)*N),1)
        END IF
        IF (i.GT.M1) THEN
          CALL DGEMV('n',N,N,-one,e1(1,1,i),N,
     &           b(1+(i-M1-1)*N),1,one,b(1+(i-1)*N),1)
        END IF
        CALL DGEMV('n',N,N,-one,e(1,1,i),N,
     &           b(1+(i-2)*N),1,one,b(1+(i-1)*N),1)
C
C        Deal with blocks in the lower-left corner coming from
C        periodic boundary conditions:
C
        IF (i+M2.GT.nblock) 
     $     CALL DGEMV('n',N,N,-one,f2(1,1,i),N,
     $       b(1+(i-nblock+M2-1)*N),1,one,b(1+(i-1)*N),1)
      ENDDO
C
C      x:=Ux
C
      CALL DCOPY(N*nblock, x, 1, work, 1)
      DO i=1,nblock-1
        CALL DGEMV('n',N,N,one,f(1,1,i),N,
     $             work(1+i*N),1,one,x(1+(i-1)*N),1)
        IF (i+M1.LE.nblock) 
     $   CALL DGEMV('n',N,N,one,f1(1,1,i),N,
     $        work(1+(i+M1-1)*N),1,one,x(1+(i-1)*N),1)
        IF (i+M2.LE.nblock) 
     $   CALL DGEMV('n',N,N,one,f2(1,1,i),N,
     $        work(1+(i+M2-1)*N),1,one,x(1+(i-1)*N),1)
C
C        Deal with blocks in the upper-right corner coming from
C        periodic boundary conditions:
C
        IF (i.LE.M2) 
     $   CALL DGEMV('n',N,N,one,e2(1,1,i),N,
     $    work(1+(i+nblock-M2-1)*N),1,one,x(1+(i-1)*N),1)
      ENDDO
      RETURN
      END

C=============================================================================

      subroutine heptapreceis(d,e,f,e1,f1,e2,f2,dm,pivot,
     &                        nblock,N,M1,M2,alf)
C
C     Auke van der Ploeg, Febuari 7 1997.
C
C     This routine performs the preprocessing necessary for constructing
C     a preconditioner for a block hepta-diagonal system. It constructs
C     an incomplete block LU-decomposition of the hepta-diagonal matrix,
C     in such a way that L+U has the same blockstructure as A.
C 
C !!! This routine also scales the original matrix A.
C
C     It can also deal with blocks in the
C     upper-right and lower-left corner that come from the treatment
C     of periodic boundary conditions.
C
C     Only the blocks on the main diagonal change, which enables the
C     efficient Eisenstat implementation.
*
*          See, for example, page 19 of the Phd thesis:
*          'Preconditioning for sparse matrices with applications.'
*          1994, University of Groningen.
*
C
C===================================================================
C         Gustafsson modification
C===================================================================
C
C     It is possible to encorporate the so-called Gustafsson
C     modification for the blocks on the main diagonal.
C     In this appoach, a splitting A=LU+R, is constructed in such a way
C     that the block row sums of R are zero. For systems of
C     linear equations coming from problems with an elliptic character
C     this can strongly improve the convergence behaviour of 
C     iterative methods. See page 22 of Phd thesis 'Preconditioning 
C     for sparse ..... ' for an illustration of this phenomenon.
C 
      IMPLICIT NONE
      INTEGER N,M1,M2,nblock,pivot(N,nblock)
      REAL*8 alf
      REAL*8 e(N,N,nblock), d(N,N,nblock), f(N,N,nblock),
     *           e1(N,N,nblock), f1(N,N,nblock), dm(N,N,nblock),
     *           e2(N,N,nblock), f2(N,N,nblock)
C
C=======================================================================
C         Description of arguments:
C=======================================================================
C
C     d, e, f, e1, f1, e2, f2,
C               The matrix L+U is assumed to be block heptadiagonal.
C               In addition, the matrix contains blocks coming from
C               the treatment of periodic boundary conditions, M2 blocks
C               in the Lower-Left corner, and M2 blocks in the 
C               Upper-Right corner. The distance from these blocks
C               to the blocks on the main diagonal is nblock-M2.
C
C               The blocks are stored as follows:
C               d(j): j=1..nblock        scaling blocks. L and U are
C                                        scaled in such a way that 
C                                        their blocks on the main
C                                        diagonal are unity blocks.
C               e(j): j=2..nblock     sub diagonal blocks.
C               f(j): j=1..nblock-1   super diagonal blocks.
C               e1(j): j=M1+1..nblock     blocks in the lower-triangular
C                                        part with distance M1 from 
C                                        the main diagonal.
C               f1(j): j=1..nblock-M1     blocks in the upper-triangular
C                                        part with distance M1 from 
C                                        the main diagonal.
C               e2(j): j=M2+1..nblock     blocks in the lower-triangular
C                                        part with distance M2 from 
C                                        the main diagonal.
C               f2(j): j=1..nblock-M2     blocks in the upper-triangular
C                                        part with distance M2 from 
C                                        the main diagonal.
C               e2(j): j=1..M2            blocks coming from periodic
C                                        boundary conditions in the 
C                                        upper-right corner, distance
C                                        nblock-M2 from blocks on the
C                                        main diagonal.
C               f2(j): j=nblock-M2+1..nblock:
C                                        blocks coming from periodic
C                                        boundary conditions in the 
C                                        lower-left corner, distance
C                                        nblock-M2 from blocks on the
C                                        main diagonal.
C
C              It is assumed that the
C              blocks are not very sparse, so the sparsity pattern
C              within the separate blocks is not exploited.
C              For example, the (i,k)-element of the j-th block on the 
C              main diagonal is stored in mat(i,k,j,3).
C
C     dm:      on entrance: main diagonal DA of the coefficient
C              matrix A.
C              on exit: DA-2I, scaled in such a way that the main
C              diagonal of L is one.
C     pivot:   integer array which contains the sequence generated
C              by partial pivoting. used in subroutine 'DGetrf'.
C     nblock:  Number of diagonal blocks.
C     N:        the size of the blocks.
C     M1:       distance of blocks to the main diagonal blocks.
C     M2:       distance of outer-most blocks to main diagonal blocks.
C               (except those coming from periodic boundary conditions)
C               1 < M1 < M2.
C               The matrix has a blockstructure corresponding with
C               a seven-point stencil on a three-dimensional, 
C               rectangular grid. The blocks corresonding with the
C               direction in which grid points are numbered
C               first, are the sub- and superdiagonal blocks.
C               The blocks corresponding with the
C               direction in which grid points are numbered
C               secondly, have distance M1 from the main diagonal
C               blocks. Finally, The blocks corresponding with the
C               direction in which grid points are numbered
C               last, have distance M2 from the main diagonal blocks.
C     alf       The parameter for Gustafsson modification. alf>=0 means 
C               no modification, 0> alf >= -1 is the valid parameter range.
C
C=======================================================================
C         Local variables:
C=======================================================================
C
      INTEGER i,j,INFO
      REAL*8 one
      PARAMETER (one=1.D0)
C
C     i,j:      integers necessary for implementing loops.
C     INFO:    (output) integer for LAPACK subroutines.
C     one:      necessary for BLAS-routine dgemm. 
C
C=======================================================================
C         External functions
C=======================================================================
C
C     dgemm,   BLAS level three Matrix-Matrix Product.
C              See 'man dgemm' for description.
C
C     DGetrf,  LAPACK routine, computes an LU factorization of a general 
C              M-by-N matrix A using partial pivoting with 
C              row interchanges. The factorization has the form
C                  A = P * L * U
C              where P is a permutation matrix, L is lower triangular 
C              with unit diagonal elements (lower trapezoidal if m > n),
C              and U is upper triangular (upper trapezoidal if m < n).
C              This is the right-looking Level 3 BLAS version of the
C              algorithm.
C
C     DGetrs,  LAPACK routine, solves a system of linear equations
C              A * X = B,  A**T * X = B,  or  A**H * X = B
C              with a general N-by-N matrix A using the LU factorization
C              computed by DGETRF.
C
      CALL dcopy(N*N*nblock,dm(1,1,1),1,d(1,1,1),1) 
      IF (alf.LT.0.D0) THEN
C
C       (Relaxed form of) Gustafsson modification:
C
        PRINT *,' '
        PRINT *,'Relaxed form of Gustafsson modification.'
        PRINT *,'Parameter alf is ',alf
        IF (alf.LT.-1.D0) THEN
          PRINT *,'Parmeter alf replaced by -1.D0'
          alf=-1.D0
        END IF
        PRINT *,' '
      ELSE
        alf=0.D0
        PRINT *,' '
        PRINT *,'No Gustafsson modification.'
        PRINT *,' '
      END IF
C
C         Compute the LU-decomposition of d(1):
C
      CALL DGETRF( N, N, d(1,1,1), N, pivot(1,1), INFO )
      CALL DGETRS('n',N,N,d(1,1,1),N,pivot(1,1),dm(1,1,1),N,INFO)
      DO i=1,N
        dm(i,i,1)=dm(i,i,1)-2.D0
      ENDDO
      CALL DGETRS('n',N,N,d(1,1,1),N,pivot(1,1),f(1,1,1),N,INFO)
      CALL DGETRS('n',N,N,d(1,1,1),N,pivot(1,1),f1(1,1,1),N,INFO)
      CALL DGETRS('n',N,N,d(1,1,1),N,pivot(1,1),f2(1,1,1),N,INFO)
      CALL DGETRS('n',N,N,d(1,1,1),N,pivot(1,1),e2(1,1,1),N,INFO)
      DO j=2,nblock
        IF (j.GT.M2) THEN
C
C         Interaction boundary conditions lower-left corner
C         (blocks stored in f2(nblock-M2+1..nblock) and
C         upper-right corner (blocks stored in e2(1..M2).
C
          IF (j+M2.GT.nblock)
     $      CALL dgemm('n','n',N,N,N,-one,f2(1,1,j),N,
     $                e2(1,1,j-nblock+M2),N,one,d(1,1,j),N)
          CALL dgemm('n','n',N,N,N,-one,e2(1,1,j),N,f2(1,1,j-M2),N,
     $                one,d(1,1,j),N)
        END IF
        IF (j.GT.M1) 
     $     CALL dgemm('n','n',N,N,N,-one,e1(1,1,j),N,f1(1,1,j-M1),N,
     $              one,d(1,1,j),N)
        CALL dgemm('n','n',N,N,N,-one,e(1,1,j),N,f(1,1,j-1),N,
     $             one,d(1,1,j),N)
        IF (alf.LT.0.D0) THEN
          IF (j.GT.M1) THEN
            IF (j+M2.GT.nblock) THEN
               CALL dgemm('n','n',N,N,N,alf,f2(1,1,j),N,
     $                f(1,1,j-nblock+M2),N,one,d(1,1,j),N)
               CALL dgemm('n','n',N,N,N,alf,f2(1,1,j),N,
     $                f1(1,1,j-nblock+M2),N,one,d(1,1,j),N)
               CALL dgemm('n','n',N,N,N,alf,f2(1,1,j),N,
     $                f2(1,1,j-nblock+M2),N,one,d(1,1,j),N)
            END IF
            IF (j.GT.M2) THEN
              CALL dgemm('n','n',N,N,N,alf,e2(1,1,j),N,f(1,1,j-M2),N,
     $                   one,d(1,1,j),N)
              CALL dgemm('n','n',N,N,N,alf,e2(1,1,j),N,f1(1,1,j-M2),N,
     $                   one,d(1,1,j),N)
              IF (j.LE.M2+M2)
     $          CALL dgemm('n','n',N,N,N,alf,e2(1,1,j),N,e2(1,1,j-M2),
     $                      N,one,d(1,1,j),N)
            END IF
            CALL dgemm('n','n',N,N,N,alf,e1(1,1,j),N,f(1,1,j-M1),N,
     $                one,d(1,1,j),N)
            CALL dgemm('n','n',N,N,N,alf,e1(1,1,j),N,f2(1,1,j-M1),N,
     $                one,d(1,1,j),N)
            IF (j.LE.M1+M2)
     $        CALL dgemm('n','n',N,N,N,alf,e1(1,1,j),N,e2(1,1,j-M1),
     $                N,one,d(1,1,j),N)
          END IF
          IF (j.LE.M2+1)
     $     CALL dgemm('n','n',N,N,N,alf,e(1,1,j),N,e2(1,1,j-1),N,
     $               one,d(1,1,j),N)
          IF (j-1+M2.LE.nblock)
     $     CALL dgemm('n','n',N,N,N,alf,e(1,1,j),N,f2(1,1,j-1),N,
     $               one,d(1,1,j),N)
          IF (j-1+M1.LE.nblock)
     $     CALL dgemm('n','n',N,N,N,alf,e(1,1,j),N,f1(1,1,j-1),N,
     $               one,d(1,1,j),N)
        END IF
C
C         Compute the LU-decomposition of d(j):
C
        CALL DGETRF( N, N, d(1,1,j), N, pivot(1,j), INFO )
        CALL DGETRS('n',N,N,d(1,1,j),N,pivot(1,j),e2(1,1,j),N,INFO)
        IF (j.GT.M1) 
     $    CALL DGETRS('n',N,N,d(1,1,j),N,pivot(1,j),e1(1,1,j),N,INFO)
        CALL DGETRS('n',N,N,d(1,1,j),N,pivot(1,j),e(1,1,j),N,INFO)
        CALL DGETRS('n',N,N,d(1,1,j),N,pivot(1,j),dm(1,1,j),N,INFO)
        DO i=1,N
          dm(i,i,j)=dm(i,i,j)-2.D0
        ENDDO
        IF (j.LT.nblock)
     $    CALL DGETRS('n',N,N,d(1,1,j),N,pivot(1,j),f(1,1,j),N,INFO)
        IF (j+M1.LE.nblock)
     $    CALL DGETRS('n',N,N,d(1,1,j),N,pivot(1,j),f1(1,1,j),N,INFO)
        CALL DGETRS('n',N,N,d(1,1,j),N,pivot(1,j),f2(1,1,j),N,INFO)
      ENDDO
      RETURN
      END

C=============================================================================
      subroutine heptasoleis(d,e,f,e1,f1,e2,f2,x,y,v,nblock,N,M1,M2)
C
C     Auke van der Ploeg, Februari 7 1997.
C
C     This routine applies the preconditonering which  
C     must have been constructed in subroutine 'heptapreceis.' 
*          Contains the matrix-vector multiplication
*
*               x := L^{-1}AU^{-1}y
*  
*          implemented according to the efficient
*          Eisenstat implementation.
*
*          See, for example, page 19 of the Phd thesis:
*          'Preconditioning for sparse matrices with applications.'
*          1994, University of Groningen.
*
C
      IMPLICIT NONE
      INTEGER N,M1,M2,nblock
      REAL*8 e(N,N,nblock),d(N,N,nblock),e1(N,N,nblock),f1(N,N,nblock),
     *           f(N,N,nblock),x(N,nblock),y(N,nblock),v(N,nblock),
     *           e2(N,N,nblock),f2(N,N,nblock)
C
C=======================================================================
C         Description of arguments:
C=======================================================================
C
C     d, e, f, e1, f1, e2, f2,
C               The matrix L+U is assumed to be block heptadiagonal.
C               In addition, the matrix contains blocks coming from
C               the treatment of periodic boundary conditions, M2 blocks
C               in the Lower-Left corner, and M2 blocks in the 
C               Upper-Right corner. The distance from these blocks
C               to the blocks on the main diagonal is nblock-M2.
C
C               The blocks are stored as follows:
C               d(j): j=1..nblock        contains D-2I 
C                                        in which D contains the blocks
C                                        on the main diagonal of A.
C               e(j): j=2..nblock     sub diagonal blocks.
C               f(j): j=1..nblock-1   super diagonal blocks.
C               e1(j): j=M1+1..nblock     blocks in the lower-triangular
C                                        part with distance M1 from 
C                                        the main diagonal.
C               f1(j): j=1..nblock-M1     blocks in the upper-triangular
C                                        part with distance M1 from 
C                                        the main diagonal.
C               e2(j): j=M2+1..nblock     blocks in the lower-triangular
C                                        part with distance M2 from 
C                                        the main diagonal.
C               f2(j): j=1..nblock-M2     blocks in the upper-triangular
C                                        part with distance M2 from 
C                                        the main diagonal.
C               e2(j): j=1..M2            blocks coming from periodic
C                                        boundary conditions in the 
C                                        upper-right corner, distance
C                                        nblock-M2 from blocks on the
C                                        main diagonal.
C               f2(j): j=nblock-M2+1..nblock:
C                                        blocks coming from periodic
C                                        boundary conditions in the 
C                                        lower-left corner, distance
C                                        nblock-M2 from blocks on the
C                                        main diagonal.
C
C              It is assumed that the
C              blocks are not very sparse, so the sparsity pattern
C              within the separate blocks is not exploited.
C              For example, the (i,k)-element of the j-th block on the 
C              main diagonal is stored in mat(i,k,j,3).
C
C     x:       On output, the solution is stored in x.
C     y:       right-hand side, unchanged on exit.
C     v:       scratch variable
C              The vector x, y and v are partitioned according to
C              the blockstructure of the matrices. For example, the 
C              (j-1)*N+i-th component is stored in x(i,j), for
C              i=0,1,...N-1 and j=1,2,...,nblock. 
C
C     nblock:  Number of diagonal blocks.
C     N:        the size of the blocks.
C     M1:       distance of blocks to the main diagonal blocks.
C     M2:       distance of outer-most blocks to main diagonal blocks.
C               (except those coming from periodic boundary conditions)
C               1 < M1 < M2.
C               The matrix has a blockstructure corresponding with
C               a seven-point stencil on a three-dimensional, 
C               rectangular grid. The blocks corresonding with the
C               direction in which grid points are numbered
C               first, are the sub- and superdiagonal blocks.
C               The blocks corresponding with the
C               direction in which grid points are numbered
C               secondly, have distance M1 from the main diagonal
C               blocks. Finally, The blocks corresponding with the
C               direction in which grid points are numbered
C               last, have distance M2 from the main diagonal blocks.
C
C=======================================================================
C         Local variables:
C=======================================================================
C
      INTEGER j, i, k, Ng
      REAL*8 one, vj0, vj1, vj2, vj3
      PARAMETER (one=1.D0, Ng=3)
C
C     Ng   :   for large values of N, each matrix-vector with a
C              subblock is implemented with a  BLAS-routine. For small 
C              values of N (N<Ng) this gives too much overhead, so it 
C              is better to write out the loops, which enables
C              more optimalizations by the compiler.
C              The user is encouraged to experiment with several
C              values for Ng, in order to determine an optimal value.
C     j,i  :   integers necessary for implementing loops.
C     INFO :   (output-)variable used for error indication by DGETRS.
C     one, :   variable necessary for BLAS-routine DGEMV. 
C
C=======================================================================
C         External function
C=======================================================================
C
C     DGEMV,   BLAS level two Matrix-Vector Product.
C              See 'man DGEMV' for description.
C
      IF (N.LT.Ng) GO TO 1000
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       Code for N >= Ng:
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C      v:=(I-U)^{-1}y
C
      CALL DCOPY( N*nblock, y(1,1), 1, v(1,1), 1)
      DO j=nblock-1,1,-1
        CALL DGEMV('n',N,N,-one,f(1,1,j),N,v(1,j+1),1,one,v(1,j),1)
        IF (j+M2.LE.nblock) 
     $   CALL DGEMV('n',N,N,-one,f2(1,1,j),N,v(1,j+M2),1,one,v(1,j),1)
        IF (j+M1.LE.nblock) 
     $   CALL DGEMV('n',N,N,-one,f1(1,1,j),N,v(1,j+M1),1,one,v(1,j),1)
C
C        Deal with blocks in the upper-right corner coming from
C        periodic boundary conditions:
C
        IF (j.LE.M2) 
     $   CALL DGEMV('n',N,N,-one,e2(1,1,j),N,v(1,j+nblock-M2),1,
     $               one,v(1,j),1)
      ENDDO
C
C      x:=y+(D-2I)v, D-2I is stored in d.
C
      CALL DCOPY( N*nblock, y(1,1), 1, x(1,1), 1)
      DO j=1,nblock
        CALL DGEMV('n',N,N,one,d(1,1,j),N,v(1,j),1,one,x(1,j),1)
      ENDDO
C
C        solve (I-L)y=x, put the solution again in x:
C
      DO j=2,nblock
        IF (j.GT.M2)
     $     CALL DGEMV('n',N,N,-one,e2(1,1,j),N,x(1,j-M2),1,one,x(1,j),1)
        IF (j.GT.M1)
     $     CALL DGEMV('n',N,N,-one,e1(1,1,j),N,x(1,j-M1),1,one,x(1,j),1)
        CALL DGEMV('n',N,N,-one,e(1,1,j),N,x(1,j-1),1,one,x(1,j),1)
C
C        Deal with blocks in the lower-left corner coming from
C        periodic boundary conditions:
C
        IF (j+M2.GT.nblock) 
     $     CALL DGEMV('n',N,N,-one,f2(1,1,j),N,x(1,j-nblock+M2),1,
     $                one,x(1,j),1)
      ENDDO
      CALL daxpy( nblock*N, one, v(1,1), 1, x(1,1), 1 )
      RETURN
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       Code for 1 < N < Ng:
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      v:=(I-U)^{-1}y
C
 1000 IF (N.EQ.1) GO TO 2000
C
      CALL DCOPY( N*nblock, y(1,1), 1, v(1,1), 1)
      DO j=nblock-1,nblock-M1+1,-1
        DO k=1,N
          vj2 = v(k,j+1)
          DO i=1,N
            v(i,j) = v(i,j) - f(i,k,j)*vj2
          ENDDO
        ENDDO
      ENDDO
      DO j=nblock-M1,nblock-M2+1,-1
        DO k=1,N
          vj2 = v(k,j+1)
          vj1 = v(k,j+M1)
          DO i=1,N
            v(i,j) = v(i,j) - f1(i,k,j)*vj1 - 
     &                     f(i,k,j)*vj2
          ENDDO
        ENDDO
      ENDDO
      DO j=nblock-M2,M2+1,-1
        DO k=1,N
          vj2 = v(k,j+1)
          vj0 = v(k,j+M2)
          vj1 = v(k,j+M1)
          DO i=1,N
            v(i,j) = v(i,j) - f2(i,k,j)*vj0 - f1(i,k,j)*vj1 - 
     &                   f(i,k,j)*vj2
          ENDDO
        ENDDO
      ENDDO
      DO j=M2,1,-1
C
C        Deal with blocks in the upper-right corner coming from
C        periodic boundary conditions:
C
        DO k=1,N
          vj2 = v(k,j+1)
          vj0 = v(k,j+M2)
          vj1 = v(k,j+M1)
          vj3 = v(k,j+Nblock-M2)
          DO i=1,N
            v(i,j) = v(i,j) - f2(i,k,j)*vj0 - f1(i,k,j)*vj1 -
     &          f(i,k,j)*vj2 - e2(i,k,j)*vj3
          ENDDO
        ENDDO
      ENDDO
C
C      x := y + (D-2I)v, D-2I is stored in d.
C
      CALL DCOPY( N*nblock, y(1,1), 1, x(1,1), 1)
      DO j=1,nblock
        DO k=1,N
          vj1 = v(k,j)
          DO i=1,N
            x(i,j) = x(i,j) + d(i,k,j)*vj1
          ENDDO
        ENDDO
      ENDDO
C
C        solve (I-L)z = x, put the solution again in x:
C
      DO j=2,M1
        DO k=1,N
          vj2 = x(k,j-1)
          DO i=1,N
            x(i,j) = x(i,j) - e(i,k,j)*vj2
          ENDDO
        ENDDO
      ENDDO
      DO j=M1+1,M2
        DO k=1,N
          vj2 = x(k,j-1)
          vj1 = x(k,j-M1)
          DO i=1,N
            x(i,j) = x(i,j) - e1(i,k,j)*vj1 - e(i,k,j)*vj2
          ENDDO
        ENDDO
      ENDDO
      DO j=M2+1,nblock-M2
        DO k=1,N
          vj2 = x(k,j-1)
          vj0 = x(k,j-M2)
          vj1 = x(k,j-M1)
          DO i=1,N
             x(i,j) = x(i,j) - e1(i,k,j)*vj1 - e2(i,k,j)*vj0 - 
     &                  e(i,k,j)*vj2 
          ENDDO
        ENDDO
      ENDDO
      DO j=nblock-M2+1,nblock
C
C        Deal with blocks in the lower-left corner coming from
C        periodic boundary conditions:
C
        DO k=1,N
          vj2 = x(k,j-1)
          vj0 = x(k,j-M2)
          vj1 = x(k,j-M1)
          vj3 = x(k,j-nblock+M2)
          DO i=1,N
            x(i,j) = x(i,j) - e1(i,k,j)*vj1 - e2(i,k,j)*vj0 - 
     &                e(i,k,j)*vj2 - f2(i,k,j)*vj3
          ENDDO
        ENDDO
      ENDDO
      CALL daxpy( nblock*N, one, v(1,1), 1, x(1,1), 1 )
      RETURN
C
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       Code for N=1:
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      v:=(I-U)^{-1}y
C
 2000 v(1,nblock)=y(1,nblock)
      DO j=nblock-1,nblock-M1+1,-1
        v(1,j) = y(1,j) - f(1,1,j)*v(1,j+1) 
      ENDDO
      DO j=nblock-M1,nblock-M2+1,-1
        v(1,j) = y(1,j) - f1(1,1,j)*v(1,j+M1) -
     &                   f(1,1,j)*v(1,j+1)
      ENDDO
      DO j=nblock-M2,M2+1,-1
        v(1,j) = y(1,j) - f1(1,1,j)*v(1,j+M1) - f(1,1,j)*v(1,j+1)-
     &                        f2(1,1,j)*v(1,j+M2)
      ENDDO
      DO j=M2,1,-1
C
C        Deal with blocks in the upper-right corner coming from
C        periodic boundary conditions:
C
        v(1,j) = y(1,j) - f2(1,1,j)*v(1,j+M2) -
     &           f1(1,1,j)*v(1,j+M1) -
     &           f(1,1,j)*v(1,j+1) - e2(1,1,j)*v(1,j+nblock-M2)
      ENDDO
C
C      x := y + (D-2I)v, D-2I is stored in d.
C
      DO j=1,nblock
        x(1,j) = y(1,j) + d(1,1,j)*v(1,j)
      ENDDO
C
C        solve (I-L)z = x, put the solution again in x:
C
      DO j=2,M1
        x(1,j) = x(1,j) - e(1,1,j)*x(1,j-1) 
      ENDDO
      DO j=M1+1,M2
        x(1,j) = x(1,j) - e1(1,1,j)*x(1,j-M1) - e(1,1,j)*x(1,j-1) 
      ENDDO
      DO j=M2+1,nblock-M2
         x(1,j) = x(1,j) - e1(1,1,j)*x(1,j-M1) - e(1,1,j)*x(1,j-1)-
     &                     e2(1,1,j)*x(1,j-M2)  
      ENDDO
      DO j=nblock-M2+1,nblock
C
C        Deal with blocks in the lower-left corner coming from
C        periodic boundary conditions:
C
         x(1,j) = x(1,j) - e2(1,1,j)*x(1,j-M2) -
     &                     e1(1,1,j)*x(1,j-M1) - 
     &            e(1,1,j)*x(1,j-1) - f2(1,1,j)*x(1,j-nblock+M2)
      ENDDO
      CALL daxpy( nblock, one, v(1,1), 1, x(1,1), 1 )
      RETURN
      END

C=============================================================================
      SUBROUTINE posthepta(f,f1,e2,f2,nblock,N,M1,M2,x)
C
C     Auke van der Ploeg, Februari 7 1997
C
C     Postprocessing heptadiagonal solver: x:=U^{-1}x
C
      IMPLICIT NONE
      INTEGER N,M1,M2,nblock
      REAL*8 f(N,N,nblock),e2(N,N,nblock),f1(N,N,nblock),
     &       f2(N,N,nblock),x(N*nblock)
C
C=======================================================================
C         Description of arguments:
C=======================================================================
C
C     f, f1, e2, f2:
C               The matrix U is the upper part of a block heptadiagonal
C               with a sparsity pattern correponding to a seven-point
C               stencil.
C               In addition, the matrix contains blocks coming from
C               the treatment of periodic boundary conditions:
C               M2 blocks in the 
C               Upper-Right corner. The distance from these blocks
C               to the blocks on the main diagonal is nblock-M2.
C
C               The blocks are stored as follows:
C               f(j): j=1..nblock-1   super diagonal blocks.
C               f1(j): j=1..nblock-M1     blocks in the upper-triangular
C                                        part with distance M1 from 
C                                        the main diagonal.
C               f2(j): j=1..nblock-M2     blocks in the upper-triangular
C                                        part with distance M2 from 
C                                        the main diagonal.
C               e2(j): j=1..M2            blocks coming from periodic
C                                        boundary conditions in the 
C                                        upper-right corner, distance
C                                        nblock-M2 from blocks on the
C                                        main diagonal.
C
C              It is assumed that the
C              blocks are not very sparse, so the sparsity pattern
C              within the separate blocks is not exploited.
C              For example, the (i,k)-element of the j-th block on the 
C              main diagonal is stored in mat(i,k,j,3).
C     nblock:  Number of diagonal blocks.
C     N:        the size of the blocks.
C     M1:       distance of blocks to the main diagonal blocks.
C     M2:       distance of outer-most blocks to main diagonal blocks.
C               (except those coming from periodic boundary conditions)
C               1 < M1 < M2.
C               The matrix has a blockstructure corresponding with
C               a seven-point stencil on a three-dimensional, 
C               rectangular grid. The blocks corresonding with the
C               direction in which grid points are numbered
C               first, are the sub- and superdiagonal blocks.
C               The blocks corresponding with the
C               direction in which grid points are numbered
C               secondly, have distance M1 from the main diagonal
C               blocks. Finally, The blocks corresponding with the
C               direction in which grid points are numbered
C               last, have distance M2 from the main diagonal blocks.
C     x:       on input:  solution of preconditioned system,
C              on output: solution of the original system.
C
C=======================================================================
C         Local variables:
C=======================================================================
C
      INTEGER i
      REAL*8 one
      PARAMETER (one=1.D0)

C     i    :   integer necessary for implementing loops.
C     one  :   variable necessary for BLAS-routine DGEMV. 
C
C=======================================================================

      DO i=nblock-1,1,-1
         CALL DGEMV('n',N,N,-one,f(1,1,i),N,
     $              x(1+i*N),1,one,x(1+(i-1)*N),1)
         IF (i+M1.LE.nblock) 
     $       CALL DGEMV('n',N,N,-one,f1(1,1,i),N,
     $                  x(1+(i+M1-1)*N),1,one,x(1+(i-1)*N),1)
         IF (i+M2.LE.nblock) 
     $       CALL DGEMV('n',N,N,-one,f2(1,1,i),N,
     $                  x(1+(i+M2-1)*N),1,one,x(1+(i-1)*N),1)
         IF (i.LE.M2) 
     $       CALL DGEMV('n',N,N,-one,e2(1,1,i),N,
     $                  x(1+(i+nblock-M2-1)*N),1,one,x(1+(i-1)*N),1)
      ENDDO

      RETURN
      END

C=============================================================================
C******************* Common subroutines **************************************
C=============================================================================
      SUBROUTINE solveb( N, A, IPIV, B)
*
*  -- modified LAPACK routine --
*     June 4, 1996
*
*     .. Scalar Arguments ..
      INTEGER            N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV(N)
      REAL*8         A(N,N), B(N)
*     ..
*
*  Purpose
*  =======
*
*  solveb solves a system of linear equations
*     A * X = B,  
*  with a general N-by-N matrix A using the LU factorization computed
*  by DGETRF.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input) REAL*8 array, dimension (N,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input/output) REAL*8 array, dimension (N)
*          On entry, the right hand side vector B.
*          On exit, the solution vector X.
*
*  =====================================================================
*
*     ..
*     .. External Subroutines ..
      EXTERNAL           DTRSV
*
*        Local variables
      INTEGER i
      REAL*8 ff
*     ..
*     .. Executable Statements ..
*
*
*        Apply row interchanges to the right hand sides.
*
      DO I=1,N
        ff=B(IPIV(I))
        B(IPIV(i))=B(i)
        B(i)=ff
      ENDDO
*
*        Solve L*X = B, overwriting B with X.
*
      CALL DTRSV('Lower','No transpose','Unit',N,A,N,B,1) 
*
*        Solve U*X = B, overwriting B with X.
*
      CALL DTRSV('Upper','No transpose','Non-unit',N,A,N,B,1) 
      RETURN
      END      
