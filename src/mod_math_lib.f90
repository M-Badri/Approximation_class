module mod_math_lib
!***************************************************************************************************
!> This module contains a bunch of low level procedures for mathematical calculations.
!> The procedures are mainly obtained from LINPACK, LAPACK and BLAS libraries (but not limited to those)
!> Some of the procedures have been modified to be more suitable for present applications.
!> Wherever is possible (and needed) the procedures are modified to be THREAD SAFE.
!
!***************************************************************************************************


   contains


   subroutine daxpy( n, da, dx, incx, dy, incy )
   !*****************************************************************************80
   !
   !! DAXPY computes constant times a vector plus a vector.
   !
   !  Discussion:
   !
   !    This routine uses double precision real arithmetic.
   !
   !    This routine uses unrolled loops for increments equal to one.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    16 May 2005
   !
   !  Author:
   !
   !    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
   !    David Kincaid, Fred Krogh.
   !    FORTRAN90 version by John Burkardt.
   !
   !  Reference:
   !
   !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
   !    LINPACK User's Guide,
   !    SIAM, 1979,
   !    ISBN13: 978-0-898711-72-1,
   !    LC: QA214.L56.
   !
   !    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
   !    Algorithm 539,
   !    Basic Linear Algebra Subprograms for Fortran Usage,
   !    ACM Transactions on Mathematical Software,
   !    Volume 5, Number 3, September 1979, pages 308-323.
   !
   !  Parameters:
   !
   !    Input, integer N, the number of elements in DX and DY.
   !
   !    Input, real ( kind = rk ) DA, the multiplier of DX.
   !
   !    Input, real ( kind = rk ) DX(*), the first vector.
   !
   !    Input, integer INCX, the increment between successive
   !    entries of DX.
   !
   !    Input/output, real ( kind = rk ) DY(*), the second vector.
   !    On output, DY(*) has been replaced by DY(*) + DA * DX(*).
   !
   !    Input, integer INCY, the increment between successive
   !    entries of DY.
   !
      implicit none

      integer, parameter :: rk = kind ( 1.0D+00 )

      real( kind = rk ) :: da
      real( kind = rk ) :: dx(*)
      real( kind = rk ) :: dy(*)
      integer :: i
      integer :: incx
      integer :: incy
      integer :: ix
      integer :: iy
      integer :: m
      integer :: n

      if ( n <= 0 ) then
         return
      end if

      if ( da == 0.0D+00 ) then
         return
      end if
      !
      !  Code for unequal increments or equal increments
      !  not equal to 1.
      !
      if ( incx /= 1 .or. incy /= 1 ) then

         if ( 0 <= incx ) then
            ix = 1
         else
            ix = ( - n + 1 ) * incx + 1
         end if

         if ( 0 <= incy ) then
            iy = 1
         else
            iy = ( - n + 1 ) * incy + 1
         end if

         do i = 1, n
            dy(iy) = dy(iy) + da * dx(ix)
            ix = ix + incx
            iy = iy + incy
         end do
      !
      !  Code for both increments equal to 1.
      !
      else

         m = mod ( n, 4 )

         dy(1:m) = dy(1:m) + da * dx(1:m)

         do i = m+1, n, 4
            dy(i  ) = dy(i  ) + da * dx(i  )
            dy(i+1) = dy(i+1) + da * dx(i+1)
            dy(i+2) = dy(i+2) + da * dx(i+2)
            dy(i+3) = dy(i+3) + da * dx(i+3)
         end do

      end if
   end subroutine daxpy



   function ddot( n, dx, incx, dy, incy )
   !*****************************************************************************80
   !
   !! DDOT forms the dot product of two vectors.
   !
   !  Discussion:
   !
   !    This routine uses double precision real arithmetic.
   !
   !    This routine uses unrolled loops for increments equal to one.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    16 May 2005
   !
   !  Author:
   !
   !    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
   !    David Kincaid, Fred Krogh.
   !    FORTRAN90 version by John Burkardt.
   !
   !  Reference:
   !
   !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
   !    LINPACK User's Guide,
   !    SIAM, 1979,
   !    ISBN13: 978-0-898711-72-1,
   !    LC: QA214.L56.
   !
   !    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
   !    Algorithm 539,
   !    Basic Linear Algebra Subprograms for Fortran Usage,
   !    ACM Transactions on Mathematical Software,
   !    Volume 5, Number 3, September 1979, pages 308-323.
   !
   !  Parameters:
   !
   !    Input, integer N, the number of entries in the vectors.
   !
   !    Input, real ( kind = rk ) DX(*), the first vector.
   !
   !    Input, integer INCX, the increment between successive
   !    entries in DX.
   !
   !    Input, real ( kind = rk ) DY(*), the second vector.
   !
   !    Input, integer INCY, the increment between successive
   !    entries in DY.
   !
   !    Output, real ( kind = rk ) DDOT, the sum of the product of the
   !    corresponding entries of DX and DY.
   !
     implicit none

     integer, parameter :: rk = kind ( 1.0D+00 )

     real ( kind = rk ) ddot
     real ( kind = rk ) dtemp
     real ( kind = rk ) dx(*)
     real ( kind = rk ) dy(*)
     integer i
     integer incx
     integer incy
     integer ix
     integer iy
     integer m
     integer n

     ddot = 0.0D+00
     dtemp = 0.0D+00

     if ( n <= 0 ) then
       return
     end if
   !
   !  Code for unequal increments or equal increments
   !  not equal to 1.
   !
     if ( incx /= 1 .or. incy /= 1 ) then

       if ( 0 <= incx ) then
         ix = 1
       else
         ix = ( - n + 1 ) * incx + 1
       end if

       if ( 0 <= incy ) then
         iy = 1
       else
         iy = ( - n + 1 ) * incy + 1
       end if

       do i = 1, n
         dtemp = dtemp + dx(ix) * dy(iy)
         ix = ix + incx
         iy = iy + incy
       end do
   !
   !  Code for both increments equal to 1.
   !
     else

       m = mod ( n, 5 )

       do i = 1, m
         dtemp = dtemp + dx(i) * dy(i)
       end do

       do i = m+1, n, 5

         dtemp = dtemp + dx(i  ) * dy(i  ) &
                       + dx(i+1) * dy(i+1) &
                       + dx(i+2) * dy(i+2) &
                       + dx(i+3) * dy(i+3) &
                       + dx(i+4) * dy(i+4)
       end do

     end if

     ddot = dtemp

   end function ddot



   function dnrm2 ( n, x, incx )
   !*****************************************************************************80
   !
   !! DNRM2 returns the euclidean norm of a vector.
   !
   !  Discussion:
   !
   !    This routine uses double precision real arithmetic.
   !
   !     DNRM2 ( X ) = sqrt ( X' * X )
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    16 May 2005
   !
   !  Author:
   !
   !    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
   !    David Kincaid, Fred Krogh.
   !    FORTRAN90 version by John Burkardt.
   !
   !  Reference:
   !
   !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
   !    LINPACK User's Guide,
   !    SIAM, 1979,
   !    ISBN13: 978-0-898711-72-1,
   !    LC: QA214.L56.
   !
   !    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
   !    Algorithm 539,
   !    Basic Linear Algebra Subprograms for Fortran Usage,
   !    ACM Transactions on Mathematical Software,
   !    Volume 5, Number 3, September 1979, pages 308-323.
   !
   !  Parameters:
   !
   !    Input, integer N, the number of entries in the vector.
   !
   !    Input, real ( kind = rk ) X(*), the vector whose norm is to be computed.
   !
   !    Input, integer INCX, the increment between successive
   !    entries of X.
   !
   !    Output, real ( kind = rk ) DNRM2, the Euclidean norm of X.
   !
     implicit none

     integer, parameter :: rk = kind ( 1.0D+00 )

     real ( kind = rk ) absxi
     real ( kind = rk ) dnrm2
     integer incx
     integer ix
     integer n
     real ( kind = rk ) norm
     real ( kind = rk ) scale
     real ( kind = rk ) ssq
     real ( kind = rk ) x(*)

     if ( n < 1 .or. incx < 1 ) then

       norm  = 0.0D+00

     else if ( n == 1 ) then

       norm  = abs ( x(1) )

     else

       scale = 0.0D+00
       ssq = 1.0D+00

       do ix = 1, 1 + ( n - 1 ) * incx, incx
         if ( x(ix) /= 0.0D+00 ) then
           absxi = abs ( x(ix) )
           if ( scale < absxi ) then
             ssq = 1.0D+00 + ssq * ( scale / absxi )**2
             scale = absxi
           else
             ssq = ssq + ( absxi / scale )**2
           end if
         end if
       end do
       norm  = scale * sqrt ( ssq )
     end if

     dnrm2 = norm

   end function dnrm2



   subroutine dqrank ( a, lda, m, n, tol, kr, jpvt, qraux, work )
   !*****************************************************************************80
   !
   !! DQRANK computes the QR factorization of a rectangular matrix.
   !
   !  Discussion:
   !
   !    This routine is used in conjunction with sqrlss to solve
   !    overdetermined, underdetermined and singular linear systems
   !    in a least squares sense.
   !
   !    DQRANK uses the LINPACK subroutine DQRDC to compute the QR
   !    factorization, with column pivoting, of an M by N matrix A.
   !    The numerical rank is determined using the tolerance TOL.
   !
   !    Note that on output, ABS ( A(1,1) ) / ABS ( A(KR,KR) ) is an estimate
   !    of the condition number of the matrix of independent columns,
   !    and of R.  This estimate will be <= 1/TOL.
   !
   !  Reference:
   !
   !    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
   !    LINPACK User's Guide,
   !    SIAM, 1979,
   !    ISBN13: 978-0-898711-72-1,
   !    LC: QA214.L56.
   !
   !  Parameters:
   !
   !    Input/output, real ( kind = rk ) A(LDA,N).  On input, the matrix whose
   !    decomposition is to be computed.  On output, the information from DQRDC.
   !    The triangular matrix R of the QR factorization is contained in the
   !    upper triangle and information needed to recover the orthogonal
   !    matrix Q is stored below the diagonal in A and in the vector QRAUX.
   !
   !    Input, integer LDA, the leading dimension of A, which must
   !    be at least M.
   !
   !    Input, integer M, the number of rows of A.
   !
   !    Input, integer N, the number of columns of A.
   !
   !    Input, real ( kind = rk ) TOL, a relative tolerance used to determine the
   !    numerical rank.  The problem should be scaled so that all the elements
   !    of A have roughly the same absolute accuracy, EPS.  Then a reasonable
   !    value for TOL is roughly EPS divided by the magnitude of the largest
   !    element.
   !
   !    Output, integer KR, the numerical rank.
   !
   !    Output, integer JPVT(N), the pivot information from DQRDC.
   !    Columns JPVT(1), ..., JPVT(KR) of the original matrix are linearly
   !    independent to within the tolerance TOL and the remaining columns
   !    are linearly dependent.
   !
   !    Output, real ( kind = rk ) QRAUX(N), will contain extra information defining
   !    the QR factorization.
   !
   !    Workspace, real ( kind = rk ) WORK(N).
   !
     implicit none

     integer, parameter :: rk = kind ( 1.0D+00 )

     integer lda
     integer n

     real ( kind = rk ) a(lda,n)
     integer j
     integer job
     integer jpvt(n)
     integer k
     integer kr
     integer m
     real ( kind = rk ) qraux(n)
     real ( kind = rk ) tol
     real ( kind = rk ) work(n)

     jpvt(1:n) = 0
     job = 1

     call dqrdc ( a, lda, m, n, qraux, jpvt, work, job )

     kr = 0
     k = min ( m, n )

     do j = 1, k
       if ( abs ( a(j,j) ) <= tol * abs ( a(1,1) ) ) then
         return
       end if
       kr = j
     end do
   end subroutine dqrank



   subroutine dqrdc ( a, lda, n, p, qraux, jpvt, work, job )
   !*****************************************************************************80
   !
   !! DQRDC computes the QR factorization of a real rectangular matrix.
   !
   !  Discussion:
   !
   !    DQRDC uses Householder transformations.
   !
   !    Column pivoting based on the 2-norms of the reduced columns may be
   !    performed at the user's option.
   !
   !  Reference:
   !
   !    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
   !    LINPACK User's Guide,
   !    SIAM, 1979,
   !    ISBN13: 978-0-898711-72-1,
   !    LC: QA214.L56.
   !
   !  Parameters:
   !
   !    Input/output, real ( kind = rk ) A(LDA,P).  On input, the N by P matrix
   !    whose decomposition is to be computed.  On output, A contains in
   !    its upper triangle the upper triangular matrix R of the QR
   !    factorization.  Below its diagonal A contains information from
   !    which the orthogonal part of the decomposition can be recovered.
   !    Note that if pivoting has been requested, the decomposition is not that
   !    of the original matrix A but that of A with its columns permuted
   !    as described by JPVT.
   !
   !    Input, integer LDA, the leading dimension of the array A.
   !    LDA must be at least N.
   !
   !    Input, integer N, the number of rows of the matrix A.
   !
   !    Input, integer P, the number of columns of the matrix A.
   !
   !    Output, real ( kind = rk ) QRAUX(P), contains further information required
   !    to recover the orthogonal part of the decomposition.
   !
   !    Input/output, integer JPVT(P).  On input, JPVT contains
   !    integers that control the selection of the pivot columns.  The K-th
   !    column A(*,K) of A is placed in one of three classes according to the
   !    value of JPVT(K).
   !      > 0, then A(K) is an initial column.
   !      = 0, then A(K) is a free column.
   !      < 0, then A(K) is a final column.
   !    Before the decomposition is computed, initial columns are moved to
   !    the beginning of the array A and final columns to the end.  Both
   !    initial and final columns are frozen in place during the computation
   !    and only free columns are moved.  At the K-th stage of the
   !    reduction, if A(*,K) is occupied by a free column it is interchanged
   !    with the free column of largest reduced norm.  JPVT is not referenced
   !    if JOB == 0.  On output, JPVT(K) contains the index of the column of the
   !    original matrix that has been interchanged into the K-th column, if
   !    pivoting was requested.
   !
   !    Workspace, real ( kind = rk ) WORK(P).  WORK is not referenced if JOB == 0.
   !
   !    Input, integer JOB, initiates column pivoting.
   !    0, no pivoting is done.
   !    nonzero, pivoting is done.
   !
     implicit none

     integer, parameter :: rk = kind ( 1.0D+00 )

     integer lda
     integer n
     integer p

     real ( kind = rk ) a(lda,p)
     integer jpvt(p)
     real ( kind = rk ) qraux(p)
     real ( kind = rk ) work(p)
     integer j
     integer job
     integer jp
     integer l
     integer lup
     integer maxj
     real ( kind = rk ) maxnrm
     real ( kind = rk ) nrmxl
     integer pl
     integer pu
   !  real ( kind = rk ) ddot
   !  real ( kind = rk ) dnrm2
     logical swapj
     real ( kind = rk ) t
     real ( kind = rk ) tt

     pl = 1
     pu = 0
   !
   !  If pivoting is requested, rearrange the columns.
   !
     if ( job /= 0 ) then

       do j = 1, p

         swapj = 0 < jpvt(j)

         if ( jpvt(j) < 0 ) then
           jpvt(j) = - j
         else
           jpvt(j) = j
         end if

         if ( swapj ) then

           if ( j /= pl ) then
             call dswap ( n, a(1,pl), 1, a(1,j), 1 )
           end if

           jpvt(j) = jpvt(pl)
           jpvt(pl) = j
           pl = pl + 1

         end if

       end do

       pu = p

       do j = p, 1, -1

         if ( jpvt(j) < 0 ) then

           jpvt(j) = - jpvt(j)

           if ( j /= pu ) then
             call dswap ( n, a(1,pu), 1, a(1,j), 1 )
             jp = jpvt(pu)
             jpvt(pu) = jpvt(j)
             jpvt(j) = jp
           end if

           pu = pu - 1

         end if

       end do

     end if
   !
   !  Compute the norms of the free columns.
   !
     do j = pl, pu
       qraux(j) = dnrm2 ( n, a(1,j), 1 )
     end do

     work(pl:pu) = qraux(pl:pu)
   !
   !  Perform the Householder reduction of A.
   !
     lup = min ( n, p )

     do l = 1, lup
   !
   !  Bring the column of largest norm into the pivot position.
   !
       if ( pl <= l .and. l < pu ) then

         maxnrm = 0.0D+00
         maxj = l
         do j = l, pu
           if ( maxnrm < qraux(j) ) then
             maxnrm = qraux(j)
             maxj = j
           end if
         end do

         if ( maxj /= l ) then
           call dswap ( n, a(1,l), 1, a(1,maxj), 1 )
           qraux(maxj) = qraux(l)
           work(maxj) = work(l)
           jp = jpvt(maxj)
           jpvt(maxj) = jpvt(l)
           jpvt(l) = jp
         end if

       end if
   !
   !  Compute the Householder transformation for column L.
   !
       qraux(l) = 0.0D+00

       if ( l /= n ) then

         nrmxl = dnrm2 ( n-l+1, a(l,l), 1 )

         if ( nrmxl /= 0.0D+00 ) then

           if ( a(l,l) /= 0.0D+00 ) then
             nrmxl = sign ( nrmxl, a(l,l) )
           end if

           call dscal ( n-l+1, 1.0D+00 / nrmxl, a(l,l), 1 )
           a(l,l) = 1.0D+00 + a(l,l)
   !
   !  Apply the transformation to the remaining columns, updating the norms.
   !
           do j = l + 1, p

             t = - ddot ( n-l+1, a(l,l), 1, a(l,j), 1 ) / a(l,l)
             call daxpy ( n-l+1, t, a(l,l), 1, a(l,j), 1 )

             if ( pl <= j .and. j <= pu ) then

               if ( qraux(j) /= 0.0D+00 ) then

                 tt = 1.0D+00 - ( abs ( a(l,j) ) / qraux(j) )**2
                 tt = max ( tt, 0.0D+00 )
                 t = tt
                 tt = 1.0D+00 + 0.05D+00 * tt * ( qraux(j) / work(j) )**2

                 if ( tt /= 1.0D+00 ) then
                   qraux(j) = qraux(j) * sqrt ( t )
                 else
                   qraux(j) = dnrm2 ( n-l, a(l+1,j), 1 )
                   work(j) = qraux(j)
                 end if

               end if

             end if

           end do
   !
   !  Save the transformation.
   !
           qraux(l) = a(l,l)
           a(l,l) = - nrmxl

         end if

       end if

     end do
   end subroutine dqrdc



   subroutine dqrls ( a, lda, m, n, tol, kr, b, x, r, jpvt, qraux, work, &
     itask, ind )
   !*****************************************************************************80
   !
   !! DQRLS factors and solves a linear system in the least squares sense.
   !
   !  Discussion:
   !
   !    The linear system may be overdetermined, underdetermined or singular.
   !    The solution is obtained using a QR factorization of the
   !    coefficient matrix.
   !
   !    DQRLS can be efficiently used to solve several least squares
   !    problems with the same matrix A.  The first system is solved
   !    with ITASK = 1.  The subsequent systems are solved with
   !    ITASK = 2, to avoid the recomputation of the matrix factors.
   !    The parameters KR, JPVT, and QRAUX must not be modified
   !    between calls to DQRLS.
   !
   !    DQRLS is used to solve in a least squares sense
   !    overdetermined, underdetermined and singular linear systems.
   !    The system is A*X approximates B where A is M by N.
   !    B is a given M-vector, and X is the N-vector to be computed.
   !    A solution X is found which minimimzes the sum of squares (2-norm)
   !    of the residual,  A*X - B.
   !
   !    The numerical rank of A is determined using the tolerance TOL.
   !
   !    DQRLS uses the LINPACK subroutine DQRDC to compute the QR
   !    factorization, with column pivoting, of an M by N matrix A.
   !
   !  Modified:
   !
   !    15 April 2012
   !
   !  Reference:
   !
   !    David Kahaner, Cleve Moler, Steven Nash,
   !    Numerical Methods and Software,
   !    Prentice Hall, 1989,
   !    ISBN: 0-13-627258-4,
   !    LC: TA345.K34.
   !
   !  Parameters:
   !
   !    Input/output, real ( kind = rk ) A(LDA,N), an M by N matrix.
   !    On input, the matrix whose decomposition is to be computed.
   !    In a least squares data fitting problem, A(I,J) is the
   !    value of the J-th basis (model) function at the I-th data point.
   !    On output, A contains the output from DQRDC.  The triangular matrix R
   !    of the QR factorization is contained in the upper triangle and
   !    information needed to recover the orthogonal matrix Q is stored
   !    below the diagonal in A and in the vector QRAUX.
   !
   !    Input, integer LDA, the leading dimension of A.
   !
   !    Input, integer M, the number of rows of A.
   !
   !    Input, integer N, the number of columns of A.
   !
   !    Input, real ( kind = rk ) TOL, a relative tolerance used to determine the
   !    numerical rank.  The problem should be scaled so that all the elements
   !    of A have roughly the same absolute accuracy EPS.  Then a reasonable
   !    value for TOL is roughly EPS divided by the magnitude of the largest
   !    element.
   !
   !    Output, integer KR, the numerical rank.
   !
   !    Input, real ( kind = rk ) B(M), the right hand side of the linear system.
   !
   !    Output, real ( kind = rk ) X(N), a least squares solution to the linear
   !    system.
   !
   !    Output, real ( kind = rk ) R(M), the residual, B - A*X.  R may
   !    overwrite B.
   !
   !    Workspace, integer JPVT(N), required if ITASK = 1.
   !    Columns JPVT(1), ..., JPVT(KR) of the original matrix are linearly
   !    independent to within the tolerance TOL and the remaining columns
   !    are linearly dependent.  ABS ( A(1,1) ) / ABS ( A(KR,KR) ) is an estimate
   !    of the condition number of the matrix of independent columns,
   !    and of R.  This estimate will be <= 1/TOL.
   !
   !    Workspace, real ( kind = rk ) QRAUX(N), required if ITASK = 1.
   !
   !    Workspace, real ( kind = rk ) WORK(N), required if ITASK = 1.
   !
   !    Input, integer ITASK.
   !    1, DQRLS factors the matrix A and solves the least squares problem.
   !    2, DQRLS assumes that the matrix A was factored with an earlier
   !       call to DQRLS, and only solves the least squares problem.
   !
   !    Output, integer IND, error code.
   !    0:  no error
   !    -1: LDA < M   (fatal error)
   !    -2: N < 1     (fatal error)
   !    -3: ITASK < 1 (fatal error)
   !
     implicit none

     integer, parameter :: rk = kind ( 1.0D+00 )

     integer lda
     integer m
     integer n

     real ( kind = rk ) a(lda,n)
     real ( kind = rk ) b(m)
     integer ind
     integer itask
     integer jpvt(n)
     integer kr
     real ( kind = rk ) qraux(n)
     real ( kind = rk ) r(m)
     real ( kind = rk ) tol
     real ( kind = rk ) work(n)
     real ( kind = rk ) x(n)

     if ( lda < m ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'DQRLS - Fatal error!'
       write ( *, '(a)' ) '  LDA < M.'
       stop
     end if

     if ( n <= 0 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'DQRLS - Fatal error!'
       write ( *, '(a)' ) '  N <= 0.'
       stop
     end if

     if ( itask < 1 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'DQRLS - Fatal error!'
       write ( *, '(a)' ) '  ITASK < 1.'
       stop
     end if

     ind = 0
   !
   !  Factor the matrix.
   !
     if ( itask == 1 ) then
       call dqrank ( a, lda, m, n, tol, kr, jpvt, qraux, work )
     end if
   !
   !  Solve the least-squares problem.
   !
     call dqrlss ( a, lda, m, n, kr, b, x, r, jpvt, qraux )

   end subroutine dqrls



   subroutine dqrlss ( a, lda, m, n, kr, b, x, r, jpvt, qraux )
   !*****************************************************************************80
   !
   !! DQRLSS solves a linear system in a least squares sense.
   !
   !  Discussion:
   !
   !    DQRLSS must be preceeded by a call to DQRANK.
   !
   !    The system is to be solved is
   !      A * X = B
   !    where
   !      A is an M by N matrix with rank KR, as determined by DQRANK,
   !      B is a given M-vector,
   !      X is the N-vector to be computed.
   !
   !    A solution X, with at most KR nonzero components, is found which
   !    minimizes the 2-norm of the residual (A*X-B).
   !
   !    Once the matrix A has been formed, DQRANK should be
   !    called once to decompose it.  Then, for each right hand
   !    side B, DQRLSS should be called once to obtain the
   !    solution and residual.
   !
   !  Modified:
   !
   !    15 April 2012
   !
   !  Parameters:
   !
   !    Input, real ( kind = rk ) A(LDA,N), the QR factorization information
   !    from DQRANK.  The triangular matrix R of the QR factorization is
   !    contained in the upper triangle and information needed to recover
   !    the orthogonal matrix Q is stored below the diagonal in A and in
   !    the vector QRAUX.
   !
   !    Input, integer LDA, the leading dimension of A, which must
   !    be at least M.
   !
   !    Input, integer M, the number of rows of A.
   !
   !    Input, integer N, the number of columns of A.
   !
   !    Input, integer KR, the rank of the matrix, as estimated
   !    by DQRANK.
   !
   !    Input, real ( kind = rk ) B(M), the right hand side of the linear system.
   !
   !    Output, real ( kind = rk ) X(N), a least squares solution to the
   !    linear system.
   !
   !    Output, real ( kind = rk ) R(M), the residual, B - A*X.  R may
   !    overwite B.
   !
   !    Input, integer JPVT(N), the pivot information from DQRANK.
   !    Columns JPVT(1), ..., JPVT(KR) of the original matrix are linearly
   !    independent to within the tolerance TOL and the remaining columns
   !    are linearly dependent.
   !
   !    Input, real ( kind = rk ) QRAUX(N), auxiliary information from DQRANK
   !    defining the QR factorization.
   !
     implicit none

     integer, parameter :: rk = kind ( 1.0D+00 )

     integer lda
     integer m
     integer n

     real ( kind = rk ) a(lda,n)
     real ( kind = rk ) b(m)
     integer info
     integer j
     integer job
     integer jpvt(n)
     integer k
     integer kr
     real ( kind = rk ) qraux(n)
     real ( kind = rk ) r(m)
     real ( kind = rk ) t
     real ( kind = rk ) x(n)

     if ( kr /= 0 ) then
       job = 110
       call dqrsl ( a, lda, m, kr, qraux, b, r, r, x, r, r, job, info )
     end if

     jpvt(1:n) = - jpvt(1:n)

     x(kr+1:n) = 0.0D+00

     do j = 1, n

       if ( jpvt(j) <= 0 ) then

         k = -jpvt(j)
         jpvt(j) = k

         do while ( k /= j )
           t = x(j)
           x(j) = x(k)
           x(k) = t
           jpvt(k) = -jpvt(k)
           k = jpvt(k)
         end do
       end if
     end do
   end subroutine dqrlss



   subroutine dqrsl ( a, lda, n, k, qraux, y, qy, qty, b, rsd, ab, job, info )
   !*****************************************************************************80
   !
   !! DQRSL computes transformations, projections, and least squares solutions.
   !
   !  Discussion:
   !
   !    DQRSL requires the output of DQRDC.
   !
   !    For K <= min(N,P), let AK be the matrix
   !
   !      AK = ( A(JPVT(1)), A(JPVT(2)), ..., A(JPVT(K)) )
   !
   !    formed from columns JPVT(1), ..., JPVT(K) of the original
   !    N by P matrix A that was input to DQRDC.  If no pivoting was
   !    done, AK consists of the first K columns of A in their
   !    original order.  DQRDC produces a factored orthogonal matrix Q
   !    and an upper triangular matrix R such that
   !
   !      AK = Q * (R)
   !               (0)
   !
   !    This information is contained in coded form in the arrays
   !    A and QRAUX.
   !
   !    The parameters QY, QTY, B, RSD, and AB are not referenced
   !    if their computation is not requested and in this case
   !    can be replaced by dummy variables in the calling program.
   !    To save storage, the user may in some cases use the same
   !    array for different parameters in the calling sequence.  A
   !    frequently occuring example is when one wishes to compute
   !    any of B, RSD, or AB and does not need Y or QTY.  In this
   !    case one may identify Y, QTY, and one of B, RSD, or AB, while
   !    providing separate arrays for anything else that is to be
   !    computed.
   !
   !    Thus the calling sequence
   !
   !      call dqrsl ( a, lda, n, k, qraux, y, dum, y, b, y, dum, 110, info )
   !
   !    will result in the computation of B and RSD, with RSD
   !    overwriting Y.  More generally, each item in the following
   !    list contains groups of permissible identifications for
   !    a single calling sequence.
   !
   !      1. (Y,QTY,B) (RSD) (AB) (QY)
   !
   !      2. (Y,QTY,RSD) (B) (AB) (QY)
   !
   !      3. (Y,QTY,AB) (B) (RSD) (QY)
   !
   !      4. (Y,QY) (QTY,B) (RSD) (AB)
   !
   !      5. (Y,QY) (QTY,RSD) (B) (AB)
   !
   !      6. (Y,QY) (QTY,AB) (B) (RSD)
   !
   !    In any group the value returned in the array allocated to
   !    the group corresponds to the last member of the group.
   !
   !  Modified:
   !
   !    15 April 2012
   !
   !  Author:
   !
   !    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart
   !
   !  Reference:
   !
   !    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
   !    LINPACK User's Guide,
   !    SIAM, 1979,
   !    ISBN13: 978-0-898711-72-1,
   !    LC: QA214.L56.
   !
   !  Parameters:
   !
   !    Input, real ( kind = rk ) A(LDA,P), contains the output of DQRDC.
   !
   !    Input, integer LDA, the leading dimension of the array A.
   !
   !    Input, integer N, the number of rows of the matrix AK.  It
   !    must have the same value as N in DQRDC.
   !
   !    Input, integer K, the number of columns of the matrix AK.  K
   !    must not be greater than min(N,P), where P is the same as in the
   !    calling sequence to DQRDC.
   !
   !    Input, real ( kind = rk ) QRAUX(P), the auxiliary output from DQRDC.
   !
   !    Input, real ( kind = rk ) Y(N), a vector to be manipulated by DQRSL.
   !
   !    Output, real ( kind = rk ) QY(N), contains Q * Y, if requested.
   !
   !    Output, real ( kind = rk ) QTY(N), contains Q' * Y, if requested.
   !
   !    Output, real ( kind = rk ) B(K), the solution of the least squares problem
   !      minimize norm2 ( Y - AK * B),
   !    if its computation has been requested.  Note that if pivoting was
   !    requested in DQRDC, the J-th component of B will be associated with
   !    column JPVT(J) of the original matrix A that was input into DQRDC.
   !
   !    Output, real ( kind = rk ) RSD(N), the least squares residual Y - AK * B,
   !    if its computation has been requested.  RSD is also the orthogonal
   !    projection of Y onto the orthogonal complement of the column space
   !    of AK.
   !
   !    Output, real ( kind = rk ) AB(N), the least squares approximation Ak * B,
   !    if its computation has been requested.  AB is also the orthogonal
   !    projection of Y onto the column space of A.
   !
   !    Input, integer JOB, specifies what is to be computed.  JOB has
   !    the decimal expansion ABCDE, with the following meaning:
   !
   !      if A /= 0, compute QY.
   !      if B /= 0, compute QTY.
   !      if C /= 0, compute QTY and B.
   !      if D /= 0, compute QTY and RSD.
   !      if E /= 0, compute QTY and AB.
   !
   !    Note that a request to compute B, RSD, or AB automatically triggers
   !    the computation of QTY, for which an array must be provided in the
   !    calling sequence.
   !
   !    Output, integer INFO, is zero unless the computation of B has
   !    been requested and R is exactly singular.  In this case, INFO is the
   !    index of the first zero diagonal element of R, and B is left unaltered.
   !
     implicit none

     integer, parameter :: rk = kind ( 1.0D+00 )

     integer k
     integer lda
     integer n

     real ( kind = rk ) a(lda,*)
     real ( kind = rk ) ab(n)
     real ( kind = rk ) b(k)
     logical cab
     logical cb
     logical cqty
     logical cqy
     logical cr
     integer info
     integer j
     integer jj
     integer job
     integer ju
     integer kp1
     real ( kind = rk ) qraux(*)
     real ( kind = rk ) qty(n)
     real ( kind = rk ) qy(n)
     real ( kind = rk ) rsd(n)
   !  real ( kind = rk ) ddot
     real ( kind = rk ) t
     real ( kind = rk ) temp
     real ( kind = rk ) y(n)
   !
   !  set info flag.
   !
     info = 0
   !
   !  Determine what is to be computed.
   !
     cqy =        job / 10000         /= 0
     cqty = mod ( job,  10000 )       /= 0
     cb =   mod ( job,   1000 ) / 100 /= 0
     cr =   mod ( job,    100 ) /  10 /= 0
     cab =  mod ( job,     10 )       /= 0

     ju = min ( k, n-1 )
   !
   !  Special action when N = 1.
   !
     if ( ju == 0 ) then

       if ( cqy ) then
         qy(1) = y(1)
       end if

       if ( cqty ) then
         qty(1) = y(1)
       end if

       if ( cab ) then
         ab(1) = y(1)
        end if

       if ( cb ) then

         if ( a(1,1) == 0.0D+00 ) then
           info = 1
         else
           b(1) = y(1) / a(1,1)
         end if

       end if

       if ( cr ) then
         rsd(1) = 0.0D+00
       end if

       return

     end if
   !
   !  Set up to compute QY or QTY.
   !
     if ( cqy ) then
       qy(1:n) = y(1:n)
     end if

     if ( cqty ) then
       qty(1:n) = y(1:n)
     end if
   !
   !  Compute QY.
   !
     if ( cqy ) then

       do jj = 1, ju

         j = ju - jj + 1

         if ( qraux(j) /= 0.0D+00 ) then
           temp = a(j,j)
           a(j,j) = qraux(j)
           t = - ddot ( n-j+1, a(j,j), 1, qy(j), 1 ) / a(j,j)
           call daxpy ( n-j+1, t, a(j,j), 1, qy(j), 1 )
           a(j,j) = temp
         end if

       end do

     end if
   !
   !  Compute Q'*Y.
   !
        if ( cqty ) then

           do j = 1, ju
              if ( qraux(j) /= 0.0D+00 ) then
                 temp = a(j,j)
                 a(j,j) = qraux(j)
                 t = - ddot ( n-j+1, a(j,j), 1, qty(j), 1 ) / a(j,j)
                 call daxpy ( n-j+1, t, a(j,j), 1, qty(j), 1 )
                 a(j,j) = temp
              end if
           end do

        end if
   !
   !  Set up to compute B, RSD, or AB.
   !
        if ( cb ) then
          b(1:k) = qty(1:k)
        end if

        kp1 = k + 1

        if ( cab ) then
          ab(1:k) = qty(1:k)
        end if

        if ( cr .and. k < n ) then
          rsd(k+1:n) = qty(k+1:n)
        end if

        if ( cab .and. k+1 <= n ) then
           ab(k+1:n) = 0.0D+00
        end if

        if ( cr ) then
           rsd(1:k) = 0.0D+00
        end if
   !
   !  Compute B.
   !
        if ( cb ) then

           do jj = 1, k

              j = k - jj + 1

              if ( a(j,j) == 0.0D+00 ) then
                 info = j
                 exit
              end if

              b(j) = b(j)/a(j,j)

              if ( j /= 1 ) then
                 t = -b(j)
                 call daxpy ( j-1, t, a(1,j), 1, b, 1 )
              end if

           end do

        end if

        if ( cr .or. cab ) then
   !
   !  Compute RSD or AB as required.
   !
           do jj = 1, ju

              j = ju - jj + 1

              if ( qraux(j) /= 0.0D+00 ) then

                 temp = a(j,j)
                 a(j,j) = qraux(j)

                 if ( cr ) then
                    t = - ddot ( n-j+1, a(j,j), 1, rsd(j), 1 ) / a(j,j)
                    call daxpy ( n-j+1, t, a(j,j), 1, rsd(j), 1 )
                 end if

                 if ( cab ) then
                    t = - ddot ( n-j+1, a(j,j), 1, ab(j), 1 ) / a(j,j)
                    call daxpy ( n-j+1, t, a(j,j), 1, ab(j), 1 )
                 end if

                 a(j,j) = temp

              end if

           end do

     end if
   end subroutine dqrsl



   subroutine drot ( n, x, incx, y, incy, c, s )
   !*****************************************************************************80
   !
   !! DROT applies a plane rotation.
   !
   !  Discussion:
   !
   !    This routine uses double precision real arithmetic.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    08 April 1999
   !
   !  Author:
   !
   !    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
   !    David Kincaid, Fred Krogh.
   !    FORTRAN90 version by John Burkardt.
   !
   !  Reference:
   !
   !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
   !    LINPACK User's Guide,
   !    SIAM, 1979,
   !    ISBN13: 978-0-898711-72-1,
   !    LC: QA214.L56.
   !
   !    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
   !    Algorithm 539,
   !    Basic Linear Algebra Subprograms for Fortran Usage,
   !    ACM Transactions on Mathematical Software,
   !    Volume 5, Number 3, September 1979, pages 308-323.
   !
   !  Parameters:
   !
   !    Input, integer N, the number of entries in the vectors.
   !
   !    Input/output, real ( kind = rk ) X(*), one of the vectors to be rotated.
   !
   !    Input, integer INCX, the increment between successive
   !    entries of X.
   !
   !    Input/output, real ( kind = rk ) Y(*), one of the vectors to be rotated.
   !
   !    Input, integer INCY, the increment between successive
   !    elements of Y.
   !
   !    Input, real ( kind = rk ) C, S, parameters (presumably the cosine and
   !    sine of some angle) that define a plane rotation.
   !
     implicit none

     integer, parameter :: rk = kind ( 1.0D+00 )

     real ( kind = rk ) c
     integer i
     integer incx
     integer incy
     integer ix
     integer iy
     integer n
     real ( kind = rk ) s
     real ( kind = rk ) stemp
     real ( kind = rk ) x(*)
     real ( kind = rk ) y(*)

     if ( n <= 0 ) then

     else if ( incx == 1 .and. incy == 1 ) then

       do i = 1, n
         stemp = c * x(i) + s * y(i)
         y(i) = c * y(i) - s * x(i)
         x(i) = stemp
       end do

     else

       if ( 0 <= incx ) then
         ix = 1
       else
         ix = ( - n + 1 ) * incx + 1
       end if

       if ( 0 <= incy ) then
         iy = 1
       else
         iy = ( - n + 1 ) * incy + 1
       end if

       do i = 1, n
         stemp = c * x(ix) + s * y(iy)
         y(iy) = c * y(iy) - s * x(ix)
         x(ix) = stemp
         ix = ix + incx
         iy = iy + incy
       end do

     end if
   end  subroutine drot



   subroutine drotg ( sa, sb, c, s )
   !*****************************************************************************80
   !
   !! DROTG constructs a Givens plane rotation.
   !
   !  Discussion:
   !
   !    This routine uses double precision real arithmetic.
   !
   !    Given values A and B, this routine computes
   !
   !    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
   !          = sign ( B ) if abs ( A ) <= abs ( B );
   !
   !    R     = SIGMA * ( A * A + B * B );
   !
   !    C = A / R if R is not 0
   !      = 1     if R is 0;
   !
   !    S = B / R if R is not 0,
   !        0     if R is 0.
   !
   !    The computed numbers then satisfy the equation
   !
   !    (  C  S ) ( A ) = ( R )
   !    ( -S  C ) ( B ) = ( 0 )
   !
   !    The routine also computes
   !
   !    Z = S     if abs ( A ) > abs ( B ),
   !      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
   !      = 1     if C is 0.
   !
   !    The single value Z encodes C and S, and hence the rotation:
   !
   !    If Z = 1, set C = 0 and S = 1;
   !    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
   !    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    15 May 2006
   !
   !  Author:
   !
   !    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
   !    David Kincaid, Fred Krogh.
   !    FORTRAN90 version by John Burkardt.
   !
   !  Reference:
   !
   !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
   !    LINPACK User's Guide,
   !    SIAM, 1979,
   !    ISBN13: 978-0-898711-72-1,
   !    LC: QA214.L56.
   !
   !    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
   !    Algorithm 539,
   !    Basic Linear Algebra Subprograms for Fortran Usage,
   !    ACM Transactions on Mathematical Software,
   !    Volume 5, Number 3, September 1979, pages 308-323.
   !
   !  Parameters:
   !
   !    Input/output, real ( kind = rk ) SA, SB.  On input, SA and SB are the values
   !    A and B.  On output, SA is overwritten with R, and SB is
   !    overwritten with Z.
   !
   !    Output, real ( kind = rk ) C, S, the cosine and sine of the
   !    Givens rotation.
   !
     implicit none

     integer, parameter :: rk = kind ( 1.0D+00 )

     real ( kind = rk ) c
     real ( kind = rk ) r
     real ( kind = rk ) roe
     real ( kind = rk ) s
     real ( kind = rk ) sa
     real ( kind = rk ) sb
     real ( kind = rk ) scale
     real ( kind = rk ) z

     if ( abs ( sb ) < abs ( sa ) ) then
       roe = sa
     else
       roe = sb
     end if

     scale = abs ( sa ) + abs ( sb )

     if ( scale == 0.0D+00 ) then
       c = 1.0D+00
       s = 0.0D+00
       r = 0.0D+00
     else
       r = scale * sqrt ( ( sa / scale )**2 + ( sb / scale )**2 )
       r = sign ( 1.0D+00, roe ) * r
       c = sa / r
       s = sb / r
     end if

     if ( 0.0D+00 < abs ( c ) .and. abs ( c ) <= s ) then
       z = 1.0D+00 / c
     else
       z = s
     end if

     sa = r
     sb = z
   end subroutine drotg


   subroutine dscal ( n, sa, x, incx )
   !*****************************************************************************80
   !
   !! DSCAL scales a vector by a constant.
   !
   !  Discussion:
   !
   !    This routine uses double precision real arithmetic.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    08 April 1999
   !
   !  Author:
   !
   !    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
   !    David Kincaid, Fred Krogh.
   !    FORTRAN90 version by John Burkardt.
   !
   !  Reference:
   !
   !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
   !    LINPACK User's Guide,
   !    SIAM, 1979,
   !    ISBN13: 978-0-898711-72-1,
   !    LC: QA214.L56.
   !
   !    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
   !    Algorithm 539,
   !    Basic Linear Algebra Subprograms for Fortran Usage,
   !    ACM Transactions on Mathematical Software,
   !    Volume 5, Number 3, September 1979, pages 308-323.
   !
   !  Parameters:
   !
   !    Input, integer N, the number of entries in the vector.
   !
   !    Input, real ( kind = rk ) SA, the multiplier.
   !
   !    Input/output, real ( kind = rk ) X(*), the vector to be scaled.
   !
   !    Input, integer INCX, the increment between successive
   !    entries of X.
   !
     implicit none

     integer, parameter :: rk = kind ( 1.0D+00 )

     integer i
     integer incx
     integer ix
     integer m
     integer n
     real ( kind = rk ) sa
     real ( kind = rk ) x(*)

     if ( n <= 0 ) then

     else if ( incx == 1 ) then

       m = mod ( n, 5 )

       x(1:m) = sa * x(1:m)

       do i = m+1, n, 5
         x(i)   = sa * x(i)
         x(i+1) = sa * x(i+1)
         x(i+2) = sa * x(i+2)
         x(i+3) = sa * x(i+3)
         x(i+4) = sa * x(i+4)
       end do

     else

       if ( 0 <= incx ) then
         ix = 1
       else
         ix = ( - n + 1 ) * incx + 1
       end if

       do i = 1, n
         x(ix) = sa * x(ix)
         ix = ix + incx
       end do

     end if
   end subroutine dscal



   subroutine dsvdc( a, lda, m, n, s, e, u, ldu, v, ldv, work, job, info )
   !*****************************************************************************80
   !
   !! DSVDC computes the singular value decomposition of a real rectangular matrix.
   !
   !  Discussion:
   !
   !    This routine reduces an M by N matrix A to diagonal form by orthogonal
   !    transformations U and V.  The diagonal elements S(I) are the singular
   !    values of A.  The columns of U are the corresponding left singular
   !    vectors, and the columns of V the right singular vectors.
   !
   !    The form of the singular value decomposition is then
   !
   !      A(MxN) = U(MxM) * S(MxN) * V(NxN)'
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    16 September 2006
   !
   !  Author:
   !
   !    FORTRAN90 version by John Burkardt.
   !
   !  Reference:
   !
   !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
   !    LINPACK User's Guide,
   !    SIAM, 1979,
   !    ISBN13: 978-0-898711-72-1,
   !    LC: QA214.L56.
   !
   !  Parameters:
   !
   !    Input/output, real ( kind = rk ) A(LDA,N).  On input, the M by N
   !    matrix whose singular value decomposition is to be computed.
   !    On output, the matrix has been destroyed.  Depending on the user's
   !    requests, the matrix may contain other useful information.
   !
   !    Input, integer LDA, the leading dimension of the array A.
   !    LDA must be at least N.
   !
   !    Input, integer M, the number of rows of the matrix.
   !
   !    Input, integer N, the number of columns of the matrix A.
   !
   !    Output, real ( kind = rk ) S(MM), where MM = max(M+1,N).  The first
   !    min(M,N) entries of S contain the singular values of A arranged in
   !    descending order of magnitude.
   !
   !    Output, real ( kind = rk ) E(MM), where MM = max(M+1,N).  Ordinarily
   !    contains zeros.  However see the discussion of INFO for exceptions.
   !
   !    Output, real ( kind = rk ) U(LDU,K).  If JOBA = 1 then K = M;
   !    if 2 <= JOBA, then K = min(M,N).  U contains the M by M matrix of
   !    left singular vectors.  U is not referenced if JOBA = 0.  If M <= N
   !    or if JOBA = 2, then U may be identified with A in the subroutine call.
   !
   !    Input, integer LDU, the leading dimension of the array U.
   !    LDU must be at least M.
   !
   !    Output, real ( kind = rk ) V(LDV,N), the N by N matrix of right singular
   !    vectors.  V is not referenced if JOB is 0.  If N <= M, then V may be
   !    identified with A in the subroutine call.
   !
   !    Input, integer LDV, the leading dimension of the array V.
   !    LDV must be at least N.
   !
   !    Workspace, real ( kind = rk ) WORK(M).
   !
   !    Input, integer JOB, controls the computation of the singular
   !    vectors.  It has the decimal expansion AB with the following meaning:
   !      A =  0, do not compute the left singular vectors.
   !      A =  1, return the M left singular vectors in U.
   !      A >= 2, return the first min(M,N) singular vectors in U.
   !      B =  0, do not compute the right singular vectors.
   !      B =  1, return the right singular vectors in V.
   !
   !    Output, integer INFO, status indicator.
   !    The singular values (and their corresponding singular vectors)
   !    S(INFO+1), S(INFO+2),...,S(MN) are correct.  Here MN = min ( M, N ).
   !    Thus if INFO is 0, all the singular values and their vectors are
   !    correct.  In any event, the matrix B = U' * A * V is the bidiagonal
   !    matrix with the elements of S on its diagonal and the elements of E on
   !    its superdiagonal.  Thus the singular values of A and B are the same.
   !
     implicit none

     integer, parameter :: rk = kind ( 1.0D+00 )

     integer lda
     integer ldu
     integer ldv
     integer m
     integer n

     real ( kind = rk ) a(lda,n)
     real ( kind = rk ) b
     real ( kind = rk ) c
     real ( kind = rk ) cs
     real ( kind = rk ) e(*)
     real ( kind = rk ) el
     real ( kind = rk ) emm1
     real ( kind = rk ) f
     real ( kind = rk ) g
     integer info
     integer iter
     integer j
     integer job
     integer jobu
     integer k
     integer kase
     integer kk
     integer l
     integer ll
     integer lls
     integer ls
     integer lu
     integer, parameter :: maxit = 30
     integer mm
     integer mm1
     integer mn
     integer nct
     integer nctp1
     integer ncu
     integer nrt
     integer nrtp1
     real ( kind = rk ) s(*)
     real ( kind = rk ) scale
   !  real ( kind = rk ) ddot
     real ( kind = rk ) shift
     real ( kind = rk ) sl
     real ( kind = rk ) sm
     real ( kind = rk ) smm1
     real ( kind = rk ) sn
   !  real ( kind = rk ) dnrm2
     real ( kind = rk ) t
     real ( kind = rk ) t1
     real ( kind = rk ) test
     real ( kind = rk ) u(ldu,m)
     real ( kind = rk ) v(ldv,n)
     logical wantu
     logical wantv
     real ( kind = rk ) work(m)
     real ( kind = rk ) ztest
   !
   !  Determine what is to be computed.
   !
     wantu = .false.
     wantv = .false.
     jobu = mod ( job, 100 ) / 10

     if ( 1 < jobu ) then
       ncu = min ( m, n )
     else
       ncu = m
     end if

     if ( jobu /= 0 ) then
       wantu = .true.
     end if

     if ( mod ( job, 10 ) /= 0 ) then
       wantv = .true.
     end if
   !
   !  Reduce A to bidiagonal form, storing the diagonal elements
   !  in S and the super-diagonal elements in E.
   !
     info = 0
     nct = min ( m-1, n )
     nrt = max ( 0, min ( m, n-2 ) )
     lu = max ( nct, nrt )

     do l = 1, lu
   !
   !  Compute the transformation for the L-th column and
   !  place the L-th diagonal in S(L).
   !
       if ( l <= nct ) then

         s(l) = dnrm2 ( m-l+1, a(l,l), 1 )

         if ( s(l) /= 0.0D+00 ) then
           if ( a(l,l) /= 0.0D+00 ) then
             s(l) = sign ( s(l), a(l,l) )
           end if
           call dscal ( m-l+1, 1.0D+00 / s(l), a(l,l), 1 )
           a(l,l) = 1.0D+00 + a(l,l)
         end if

         s(l) = -s(l)

       end if

       do j = l+1, n
   !
   !  Apply the transformation.
   !
         if ( l <= nct .and. s(l) /= 0.0D+00 ) then
           t = -ddot ( m-l+1, a(l,l), 1, a(l,j), 1 ) / a(l,l)
           call daxpy ( m-l+1, t, a(l,l), 1, a(l,j), 1 )
         end if
   !
   !  Place the L-th row of A into E for the
   !  subsequent calculation of the row transformation.
   !
         e(j) = a(l,j)

       end do
   !
   !  Place the transformation in U for subsequent back multiplication.
   !
       if ( wantu .and. l <= nct ) then
         u(l:m,l) = a(l:m,l)
       end if
   !
   !  Compute the L-th row transformation and place the
   !  L-th superdiagonal in E(L).
   !
       if ( l <= nrt ) then

         e(l) = dnrm2 ( n-l, e(l+1), 1 )

         if ( e(l) /= 0.0D+00 ) then
           if ( e(l+1) /= 0.0D+00 ) then
             e(l) = sign ( e(l), e(l+1) )
           end if
           call dscal ( n-l, 1.0D+00 / e(l), e(l+1), 1 )
           e(l+1) = 1.0D+00 + e(l+1)
         end if

         e(l) = -e(l)
   !
   !  Apply the transformation.
   !
         if ( l + 1 <= m .and. e(l) /= 0.0D+00 ) then

           work(l+1:m) = 0.0D+00

           do j = l+1, n
             call daxpy ( m-l, e(j), a(l+1,j), 1, work(l+1), 1 )
           end do

           do j = l+1, n
             call daxpy ( m-l, -e(j)/e(l+1), work(l+1), 1, a(l+1,j), 1 )
           end do

         end if
   !
   !  Place the transformation in V for subsequent back multiplication.
   !
         if ( wantv ) then
           v(l+1:n,l) = e(l+1:n)
         end if

       end if

     end do
   !
   !  Set up the final bidiagonal matrix of order MN.
   !
     mn = min ( m + 1, n )
     nctp1 = nct + 1
     nrtp1 = nrt + 1

     if ( nct < n ) then
       s(nctp1) = a(nctp1,nctp1)
     end if

     if ( m < mn ) then
       s(mn) = 0.0D+00
     end if

     if ( nrtp1 < mn ) then
       e(nrtp1) = a(nrtp1,mn)
     end if

     e(mn) = 0.0D+00
   !
   !  If required, generate U.
   !
     if ( wantu ) then

       u(1:m,nctp1:ncu) = 0.0D+00

       do j = nctp1, ncu
         u(j,j) = 1.0D+00
       end do

       do ll = 1, nct

         l = nct - ll + 1

         if ( s(l) /= 0.0D+00 ) then

           do j = l+1, ncu
             t = -ddot ( m-l+1, u(l,l), 1, u(l,j), 1 ) / u(l,l)
             call daxpy ( m-l+1, t, u(l,l), 1, u(l,j), 1 )
           end do

           u(l:m,l) = -u(l:m,l)
           u(l,l) = 1.0D+00 + u(l,l)
           u(1:l-1,l) = 0.0D+00

         else

           u(1:m,l) = 0.0D+00
           u(l,l) = 1.0D+00

         end if

       end do

     end if
   !
   !  If it is required, generate V.
   !
     if ( wantv ) then

       do ll = 1, n

         l = n - ll + 1

         if ( l <= nrt .and. e(l) /= 0.0D+00 ) then

           do j = l + 1, n
             t = -ddot ( n-l, v(l+1,l), 1, v(l+1,j), 1 ) / v(l+1,l)
             call daxpy ( n-l, t, v(l+1,l), 1, v(l+1,j), 1 )
           end do

         end if

         v(1:n,l) = 0.0D+00
         v(l,l) = 1.0D+00

       end do

     end if
   !
   !  Main iteration loop for the singular values.
   !
     mm = mn
     iter = 0

     do while ( 0 < mn )
   !
   !  If too many iterations have been performed, set flag and return.
   !
       if ( maxit <= iter ) then
         info = mn
         return
       end if
   !
   !  This section of the program inspects for
   !  negligible elements in the S and E arrays.
   !
   !  On completion the variables KASE and L are set as follows:
   !
   !  KASE = 1     if S(MN) and E(L-1) are negligible and L < MN
   !  KASE = 2     if S(L) is negligible and L < MN
   !  KASE = 3     if E(L-1) is negligible, L < MN, and
   !               S(L), ..., S(MN) are not negligible (QR step).
   !  KASE = 4     if E(MN-1) is negligible (convergence).
   !
       do ll = 1, mn

         l = mn - ll

         if ( l == 0 ) then
           exit
         end if

         test = abs ( s(l) ) + abs ( s(l+1) )
         ztest = test + abs ( e(l) )

         if ( ztest == test ) then
           e(l) = 0.0D+00
           exit
         end if

       end do

       if ( l == mn - 1 ) then

         kase = 4

       else

         do lls = l + 1, mn + 1

           ls = mn - lls + l + 1

           if ( ls == l ) then
             exit
           end if

           test = 0.0D+00
           if ( ls /= mn ) then
             test = test + abs ( e(ls) )
           end if

           if ( ls /= l + 1 ) then
             test = test + abs ( e(ls-1) )
           end if

           ztest = test + abs ( s(ls) )

           if ( ztest == test ) then
             s(ls) = 0.0D+00
             exit
           end if

         end do

         if ( ls == l ) then
           kase = 3
         else if ( ls == mn ) then
           kase = 1
         else
           kase = 2
           l = ls
         end if

       end if

       l = l + 1
   !
   !  Deflate negligible S(MN).
   !
       if ( kase == 1 ) then

         mm1 = mn - 1
         f = e(mn-1)
         e(mn-1) = 0.0D+00

         do kk = l, mm1

           k = mm1 - kk + l
           t1 = s(k)
           call drotg ( t1, f, cs, sn )
           s(k) = t1

           if ( k /= l ) then
             f = -sn * e(k-1)
             e(k-1) = cs * e(k-1)
           end if

           if ( wantv ) then
             call drot ( n, v(1,k), 1, v(1,mn), 1, cs, sn )
           end if

         end do
   !
   !  Split at negligible S(L).
   !
       else if ( kase == 2 ) then

         f = e(l-1)
         e(l-1) = 0.0D+00

         do k = l, mn

           t1 = s(k)
           call drotg ( t1, f, cs, sn )
           s(k) = t1
           f = -sn * e(k)
           e(k) = cs * e(k)
           if ( wantu ) then
             call drot ( m, u(1,k), 1, u(1,l-1), 1, cs, sn )
           end if

         end do
   !
   !  Perform one QR step.
   !
       else if ( kase == 3 ) then
   !
   !  Calculate the shift.
   !
         scale = max ( abs ( s(mn) ), abs ( s(mn-1) ), abs ( e(mn-1) ), &
                       abs ( s(l) ), abs ( e(l) ) )

         sm = s(mn) / scale
         smm1 = s(mn-1) / scale
         emm1 = e(mn-1) / scale
         sl = s(l) / scale
         el = e(l) / scale
         b = ( ( smm1 + sm ) * ( smm1 - sm ) + emm1 * emm1 ) / 2.0D+00
         c = sm  * sm * emm1 * emm1
         shift = 0.0D+00

         if ( b /= 0.0D+00 .or. c /= 0.0D+00 ) then
           shift = sqrt ( b * b + c )
           if ( b < 0.0D+00 ) then
             shift = -shift
           end if
           shift = c / ( b + shift )
         end if

         f = ( sl + sm ) * ( sl - sm ) + shift
         g = sl * el
   !
   !  Chase zeros.
   !
         mm1 = mn - 1

         do k = l, mm1

           call drotg ( f, g, cs, sn )

           if ( k /= l ) then
             e(k-1) = f
           end if

           f = cs * s(k) + sn * e(k)
           e(k) = cs * e(k) - sn * s(k)
           g = sn * s(k+1)
           s(k+1) = cs * s(k+1)

           if ( wantv ) then
             call drot ( n, v(1,k), 1, v(1,k+1), 1, cs, sn )
           end if

           call drotg ( f, g, cs, sn )
           s(k) = f
           f = cs * e(k) + sn * s(k+1)
           s(k+1) = -sn * e(k) + cs * s(k+1)
           g = sn * e(k+1)
           e(k+1) = cs * e(k+1)

           if ( wantu .and. k < m ) then
             call drot ( m, u(1,k), 1, u(1,k+1), 1, cs, sn )
           end if

         end do

         e(mn-1) = f
         iter = iter + 1
   !
   !  Convergence.
   !
       else if ( kase == 4 ) then
   !
   !  Make the singular value nonnegative.
   !
         if ( s(l) < 0.0D+00 ) then
           s(l) = -s(l)
           if ( wantv ) then
             v(1:n,l) = -v(1:n,l)
           end if
         end if
   !
   !  Order the singular value.
   !
         do

           if ( l == mm ) then
             exit
           end if

           if ( s(l+1) <= s(l) ) then
             exit
           end if

           t = s(l)
           s(l) = s(l+1)
           s(l+1) = t

           if ( wantv .and. l < n ) then
             call dswap ( n, v(1,l), 1, v(1,l+1), 1 )
           end if

           if ( wantu .and. l < m ) then
             call dswap ( m, u(1,l), 1, u(1,l+1), 1 )
           end if

           l = l + 1

         end do

         iter = 0
         mn = mn - 1
       end if
     end do
   end subroutine dsvdc



   subroutine dswap( n, x, incx, y, incy )
   !*****************************************************************************80
   !
   !! DSWAP interchanges two vectors.
   !
   !  Discussion:
   !
   !    This routine uses double precision real arithmetic.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    08 April 1999
   !
   !  Author:
   !
   !    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
   !    David Kincaid, Fred Krogh.
   !    FORTRAN90 version by John Burkardt.
   !
   !  Reference:
   !
   !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
   !    LINPACK User's Guide,
   !    SIAM, 1979,
   !    ISBN13: 978-0-898711-72-1,
   !    LC: QA214.L56.
   !
   !    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
   !    Algorithm 539,
   !    Basic Linear Algebra Subprograms for Fortran Usage,
   !    ACM Transactions on Mathematical Software,
   !    Volume 5, Number 3, September 1979, pages 308-323.
   !
   !  Parameters:
   !
   !    Input, integer N, the number of entries in the vectors.
   !
   !    Input/output, real ( kind = rk ) X(*), one of the vectors to swap.
   !
   !    Input, integer INCX, the increment between successive
   !    entries of X.
   !
   !    Input/output, real ( kind = rk ) Y(*), one of the vectors to swap.
   !
   !    Input, integer INCY, the increment between successive
   !    elements of Y.
   !
      implicit none
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer :: i
      integer :: incx
      integer :: incy
      integer :: ix
      integer :: iy
      integer :: m
      integer :: n
      real ( kind = rk ) :: temp
      real ( kind = rk ) :: x(*)
      real ( kind = rk ) :: y(*)

      if ( n <= 0 ) then

      else if ( incx == 1 .and. incy == 1 ) then
         m = mod ( n, 3 )
         do i = 1, m
            temp = x(i)
            x(i) = y(i)
            y(i) = temp
         end do

         do i = m + 1, n, 3
            temp = x(i)
            x(i) = y(i)
            y(i) = temp

            temp = x(i+1)
            x(i+1) = y(i+1)
            y(i+1) = temp

            temp = x(i+2)
            x(i+2) = y(i+2)
            y(i+2) = temp
         end do
      else
         if ( 0 <= incx ) then
            ix = 1
         else
            ix = ( - n + 1 ) * incx + 1
         end if

         if ( 0 <= incy ) then
            iy = 1
         else
            iy = ( - n + 1 ) * incy + 1
         end if

         do i = 1, n
            temp = x(ix)
            x(ix) = y(iy)
            y(iy) = temp
            ix = ix + incx
            iy = iy + incy
         end do
      end if
   end subroutine dswap


   subroutine normal_solve( m, n, a, b, x, flag )
   !*****************************************************************************80
   !
   !! NORMAL_SOLVE solves a linear system using the normal equations.
   !
   !  Discussion:
   !
   !    Given a presumably rectangular MxN system of equations A*x=b, this routine
   !    sets up the NxN system A'*A*x=A'b.  Assuming N <= M, and that A has full
   !    column rank, the system will be solvable, and the vector x that is returned
   !    will minimize the Euclidean norm of the residual.
   !
   !    One drawback to this approach is that the condition number of the linear
   !    system A'*A is effectively the square of the condition number of A,
   !    meaning that there is a substantial loss of accuracy.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    17 April 2012
   !
   !  Author:
   !
   !    John Burkardt
   !
   !  Reference:
   !
   !    David Kahaner, Cleve Moler, Steven Nash,
   !    Numerical Methods and Software,
   !    Prentice Hall, 1989,
   !    ISBN: 0-13-627258-4,
   !    LC: TA345.K34.
   !
   !  Parameters:
   !
   !    Input, integer M, the number of rows of A.
   !
   !    Input, integer N, the number of columns of A.
   !    It must be the case that N <= M.
   !
   !    Input, real ( kind = rk ) A(M,N), the matrix.
   !    The matrix must have full column rank.
   !
   !    Input, real ( kind = rk ) B(M), the right hand side.
   !
   !    Output, real ( kind = rk ) X(N), the least squares solution.
   !
   !    Output, integer FLAG,
   !    0, no error was detected.
   !    1, an error occurred.
   !
      implicit none
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer :: m
      integer :: n
      real ( kind = rk ) :: a(m,n)
      real ( kind = rk ) :: ata(n,n)
      real ( kind = rk ) :: ata_c(n,n)
      real ( kind = rk ) :: atb(n)
      real ( kind = rk ) :: b(m)
      integer :: flag
      real ( kind = rk ) :: x(n)

      flag = 0
      if ( m < n ) then
         flag = 1
         return
      end if

      ata(1:n,1:n) = matmul( transpose( a(1:m,1:n) ), a(1:m,1:n) )
      atb(1:n) = matmul( transpose( a(1:m,1:n) ), b(1:m) )

      call r8mat_cholesky_factor ( n, ata, ata_c, flag )

      if ( flag /= 0 ) then
         return
      end if

      call r8mat_cholesky_solve ( n, ata_c, atb, x )
   end subroutine normal_solve


   subroutine Get_qr_solve( m, n, a, b, x )
   !*****************************************************************************80
   !
   !! QR_SOLVE solves a linear system in the least squares sense.
   !
   !  Discussion:
   !
   !    If the matrix A has full column rank, then the solution X should be the
   !    unique vector that minimizes the Euclidean norm of the residual.
   !
   !    If the matrix A does not have full column rank, then the solution is
   !    not unique; the vector X will minimize the residual norm, but so will
   !    various other vectors.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    15 April 2012
   !
   !  Author:
   !
   !    John Burkardt
   !
   !  Reference:
   !
   !    David Kahaner, Cleve Moler, Steven Nash,
   !    Numerical Methods and Software,
   !    Prentice Hall, 1989,
   !    ISBN: 0-13-627258-4,
   !    LC: TA345.K34.
   !
   !  Parameters:
   !
   !    Input, integer M, the number of rows of A.
   !
   !    Input, integer N, the number of columns of A.
   !
   !    Input, real ( kind = rk ) A(M,N), the matrix.
   !
   !    Input, real ( kind = rk ) B(M), the right hand side.
   !
   !    Output, real ( kind = rk ) X(N), the least squares solution.
   !
     implicit none
     integer, parameter :: rk = kind ( 1.0D+00 )
     integer :: m
     integer :: n
     real ( kind = rk ) :: a(m,n)
     real ( kind = rk ) :: a_qr(m,n)
     real ( kind = rk ) :: b(m)
     integer :: ind
     integer :: itask
     integer :: jpvt(n)
     integer :: kr
     integer :: lda
     real ( kind = rk ) :: qraux(n)
     real ( kind = rk ) :: r(m)
     real ( kind = rk ) :: tol
     real ( kind = rk ) :: work(n)
     real ( kind = rk ) :: x(n)

     a_qr(1:m,1:n) = a(1:m,1:n)
     lda = m
     tol = epsilon ( tol ) / maxval ( abs ( a_qr(1:m,1:n) ) )
     itask = 1

     call dqrls ( a_qr, lda, m, n, tol, kr, b, x, r, jpvt, qraux, work, &
       itask, ind )
   end subroutine Get_qr_solve



   subroutine Get_svd_solve( m, n, a, b, x )
   !*****************************************************************************80
   !
   !! SVD_SOLVE solves a linear system in the least squares sense.
   !
   !  Discussion:
   !
   !    The vector X returned by this routine should always minimize the
   !    Euclidean norm of the residual ||A*x-b||.
   !
   !    If the matrix A does not have full column rank, then there are multiple
   !    vectors that attain the minimum residual.  In that case, the vector
   !    X returned by this routine is the unique such minimizer that has the
   !    the minimum possible Euclidean norm, that is, ||A*x-b|| and ||x||
   !    are both minimized.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    07 September 2012
   !
   !  Reference:
   !
   !    David Kahaner, Cleve Moler, Steven Nash,
   !    Numerical Methods and Software,
   !    Prentice Hall, 1989,
   !    ISBN: 0-13-627258-4,
   !    LC: TA345.K34.
   !
   !  Parameters:
   !
   !    Input, integer M, the number of rows of A.
   !
   !    Input, integer N, the number of columns of A.
   !
   !    Input, real ( kind = rk ) A(M,N), the matrix.
   !
   !    Input, real ( kind = rk ) B(M), the right hand side.
   !
   !    Output, real ( kind = rk ) X(N), the least squares solution.
   !
      implicit none
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer :: m
      integer :: n
      real ( kind = rk ) :: a(m,n)
      real ( kind = rk ) :: a_copy(m,n)
      real ( kind = rk ) :: b(m)
      real ( kind = rk ) :: e(max(m+1,n))
      integer :: i
      integer :: info
      integer :: lda
      integer :: ldu
      integer :: ldv
      integer :: job
      real ( kind = rk ) :: sdiag(max(m+1,n))
      real ( kind = rk ) :: smax
      real ( kind = rk ) :: stol
      real ( kind = rk ) :: sub(n)
      real ( kind = rk ) :: u(m,m)
      real ( kind = rk ) :: ub(m)
      real ( kind = rk ) :: v(n,n)
      real ( kind = rk ) :: work(m)
      real ( kind = rk ) :: x(n)

      !  Get the SVD.

      a_copy(1:m,1:n) = a(1:m,1:n)
      lda = m
      ldu = m
      ldv = n
      job = 11

      call dsvdc( a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, job, info )

      if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SVD_SOLVE - Failure!'
      write ( *, '(a)' ) '  The SVD could not be calculated.'
      write ( *, '(a)' ) '  LINPACK routine DSVDC returned a nonzero'
      write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info
      stop
      end if

      ub(1:m) = matmul ( transpose ( u(1:m,1:m) ), b(1:m) )
      sub(1:n) = 0.0D+00
      !
      !  For singular problems, there may be tiny but nonzero singular values
      !  that should be ignored.  This is a reasonable attempt to avoid such
      !  problems, although in general, the user might wish to control the tolerance.
      !
      smax = maxval ( sdiag(1:n) )
      if ( smax <= epsilon ( smax ) ) then
         smax = 1.0D+00
      end if

      stol = epsilon ( smax ) * smax

      do i = 1, n
         if ( i <= m ) then
            if ( stol <= sdiag(i) ) then
               sub(i) = ub(i) / sdiag(i)
            end if
         end if
      end do

     x(1:n) = matmul ( v(1:n,1:n), sub(1:n) )
   end subroutine Get_svd_solve



   subroutine r8mat_cholesky_factor( n, a, c, flag )
   !*****************************************************************************80
   !
   !! r8mat_cholesky_factor computes the Cholesky factor of a symmetric matrix.
   !
   !  Discussion:
   !
   !    The matrix must be symmetric and positive semidefinite.
   !
   !    For a positive semidefinite symmetric matrix A, the Cholesky factorization
   !    is a lower triangular matrix L such that:
   !
   !      A = L * L'
   !
   !    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    08 April 2009
   !
   !  Author:
   !
   !    John Burkardt
   !
   !  Input:
   !
   !    integer N, the number of rows and columns.
   !
   !    real ( kind = rk ) A(N,N), the N by N matrix.
   !
   !  Output:
   !
   !    real ( kind = rk ) C(N,N), the lower triangular Cholesky factor.
   !
   !    integer FLAG:
   !    0, no error occurred.
   !    1, the matrix is not positive definite.
   !
      implicit none

      integer, parameter :: rk = kind ( 1.0D+00 )
      integer :: n
      real ( kind = rk ) :: a(n,n)
      real ( kind = rk ) :: c(n,n)
      integer :: flag
      integer :: i
      integer :: j
      real ( kind = rk ) :: sum2

      flag = 0
      c(1:n,1:n) = a(1:n,1:n)

      do j = 1, n
         c(1:j-1,j) = 0.0D+00
         do i = j, n
            sum2 = c(j,i) - dot_product ( c(j,1:j-1), c(i,1:j-1) )
            if ( i == j ) then
               if ( sum2 <= 0.0D+00 ) then
                  flag = 1
               else
                  c(i,j) = sqrt ( sum2 )
               end if
            else
               if ( c(j,j) /= 0.0D+00 ) then
                  c(i,j) = sum2 / c(j,j)
               else
                  c(i,j) = 0.0D+00
               end if
            end if
         end do
      end do
   end subroutine r8mat_cholesky_factor



   subroutine r8mat_cholesky_solve( n, l, b, x )
   !*****************************************************************************80
   !
   !! r8mat_cholesky_solve solves a Cholesky factored linear system A * x = b.
   !
   !  Discussion:
   !
   !    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
   !
   !    This routine works with the lower triangular Cholesky factor A = L * L'.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    07 December 2004
   !
   !  Author:
   !
   !    John Burkardt
   !
   !  Input:
   !
   !    integer N, the number of rows and columns.
   !
   !    real ( kind = rk ) L(N,N), the N by N lower Cholesky factor of the
   !    system matrix A.
   !
   !    real ( kind = rk ) B(N), the right hand side of the linear system.
   !
   !  Output:
   !
   !    real ( kind = rk ) X(N), the solution of the linear system.
   !
      implicit none
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer :: n
      real ( kind = rk ) :: b(n)
      real ( kind = rk ) :: l(n,n)
      real ( kind = rk ) :: x(n)
      !
      !  Solve L * y = b.
      !
      call r8mat_l_solve ( n, l, b, x )
      !
      !  Solve L' * x = y.
      !
      call r8mat_lt_solve ( n, l, x, x )
   end subroutine r8mat_cholesky_solve



   subroutine r8mat_l_solve( n, a, b, x )
   !*****************************************************************************80
   !
   !! r8mat_l_solve solves a lower triangular linear system.
   !
   !  Discussion:
   !
   !    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    07 December 2004
   !
   !  Author:
   !
   !    John Burkardt
   !
   !  Input:
   !
   !    integer N, the number of rows and columns.
   !
   !    real ( kind = rk ) A(N,N), the N by N lower triangular matrix.
   !
   !    real ( kind = rk ) B(N), the right hand side of the linear system.
   !
   !  Output:
   !
   !    real ( kind = rk ) X(N), the solution of the linear system.
   !
      implicit none
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer :: n
      real ( kind = rk ) :: a(n,n)
      real ( kind = rk ) :: b(n)
      integer :: i
      real ( kind = rk ) :: x(n)
      !
      !  Solve L * x = b.
      !
      do i = 1, n
         x(i) = ( b(i) - dot_product ( a(i,1:i-1), x(1:i-1) ) ) / a(i,i)
      end do
   end subroutine r8mat_l_solve



   subroutine r8mat_lt_solve( n, a, b, x )
   !*****************************************************************************80
   !
   !! r8mat_lt_solve solves a transposed lower triangular linear system.
   !
   !  Discussion:
   !
   !    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
   !
   !    Given the lower triangular matrix A, the linear system to be solved is:
   !
   !      A' * x = b
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    07 December 2004
   !
   !  Author:
   !
   !    John Burkardt
   !
   !  Input:
   !
   !    integer N, the number of rows and columns.
   !
   !    real ( kind = rk ) A(N,N), the N by N lower triangular matrix.
   !
   !    real ( kind = rk ) B(N), the right hand side of the linear system.
   !
   !  Output:
   !
   !    real ( kind = rk ) X(N), the solution of the linear system.
   !
      implicit none
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer :: n
      real ( kind = rk ) :: a(n,n)
      real ( kind = rk ) :: b(n)
      integer :: i
      real ( kind = rk ) :: x(n)
      !
      !  Solve L'*x = b.
      !
      do i = n, 1, -1
         x(i) = ( b(i) - dot_product( x(i+1:n), a(i+1:n,i) ) ) / a(i,i)
      end do
   end subroutine r8mat_lt_solve



end module mod_math_lib
