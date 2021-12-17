module mod_ac_linear
!***************************************************************************************************
!> This module contains the required type and procedures for linear interpolation in 2D and 3D space.
!> The procedures mainly do Least Square approximations.
!> Currently it contains the methods and required types just for scalar-value fields interpolation.
!
!***************************************************************************************************

   use mod_aux_precision
   use mod_math_lib

   implicit none

   type ac_linear_scalar_t
      real(WP), allocatable :: a(:,:)
      real(WP), allocatable :: b(:)
      real(WP), allocatable :: u(:)
      integer               :: ndim
      integer               :: order
      integer               :: nbasis
      integer               :: nd
   end type ac_linear_scalar_t


   interface ac_build
         module procedure Build_linear_scalar
   end interface


   interface ac_set_linear_solver_scalar
         module procedure Set_linear_scalar_2d
         module procedure Set_linear_scalar_3d
   end interface


   interface ac_solve_sys_scalar
         module procedure Solve_sys_2d
         module procedure Solve_sys_3d
   end interface


contains


   subroutine Build_linear_scalar( sys, nd, ndim, order )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      type(ac_linear_scalar_t), intent(inout) :: sys
      integer                 , intent(in   ) :: nd
      integer                 , intent(in   ) :: ndim
      integer                 , intent(in   ) :: order

      sys%order = order
      sys%ndim  = ndim
      sys%nd    = nd

      call Get_num_monomial_basis( sys%ndim, sys%order, sys%nbasis )

      allocate( sys%a(1:nd, 1:sys%nbasis) )
      allocate( sys%u(1:sys%nbasis) )
      allocate( sys%b(1:nd) )
   end subroutine  Build_linear_scalar
   !************************************************************************************************



   pure subroutine Set_linear_scalar_2d( sys, xc, yc )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      ! Arguments
      type(ac_linear_scalar_t), intent(inout) :: sys
      real(WP),                 intent(in   ) :: xc(:)
      real(WP),                 intent(in   ) :: yc(:)

      call Set_vandermond_matrix_2d( a=sys%a, ncol=sys%nbasis, nrow=sys%nd, order=sys%order, &
                                     xc=xc, yc=yc )
   end subroutine Set_linear_scalar_2d
   !************************************************************************************************



   pure subroutine Set_vandermond_matrix_2d( a, ncol, nrow, order,  xc, yc )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      real(WP), intent(inout) :: a(:,:)
      real(WP), intent(in   ) :: xc(:)
      real(WP), intent(in   ) :: yc(:)
      integer,  intent(in   ) :: ncol
      integer,  intent(in   ) :: nrow
      integer,  intent(in   ) :: order
      !> Locals
      integer                 :: ex, ey
      integer                 :: s, j

      j = 0
      do s = 0, order
         do ex = s, 0, -1
            ey = s - ex
            j = j + 1
            a(1:nrow, j) = (xc(1:nrow) ** ex) * (yc(1:nrow) ** ey)
         end do
      end do
   end subroutine Set_vandermond_matrix_2d
   !************************************************************************************************



   pure subroutine Set_vandermond_matrix_3d( a, ncol, nrow, order, xc, yc, zc )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      real(WP), intent(inout) :: a(:,:)
      real(WP), intent(in   ) :: xc(:)
      real(WP), intent(in   ) :: yc(:)
      real(WP), intent(in   ) :: zc(:)
      integer,  intent(in   ) :: ncol
      integer,  intent(in   ) :: nrow
      integer,  intent(in   ) :: order
      !> Locals
      integer                 :: ex, ey, ez
      integer                 :: s, j

      j = 0
      do s = 0, order
         do ex = s, 0, -1
            do ey = s - ex, 0, -1
               ez = s - ex - ey
               j = j + 1
               a(1:nrow, j) = (xc(1:nrow) ** ex) * (yc(1:nrow) ** ey) * (zc(1:nrow) ** ez)
            end do
         end do
      end do
   end subroutine Set_vandermond_matrix_3d
   !************************************************************************************************



   pure subroutine Set_linear_scalar_3d( sys, xc, yc, zc )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      ! Arguments
      type(ac_linear_scalar_t), intent(inout) :: sys
      real(WP),                 intent(in   ) :: xc(:)
      real(WP),                 intent(in   ) :: yc(:)
      real(WP),                 intent(in   ) :: zc(:)

      call Set_vandermond_matrix_3d( a=sys%a, ncol=sys%nbasis, nrow=sys%nd, order=sys%order, &
                                     xc=xc, yc=yc, zc=zc )
   end subroutine Set_linear_scalar_3d
   !************************************************************************************************



   subroutine Solve_sys_2d( sys, xi, yi, fd, fi )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      type(ac_linear_scalar_t), intent(in   ) :: sys
      real(WP),                 intent(in   ) :: fd(:)
      real(WP),                 intent(in   ) :: xi(:)
      real(WP),                 intent(in   ) :: yi(:)
      real(WP),                 intent(  out) :: fi(:)

!      call Get_qr_solve( sys%nd, sys%nbasis, sys%a, fd, sys%u )
      call Get_svd_solve( sys%nd, sys%nbasis, sys%a, fd, sys%u )
      call Get_poly_value_2d( sys%order, sys%u, xi, yi, fi )
   end subroutine Solve_sys_2d



   subroutine Solve_sys_3d( sys, xi, yi, zi, fd, fi )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      type(ac_linear_scalar_t), intent(in   ) :: sys
      real(WP),                 intent(in   ) :: fd(:)
      real(WP),                 intent(in   ) :: xi(:)
      real(WP),                 intent(in   ) :: yi(:)
      real(WP),                 intent(in   ) :: zi(:)
      real(WP),                 intent(  out) :: fi(:)

!      call Get_qr_solve( sys%nd, sys%nbasis, sys%a, fd, sys%u )
      call Get_svd_solve( sys%nd, sys%nbasis, sys%a, fd, sys%u )
      call Get_poly_value_3d( sys%order, sys%u, xi, yi, zi, fi )
   end subroutine Solve_sys_3d
   !************************************************************************************************



   subroutine Get_poly_value_2d( m, c, x, y, p )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      real(kind=WP), intent(in   ) :: c(:)
      real(kind=WP), intent(in   ) :: x(:)
      real(kind=WP), intent(in   ) :: y(:)
      real(kind=WP), intent(  out) :: p(:)
      integer,       intent(in   ) :: m
      !> Locals
      integer                      :: ex
      integer                      :: ey
      integer                      :: j
      integer                      :: s

      p(:) = 0.0D+00

      j = 0
      do s = 0, m
         do ex = s, 0, -1
            ey = s - ex
            j = j + 1
            p(:) = p(:) + c(j) * x(:) ** ex * y(:) ** ey
         end do
      end do
   end subroutine Get_poly_value_2d
   !************************************************************************************************



   subroutine Get_poly_value_3d ( m, c, x, y, z, p )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      real(kind=WP), intent(in   ) :: c(:)
      real(kind=WP), intent(in   ) :: x(:)
      real(kind=WP), intent(in   ) :: y(:)
      real(kind=WP), intent(in   ) :: z(:)
      real(kind=WP), intent(  out) :: p(:)
      integer,       intent(in   ) :: m
      !> Locals
      integer                      :: ex
      integer                      :: ey
      integer                      :: ez
      integer                      :: j
      integer                      :: s

      p(:) = 0.0D+00

      j = 0
      do s = 0, m
         do ex = s, 0, -1
            do ey = s - ex, 0, -1
               ez = s - ex - ey
               j = j + 1
               p(:) = p(:) + c(j) * x(:) ** ex * y(:) ** ey * z(:) ** ez
            end do
         end do
      end do
   end subroutine Get_poly_value_3d
   !************************************************************************************************



   pure function Factor( n ) result (res)
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      integer, intent(in) :: n
      !> Locals
      real(WP)            :: res
      integer             :: i

      res = product( [ (i, i=1, n) ] )
   end function Factor
   !************************************************************************************************



   pure function Conbination( n, k ) result (res)
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      integer, intent(in) :: n, k
      !> Locals
      real(WP)            :: res

      res = Factor( n ) / ( Factor( k ) * Factor( n-k ) )
   end function Conbination
   !************************************************************************************************



   pure function num_monomial_terms( ndim, order ) result(res)
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      integer, intent(in) :: ndim, order
      !> Locals!************************************************************************************************
      integer             :: res

      res  = Conbination ( ndim + order, order )
   end function num_monomial_terms
   !************************************************************************************************



   subroutine Get_num_monomial_basis( ndim, order, nbasis )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      integer, intent(in) :: ndim, order
      !> Locals
      integer             :: nbasis

      nbasis = Conbination( ndim + order, order )
   end subroutine Get_num_monomial_basis
   !************************************************************************************************



   subroutine ac_destroy_sys (sys)
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      type(ac_linear_scalar_t), intent(inout) :: sys

      if ( allocated( sys%a ) ) deallocate( sys%a )
      if ( allocated( sys%u ) ) deallocate( sys%u )
      if ( allocated( sys%b ) ) deallocate( sys%b )
   end subroutine ac_destroy_sys
   !************************************************************************************************

end module mod_ac_linear

