module mod_aux_procs_for_tests
!***************************************************************************************************
!> This module contains some auxiliary procedures for making 2d and 3d examples to check the
!> interpolations out.
!
!***************************************************************************************************

   use mod_aux_precision
   use mod_feqparse

   implicit none

   private
   public :: ac_set_sample_nodes
   public :: ac_get_sample_data
   public :: ac_get_conclusion_message

   !> Module-wise variables
   real(WP)              :: lx(2), ly(2), lz(2)
   real(WP), allocatable :: f_exact(:)


   interface ac_get_sample_data
      module procedure ac_get_sample_data_2d
      module procedure ac_get_sample_data_3d
   end interface


   interface ac_set_sample_nodes
      module procedure ac_set_sample_nodes_2d
      module procedure ac_set_sample_nodes_3d
   end interface


contains



   subroutine ac_get_sample_data_2d( name, xd, yd, fd )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      character(len=*), intent(in   )     :: name
      real(wp),allocatable, intent(  out) :: xd(:)
      real(wp),allocatable, intent(  out) :: yd(:)
      real(wp),allocatable, intent(  out) :: fd(:)
      !> Locals
      character(len=:), allocatable   :: errmsg

      select case ( name )

         case ( "uni_cart_1")
!            print*, "uni_cart_1"
            call Get_uni_cart_2d_1( xd, yd, fd )

         case default
            errmsg = "The sample case << " // name // " >> is not defined!"
            stop errmsg

      end select
   end subroutine ac_get_sample_data_2d
   !************************************************************************************************


   subroutine ac_get_sample_data_3d( name, xd, yd, zd, fd )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      character(len=*), intent(in   )     :: name
      real(wp),allocatable, intent(  out) :: xd(:)
      real(wp),allocatable, intent(  out) :: yd(:)
      real(wp),allocatable, intent(  out) :: zd(:)
      real(wp),allocatable, intent(  out) :: fd(:)
      !> Locals
      character(len=:), allocatable   :: errmsg

      select case ( name )

         case ( "uni_cart_1")
!            print*, "uni_cart_1"
            call Get_uni_cart_3d_1( xd, yd, zd, fd )

         case default
            errmsg = "The sample case << " // name // " >> is not defined!"
            stop errmsg

      end select
   end subroutine ac_get_sample_data_3d
   !************************************************************************************************



   subroutine ac_set_sample_nodes_2d( name, xi, yi, fi)
   !************************************************************************************************
   !
   !************************************************************************************************
      character(len=*), intent(in   )      :: name
      real(WP), allocatable, intent(  out) :: xi(:)
      real(WP), allocatable, intent(  out) :: yi(:)
      real(WP), allocatable, intent(  out) :: fi(:)
      !> Locals
      character(len=:), allocatable   :: errmsg

      select case ( name )

         case ( "sample_nodes_1")
            call Get_sample_node_2d_1( xi, yi, fi)

         case default
            errmsg = "The sample case << " // name // " >> is not defined!"
            stop errmsg

      end select
   end subroutine ac_set_sample_nodes_2d
   !************************************************************************************************



   subroutine ac_set_sample_nodes_3d( name, xi, yi, zi, fi)
   !************************************************************************************************
   !
   !************************************************************************************************
      character(len=*), intent(in   )      :: name
      real(WP), allocatable, intent(  out) :: xi(:)
      real(WP), allocatable, intent(  out) :: yi(:)
      real(WP), allocatable, intent(  out) :: zi(:)
      real(WP), allocatable, intent(  out) :: fi(:)
      !> Locals
      character(len=:), allocatable   :: errmsg

      select case ( name )

         case ( "sample_nodes_1")
            call Get_sample_node_3d_1( xi, yi, zi, fi)

         case default
            errmsg = "The sample case << " // name // " >> is not defined!"
            stop errmsg

      end select
   end subroutine ac_set_sample_nodes_3d
   !************************************************************************************************



   subroutine Get_uni_cart_2d_1( xd, yd, fd )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      real(WP),allocatable, intent(  out) :: xd(:)
      real(WP),allocatable, intent(  out) :: yd(:)
      real(WP),allocatable, intent(  out) :: fd(:)
      !> Locals
      real(Wp), allocatable           :: xyd(:, :)
      integer, parameter              :: nx = 9
      integer, parameter              :: ny = 5
      integer, parameter              :: nd = nx * ny
      integer                         :: inode
      character(len=1)                :: independentvars(2)
      character(len=30)               :: eqChar
      type(equationparser)            :: eparser

      !> Specify the independent variables
      independentVars = ['x', 'y']

      !> Specify an equation string that we want to evaluate
      eqChar = 'f = (cos(x) + sin(y))'

      !> Create the EquationParser object
      eparser = EquationParser( eqChar, independentVars )

      !> Evaluate the equation
      lx(:) = [0.0d0, 1.0d0]
      ly(:) = [0.0d0, 1.0d0]
      lz(:) = [0.0d0, 1.0d0]
      call Get_rect_mesh_2d( nx, ny, lx, ly, xd, yd )
      call Get_rect_coordiate_vector( xd, yd, xyd )

      allocate( fd(nd) )
      do inode = 1, nd
         fd(inode) = eparser%evaluate( xyd(:, inode) )
      end do

      if ( allocated( xd ) ) deallocate( xd )
      if ( allocated( yd ) ) deallocate( yd )

      allocate( xd, source=xyd(1,:) )
      allocate( yd, source=xyd(2,:) )
   end subroutine Get_uni_cart_2d_1
   !************************************************************************************************



   subroutine Get_uni_cart_3d_1( xd, yd, zd, fd )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      real(WP),allocatable, intent(  out) :: xd(:)
      real(WP),allocatable, intent(  out) :: yd(:)
      real(WP),allocatable, intent(  out) :: zd(:)
      real(WP),allocatable, intent(  out) :: fd(:)
      !> Locals
      real(Wp), allocatable           :: xyzd(:, :)
      integer, parameter              :: nx = 9
      integer, parameter              :: ny = 5
      integer, parameter              :: nz = 2
      integer, parameter              :: nd = nx * ny * nz
      integer                         :: inode
      character(len=1)                :: independentvars(3)
      character(len=30)               :: eqChar
      type(equationparser)            :: eparser

      !> Specify the independent variables
      independentVars = ['x', 'y', 'z']

      !> Specify an equation string that we want to evaluate
      eqChar = 'f = (cos(x) + sin(y) + cos(z))'

      !> Create the EquationParser object
      eparser = EquationParser( eqChar, independentVars )

      !> Evaluate the equation
      lx(:) = [0.0d0, 1.0d0]
      ly(:) = [0.0d0, 1.0d0]
      lz(:) = [0.0d0, 1.0d0]
      call Get_rect_mesh_2d( nx, ny, lx, ly, xd, yd )
      call Get_rect_mesh_3d( nx, ny, nz, lx, ly, lz, xd, yd, zd )
      call Get_rect_coordiate_vector_3d( xd, yd, zd, xyzd )

      allocate( fd(nd) )
      do inode = 1, nd
         fd(inode) = eparser%evaluate( xyzd(:, inode) )
      end do

      if ( allocated( xd ) ) deallocate( xd )
      if ( allocated( yd ) ) deallocate( yd )
      if ( allocated( zd ) ) deallocate( zd )

      allocate( xd, source=xyzd(1,:) )
      allocate( yd, source=xyzd(2,:) )
      allocate( zd, source=xyzd(3,:) )
   end subroutine Get_uni_cart_3d_1
   !************************************************************************************************



   subroutine Get_sample_node_2d_1( xi, yi, fi )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      real(WP), allocatable, intent(  out) :: xi(:)
      real(WP), allocatable, intent(  out) :: yi(:)
      real(WP), allocatable, intent(  out) :: fi(:)
      integer                              :: ni

      ni = 1
      allocate( fi(1:ni) )
      allocate( xi(1:ni) )
      allocate( yi(1:ni) )
      xi = 0.5d0 * sum( lx(:) ) + 0.2d0
      yi = 0.5d0 * sum( ly(:) ) - 0.1d0
      fi = 0.0d0

      allocate( f_exact(1:ni) )
      f_exact(:) = cos( xi(:) ) + sin( yi(:) )
   end subroutine Get_sample_node_2d_1
   !************************************************************************************************



   subroutine Get_sample_node_3d_1( xi, yi, zi, fi )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      real(WP), allocatable, intent(  out) :: xi(:)
      real(WP), allocatable, intent(  out) :: yi(:)
      real(WP), allocatable, intent(  out) :: zi(:)
      real(WP), allocatable, intent(  out) :: fi(:)
      integer                              :: ni

      ni = 1
      allocate( fi(1:ni) )
      allocate( xi(1:ni) )
      allocate( yi(1:ni) )
      allocate( zi(1:ni) )
      xi = 0.5d0 * sum( lx(:) ) + 0.2d0
      yi = 0.5d0 * sum( ly(:) ) - 0.1d0
      zi = 0.5d0 * sum( lz(:) ) + 0.3d0
      fi = 0.0d0

      allocate( f_exact(1:ni) )
      f_exact(:) = cos( xi(:) ) + sin( yi(:) ) + cos( zi(:) )
   end subroutine Get_sample_node_3d_1
   !************************************************************************************************



   pure subroutine Get_rect_mesh_2d ( nx, ny, lx, ly, xs, ys )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      integer,               intent(in   ) :: nx, ny
      real(WP),              intent(in   ) :: lx(:), ly(:)
      real(WP), allocatable, intent(  out) :: xs(:), ys(:)

      if ( allocated( xs ) ) deallocate( xs )
      if ( allocated( ys ) ) deallocate( ys )
      allocate( xs(nx), ys(ny) )

      xs(:) = Uniform_linespace( nx, lx(1), lx(2) )
      ys(:) = Uniform_linespace( ny, ly(1), ly(2) )
   end subroutine Get_rect_mesh_2d
   !************************************************************************************************



   pure subroutine Get_rect_mesh_3d ( nx, ny, nz, lx, ly, lz, xs, ys, zs )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      integer,               intent(in   ) :: nx, ny, nz
      real(WP),              intent(in   ) :: lx(:), ly(:), lz(:)
      real(WP), allocatable, intent(  out) :: xs(:), ys(:), zs(:)

      if ( allocated( xs ) ) deallocate( xs )
      if ( allocated( ys ) ) deallocate( ys )
      if ( allocated( zs ) ) deallocate( zs )
      allocate( xs(nx), ys(ny), zs(nz) )

      xs(:) = Uniform_linespace( nx, lx(1), lx(2) )
      ys(:) = Uniform_linespace( ny, ly(1), ly(2) )
      zs(:) = Uniform_linespace( nz, lz(1), lz(2) )
   end subroutine Get_rect_mesh_3d
   !************************************************************************************************



   pure function Uniform_linespace ( n, lmin, lmax ) result(res)
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      real(WP), intent(in) :: lmin, lmax
      integer,  intent(in) :: n
      !> Locals
      real(WP)             :: res(n)
      real(WP)             :: dl
      integer              :: i

      dl = (lmax - lmin) / real( n-1, kind=WP )
      res(:) = [ ((lmin + real( i, kind=WP ) * dl), i=0, n-1) ]
   end function Uniform_linespace
   !************************************************************************************************



   subroutine Get_rect_coordiate_vector( xs, ys, cords )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      real(WP),              intent(in)    :: xs(:)
      real(WP),              intent(in)    :: ys(:)
      real(WP), allocatable, intent(  out) :: cords(:,:)
      !> Locals
      integer                              :: i, j
      integer                              :: inode, nnodes
      integer                              :: nx, ny

      nx = size(xs)
      ny = size(ys)
      nnodes = nx * ny

      if ( allocated( cords ) ) deallocate (cords)
      allocate( cords(2, nnodes) )

      inode = 0
      do j = 1, ny
         do i = 1, nx
            inode = inode + 1
            cords(1, inode) = xs(i)
            cords(2, inode) = ys(j)
         end do
      end do
   end subroutine Get_rect_coordiate_vector
   !************************************************************************************************



   subroutine Get_rect_coordiate_vector_3d( xs, ys, zs, cords )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      !> Arguments
      real(WP),              intent(in)    :: xs(:)
      real(WP),              intent(in)    :: ys(:)
      real(WP),              intent(in)    :: zs(:)
      real(WP), allocatable, intent(  out) :: cords(:,:)
      !> Locals
      integer                              :: i, j, k
      integer                              :: inode, nnodes
      integer                              :: nx, ny, nz

      nx = size(xs)
      ny = size(ys)
      nz = size(zs)
      nnodes = nx * ny * nz

      if ( allocated( cords ) ) deallocate (cords)
      allocate( cords(3, nnodes) )

      inode = 0
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               inode = inode + 1
               cords(1, inode) = xs(i)
               cords(2, inode) = ys(j)
               cords(3, inode) = zs(k)
            end do
         end do
      end do
   end subroutine Get_rect_coordiate_vector_3d
   !************************************************************************************************



   subroutine ac_get_conclusion_message( fi, case_name )
   !************************************************************************************************
   !
   !************************************************************************************************
      implicit none
      real(WP),         intent(in) :: fi(:)
      character(len=*), intent(in) :: case_name

      write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++++++"
      write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++++++"
      write(*,*) " "
      write(*,*) "The case:           << ", case_name, " >> "
      write(*,*) " "
      write(*,*) "The calculated value(s): ", fi
      write(*,*) " "
      write(*,*) "The exact values(s):     ", f_exact
      write(*,*) " "
      write(*,*) "The absolute error(s):    ", abs( fi - f_exact )
      write(*,*) " "
      write(*,*) "The relative error(s):    ", abs( fi - f_exact ) / abs( f_exact )
      write(*,*) " "
      write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++++++"
      write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++++++"
      write(*,*) " "
      write(*,*) " "
   end subroutine ac_get_conclusion_message
   !************************************************************************************************


end module mod_aux_procs_for_tests
