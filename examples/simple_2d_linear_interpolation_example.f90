
program simple_2d_linear_interpolation_example

   use mod_ac_linear
   use mod_aux_procs_for_tests

   implicit none
   type(ac_linear_scalar_t) :: las
   real(dp), allocatable    :: xd(:), yd(:)
   real(dp), allocatable    :: xi(:), yi(:)
   real(dp), allocatable    :: fd(:)
   real(dp), allocatable    :: fi(:)

   !> According to the given name, it makes the nodes and their values form which we do interpolate
   call ac_get_sample_data( name="uni_cart_1", xd=xd, yd=yd, fd=fd )

   !> According to the given name, it gives the location(s) we will interpolate on.
   call ac_set_sample_nodes( name="sample_nodes_1", xi=xi, yi=yi, fi=fi )

   !> Initialize the system for interpolation.
   !> This step must be taken BEFORE the parallel region
   call ac_build( sys=las, nd=size(fd), ndim=2, order=3 )

   !> Sets the system. This step can be taken WITHIN the parallel region
   call ac_set_linear_solver_scalar( las, xd, yd )

   !> Solves the system to get the interpolated value(s)
   !> This step can be taken WITHIN the parallel region
   call ac_solve_sys_scalar( las, xi, yi, fd, fi  )


   call ac_get_conclusion_message( fi, case_name="'2D' scalar function on a uniform Cartesian mesh" )

   !> Frees the taken memory
   !> This step must be taken AFTER the parallel region
   call ac_destroy_sys(las)



end program simple_2d_linear_interpolation_example
