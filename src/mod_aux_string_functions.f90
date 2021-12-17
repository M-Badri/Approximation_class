module mod_aux_string_functions
!***************************************************************************************************
!> This module contains some auxiliary procedures and parameters for working with string variables.
!
!***************************************************************************************************

   use mod_aux_precision, only: WP

   implicit none

   !> Module parameters
   character(len=1), parameter :: CBLANK  = " "
   integer,          parameter :: CLENGTH = 30


contains


   pure function itoa( i ) result( res )
   !************************************************************************************************
   !> Returns a string of an integer number
   !
   !************************************************************************************************
      implicit none
      integer,                  intent(in) :: i
      character(:), allocatable            :: res

      !> Local
      character(range( i ) + 2) :: tmp

      write( tmp,'(i0)' ) i
      res = trim( tmp )
   end function itoa
   !************************************************************************************************



   pure function rtoa( r ) result( res )
   !************************************************************************************************
   !> Returns a string of a real number
   !
   !************************************************************************************************
      implicit none
      real(WP),                 intent (in) :: r
      character(:), allocatable             :: res

      !> Local
      character(range( r ) + 2) :: tmp

      write( tmp,'(g0)' ) r
      res = trim( tmp )
   end function rtoa
   !************************************************************************************************

end module mod_aux_string_functions
