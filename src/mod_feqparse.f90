! copyright 2020 fluid numerics llc
! all rights reserved.
!
! author : joe schoonover ( joe@fluidnumerics.com )
!
! equationparser defines a public class that can be used to parse and evaluate strings
! representative of equations. an equation, written in infix form, is converted to
! postfix form and evaluated using a postfix calculator.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

module mod_feqparse

use iso_fortran_env

implicit none

  real(real64), parameter :: pi   = 4.0_real64*atan(1.0_real64)

  integer, parameter, private :: error_message_length = 256
  integer, parameter, private :: max_equation_length  = 1024
  integer, parameter, private :: max_function_length  = 5
  integer, parameter, private :: max_variable_length  = 12
  integer, parameter, private :: token_length         = 48
  integer, parameter, private :: stack_length         = 128

  ! token types
  integer, parameter, private :: none_token               = 0
  integer, parameter, private :: number_token             = 1
  integer, parameter, private :: variable_token           = 2
  integer, parameter, private :: operator_token           = 3
  integer, parameter, private :: function_token           = 4
  integer, parameter, private :: openingparentheses_token = 5
  integer, parameter, private :: closingparentheses_token = 6
  integer, parameter, private :: monadic_token            = 7

  integer, parameter, private :: nfunctions = 14
  integer, parameter, private :: nseparators = 7

  type string
    character(10) :: str
  end type string

  character(1), dimension(7), private  :: separators = (/ "+", "-", "*", "/", "(", ")", "^" /)
  character(1), dimension(5), private  :: operators  = (/ "+", "-", "*", "/", "^" /)
  character(1), dimension(10), private :: numbers    = (/ "0", "1", "2", "3", "4", "5", "6", "7", "8", "9" /)
  type(string), dimension(14), private :: functions

  ! private types !

  type token
    character(token_length) :: tokenstring
    integer                 :: tokentype

    contains
      procedure :: equals_token

  end type token


  type tokenstack
    type(token), allocatable :: tokens(:)
    integer                  :: top_index = 0

    contains

      procedure :: construct => construct_tokenstack
      procedure :: destruct  => destruct_tokenstack

      procedure :: push      => push_tokenstack
      procedure :: pop       => pop_tokenstack
      procedure :: peek      => peek_tokenstack

      procedure :: isempty   => isempty_tokenstack
      procedure :: toptoken

  end type tokenstack

  type numberstack
    real(real64), allocatable :: tokens(:)
    integer                   :: top_index

    contains

      procedure :: construct => construct_numberstack
      procedure :: destruct  => destruct_numberstack

      procedure :: push      => push_numberstack
      procedure :: pop       => pop_numberstack
      procedure :: peek      => peek_numberstack

      procedure :: isempty   => isempty_numberstack

  end type numberstack

  private :: token, tokenstack, numberstack
  private :: construct_tokenstack, destruct_tokenstack, push_tokenstack, pop_tokenstack, peek_tokenstack, isempty_tokenstack
  private :: construct_numberstack, destruct_numberstack, push_numberstack, pop_numberstack, peek_numberstack, isempty_numberstack
  private :: isnumber, isvariable, isfunction, isoperator, isseparator, findlastfunctionindex, f_of_x, priority


  type equationparser
    character(max_equation_length)     :: equation
    character(max_variable_length)     :: variablename
    character(max_equation_length)     :: infixformula
    integer                            :: nindepvars
    character(len=1), allocatable      :: indepvars(:)
    type( tokenstack )                 :: infix
    type( tokenstack )                 :: postfix

    contains

      procedure :: destruct => destruct_equationparser
      procedure :: cleanequation
      procedure :: tokenize
      procedure :: converttopostfix

      generic :: evaluate => evaluate_real32, evaluate_real64
      procedure, private :: evaluate_real32, evaluate_real64

      procedure :: print_infixtokens
      procedure :: print_postfixtokens

  end type equationparser

  interface equationparser
    procedure construct_equationparser
  end interface equationparser


contains

   function construct_equationparser( equation, indepvars ) result( parser )
    type( equationparser ) :: parser
    character(*)           :: equation
    character(1)           :: indepvars(1:)
    ! local
    integer :: i
    character(error_message_length) :: errormsg
    logical                         :: equationisclean, tokenized, success
    integer                         :: nindepvars

      functions(1) % str = "cos"
      functions(2) % str = "sin"
      functions(3) % str = "tan"
      functions(4) % str = "tanh"
      functions(5) % str = "sqrt"
      functions(6) % str = "abs"
      functions(7) % str = "exp"
      functions(8) % str = "ln"
      functions(9) % str = "log"
      functions(10) % str = "acos"
      functions(11) % str = "asin"
      functions(12) % str = "atan"
      functions(13) % str = "sech"
      functions(14) % str = "rand"

      nindepvars = size(indepvars)
      allocate( parser % indepvars(1:nindepvars) )
      parser % nindepvars = nindepvars
      do i = 1, nindepvars
        parser % indepvars(i) = indepvars(i)
      enddo

      parser % infixformula = " "
      parser % equation = equation
      errormsg = " "

      call parser % cleanequation( equationisclean, errormsg )

      if( equationisclean )then

        call parser % tokenize( tokenized, errormsg )

        if( tokenized )then

          call parser % converttopostfix( )

        else

           print*, trim( errormsg )
           success = .false.

        endif

      end if

  end function construct_equationparser

  subroutine destruct_equationparser( parser )
    implicit none
    class(equationparser), intent(inout) :: parser

      deallocate( parser % indepvars )
      parser % infixformula = ""
      parser % equation = ""

      call parser % infix % destruct()
      call parser % postfix % destruct()

  end subroutine destruct_equationparser

  subroutine cleanequation( parser, equationcleaned, errormsg )
    class( equationparser ), intent(inout)       :: parser
    logical, intent(out)                         :: equationcleaned
    character(error_message_length), intent(out) :: errormsg
    ! local
    integer :: nchar, equalsignloc, j, i


      equationcleaned = .false.
      parser % variablename    = '#noname'

      nchar = len_trim( parser % equation )
      equalsignloc = index( parser % equation, "=" )

      if( equalsignloc == 0 )then
        errormsg = "no equal sign found"
        return
      endif

      parser % variablename = trim( parser % equation(1:equalsignloc-1) )

      ! grab the formula to the right of the equal sign and left adjust the formula
      parser % infixformula = parser % equation(equalsignloc+1:)
      parser % infixformula = adjustl(parser % infixformula)


      ! remove any spaces
      j = 1
      do i = 1, len_trim(parser % infixformula)
        if( parser % infixformula(i:i) /= " " )then
          parser % infixformula(j:j) = parser % infixformula(i:i)
          j = j + 1
        endif
      enddo

      parser % infixformula(j:max_equation_length) = " "

      equationcleaned = .true.

  end subroutine cleanequation

  subroutine tokenize( parser, tokenized, errormsg )
    class( equationparser ), intent(inout) :: parser
    logical, intent(out)                   :: tokenized
    character(error_message_length)        :: errormsg
    ! local
    integer :: i, j


      tokenized = .false.
      errormsg  = " "

      call parser % infix % construct( stack_length )

      i = 1
      do while( parser % infixformula(i:i) /= " " )

        if( isvariable( parser % infixformula(i:i), parser % indepvars, parser % nindepvars ) )then

          parser % infix % top_index = parser % infix % top_index + 1
          parser % infix % tokens( parser % infix % top_index ) % tokenstring = parser % infixformula(i:i)
          parser % infix % tokens( parser % infix % top_index ) % tokentype   = variable_token
          i = i+1

          ! next item must be an operator, closing parentheses, or end of equation

          if( .not. isoperator( parser % infixformula(i:i) ) .and. &
              parser % infixformula(i:i) /= ")" .and. parser % infixformula(i:i) /= " "  )then

            errormsg = "missing operator or closing parentheses after token : "//&
                       trim( parser % infix % tokens( parser % infix % top_index ) % tokenstring )
            return

          endif

        elseif( isnumber( parser % infixformula(i:i) ) )then

          parser % infix % top_index = parser % infix % top_index + 1
          parser % infix % tokens( parser % infix % top_index ) % tokenstring = ''


          if( parser % infixformula(i:i) == 'p' .or. parser % infixformula(i:i) == 'p' )then

            ! conditional for using built in "pi" definition
            parser % infix % tokens( parser % infix % top_index ) % tokenstring(1:2) = parser % infixformula(i:i+1)
            j = 2

          else

            j = 0
            do while( isnumber( parser % infixformula(i+j:i+j) ) )

              parser % infix % tokens( parser % infix % top_index ) % tokenstring(j+1:j+1) = parser % infixformula(i+j:i+j)
              j = j+1

            enddo

          endif

          parser % infix % tokens( parser % infix % top_index ) % tokentype = number_token

          i = i + j

          ! next item must be an operator or a closing parentheses
          if( .not. isoperator( parser % infixformula(i:i) ) .and. &
              parser % infixformula(i:i) /= ")" .and. parser % infixformula(i:i) /= " " )then

            errormsg = "missing operator or closing parentheses after token : "//&
                       trim( parser % infix % tokens( parser % infix % top_index ) % tokenstring )
            return

          endif

        elseif( isseparator( parser % infixformula(i:i) ) )then


          parser % infix % top_index = parser % infix % top_index + 1
          parser % infix % tokens( parser % infix % top_index ) % tokenstring = parser % infixformula(i:i)

          if( parser % infixformula(i:i) == "(" )then
            parser % infix % tokens( parser % infix % top_index ) % tokentype   = openingparentheses_token
          elseif( parser % infixformula(i:i) == ")" )then
            parser % infix % tokens( parser % infix % top_index ) % tokentype   = closingparentheses_token
          else
            parser % infix % tokens( parser % infix % top_index ) % tokentype   = operator_token
          endif

          i = i + 1


        elseif( isfunction( parser % infixformula(i:i) ) )then

          parser % infix % top_index = parser % infix % top_index + 1
          parser % infix % tokens( parser % infix % top_index ) % tokenstring = ''

          j = findlastfunctionindex( parser % infixformula(i:i+max_function_length-1) )

          parser % infix % tokens( parser % infix % top_index ) % tokenstring = parser % infixformula(i:i+j)
          parser % infix % tokens( parser % infix % top_index ) % tokentype   = function_token
          i = i+j+1

          ! check to see if the next string
          if( parser % infixformula(i:i) /= "(" )then
            errormsg = "missing opening parentheses after token : "//&
                       trim( parser % infix % tokens( parser % infix % top_index ) % tokenstring )

            return
          endif


        else

          errormsg = "invalid token : "//&
                     trim( parser % infixformula(i:i) )

          return

        endif

      enddo


      if( parser % infix % tokens(1) % tokentype == operator_token )then
         if( trim( parser % infix % tokens(1) % tokenstring ) == "+" .or. &
              trim( parser % infix % tokens(1) % tokenstring ) == "-" ) then
            parser % infix % tokens(1) % tokentype = monadic_token
         end if
      end if

      do i = 2, parser % infix % top_index
         if( parser % infix % tokens(i) % tokentype == operator_token .and. &
              parser % infix % tokens(i-1) % tokentype == openingparentheses_token ) then
            parser % infix % tokens(i) % tokentype = monadic_token
         end if
      end do


      tokenized = .true.

  end subroutine tokenize


  subroutine converttopostfix( parser )
    class( equationparser ), intent(inout) :: parser
    ! local
    character(error_message_length) :: errormsg
    type( tokenstack )              :: operator_stack
    type( token )                   :: tok
    integer                         :: i

      !success = .false.

      call parser % postfix % construct( stack_length )
      call operator_stack % construct( stack_length )

      do i = 1, parser % infix % top_index

        if( parser % infix % tokens(i) % tokentype == variable_token .or. &
            parser % infix % tokens(i) % tokentype == number_token )then


          call parser % postfix % push( parser % infix % tokens(i) )


        elseif( parser % infix % tokens(i) % tokentype == function_token )then

          call operator_stack % push( parser % infix % tokens(i) )

        elseif( parser % infix % tokens(i) % tokentype == operator_token &
                .or. parser % infix % tokens(i) % tokentype == monadic_token )then


          if( .not. operator_stack % isempty( ) )then

            tok = operator_stack % toptoken( )

            do while( trim(tok % tokenstring) /= "(" .and. &
                      priority( trim(tok % tokenstring) ) > &
                       priority( trim(parser % infix % tokens(i) % tokenstring) ) )

              call parser % postfix % push( tok )
              call operator_stack % pop( tok )
              tok = operator_stack % toptoken( )

            enddo

          endif

          call operator_stack % push( parser % infix % tokens(i) )

        elseif( parser % infix % tokens(i) % tokentype == openingparentheses_token )then

          call operator_stack % push( parser % infix % tokens(i) )


        elseif( parser % infix % tokens(i) % tokentype == closingparentheses_token )then

          tok = operator_stack % toptoken( )

          do while( .not.( operator_stack % isempty( ) ) .and. trim(tok % tokenstring) /= "(" )

            call parser % postfix % push( tok )
            call operator_stack % pop( tok )
            tok = operator_stack % toptoken( )

          enddo

          ! pop the opening parenthesis
          call operator_stack % pop( tok )

        endif

      enddo

      ! pop the remaining operators
      do while( .not.( operator_stack % isempty( ) ) )

        tok = operator_stack % toptoken( )
        call parser % postfix % push( tok )
        call operator_stack % pop( tok )

      enddo

  end subroutine converttopostfix

  function evaluate_real32( parser, x ) result( f )
    class(equationparser) :: parser
    real(real32) :: x(1:parser % nindepvars)
    real(real32) :: f
    ! local
    integer :: i, k
    type(token) :: t
    type(numberstack) :: stack
    real(real64) :: v, a, b, c

      call stack % construct( stack_length )

      if( .not.( allocated( parser % postfix % tokens ) ) )then

        f = 0.0_real32

      else

        do k = 1, parser % postfix % top_index

          t = parser % postfix % tokens(k) % equals_token( )

          select case ( t % tokentype )

            case( number_token )

              if( t % tokenstring == "pi" .or. t % tokenstring == "pi" )     then
                 v = pi
              else
                read( t % tokenstring, * ) v
              end if

              call stack % push( v )

            case ( variable_token )

              do i = 1, parser % nindepvars
                if( trim( t % tokenstring ) == parser % indepvars(i) )then
                  call stack % push( real(x(i),real64) )
                  exit
                endif
              enddo

              !if( trim( t % tokenstring ) == "x" )then

              !   call stack % push( x(1) )

              !elseif( trim( t % tokenstring ) == "y" )then

              !   call stack % push( x(2) )

              !elseif( trim( t % tokenstring ) == "z" )then

              !   call stack % push( x(3) )

              !endif

            case ( operator_token )

              call stack % pop( a )
              call stack % pop( b )

              select case ( trim(t % tokenstring) )

                 case ( "+" )

                    c = a + b

                 case ( "-" )

                    c = b - a

                 case ( "*" )

                    c = a*b

                 case ( "/" )

                    c = b/a

                 case ( "^" )

                    c = b**a
                 case default

              end select

              call stack % push( c )

           case ( function_token )

              call stack % pop( a )

              b = f_of_x( trim(t % tokenstring), a )

              call stack % push( b )

           case ( monadic_token )

             if( trim(t % tokenstring) == "-" )     then

                call stack % pop( a )
                a = -a
                call stack % push( a )

             end if

           case default

         end select

       end do

       call stack % pop( a )
       f = a

       call stack % destruct( )

     endif

  end function evaluate_real32

  function evaluate_real64( parser, x ) result( f )
    class(equationparser) :: parser
    real(real64) :: x(1:parser % nindepvars)
    real(real64) :: f
    ! local
    integer :: i, k
    type(token) :: t
    type(numberstack) :: stack
    real(real64) :: v, a, b, c

      call stack % construct( stack_length )

      if( .not.( allocated( parser % postfix % tokens ) ) )then

        f = 0.0_real64

      else

        do k = 1, parser % postfix % top_index

          t = parser % postfix % tokens(k) % equals_token( )

          select case ( t % tokentype )

            case( number_token )

              if( t % tokenstring == "pi" .or. t % tokenstring == "pi" )     then
                 v = pi
              else
                read( t % tokenstring, * ) v
              end if

              call stack % push( v )

            case ( variable_token )

              do i = 1, parser % nindepvars
                if( trim( t % tokenstring ) == parser % indepvars(i) )then
                  call stack % push( x(i) )
                  exit
                endif
              enddo

              !if( trim( t % tokenstring ) == "x" )then

              !   call stack % push( x(1) )

              !elseif( trim( t % tokenstring ) == "y" )then

              !   call stack % push( x(2) )

              !elseif( trim( t % tokenstring ) == "z" )then

              !   call stack % push( x(3) )

              !endif

            case ( operator_token )

              call stack % pop( a )
              call stack % pop( b )

              select case ( trim(t % tokenstring) )

                 case ( "+" )

                    c = a + b

                 case ( "-" )

                    c = b - a

                 case ( "*" )

                    c = a*b

                 case ( "/" )

                    c = b/a

                 case ( "^" )

                    c = b**a
                 case default

              end select

              call stack % push( c )

           case ( function_token )

              call stack % pop( a )

              b = f_of_x( trim(t % tokenstring), a )

              call stack % push( b )

           case ( monadic_token )

             if( trim(t % tokenstring) == "-" )     then

                call stack % pop( a )
                a = -a
                call stack % push( a )

             end if

           case default

         end select

       end do

       call stack % pop( a )
       f = a

       call stack % destruct( )

     endif

  end function evaluate_real64

  subroutine print_infixtokens( parser )
    class( equationparser ), intent(in) :: parser
    ! local
    integer :: i

      do i = 1, parser % infix % top_index
        print*, trim( parser % infix % tokens(i) % tokenstring )
      enddo

  end subroutine print_infixtokens

  subroutine print_postfixtokens( parser )
    class( equationparser ), intent(in) :: parser
    ! local
    integer :: i

      do i = 1, parser % postfix % top_index
        print*, trim( parser % postfix % tokens(i) % tokenstring )
      enddo

  end subroutine print_postfixtokens

  ! tokenstack and numberstack

  subroutine construct_tokenstack( stack, n )
   class(tokenstack), intent(out) :: stack
   integer, intent(in)            :: n

     allocate( stack % tokens(1:n) )
     stack % top_index = 0

  end subroutine construct_tokenstack

  subroutine destruct_tokenstack( stack )
    class(tokenstack), intent(inout) :: stack

      if( allocated( stack % tokens ) ) deallocate( stack % tokens )
      stack % top_index = 0

  end subroutine destruct_tokenstack

  subroutine push_tokenstack( stack, tok )
    class(tokenstack), intent(inout) :: stack
    type(token), intent(in)         :: tok

      stack % top_index                  = stack % top_index + 1
      stack % tokens(stack % top_index)  % tokenstring = tok % tokenstring
      stack % tokens(stack % top_index)  % tokentype   = tok % tokentype

  end subroutine push_tokenstack

  subroutine pop_tokenstack( stack, tok )
    class(tokenstack), intent(inout) :: stack
    type(token), intent(out)        :: tok

      if( stack % top_index <= 0 ) then
        print *, "attempt to pop from empty token stack"
      else
        tok % tokenstring         = stack % tokens( stack % top_index ) % tokenstring
        tok % tokentype           = stack % tokens( stack % top_index ) % tokentype
        stack % top_index = stack % top_index - 1
      end if


  end subroutine pop_tokenstack

  subroutine peek_tokenstack( stack, tok )
    class(tokenstack), intent(in) :: stack
    type(token), intent(out)     :: tok

      if( stack % top_index <= 0 ) then
        print *, "attempt to peek from empty token stack"
      else
        tok % tokenstring = stack % tokens( stack % top_index ) % tokenstring
        tok % tokentype   = stack % tokens( stack % top_index ) % tokentype
      end if
  end subroutine peek_tokenstack

  logical function isempty_tokenstack( stack )
    class( tokenstack ) :: stack

      isempty_tokenstack = .false.

      if( stack % top_index <= 0 )then
        isempty_tokenstack = .true.
      endif

  end function isempty_tokenstack

  type( token ) function toptoken( stack )
    class( tokenstack ) :: stack

      if( stack % top_index > 0 )then
        toptoken % tokenstring = stack % tokens( stack % top_index ) % tokenstring
        toptoken % tokentype   = stack % tokens( stack % top_index ) % tokentype
      else
        toptoken % tokenstring = ''
      endif

  end function toptoken

  function equals_token( tok1 ) result( tok2 )
    class(token) :: tok1
    type(token)  :: tok2

      tok2 % tokenstring = tok1 % tokenstring
      tok2 % tokentype   = tok1 % tokentype

  end function equals_token

  subroutine construct_numberstack( stack, n )
   class(numberstack), intent(out) :: stack
   integer, intent(in)            :: n

     allocate( stack % tokens(1:n) )
     stack % top_index = 0

  end subroutine construct_numberstack

  subroutine destruct_numberstack( stack )
    class(numberstack), intent(inout) :: stack

      if( allocated( stack % tokens) ) deallocate( stack % tokens )
      stack % top_index = 0

  end subroutine destruct_numberstack

  subroutine push_numberstack( stack, tok )
    class(numberstack), intent(inout) :: stack
    real(real64), intent(in)         :: tok

      stack % top_index                  = stack % top_index + 1
      stack % tokens(stack % top_index) = tok

  end subroutine push_numberstack

  subroutine pop_numberstack( stack, tok )
    class(numberstack), intent(inout) :: stack
    real(real64), intent(out)        :: tok

      if( stack % top_index <= 0 ) then
        print *, "attempt to pop from empty token stack"
      else
        tok               = stack % tokens( stack % top_index )
        stack % top_index = stack % top_index - 1
      end if


  end subroutine pop_numberstack

  subroutine peek_numberstack( stack, tok )
    class(numberstack), intent(in) :: stack
    real(real64), intent(out)        :: tok

      if( stack % top_index <= 0 ) then
        print *, "attempt to peek from empty token stack"
      else
        tok = stack % tokens( stack % top_index )
      end if
  end subroutine peek_numberstack

  logical function isempty_numberstack( stack )
    class( numberstack ) :: stack

      isempty_numberstack = .false.

      if( stack % top_index <= 0 )then
        isempty_numberstack = .true.
      endif

  end function isempty_numberstack

  ! support functions !

  logical function isseparator( eqchar )
    character(1) :: eqchar
    ! local
    integer :: i

      isseparator = .false.
      do i = 1, nseparators

        if( eqchar == separators(i) )then
          isseparator = .true.
        endif

      enddo

  end function isseparator

  logical function isnumber( eqchar )
    character(1) :: eqchar
    ! local
    integer :: i

      isnumber = .false.

      if( eqchar == '.' .or. eqchar == 'p' .or. eqchar == 'p' )then
        isnumber = .true.
        return
      endif

      do i = 1, 10

        if( eqchar == numbers(i) )then
          isnumber = .true.
        endif

      enddo

  end function isnumber

  logical function isvariable( eqchar, variables, nvariables )
    character(1) :: eqchar
    integer      :: nvariables
    character(1) :: variables(1:nvariables)
    ! local
    integer :: i

      isvariable = .false.
      do i = 1, nvariables

        if( eqchar == variables(i) )then
          isvariable = .true.
        endif

      enddo

  end function isvariable

  logical function isoperator( eqchar )
    character(1) :: eqchar
    ! local
    integer :: i

      isoperator = .false.
      do i = 1, 5

        if( eqchar == operators(i) )then
          isoperator = .true.
        endif

      enddo

  end function isoperator

  logical function isfunction( eqchar )
    character(1) :: eqchar
    ! local
    integer :: i

      isfunction = .false.
      do i = 1, nfunctions

        if( eqchar == functions(i) % str(1:1) )then
          isfunction = .true.
        endif

      enddo

  end function isfunction

  function findlastfunctionindex( eqchar ) result( j )
    character(max_function_length) :: eqchar
    integer                        :: i, j

      do i = 1, max_function_length

        if( eqchar(i:i) == "(" )then
          j = i-2
          exit
        endif

      enddo

  end function findlastfunctionindex

  real(real64) function f_of_x( func, x )
    character(*) :: func
    real(real64)   :: x
    ! local
    real(real64)   :: r

      if( trim( func ) == "cos" .or. trim( func ) == "cos" )then

        f_of_x = cos( x )

      elseif( trim( func ) == "sin" .or. trim( func ) == "sin" )then

        f_of_x = sin( x )

      elseif( trim( func ) == "tan" .or. trim( func ) == "tan" )then

        f_of_x = tan( x )

      elseif( trim( func ) == "tanh" .or. trim( func ) == "tanh" )then

        f_of_x = tanh( x )

      elseif( trim( func ) == "sech" .or. trim( func ) == "sech" )then

        f_of_x = 2.0_real64/( exp(x) + exp(-x) )

      elseif( trim( func ) == "sqrt" .or. trim( func ) == "sqrt" )then

        f_of_x = sqrt( x )

      elseif( trim( func ) == "abs" .or. trim( func ) == "abs" )then

        f_of_x = abs( x )

      elseif( trim( func ) == "exp" .or. trim( func ) == "exp" )then

        f_of_x = exp( x )

      elseif( trim( func ) == "ln" .or. trim( func ) == "ln" )then

        f_of_x = log( x )

      elseif( trim( func ) == "log" .or. trim( func ) == "log" )then

        f_of_x = log10( x )

      elseif( trim( func ) == "acos" .or. trim( func ) == "acos" )then

        f_of_x = acos( x )

      elseif( trim( func ) == "asin" .or. trim( func ) == "asin" )then

        f_of_x = asin( x )

      elseif( trim( func ) == "atan" .or. trim( func ) == "atan" )then

        f_of_x = atan( x )

      elseif( trim( func ) == "rand" .or. trim( func ) == "rand" )then

        call random_number( r )
        f_of_x = r*x

      else

        f_of_x = 0.0_real64

      endif


  end function f_of_x

  integer function priority( operatorstring )
    character(1) :: operatorstring


      if( isfunction( operatorstring ) )then

        priority = 4

      elseif( operatorstring == '^' )then

        priority = 3

      elseif( operatorstring == '*' .or. operatorstring == '/' )then

        priority = 2

      elseif( operatorstring == '+' .or. operatorstring == '-' )then

        priority = 1

      else

        priority = 0

      endif

  end function priority

end module mod_feqparse
