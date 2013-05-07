module aspicmodels
  use infprec, only : kp
  implicit none

  private

  include 'aspicmodels.h'

  logical, save :: CheckAspicPtrs = .false.

  procedure(onep_x_rrad), pointer :: ptronep_x_rrad => null()
  procedure(onep_x_rreh), pointer :: ptronep_x_rreh => null()
  procedure(onep_norm_potential), pointer :: ptronep_norm_potential => null()
  procedure(onep_epsilon_one), pointer :: ptronep_epsilon_one => null()
  procedure(onep_epsilon_two), pointer :: ptronep_epsilon_two => null()
  procedure(onep_epsilon_three), pointer :: ptronep_epsilon_three => null()

  procedure(twop_x_rrad), pointer :: ptrtwop_x_rrad => null()
  procedure(twop_x_rreh), pointer :: ptrtwop_x_rreh => null()
  procedure(twop_norm_potential), pointer :: ptrtwop_norm_potential => null()
  procedure(twop_epsilon_one), pointer :: ptrtwop_epsilon_one => null()
  procedure(twop_epsilon_two), pointer :: ptrtwop_epsilon_two => null()
  procedure(twop_epsilon_three), pointer :: ptrtwop_epsilon_three => null()


  public check_aspic_ptrs, free_aspic_ptrs

  public aspic_x_rrad, aspic_x_rreh
  public aspic_norm_potential
  public aspic_epsilon_one, aspic_epsilon_two, aspic_epsilon_three
  public aspic_x_endinf

contains


  
  subroutine initialize_aspic_ptrs(name)
    include 'aspic.h'
    implicit none
    character(len=*), intent(in) :: name
    
    if (CheckAspicPtrs) return

    select case (trim(name))

!reminder: cpp is only used in "traditional mode" for fortran
#define ONEPNAME(fooi) \
    case ('fooi'); \
       ptronep_x_rrad => fooi/**/_x_rrad ; \
       ptronep_x_rreh => fooi/**/_x_rreh ; \
       ptronep_norm_potential => fooi/**/_norm_potential ; \
       ptronep_epsilon_one => fooi/**/_epsilon_one ; \
       ptronep_epsilon_two => fooi/**/_epsilon_two ; \
       ptronep_epsilon_three => fooi/**/_epsilon_three

#define TWOPNAME(fooi) \
    case ('fooi') ; \
       ptrtwop_x_rrad => fooi/**/_x_rrad ; \
       ptrtwop_x_rreh => fooi/**/_x_rreh ; \
       ptrtwop_norm_potential => fooi/**/_norm_potential ; \
       ptrtwop_epsilon_one => fooi/**/_epsilon_one ; \
       ptrtwop_epsilon_two => fooi/**/_epsilon_two ; \
       ptrtwop_epsilon_three => fooi/**/_epsilon_three

#include "aspicnames.pp"
#undef ONEPNAME
#undef TWOPNAME

    case default
       stop 'initialize_aspic_ptrs: model not found!'

    end select

    CheckAspicPtrs = .true.

  end subroutine initialize_aspic_ptrs

   
  function check_aspic_ptrs()
    implicit none
    logical :: check_aspic_ptrs

    check_aspic_ptrs = CheckAspicPtrs

  end function check_aspic_ptrs


  subroutine free_aspic_ptrs()
    implicit none

    if (.not.CheckAspicPtrs) return

    if (associated(ptronep_x_rrad)) ptronep_x_rrad => null()    
    if (associated(ptronep_x_rreh)) ptronep_x_rreh => null()
    if (associated(ptronep_norm_potential)) ptronep_norm_potential => null()
    if (associated(ptronep_epsilon_one)) ptronep_epsilon_one => null()
    if (associated(ptronep_epsilon_two)) ptronep_epsilon_two => null()
    if (associated(ptronep_epsilon_three)) ptronep_epsilon_two => null()

    if (associated(ptrtwop_x_rrad)) ptrtwop_x_rrad => null()
    if (associated(ptrtwop_x_rreh)) ptrtwop_x_rreh => null()
    if (associated(ptrtwop_norm_potential)) ptrtwop_norm_potential => null()
    if (associated(ptrtwop_epsilon_one)) ptrtwop_epsilon_one => null()
    if (associated(ptrtwop_epsilon_two)) ptrtwop_epsilon_two => null()
    if (associated(ptrtwop_epsilon_three)) ptrtwop_epsilon_three => null()
    CheckAspicPtrs = .false.

  end subroutine free_aspic_ptrs




  function aspic_x_rrad(name,params,lnRrad,Pstar,bfoldstar)
    implicit none
    real(kp) :: aspic_x_rrad
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params
    real(kp), intent(in) :: lnRrad, Pstar
    real(kp), optional, intent(out) :: bfoldstar

    call initialize_aspic_ptrs(name)
   
    select case (size(params,1))

    case (1)
       aspic_x_rrad = ptronep_x_rrad(params(1),lnRrad,Pstar,bfoldstar)

    case (2)
       aspic_x_rrad = ptrtwop_x_rrad(params(1),params(2) &
            ,lnRrad,Pstar,bfoldstar)

    case default
       stop 'aspic_x_rrad: wrong nparams!'
    end select

  end function aspic_x_rrad



  function aspic_x_rreh(name,params,lnRreh,bfoldstar)
    implicit none
    real(kp) :: aspic_x_rreh
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params
    real(kp), intent(in) :: lnRreh
    real(kp), intent(out), optional :: bfoldstar

    call initialize_aspic_ptrs(name)

    select case (size(params,1))

    case (1)
       aspic_x_rreh = ptronep_x_rreh(params(1),lnRreh,bfoldstar)

    case (2)
       aspic_x_rreh = ptrtwop_x_rreh(params(1),params(2) &
            ,lnRreh,bfoldstar)

    case default
       stop 'aspic_x_rreh: wrong nparams!'
    end select

  end function aspic_x_rreh


  function aspic_norm_potential(name,x,params)
    implicit none
    real(kp) :: aspic_norm_potential
    character(len=*), intent(in) :: name
    real(kp), intent(in) :: x
    real(kp), dimension(:), intent(in) :: params

    call initialize_aspic_ptrs(name)

    select case (size(params,1))

    case (1)
       aspic_norm_potential = ptronep_norm_potential(x &
            ,params(1))

    case (2)
       aspic_norm_potential = ptrtwop_norm_potential(x &
            ,params(1),params(2))

    case default
       stop 'aspic_norm_potential: wrong nparams!'
    end select

  end function aspic_norm_potential


  function aspic_epsilon_one(name,x,params)
    implicit none
    real(kp) :: aspic_epsilon_one
    character(len=*), intent(in) :: name
    real(kp), intent(in) :: x
    real(kp), dimension(:), intent(in) :: params

    call initialize_aspic_ptrs(name)

    select case (size(params,1))

    case (1)
       aspic_epsilon_one = ptronep_epsilon_one(x &
            ,params(1))

    case (2)
       aspic_epsilon_one = ptrtwop_epsilon_one(x &
            ,params(1),params(2))

    case default
       stop 'aspic_epsilon_one: wrong nparams!'
    end select

  end function aspic_epsilon_one


  function aspic_epsilon_two(name,x,params)
    implicit none
    real(kp) :: aspic_epsilon_two
    character(len=*), intent(in) :: name
    real(kp), intent(in) :: x
    real(kp), dimension(:), intent(in) :: params

    call initialize_aspic_ptrs(name)

    select case (size(params,1))

    case (1)
       aspic_epsilon_two = ptronep_epsilon_two(x &
            ,params(1))

    case (2)
       aspic_epsilon_two = ptrtwop_epsilon_two(x &
            ,params(1),params(2))

    case default
       stop 'aspic_epsilon_two: wrong nparams!'
    end select

  end function aspic_epsilon_two



  function aspic_epsilon_three(name,x,params)
    implicit none
    real(kp) :: aspic_epsilon_three
    character(len=*), intent(in) :: name
    real(kp), intent(in) :: x
    real(kp), dimension(:), intent(in) :: params

    call initialize_aspic_ptrs(name)

    select case (size(params,1))

    case (1)
       aspic_epsilon_three = ptronep_epsilon_three(x &
            ,params(1))

    case (2)
       aspic_epsilon_three = ptrtwop_epsilon_three(x &
            ,params(1),params(2))

    case default
       stop 'aspic_epsilon_three: wrong nparams!'
    end select

  end function aspic_epsilon_three



  function aspic_x_endinf(name,params)
    implicit none
    real(kp) :: aspic_x_endinf
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params

    aspic_x_endinf = 0._kp

  end function aspic_x_endinf



end module aspicmodels
