module aspicmodels
  use infprec, only : kp
  implicit none

  private

  include 'aspicmodels.h'

!if pointers have been associated to their respective routines
  logical, save :: CheckAspicPtrs = .false.

!if xend is the last model parameter
  logical, save :: xEndIsParam = .false.
!if xini or xuv is a model parameter
  logical, save :: xIniIsParam = .false.

!number of aspic parameters + xend if required
  integer, save :: numAspicParams = 0

  procedure(zerop_x_rrad), pointer :: ptrzerop_x_rrad => null()
  procedure(zerop_x_rreh), pointer :: ptrzerop_x_rreh => null()
  procedure(zerop_x_rhow), pointer :: ptrzerop_x_rhow => null()
  procedure(onep_x_rrad), pointer :: ptronep_x_rrad => null()
  procedure(onep_x_rreh), pointer :: ptronep_x_rreh => null()
  procedure(onep_x_rhow), pointer :: ptronep_x_rhow => null()
  procedure(twop_x_rrad), pointer :: ptrtwop_x_rrad => null()
  procedure(twop_x_rreh), pointer :: ptrtwop_x_rreh => null()
  procedure(twop_x_rhow), pointer :: ptrtwop_x_rhow => null()
  procedure(threep_x_rrad), pointer :: ptrthreep_x_rrad => null()
  procedure(threep_x_rreh), pointer :: ptrthreep_x_rreh => null()
  procedure(threep_x_rhow), pointer :: ptrthreep_x_rhow => null()
  procedure(fourp_x_rrad), pointer :: ptrfourp_x_rrad => null()
  procedure(fourp_x_rreh), pointer :: ptrfourp_x_rreh => null()
  procedure(fourp_x_rhow), pointer :: ptrfourp_x_rhow => null()

  procedure(zerop_norm_potential), pointer :: ptrzerop_norm_potential => null()
  procedure(zerop_epsilon_one), pointer :: ptrzerop_epsilon_one => null()
  procedure(zerop_epsilon_two), pointer :: ptrzerop_epsilon_two => null()
  procedure(zerop_epsilon_three), pointer :: ptrzerop_epsilon_three => null()
  procedure(zerop_x_endinf), pointer :: ptrzerop_x_endinf => null()

  procedure(onep_norm_potential), pointer :: ptronep_norm_potential => null()
  procedure(onep_epsilon_one), pointer :: ptronep_epsilon_one => null()
  procedure(onep_epsilon_two), pointer :: ptronep_epsilon_two => null()
  procedure(onep_epsilon_three), pointer :: ptronep_epsilon_three => null()
  procedure(onep_x_endinf), pointer :: ptronep_x_endinf => null()
  
  procedure(twop_norm_potential), pointer :: ptrtwop_norm_potential => null()
  procedure(twop_epsilon_one), pointer :: ptrtwop_epsilon_one => null()
  procedure(twop_epsilon_two), pointer :: ptrtwop_epsilon_two => null()
  procedure(twop_epsilon_three), pointer :: ptrtwop_epsilon_three => null()
  procedure(twop_x_endinf), pointer :: ptrtwop_x_endinf => null()

  procedure(threep_norm_potential), pointer :: ptrthreep_norm_potential => null()
  procedure(threep_epsilon_one), pointer :: ptrthreep_epsilon_one => null()
  procedure(threep_epsilon_two), pointer :: ptrthreep_epsilon_two => null()
  procedure(threep_epsilon_three), pointer :: ptrthreep_epsilon_three => null()
  procedure(threep_x_endinf), pointer :: ptrthreep_x_endinf => null()


  public check_aspic_ptrs, free_aspic_ptrs, get_aspic_numparams
  public check_aspic_xend, initialize_aspic_ptrs

  public aspic_x_rrad, aspic_x_rreh, aspic_x_rhow
  public aspic_norm_potential
  public aspic_epsilon_one, aspic_epsilon_two, aspic_epsilon_three
  public aspic_x_endinf

contains


  
  subroutine initialize_aspic_ptrs(name)
#ifdef ASPIC
    include 'aspic.h'
#endif
#ifdef ASPICQ
    include 'aspicq.h'
#endif
    implicit none
    character(len=*), intent(in) :: name
    
    if (CheckAspicPtrs) return

    select case (trim(name))

!reminder: cpp is only used in "traditional mode" for fortran
#define ZEROPP_PURE(fooi) \
    case ('fooi'); \
       ptrzerop_x_rrad => fooi/**/_x_rrad ; \
       ptrzerop_x_rreh => fooi/**/_x_rreh ; \
       ptrzerop_x_rhow => fooi/**/_x_star ; \
       ptrzerop_norm_potential => fooi/**/_norm_potential ; \
       ptrzerop_epsilon_one => fooi/**/_epsilon_one ; \
       ptrzerop_epsilon_two => fooi/**/_epsilon_two ; \
       ptrzerop_epsilon_three => fooi/**/_epsilon_three ; \
       ptrzerop_x_endinf => fooi/**/_x_endinf ; \
       numAspicParams = 0

#define ZEROPP_XEND(fooi) \
    case ('fooi'); \
       ptronep_x_rrad => fooi/**/_x_rrad ; \
       ptronep_x_rreh => fooi/**/_x_rreh ; \
       ptronep_x_rhow => fooi/**/_x_star ; \
       ptrzerop_norm_potential => fooi/**/_norm_potential ; \
       ptrzerop_epsilon_one => fooi/**/_epsilon_one ; \
       ptrzerop_epsilon_two => fooi/**/_epsilon_two ; \
       ptrzerop_epsilon_three => fooi/**/_epsilon_three ; \
       ptrzerop_x_endinf => null() ; \
       xEndIsParam = .true. ; \
       numAspicParams = 1

#define ONEPP_PURE(fooi) \
    case ('fooi'); \
       ptronep_x_rrad => fooi/**/_x_rrad ; \
       ptronep_x_rreh => fooi/**/_x_rreh ; \
       ptronep_x_rhow => fooi/**/_x_star ; \
       ptronep_norm_potential => fooi/**/_norm_potential ; \
       ptronep_epsilon_one => fooi/**/_epsilon_one ; \
       ptronep_epsilon_two => fooi/**/_epsilon_two ; \
       ptronep_epsilon_three => fooi/**/_epsilon_three ; \
       ptronep_x_endinf => fooi/**/_x_endinf ; \
       numAspicParams = 1


#define ONEPP_XEND(fooi) \
    case ('fooi'); \
       ptrtwop_x_rrad => fooi/**/_x_rrad ; \
       ptrtwop_x_rreh => fooi/**/_x_rreh ; \
       ptrtwop_x_rhow => fooi/**/_x_star ; \
       ptronep_norm_potential => fooi/**/_norm_potential ; \
       ptronep_epsilon_one => fooi/**/_epsilon_one ; \
       ptronep_epsilon_two => fooi/**/_epsilon_two ; \
       ptronep_epsilon_three => fooi/**/_epsilon_three ;\
       ptronep_x_endinf => null() ; \
       xEndIsParam = .true. ;\
       numAspicParams = 2

#define TWOPP_PURE(fooi) \
    case ('fooi') ; \
       ptrtwop_x_rrad => fooi/**/_x_rrad ; \
       ptrtwop_x_rreh => fooi/**/_x_rreh ; \
       ptrtwop_x_rhow => fooi/**/_x_star ; \
       ptrtwop_norm_potential => fooi/**/_norm_potential ; \
       ptrtwop_epsilon_one => fooi/**/_epsilon_one ; \
       ptrtwop_epsilon_two => fooi/**/_epsilon_two ; \
       ptrtwop_epsilon_three => fooi/**/_epsilon_three; \
       ptrtwop_x_endinf => fooi/**/_x_endinf; \
       xEndIsParam = .false. ;\
       numAspicParams = 2

#define TWOPP_XEND(fooi) \
    case ('fooi') ; \
       ptrthreep_x_rrad => fooi/**/_x_rrad ; \
       ptrthreep_x_rreh => fooi/**/_x_rreh ; \
       ptrthreep_x_rhow => fooi/**/_x_star ; \
       ptrtwop_norm_potential => fooi/**/_norm_potential ; \
       ptrtwop_epsilon_one => fooi/**/_epsilon_one ; \
       ptrtwop_epsilon_two => fooi/**/_epsilon_two ; \
       ptrtwop_epsilon_three => fooi/**/_epsilon_three ;\
       ptrtwop_x_endinf => null() ; \
       xEndIsParam = .true. ;\
       numAspicParams = 3

#define TWOPP_XEND_XINI(fooi) \
    case ('fooi') ; \
       ptrthreep_x_rrad => fooi/**/_x_rrad ; \
       ptrthreep_x_rreh => fooi/**/_x_rreh ; \
       ptrthreep_x_rhow => fooi/**/_x_star ; \
       ptrtwop_norm_potential => fooi/**/_norm_potential ; \
       ptrtwop_epsilon_one => fooi/**/_epsilon_one ; \
       ptrtwop_epsilon_two => fooi/**/_epsilon_two ; \
       ptrtwop_epsilon_three => fooi/**/_epsilon_three ;\
       ptrtwop_x_endinf => null() ; \
       xEndIsParam = .true. ;\
       xIniIsParam = .true. ;\
       numAspicParams = 4


#define THREEPP_PURE(fooi) \
    case ('fooi') ; \
       ptrthreep_x_rrad => fooi/**/_x_rrad ; \
       ptrthreep_x_rreh => fooi/**/_x_rreh ; \
       ptrthreep_x_rhow => fooi/**/_x_star ; \
       ptrthreep_norm_potential => fooi/**/_norm_potential ; \
       ptrthreep_epsilon_one => fooi/**/_epsilon_one ; \
       ptrthreep_epsilon_two => fooi/**/_epsilon_two ; \
       ptrthreep_epsilon_three => fooi/**/_epsilon_three; \
       ptrthreep_x_endinf => fooi/**/_x_endinf; \
       xEndIsParam = .false. ;\
       numAspicParams = 3

#define THREEPP_XEND(fooi) \
    case ('fooi') ; \
       ptrfourp_x_rrad => fooi/**/_x_rrad ; \
       ptrfourp_x_rreh => fooi/**/_x_rreh ; \
       ptrfourp_x_rhow => fooi/**/_x_star ; \
       ptrfourp_norm_potential => fooi/**/_norm_potential ; \
       ptrthreep_epsilon_one => fooi/**/_epsilon_one ; \
       ptrthreep_epsilon_two => fooi/**/_epsilon_two ; \
       ptrthreep_epsilon_three => fooi/**/_epsilon_three ;\
       ptrthreep_x_endinf => null() ; \
       xEndIsParam = .true. ;\
       numAspicParams = 4


#include "aspicmodels.pp"
#undef ZEROPP_PURE
#undef ZEROPP_XEND
#undef ONEPP_PURE
#undef ONEPP_XEND
#undef TWOPP_PURE
#undef TWOPP_XEND
#undef TWOPP_XEND_XINI
#undef THREEPP_PURE
#undef THREEPP_XEND


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


  function check_aspic_xend()
    implicit none
    logical :: check_aspic_xend

    check_aspic_xend = xEndIsParam

  end function check_aspic_xend


  function get_aspic_numparams()
    implicit none
    integer :: get_aspic_numparams

    get_aspic_numparams = numAspicParams

  end function get_aspic_numparams


  subroutine free_aspic_ptrs()
    implicit none

    if (.not.CheckAspicPtrs) return

    if (associated(ptrzerop_norm_potential)) ptrzerop_norm_potential => null()
    if (associated(ptrzerop_epsilon_one)) ptrzerop_epsilon_one => null()
    if (associated(ptrzerop_epsilon_two)) ptrzerop_epsilon_two => null()
    if (associated(ptrzerop_epsilon_three)) ptrzerop_epsilon_three => null()
    if (associated(ptrzerop_x_endinf)) ptrzerop_x_endinf => null()

    if (associated(ptronep_norm_potential)) ptronep_norm_potential => null()
    if (associated(ptronep_epsilon_one)) ptronep_epsilon_one => null()
    if (associated(ptronep_epsilon_two)) ptronep_epsilon_two => null()
    if (associated(ptronep_epsilon_three)) ptronep_epsilon_three => null()
    if (associated(ptronep_x_endinf)) ptronep_x_endinf => null()

    if (associated(ptrtwop_norm_potential)) ptrtwop_norm_potential => null()
    if (associated(ptrtwop_epsilon_one)) ptrtwop_epsilon_one => null()
    if (associated(ptrtwop_epsilon_two)) ptrtwop_epsilon_two => null()
    if (associated(ptrtwop_epsilon_three)) ptrtwop_epsilon_three => null()
    if (associated(ptrtwop_x_endinf)) ptrtwop_x_endinf => null()

    if (associated(ptrthreep_norm_potential)) ptrthreep_norm_potential => null()
    if (associated(ptrthreep_epsilon_one)) ptrthreep_epsilon_one => null()
    if (associated(ptrthreep_epsilon_two)) ptrthreep_epsilon_two => null()
    if (associated(ptrthreep_epsilon_three)) ptrthreep_epsilon_three => null()
    if (associated(ptrthreep_x_endinf)) ptrthreep_x_endinf => null()

    if (associated(ptrzerop_x_rrad)) ptrzerop_x_rrad => null()    
    if (associated(ptrzerop_x_rreh)) ptrzerop_x_rreh => null()
    if (associated(ptrzerop_x_rhow)) ptrzerop_x_rhow => null()
    if (associated(ptronep_x_rrad)) ptronep_x_rrad => null()    
    if (associated(ptronep_x_rreh)) ptronep_x_rreh => null()
    if (associated(ptronep_x_rhow)) ptronep_x_rhow => null()
    if (associated(ptrtwop_x_rrad)) ptrtwop_x_rrad => null()
    if (associated(ptrtwop_x_rreh)) ptrtwop_x_rreh => null()
    if (associated(ptrtwop_x_rhow)) ptrtwop_x_rhow => null()
    if (associated(ptrthreep_x_rrad)) ptrthreep_x_rrad => null()
    if (associated(ptrthreep_x_rreh)) ptrthreep_x_rreh => null()
    if (associated(ptrthreep_x_rhow)) ptrthreep_x_rhow => null()
    if (associated(ptrfourp_x_rrad)) ptrfourp_x_rrad => null()
    if (associated(ptrfourp_x_rreh)) ptrfourp_x_rreh => null()
    if (associated(ptrfourp_x_rhow)) ptrfourp_x_rhow => null()

    CheckAspicPtrs = .false.

  end subroutine free_aspic_ptrs



  function aspic_x_rrad(name,params,lnRrad,Pstar,bfoldstar)
    implicit none
    real(kp) :: aspic_x_rrad
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params
    real(kp), intent(in) :: lnRrad, Pstar
    real(kp), optional, intent(out) :: bfoldstar

    if (.not.CheckAspicPtrs) stop 'aspic not initialized!'
   
    if (associated(ptrzerop_x_rrad)) then
       aspic_x_rrad = ptrzerop_x_rrad(lnRrad,Pstar,bfoldstar)

    elseif (associated(ptronep_x_rrad)) then
       aspic_x_rrad = ptronep_x_rrad(params(1),lnRrad,Pstar,bfoldstar)

    elseif (associated(ptrtwop_x_rrad)) then
       aspic_x_rrad = ptrtwop_x_rrad(params(1),params(2) &
            ,lnRrad,Pstar,bfoldstar)

    elseif (associated(ptrthreep_x_rrad)) then
       aspic_x_rrad = ptrthreep_x_rrad(params(1),params(2) &
            ,params(3),lnRrad,Pstar,bfoldstar)
    
    elseif (associated(ptrfourp_x_rrad)) then
       aspic_x_rrad = ptrfourp_x_rrad(params(1),params(2) &
            ,params(3),params(4),lnRrad,Pstar,bfoldstar)

    else
       stop 'aspic_x_rrad: wrong nparams!'
    end if

  end function aspic_x_rrad



  function aspic_x_rreh(name,params,lnRreh,bfoldstar)
    implicit none
    real(kp) :: aspic_x_rreh
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params
    real(kp), intent(in) :: lnRreh
    real(kp), intent(out), optional :: bfoldstar

    if (.not.CheckAspicPtrs) stop 'aspic not initialized!'

    if (associated(ptrzerop_x_rreh)) then
       aspic_x_rreh = ptrzerop_x_rreh(lnRreh,bfoldstar)

    elseif (associated(ptronep_x_rreh)) then
       aspic_x_rreh = ptronep_x_rreh(params(1),lnRreh,bfoldstar)

    elseif (associated(ptrtwop_x_rreh)) then
       aspic_x_rreh = ptrtwop_x_rreh(params(1),params(2) &
            ,lnRreh,bfoldstar)
    
    elseif (associated(ptrthreep_x_rreh)) then
       aspic_x_rreh = ptrthreep_x_rreh(params(1),params(2) &
            ,params(3),lnRreh,bfoldstar)

    elseif (associated(ptrfourp_x_rreh)) then
       aspic_x_rreh = ptrfourp_x_rreh(params(1),params(2) &
            ,params(3),params(4),lnRreh,bfoldstar)

    else
       stop 'aspic_x_rreh: wrong nparams!'
    end if

  end function aspic_x_rreh



  function aspic_x_rhow(name,params,lnRhoReh,w,Pstar,bfoldstar)
    implicit none
    real(kp) :: aspic_x_rhow
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params
    real(kp), intent(in) :: lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    if (.not.CheckAspicPtrs) stop 'aspic not initialized!'

    if (associated(ptrzerop_x_rhow)) then
       aspic_x_rhow = ptrzerop_x_rhow(w,lnrhoReh,Pstar,bfoldstar)

    elseif (associated(ptronep_x_rhow)) then
       aspic_x_rhow = ptronep_x_rhow(params(1),w,lnRhoReh,Pstar,bfoldstar)

    elseif (associated(ptrtwop_x_rhow)) then
       aspic_x_rhow = ptrtwop_x_rhow(params(1),params(2) &
            ,w,lnRhoReh,Pstar,bfoldstar)
    
    elseif (associated(ptrthreep_x_rhow)) then
       aspic_x_rhow = ptrthreep_x_rhow(params(1),params(2) &
            ,params(3),w,lnRhoReh,Pstar,bfoldstar)

    elseif (associated(ptrfourp_x_rhow)) then
       aspic_x_rhow = ptrfourp_x_rhow(params(1),params(2) &
            ,params(3),params(4),w,lnRhoReh,Pstar,bfoldstar)

    else
       stop 'aspic_x_rhow: wrong nparams!'
    end if

  end function aspic_x_rhow



  function aspic_norm_potential(name,x,params)
    implicit none
    real(kp) :: aspic_norm_potential
    character(len=*), intent(in) :: name
    real(kp), intent(in) :: x
    real(kp), dimension(:), intent(in) :: params

    if (.not.CheckAspicPtrs) stop 'aspic not initialized!'

    if (associated(ptrzerop_norm_potential)) then

       aspic_norm_potential = ptrzerop_norm_potential(x)

    elseif (associated(ptronep_norm_potential)) then

       aspic_norm_potential = ptronep_norm_potential(x &
            ,params(1))
       
    elseif (associated(ptrtwop_norm_potential)) then
       
       aspic_norm_potential = ptrtwop_norm_potential(x &
            ,params(1),params(2))

    elseif (associated(ptrthreep_norm_potential)) then
       
       aspic_norm_potential = ptrthreep_norm_potential(x &
            ,params(1),params(2),params(3))

    else
       stop 'aspic_norm_potential: wrong nparams!'
    end if

  end function aspic_norm_potential



  function aspic_epsilon_one(name,x,params)
    implicit none
    real(kp) :: aspic_epsilon_one
    character(len=*), intent(in) :: name
    real(kp), intent(in) :: x
    real(kp), dimension(:), intent(in) :: params

    if (.not.CheckAspicPtrs) stop 'aspic not initialized!'

    if (associated(ptrzerop_epsilon_one)) then
       
       aspic_epsilon_one = ptrzerop_epsilon_one(x)

    elseif (associated(ptronep_epsilon_one)) then
       
       aspic_epsilon_one = ptronep_epsilon_one(x &
            ,params(1))

    elseif (associated(ptrtwop_epsilon_one)) then
       aspic_epsilon_one = ptrtwop_epsilon_one(x &
            ,params(1),params(2))

    elseif (associated(ptrthreep_epsilon_one)) then
       aspic_epsilon_one = ptrthreep_epsilon_one(x &
            ,params(1),params(2),params(3))

    else
       stop 'aspic_epsilon_one: wrong nparams!'
    end if

  end function aspic_epsilon_one


  function aspic_epsilon_two(name,x,params)
    implicit none
    real(kp) :: aspic_epsilon_two
    character(len=*), intent(in) :: name
    real(kp), intent(in) :: x
    real(kp), dimension(:), intent(in) :: params

    if (.not.CheckAspicPtrs) stop 'aspic not initialized!'

    if (associated(ptrzerop_epsilon_two)) then
       
       aspic_epsilon_two = ptrzerop_epsilon_two(x)

    elseif (associated(ptronep_epsilon_two)) then
       
       aspic_epsilon_two = ptronep_epsilon_two(x &
            ,params(1))

    elseif (associated(ptrtwop_epsilon_two)) then
       aspic_epsilon_two = ptrtwop_epsilon_two(x &
            ,params(1),params(2))

    elseif (associated(ptrthreep_epsilon_two)) then
       aspic_epsilon_two = ptrthreep_epsilon_two(x &
            ,params(1),params(2),params(3))

    else
       stop 'aspic_epsilon_two: wrong nparams!'
    end if


  end function aspic_epsilon_two



  function aspic_epsilon_three(name,x,params)
    implicit none
    real(kp) :: aspic_epsilon_three
    character(len=*), intent(in) :: name
    real(kp), intent(in) :: x
    real(kp), dimension(:), intent(in) :: params
    
    if (.not.CheckAspicPtrs) stop 'aspic not initialized!'

    if (associated(ptrzerop_epsilon_three)) then
       
       aspic_epsilon_three = ptrzerop_epsilon_three(x)

    elseif (associated(ptronep_epsilon_three)) then
       
       aspic_epsilon_three = ptronep_epsilon_three(x &
            ,params(1))

    elseif (associated(ptrtwop_epsilon_three)) then
       aspic_epsilon_three = ptrtwop_epsilon_three(x &
            ,params(1),params(2))

    elseif (associated(ptrthreep_epsilon_three)) then
       aspic_epsilon_three = ptrthreep_epsilon_three(x &
            ,params(1),params(2),params(3))

    else
       stop 'aspic_epsilon_three: wrong nparams!'
    end if

   
  end function aspic_epsilon_three



  function aspic_x_endinf(name,params)
    implicit none
    real(kp) :: aspic_x_endinf
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params
    
    if (.not.CheckAspicPtrs) stop 'aspic not initialized!'

!if xini is a parameter, it is the last one and xend the before to
!last. Otherwise, xend is the last one.
       
    if (xEndIsParam) then
       if (xIniIsParam) then
          aspic_x_endinf = params(size(params,1)-1)
       else
          aspic_x_endinf = params(size(params,1))
       endif
       return
    endif

    if (associated(ptrzerop_x_endinf)) then
       
       aspic_x_endinf = ptrzerop_x_endinf()

    elseif (associated(ptronep_x_endinf)) then
       
       aspic_x_endinf = ptronep_x_endinf(params(1))

    elseif (associated(ptrtwop_x_endinf)) then
       
       aspic_x_endinf = ptrtwop_x_endinf(params(1),params(2))

    elseif (associated(ptrthreep_x_endinf)) then

       aspic_x_endinf = ptrthreep_x_endinf(params(1),params(2),params(3))

    else
       stop 'aspic_x_endinf: wrong nparams!'
    end if
    
  end function aspic_x_endinf



end module aspicmodels
