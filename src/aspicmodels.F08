!   This file is part of bayaspic
!
!   Copyright (C) 2013-2023 C. Ringeval
!   
!   bayaspic is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   bayaspic is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with bayaspic.  If not, see <https://www.gnu.org/licenses/>.

module aspicmodels
  use infprec, only : kp
  implicit none

  private

  include 'aspicmodels.h'

!if pointers have been associated to their respective routines
  logical, save :: CheckAspicPtrs = .false.

!if xend is the last model parameter
  logical, save :: xEndIsParam = .false.

!number of aspic parameters including xend
  integer, save :: numAspicParams = 0

!number of derived aspic parameters
  integer, save :: numAspicDerived = 0
  
!non-standard reheating history (see aspicnonstd.F08 module)
  logical, save :: nonStandardReheat = .false.

  procedure(onep_ln_omega4), pointer :: ptronep_ln_omega4 => null()
  
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
  public get_aspic_numderived, check_aspic_xend, initialize_aspic_ptrs

  public aspic_x_rrad, aspic_x_rreh, aspic_x_rhow
  public aspic_norm_potential
  public aspic_epsilon_one, aspic_epsilon_two, aspic_epsilon_three
  public aspic_x_endinf, aspic_ln_rho_endinf

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

    numAspicParams = 0
    numAspicDerived = 0
    nonStandardReheat = .false.

    
!initialize pointers to conformal factor function for non-minimal
!models
    select case (trim(name))

#define NONMIN(fooi) \
    case ('fooi'); \
       ptronep_ln_omega4 => fooi/**/_ln_omega4

#include "aspicnonmin.pp"
#undef NONMIN

    case default
       ptronep_ln_omega4 => null()
       
    end select


!for some models, we could keep trace of some extra parameter in derived params
!    select case (trim(name))
!
!    case ('rclfi1')
!       numAspicDerived = 1
!       
!    end select


    
    
!initialize pointers to Hubble flow and reheating functions for all
!models
    select case (trim(name))

!reminder: cpp is only used in "traditional mode" for fortran
#define ZERO(fooi) \
    case ('fooi'); \
       ptrzerop_x_rrad => fooi/**/_x_rrad ; \
       ptrzerop_x_rreh => fooi/**/_x_rreh ; \
       ptrzerop_x_rhow => fooi/**/_x_star ; \
       ptrzerop_norm_potential => fooi/**/_norm_potential ; \
       ptrzerop_epsilon_one => fooi/**/_epsilon_one ; \
       ptrzerop_epsilon_two => fooi/**/_epsilon_two ; \
       ptrzerop_epsilon_three => fooi/**/_epsilon_three ; \
       ptrzerop_x_endinf => fooi/**/_x_endinf ; \
       xEndIsParam = .false. ; \
       numAspicParams = 0
          
#define ONE(fooi) \
    case ('fooi'); \
       ptronep_x_rrad => fooi/**/_x_rrad ; \
       ptronep_x_rreh => fooi/**/_x_rreh ; \
       ptronep_x_rhow => fooi/**/_x_star ; \
       ptronep_norm_potential => fooi/**/_norm_potential ; \
       ptronep_epsilon_one => fooi/**/_epsilon_one ; \
       ptronep_epsilon_two => fooi/**/_epsilon_two ; \
       ptronep_epsilon_three => fooi/**/_epsilon_three ; \
       xEndIsParam = .false. ; \
       ptronep_x_endinf => fooi/**/_x_endinf ; \
       numAspicParams = 1


#define ONE_XEND(fooi) \
    case ('fooi'); \
       ptrzerop_x_rrad => fooi/**/_x_rrad ; \
       ptrzerop_x_rreh => fooi/**/_x_rreh ; \
       ptrzerop_x_rhow => fooi/**/_x_star ; \
       ptrzerop_norm_potential => fooi/**/_norm_potential ; \
       ptrzerop_epsilon_one => fooi/**/_epsilon_one ; \
       ptrzerop_epsilon_two => fooi/**/_epsilon_two ; \
       ptrzerop_epsilon_three => fooi/**/_epsilon_three ; \
       xEndIsParam = .true. ; \
       ptrzerop_x_endinf => null() ; \
       numAspicParams = 1


#define TWO(fooi) \
    case ('fooi') ; \
       ptrtwop_x_rrad => fooi/**/_x_rrad ; \
       ptrtwop_x_rreh => fooi/**/_x_rreh ; \
       ptrtwop_x_rhow => fooi/**/_x_star ; \
       ptrtwop_norm_potential => fooi/**/_norm_potential ; \
       ptrtwop_epsilon_one => fooi/**/_epsilon_one ; \
       ptrtwop_epsilon_two => fooi/**/_epsilon_two ; \
       ptrtwop_epsilon_three => fooi/**/_epsilon_three ; \
       xEndIsParam = .false. ; \
       ptrtwop_x_endinf => fooi/**/_x_endinf ; \
       numAspicParams = 2

       
#define TWO_XEND(fooi) \
    case ('fooi') ; \
       ptronep_x_rrad => fooi/**/_x_rrad ; \
       ptronep_x_rreh => fooi/**/_x_rreh ; \
       ptronep_x_rhow => fooi/**/_x_star ; \
       ptronep_norm_potential => fooi/**/_norm_potential ; \
       ptronep_epsilon_one => fooi/**/_epsilon_one ; \
       ptronep_epsilon_two => fooi/**/_epsilon_two ; \
       ptronep_epsilon_three => fooi/**/_epsilon_three ; \
       xEndIsParam = .true. ; \
       ptronep_x_endinf => null() ; \
       numAspicParams = 2       


#define THREE(fooi) \
    case ('fooi') ; \
       ptrthreep_x_rrad => fooi/**/_x_rrad ; \
       ptrthreep_x_rreh => fooi/**/_x_rreh ; \
       ptrthreep_x_rhow => fooi/**/_x_star ; \
       ptrthreep_norm_potential => fooi/**/_norm_potential ; \
       ptrthreep_epsilon_one => fooi/**/_epsilon_one ; \
       ptrthreep_epsilon_two => fooi/**/_epsilon_two ; \
       ptrthreep_epsilon_three => fooi/**/_epsilon_three ; \
       xEndIsParam = .false. ; \
       ptrthreep_x_endinf => fooi/**/_x_endinf ; \
       numAspicParams = 3

       
#define THREE_XEND(fooi) \
    case ('fooi') ; \
       ptrtwop_x_rrad => fooi/**/_x_rrad ; \
       ptrtwop_x_rreh => fooi/**/_x_rreh ; \
       ptrtwop_x_rhow => fooi/**/_x_star ; \
       ptrtwop_norm_potential => fooi/**/_norm_potential ; \
       ptrtwop_epsilon_one => fooi/**/_epsilon_one ; \
       ptrtwop_epsilon_two => fooi/**/_epsilon_two ; \
       ptrtwop_epsilon_three => fooi/**/_epsilon_three ; \
       xEndIsParam = .true. ; \
       ptrtwop_x_endinf => null() ; \
       numAspicParams = 3
       
#include "aspicmodels.pp"
#undef ONE
#undef ONE_XEND       
#undef TWO
#undef TWO_XEND
#undef THREE
#undef THREE_XEND

!the remaining case are non-standard reheating histories and require a
!case-by-case treatment
    case ('hi')
       nonStandardReheat = .true.
       xEndIsParam = .false.
       numAspicParams = 0
!we solve for xistar matching the CMB amplitude in aspic
       numAspicDerived = 1
       ptronep_ln_omega4 => hi_parametric_ln_omega4
       
    case ('di')
       nonStandardReheat = .true.
       xEndIsParam = .false.
!this is the dimensionless parameter f
       numAspicParams = 1
       numAspicDerived = 1
       ptronep_ln_omega4 => null()

    case ('nmlfi1')
       nonStandardReheat = .true.
       xEndIsParam = .false.
       numAspicParams = 2
       numAspicDerived = 0
       ptronep_ln_omega4 => nmlfi_parametric_ln_omega4

    case ('nmlfi3')
       nonStandardReheat = .true.
       xEndIsParam = .true.
       numAspicParams = 3
       numAspicDerived = 0
       ptronep_ln_omega4 => nmlfi_parametric_ln_omega4
       
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


  function get_aspic_numderived()
    implicit none
    integer :: get_aspic_numderived

    get_aspic_numderived = numAspicDerived

  end function get_aspic_numderived
  

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

    if (associated(ptronep_ln_omega4)) ptronep_ln_omega4 => null()
    
    CheckAspicPtrs = .false.

  end subroutine free_aspic_ptrs



  function aspic_x_rrad(name,params,xend,lnRrad,Pstar,bfoldstar)
    use aspicnonstd, only : aspic_nonstd_rrad
    implicit none
    real(kp) :: aspic_x_rrad
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params
    real(kp), intent(in) :: lnRrad, Pstar
    real(kp), intent(inout) :: xend
    real(kp), optional, intent(out) :: bfoldstar

    if (.not.CheckAspicPtrs) stop 'aspic not initialized!'
   
    if (associated(ptrzerop_x_rrad)) then
       aspic_x_rrad = ptrzerop_x_rrad(xend,lnRrad,Pstar,bfoldstar)

    elseif (associated(ptronep_x_rrad)) then
       aspic_x_rrad = ptronep_x_rrad(params(1),xend,lnRrad,Pstar,bfoldstar)

    elseif (associated(ptrtwop_x_rrad)) then
       aspic_x_rrad = ptrtwop_x_rrad(params(1),params(2) &
            ,xend,lnRrad,Pstar,bfoldstar)

    elseif (associated(ptrthreep_x_rrad)) then
       aspic_x_rrad = ptrthreep_x_rrad(params(1),params(2) &
            ,params(3),xend,lnRrad,Pstar,bfoldstar)

    elseif (nonStandardReheat) then
       aspic_x_rrad = aspic_nonstd_rrad(name,params,xend,lnRrad,Pstar,bfoldstar)
       
    else
       stop 'aspic_x_rrad: wrong nparams!'
    end if

  end function aspic_x_rrad



  function aspic_x_rreh(name,params,xend,lnRreh,Pstar,bfoldstar)
    use aspicnonstd, only : aspic_nonstd_rreh
    implicit none
    real(kp) :: aspic_x_rreh
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params   
    real(kp), intent(in) :: lnRreh,Pstar
    real(kp), intent(inout) :: xend
    real(kp), intent(out), optional :: bfoldstar

    if (.not.CheckAspicPtrs) stop 'aspic not initialized!'

    if (associated(ptrzerop_x_rreh)) then
       aspic_x_rreh = ptrzerop_x_rreh(xend,lnRreh,bfoldstar)

    elseif (associated(ptronep_x_rreh)) then
       aspic_x_rreh = ptronep_x_rreh(params(1),xend,lnRreh,bfoldstar)

    elseif (associated(ptrtwop_x_rreh)) then
       aspic_x_rreh = ptrtwop_x_rreh(params(1),params(2) &
            ,xend,lnRreh,bfoldstar)
    
    elseif (associated(ptrthreep_x_rreh)) then
       aspic_x_rreh = ptrthreep_x_rreh(params(1),params(2) &
            ,params(3),xend,lnRreh,bfoldstar)

    elseif (nonStandardReheat) then
       aspic_x_rreh = aspic_nonstd_rreh(name,params,xend,lnRreh,Pstar,bfoldstar)
       
    else
       stop 'aspic_x_rreh: wrong nparams!'
    end if

  end function aspic_x_rreh



  function aspic_x_rhow(name,params,xend,w,lnRhoReh,Pstar,bfoldstar)
    use aspicnonstd, only : aspic_nonstd_rhow
    implicit none
    real(kp) :: aspic_x_rhow
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params
    real(kp), intent(in) :: w,lnRhoReh,Pstar
    real(kp), intent(inout) :: xend
    real(kp), intent(out), optional :: bfoldstar

    if (.not.CheckAspicPtrs) stop 'aspic not initialized!'

    if (associated(ptrzerop_x_rhow)) then
       aspic_x_rhow = ptrzerop_x_rhow(xend,w,lnrhoReh,Pstar,bfoldstar)

    elseif (associated(ptronep_x_rhow)) then
       aspic_x_rhow = ptronep_x_rhow(params(1),xend,w,lnRhoReh,Pstar,bfoldstar)

    elseif (associated(ptrtwop_x_rhow)) then
       aspic_x_rhow = ptrtwop_x_rhow(params(1),params(2) &
            ,xend,w,lnRhoReh,Pstar,bfoldstar)
    
    elseif (associated(ptrthreep_x_rhow)) then
       aspic_x_rhow = ptrthreep_x_rhow(params(1),params(2) &
            ,params(3),xend,w,lnRhoReh,Pstar,bfoldstar)

    elseif (nonStandardReheat) then
       aspic_x_rhow = aspic_nonstd_rhow(name,params,xend,w,lnRhoReh,Pstar,bfoldstar)
       
    else
       stop 'aspic_x_rhow: wrong nparams!'
    end if

  end function aspic_x_rhow



  function aspic_norm_potential(name,x,params)
    use aspicnonstd, only : aspic_nonstd_potential
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

    elseif (nonStandardReheat) then
       aspic_norm_potential = aspic_nonstd_potential(name,x,params)
       
    else
       stop 'aspic_norm_potential: wrong nparams!'
    end if

  end function aspic_norm_potential



  function aspic_epsilon_one(name,x,params)
    use aspicnonstd, only : aspic_nonstd_epsilon_one
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

    elseif (nonStandardReheat) then
       aspic_epsilon_one = aspic_nonstd_epsilon_one(name,x,params)
       
    else
       stop 'aspic_epsilon_one: wrong nparams!'
    end if

  end function aspic_epsilon_one


  function aspic_epsilon_two(name,x,params)
    use aspicnonstd, only : aspic_nonstd_epsilon_two
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

    elseif (nonStandardReheat) then
       aspic_epsilon_two = aspic_nonstd_epsilon_two(name,x,params)
              
    else
       stop 'aspic_epsilon_two: wrong nparams!'
    end if


  end function aspic_epsilon_two



  function aspic_epsilon_three(name,x,params)
    use aspicnonstd, only : aspic_nonstd_epsilon_three
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

    elseif (nonStandardReheat) then
       aspic_epsilon_three = aspic_nonstd_epsilon_three(name,x,params)
       
    else
       stop 'aspic_epsilon_three: wrong nparams!'
    end if

   
  end function aspic_epsilon_three



  function aspic_x_endinf(name,params)
    use aspicnonstd, only : aspic_nonstd_endinf
    implicit none
    real(kp) :: aspic_x_endinf
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params
    
    if (.not.CheckAspicPtrs) stop 'aspic not initialized!'

!when xend is present it is the the last one.
       
    if (xEndIsParam) then
       if (size(params,1).lt.numAspicParams) then
          stop 'aspic_x_endinf: time to debug...!'
       endif
       aspic_x_endinf = params(numAspicParams)
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

    elseif (nonStandardReheat) then
       aspic_x_endinf = aspic_nonstd_endinf(name,params)
       
    else
       stop 'aspic_x_endinf: wrong nparams!'
    end if
    
  end function aspic_x_endinf


  
  function aspic_ln_rho_endinf(name,Pstar,xstar,epsStar,Vstar,xend,epsEnd,Vend)
    use srreheat, only : ln_rho_endinf
    implicit none
    real(kp) :: aspic_ln_rho_endinf
    character(len=*), intent(in) :: name
    real(kp), intent(in) :: Pstar,xstar,epsStar,Vstar,xend,epsEnd,Vend
    real(kp) :: lnOmega4End

    if (associated(ptronep_ln_omega4)) then

       lnOmega4End = ptronep_ln_omega4(xend)
       aspic_ln_rho_endinf = ln_rho_endinf(Pstar,epsStar,epsEnd,Vend/Vstar,lnOmega4End)

    else

       aspic_ln_rho_endinf = ln_rho_endinf(Pstar,epsStar,epsEnd,Vend/Vstar)

    end if
       
  end function aspic_ln_rho_endinf

  

end module aspicmodels
