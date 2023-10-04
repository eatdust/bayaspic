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

module aspicnonstd
  use infprec, only : kp
  use aspicvars, only : AspicModel
#ifdef ASPIC
    include 'aspic.h'
#endif
#ifdef ASPICQ
    include 'aspicq.h'
#endif
    
  implicit none

  private

  public aspic_nonstd_rrad, aspic_nonstd_rreh, aspic_nonstd_rhow
  public aspic_nonstd_endinf, aspic_nonstd_ln_rho_endinf, aspic_nonstd_potential
  public aspic_nonstd_epsilon_one, aspic_nonstd_epsilon_two, aspic_nonstd_epsilon_three
  
contains

  
  function aspic_nonstd_rrad(name,params,xend,lnRrad,Pstar,bfoldstar)
    implicit none
    real(kp) :: aspic_nonstd_rrad
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params
    real(kp), intent(in) :: lnRrad, Pstar
    real(kp), intent(inout) :: xend
    real(kp), optional, intent(out) :: bfoldstar

    select case(trim(name))

    case('hi')
!we output the parametric field hbarstar and the first derived aspic
!parameter is updated to xistar. We also update xend because it depends on xistar
       aspic_nonstd_rrad = hi_hbar_rrad(lnRrad,Pstar,bfoldstar,AspicModel%derived(1))
       xend = hi_hbar_endinf(AspicModel%derived(1))

    case default
       stop 'aspic_nonstd_rrad: name not found!'
       
    end select
    
    

  end function aspic_nonstd_rrad

  

  function aspic_nonstd_rreh(name,params,xend,lnRreh,Pstar,bfoldstar)
    implicit none
    real(kp) :: aspic_nonstd_rreh
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params   
    real(kp), intent(inout) :: xend
    real(kp), intent(in) :: lnRreh,Pstar
    real(kp), intent(out), optional :: bfoldstar

    select case(trim(name))

    case('hi')
!we output the parametric field hbarstar and the first derived aspic
!parameter is updated to xistar. We also update xend because it depends on xistar
       aspic_nonstd_rreh = hi_hbar_rreh(lnRreh,Pstar,bfoldstar,AspicModel%derived(1))
       xend = hi_hbar_endinf(AspicModel%derived(1))

    case default
       stop 'aspic_nonstd_rreh: name not found!'
       
    end select
    
    
  end function aspic_nonstd_rreh


  function aspic_nonstd_rhow(name,params,xend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: aspic_nonstd_rhow
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params
    real(kp), intent(in) :: w,lnRhoReh,Pstar
    real(kp), intent(inout) :: xend
    real(kp), intent(out), optional :: bfoldstar


    select case(trim(name))

    case('hi')
!we output the parametric field hbarstar and the first derived aspic
!parameter is updated to xistar. We also update xend because it depends on xistar
       aspic_nonstd_rhow = hi_hbar_star(w,lnRhoReh,Pstar,bfoldstar,AspicModel%derived(1))
       xend = hi_hbar_endinf(AspicModel%derived(1))

    case default
       stop 'aspic_nonstd_rhow: name not found!'
       
    end select

    
  end function aspic_nonstd_rhow


  
  function aspic_nonstd_endinf(name,params)
    implicit none
    real(kp) :: aspic_nonstd_endinf
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params

    select case(trim(name))

    case('hi')
!cannot be determined without solving the reheating equation. Let's return junk inside
       aspic_nonstd_endinf = 0._kp

    case default
       stop 'aspic_nonstd_rreh: name not found!'
       
    end select
       
       
  end function aspic_nonstd_endinf
    
  

   function aspic_nonstd_potential(name,x,params)
    implicit none
    real(kp) :: aspic_nonstd_potential
    character(len=*), intent(in) :: name
    real(kp), intent(in) :: x
    real(kp), dimension(:), intent(in) :: params
  
    select case(trim(name))

    case('hi')
!we are using the parametric field hbar in the computations (x should
!be hbar on input and the first derived parameter should have been
!computed to xistar before calling this function)
       aspic_nonstd_potential = hi_norm_parametric_potential(x,AspicModel%derived(1))

    case default
       stop 'aspic_nonstd_rreh: name not found!'
       
    end select

    
  end function aspic_nonstd_potential

  

  function aspic_nonstd_epsilon_one(name,x,params)
    implicit none
    real(kp) :: aspic_nonstd_epsilon_one
    character(len=*), intent(in) :: name
    real(kp), intent(in) :: x
    real(kp), dimension(:), intent(in) :: params

    select case(trim(name))

    case('hi')
!we are using the parametric field hbar in the computations (x should
!be hbar on input and the first derived parameter should have been
!computed to xistar before calling this function)
       aspic_nonstd_epsilon_one = hi_parametric_epsilon_one(x,AspicModel%derived(1))

    case default
       stop 'aspic_nonstd_rreh: name not found!'
       
    end select


    
  end function aspic_nonstd_epsilon_one

  

  function aspic_nonstd_epsilon_two(name,x,params)
    implicit none
    real(kp) :: aspic_nonstd_epsilon_two
    character(len=*), intent(in) :: name
    real(kp), intent(in) :: x
    real(kp), dimension(:), intent(in) :: params

    select case(trim(name))
    case('hi')
!we are using the parametric field hbar in the computations (x should
!be hbar on input and the first derived parameter should have been
!computed to xistar before calling this function)
       aspic_nonstd_epsilon_two = hi_parametric_epsilon_two(x,AspicModel%derived(1))

    case default
       stop 'aspic_nonstd_rreh: name not found!'
       
    end select

    
  end function aspic_nonstd_epsilon_two



  function aspic_nonstd_epsilon_three(name,x,params)
    implicit none
    real(kp) :: aspic_nonstd_epsilon_three
    character(len=*), intent(in) :: name
    real(kp), intent(in) :: x
    real(kp), dimension(:), intent(in) :: params

    select case(trim(name))
    case('hi')
!we are using the parametric field hbar in the computations (x should
!be hbar on input and the first derived parameter should have been
!computed to xistar before calling this function)
       aspic_nonstd_epsilon_three = hi_parametric_epsilon_three(x,AspicModel%derived(1))

    case default
       stop 'aspic_nonstd_rreh: name not found!'
       
    end select    
    
  end function aspic_nonstd_epsilon_three


  
  function aspic_nonstd_ln_rho_endinf(name,Pstar,xstar,epsStar,Vstar,xend,epsEnd,Vend)
    use srreheat, only : ln_rho_endinf
    implicit none
    real(kp) :: aspic_nonstd_ln_rho_endinf
    character(len=*), intent(in) :: name
    real(kp), intent(in) :: Pstar,xstar,epsStar,Vstar
    real(kp), intent(in) :: xend,epsEnd,Vend

    real(kp) :: lnOmega4End
    
    select case(trim(name))

    case('hi')
!we are using the parametric field hbar in the computations (x should
!be hbar on input and the first derived parameter should have been
!computed to xistar before calling this function)

       lnOmega4End = 2._kp*log(1._kp + xend*xend)
       
       aspic_nonstd_ln_rho_endinf = ln_rho_endinf(Pstar,epsStar,epsEnd,Vend/Vstar,lnOmega4End)

    case default
       stop 'aspic_nonstd_rreh: name not found!'
       
    end select    
    
  end function aspic_nonstd_ln_rho_endinf


  
end module aspicnonstd
