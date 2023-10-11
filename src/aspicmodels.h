!   This file is part of bayaspic
!
!   Copyright (C) 2013-2021 C. Ringeval
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

  abstract interface

!ln[1/A^4] for non-minimal models
     function onep_ln_omega4(x)
       use infprec, only : kp
       implicit none
       real(kp) :: onep_ln_omega4
       real(kp), intent(in) :: x
     end function onep_ln_omega4
     
!for zero potential parameter models
     function zerop_norm_potential(x)
       use infprec, only : kp
       implicit none
       real(kp) :: zerop_norm_potential
       real(kp), intent(in) :: x
     end function zerop_norm_potential

     function zerop_epsilon_one(x)    
       use infprec, only : kp
       implicit none
       real(kp) :: zerop_epsilon_one
       real(kp), intent(in) :: x
     end function zerop_epsilon_one

     function zerop_epsilon_two(x)    
       use infprec, only : kp
       implicit none
       real(kp) :: zerop_epsilon_two
       real(kp), intent(in) :: x
     end function zerop_epsilon_two

     function zerop_epsilon_three(x)    
       use infprec, only : kp
       implicit none
       real(kp) :: zerop_epsilon_three
       real(kp), intent(in) :: x
     end function zerop_epsilon_three

     function zerop_x_endinf()
       use infprec, only : kp
       implicit none
       real(kp) :: zerop_x_endinf
     end function zerop_x_endinf


!for one potential parameter models
     function onep_norm_potential(x,p1)
       use infprec, only : kp
       implicit none
       real(kp) :: onep_norm_potential
       real(kp), intent(in) :: x,p1
     end function onep_norm_potential

     function onep_epsilon_one(x,p1)    
       use infprec, only : kp
       implicit none
       real(kp) :: onep_epsilon_one
       real(kp), intent(in) :: x,p1
     end function onep_epsilon_one

     function onep_epsilon_two(x,p1)    
       use infprec, only : kp
       implicit none
       real(kp) :: onep_epsilon_two
       real(kp), intent(in) :: x,p1
     end function onep_epsilon_two

     function onep_epsilon_three(x,p1)    
       use infprec, only : kp
       implicit none
       real(kp) :: onep_epsilon_three
       real(kp), intent(in) :: x,p1
     end function onep_epsilon_three

     function onep_x_endinf(p1)
       use infprec, only : kp
       implicit none
       real(kp), intent(in) :: p1
       real(kp) :: onep_x_endinf
     end function onep_x_endinf


!for two potential parameters models
     function twop_norm_potential(x,p1,p2)
       use infprec, only : kp
       implicit none
       real(kp) :: twop_norm_potential
       real(kp), intent(in) :: x,p1,p2
     end function twop_norm_potential
     
     function twop_epsilon_one(x,p1,p2)    
       use infprec, only : kp
       implicit none
       real(kp) :: twop_epsilon_one
       real(kp), intent(in) :: x,p1,p2
     end function twop_epsilon_one

     function twop_epsilon_two(x,p1,p2)    
       use infprec, only : kp
       implicit none
       real(kp) :: twop_epsilon_two
       real(kp), intent(in) :: x,p1,p2
     end function twop_epsilon_two

     function twop_epsilon_three(x,p1,p2)    
       use infprec, only : kp
       implicit none
       real(kp) :: twop_epsilon_three
       real(kp), intent(in) :: x,p1,p2
     end function twop_epsilon_three

     function twop_x_endinf(p1,p2)
       use infprec, only : kp
       implicit none
       real(kp), intent(in) :: p1,p2
       real(kp) :: twop_x_endinf
     end function twop_x_endinf

!three potential parameters
     function threep_norm_potential(x,p1,p2,p3)
       use infprec, only : kp
       implicit none
       real(kp) :: threep_norm_potential
       real(kp), intent(in) :: x,p1,p2,p3
     end function threep_norm_potential
     
     function threep_epsilon_one(x,p1,p2,p3)    
       use infprec, only : kp
       implicit none
       real(kp) :: threep_epsilon_one
       real(kp), intent(in) :: x,p1,p2,p3
     end function threep_epsilon_one

     function threep_epsilon_two(x,p1,p2,p3)    
       use infprec, only : kp
       implicit none
       real(kp) :: threep_epsilon_two
       real(kp), intent(in) :: x,p1,p2,p3
     end function threep_epsilon_two

     function threep_epsilon_three(x,p1,p2,p3)    
       use infprec, only : kp
       implicit none
       real(kp) :: threep_epsilon_three
       real(kp), intent(in) :: x,p1,p2,p3
     end function threep_epsilon_three

     function threep_x_endinf(p1,p2,p3)
       use infprec, only : kp
       implicit none
       real(kp), intent(in) :: p1,p2,p3
       real(kp) :: threep_x_endinf
     end function threep_x_endinf


!for zero param reheat
     function zerop_x_rrad(xend,lnRrad,Pstar,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: zerop_x_rrad
       real(kp), intent(in) :: xend, lnRrad,Pstar
       real(kp), intent(out), optional :: bfoldstar
     end function zerop_x_rrad

     function zerop_x_rreh(xend,lnRreh,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: zerop_x_rreh
       real(kp), intent(in) :: xend, lnRreh
       real(kp), intent(out), optional :: bfoldstar
     end function zerop_x_rreh

     function zerop_x_rhow(xend,w,lnRhoReh,Pstar,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: zerop_x_rhow
       real(kp), intent(in) :: xend, w,lnRhoReh,Pstar
       real(kp), intent(out), optional :: bfoldstar
     end function zerop_x_rhow


!for one param reheat
     function onep_x_rrad(p1,xend,lnRrad,Pstar,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: onep_x_rrad
       real(kp), intent(in) :: p1,xend,lnRrad,Pstar
       real(kp), intent(out), optional :: bfoldstar
     end function onep_x_rrad

     function onep_x_rreh(p1,xend,lnRreh,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: onep_x_rreh
       real(kp), intent(in) :: p1,xend,lnRreh
       real(kp), intent(out), optional :: bfoldstar
     end function onep_x_rreh

     function onep_x_rhow(p1,xend,w,lnRhoReh,Pstar,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: onep_x_rhow
       real(kp), intent(in) :: p1,xend,w,lnRhoReh,Pstar
       real(kp), intent(out), optional :: bfoldstar
     end function onep_x_rhow

!for two params reheat
    function twop_x_rrad(p1,p2,xend,lnRrad,Pstar,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: twop_x_rrad
       real(kp), intent(in) :: p1,p2,xend,lnRrad,Pstar
       real(kp), intent(out), optional :: bfoldstar
     end function twop_x_rrad

     function twop_x_rreh(p1,p2,xend,lnRreh,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: twop_x_rreh
       real(kp), intent(in) :: p1,p2,xend,lnRreh
       real(kp), intent(out), optional :: bfoldstar
     end function twop_x_rreh

     function twop_x_rhow(p1,p2,xend,w,lnRhoReh,Pstar,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: twop_x_rhow
       real(kp), intent(in) :: p1,p2,xend,w,lnRhoReh,Pstar
       real(kp), intent(out), optional :: bfoldstar
     end function twop_x_rhow

!for three params reheat
    function threep_x_rrad(p1,p2,p3,xend,lnRrad,Pstar,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: threep_x_rrad
       real(kp), intent(in) :: p1,p2,p3,xend,lnRrad,Pstar
       real(kp), intent(out), optional :: bfoldstar
     end function threep_x_rrad

     function threep_x_rreh(p1,p2,p3,xend,lnRreh,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: threep_x_rreh
       real(kp), intent(in) :: p1,p2,p3,xend,lnRreh
       real(kp), intent(out), optional :: bfoldstar
     end function threep_x_rreh

     function threep_x_rhow(p1,p2,p3,xend,w,lnRhoReh,Pstar,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: threep_x_rhow
       real(kp), intent(in) :: p1,p2,p3,xend,w,lnRhoReh,Pstar
       real(kp), intent(out), optional :: bfoldstar
     end function threep_x_rhow

     
  end interface
