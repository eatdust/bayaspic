  abstract interface

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
     function zerop_x_rrad(lnRrad,Pstar,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: zerop_x_rrad
       real(kp), intent(in) :: lnRrad,Pstar
       real(kp), intent(out), optional :: bfoldstar
     end function zerop_x_rrad

     function zerop_x_rreh(lnRreh,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: zerop_x_rreh
       real(kp), intent(in) :: lnRreh
       real(kp), intent(out), optional :: bfoldstar
     end function zerop_x_rreh


!for one param reheat
     function onep_x_rrad(p1,lnRrad,Pstar,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: onep_x_rrad
       real(kp), intent(in) :: p1,lnRrad,Pstar
       real(kp), intent(out), optional :: bfoldstar
     end function onep_x_rrad

     function onep_x_rreh(p1,lnRreh,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: onep_x_rreh
       real(kp), intent(in) :: p1,lnRreh
       real(kp), intent(out), optional :: bfoldstar
     end function onep_x_rreh

!for two params reheat
    function twop_x_rrad(p1,p2,lnRrad,Pstar,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: twop_x_rrad
       real(kp), intent(in) :: p1,p2,lnRrad,Pstar
       real(kp), intent(out), optional :: bfoldstar
     end function twop_x_rrad

     function twop_x_rreh(p1,p2,lnRreh,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: twop_x_rreh
       real(kp), intent(in) :: p1,p2,lnRreh
       real(kp), intent(out), optional :: bfoldstar
     end function twop_x_rreh

!for three params reheat
    function threep_x_rrad(p1,p2,p3,lnRrad,Pstar,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: threep_x_rrad
       real(kp), intent(in) :: p1,p2,p3,lnRrad,Pstar
       real(kp), intent(out), optional :: bfoldstar
     end function threep_x_rrad

     function threep_x_rreh(p1,p2,p3,lnRreh,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: threep_x_rreh
       real(kp), intent(in) :: p1,p2,p3,lnRreh
       real(kp), intent(out), optional :: bfoldstar
     end function threep_x_rreh

!for four params reheat
    function fourp_x_rrad(p1,p2,p3,p4,lnRrad,Pstar,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: fourp_x_rrad
       real(kp), intent(in) :: p1,p2,p3,p4,lnRrad,Pstar
       real(kp), intent(out), optional :: bfoldstar
     end function fourp_x_rrad

     function fourp_x_rreh(p1,p2,p3,p4,lnRreh,bfoldstar)
       use infprec, only : kp
       implicit none
       real(kp) :: fourp_x_rreh
       real(kp), intent(in) :: p1,p2,p3,p4,lnRreh
       real(kp), intent(out), optional :: bfoldstar
     end function fourp_x_rreh

  end interface
