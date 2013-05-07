  abstract interface
!for one potential parameter models
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

  end interface
