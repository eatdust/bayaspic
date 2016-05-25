module aspicvars  
  use infprec, only : kp, pi
  implicit none

  private

!  integer, parameter :: kp = kind(1._8)

  integer, parameter :: lname = 12

  type infaspic     
     character(len=lname) :: name     
     character(len=lname) :: extname
     integer :: nasp
     real(kp), dimension(:), pointer :: params
     character(len=lname), dimension(:), pointer :: cmaps
     real(kp) :: Pstar, lnRrad, lnRreh
     real(kp) :: lnM, logeps, eps2, eps3
     real(kp) :: lnRhoEnd, bfold
     real(kp) :: ns, logr, alpha
  end type infaspic

  public infaspic, kp, pi, lname


end module aspicvars
