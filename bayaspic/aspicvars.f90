module aspicvars  
  implicit none

  private

  integer, parameter :: kp = kind(1._8)

  integer, parameter :: lname = 6

  type infaspic     
     character(len=lname) :: name     
     integer :: nasp
     real(kp), dimension(:), pointer :: params
     real(kp) :: Pstar, lnRrad, lnRreh
     real(kp) :: lnM, lnRhoEnd, bfold
     real(kp) :: ns, r, alpha
  end type infaspic

  public infaspic, kp, lname


end module aspicvars
