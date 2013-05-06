module aspicvars  
  implicit none

  private

  integer, parameter :: kp = kind(1._8)

  integer, parameter :: lname = 6
  integer, parameter :: neps = 3

  type infaspic     
     character(len=lname) :: name     
     integer :: nasp
     real(kp), dimension(:), pointer :: params
     real(kp) :: Pstar, lnRrad
     real(kp) :: lnM, lnRhoEnd, lnR, bfold
     real(kp), dimension(neps) :: eps
     real(kp) :: ns, r, alpha
  end type infaspic

  public infaspic, kp, neps


end module aspicvars
