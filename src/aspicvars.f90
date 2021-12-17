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



module aspicvars  
  use infprec, only : kp, pi
  implicit none

  private

!  integer, parameter :: kp = kind(1._8)

  integer, parameter :: lname = 12

  type infaspic     
     character(len=lname) :: name     
     character(len=lname) :: extname
     integer :: nasp, nhid
     real(kp), dimension(:), pointer :: params
     character(len=lname), dimension(:), pointer :: cmaps
     real(kp) :: Pstar, lnRrad, lnRreh
     real(kp) :: lnM, logeps, eps2, eps3
     real(kp) :: lnRhoEnd, bfold
     real(kp) :: ns, logr, alpha
  end type infaspic

  
  
  public infaspic, kp, pi, lname



end module aspicvars
