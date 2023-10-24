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



module aspicvars  
  use infprec, only : kp, pi
  implicit none

  private

!  integer, parameter :: kp = kind(1._8)

  integer, parameter :: lname = 12

  type infaspic     
     character(len=lname) :: name     
     character(len=lname) :: extname
     integer :: nasp, nhid, nder
     real(kp), dimension(:), pointer :: params
     character(len=lname), dimension(:), pointer :: cmaps
     real(kp) :: Pstar, lnRrad, lnRreh
     real(kp) :: lnM, logeps, eps2, eps3
     real(kp) :: lnRhoEnd, bfold
     real(kp) :: ns, logr, alpha
     real(kp), dimension(:), pointer :: derived
   contains
     procedure :: print => print_infaspic
  end type infaspic

  type(infaspic), save :: AspicModel

  
  public infaspic, kp, pi, lname
  public AspicModel

contains


  subroutine print_infaspic(this)
    implicit none
    class(infaspic), intent(in) :: this

    write(*,*)'--------------------------------------------------------'
    write(*,*)'name:                ',trim(this%name)
    write(*,*)'extname:             ',trim(this%extname)
    write(*,*)'Pstar=               ',this%Pstar
    write(*,*)'lnRrad= lnRreh=      ',this%lnRrad, this%lnRreh
    write(*,*)'lnM=                 ',this%lnM
    write(*,*)'log(epsV1)= epsV23=  ',this%logeps, this%eps2, this%eps3
    write(*,*)'ln(rhoend)= N*-Nend= ',this%lnRhoEnd,this%bfold
    write(*,*)'ns= log(r)= alpha=   ',this%ns,this%logr,this%alpha
    write(*,*)'cmaps=               ',this%cmaps
    write(*,*)'params=              ',this%params
    write(*,*)'derived=             ',this%derived
    write(*,*)'--------------------------------------------------------'

    
  end subroutine print_infaspic
  
  

end module aspicvars
