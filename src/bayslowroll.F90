!   This file is part of bayaspic
!
!   Copyright (C) 2013-2025 C. Ringeval
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



program bayslowroll
  use samplsr, only : nest_init_slowroll, nest_sample_slowroll
  use samplsr, only : chord_init_slowroll, chord_sample_slowroll
#if defined MPISCHED
  use mpi
#endif
  implicit none

#if defined MPISCHED
  integer :: mpiCode
#endif
  
  character(len=*), parameter :: sampler = 'multinest'
!  character(len=*), parameter :: sampler = 'polychord'
  
  write(*,*)
  write(*,*)'Sampling slow-roll parameter space with: ', sampler
  write(*,*)


#ifdef MPISCHED
  call MPI_INIT(mpiCode)
#endif
  
  select case (sampler)

  case ('multinest')
     call nest_init_slowroll()
     call nest_sample_slowroll()

  case ('polychord')
     call chord_init_slowroll()
     call chord_sample_slowroll()

  case default
     stop 'bayslowroll: sampler not found!'

  end select


#ifdef MPISCHED
  call MPI_BARRIER(MPI_COMM_WORLD,mpiCode)
  call MPI_FINALIZE(mpiCode)
#endif

  
end program bayslowroll
