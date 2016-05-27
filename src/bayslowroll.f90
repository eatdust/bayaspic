program bayslowroll
  use samplsr, only : nest_init_slowroll, nest_sample_slowroll
  use samplsr, only : chord_init_slowroll, chord_sample_slowroll
  implicit none

  character(len=*), parameter :: sampler = 'polychord'
  
  write(*,*)
  write(*,*)'Sampling slow-roll parameter space with: ', sampler
  write(*,*)

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

end program bayslowroll
