program bayaspic
  use nestprec
  use nestwrap, only : nest_init_slowroll, nest_sample_slowroll
  use nestwrap, only : nest_init_aspic, nest_sample_aspic
  implicit none

!  call nest_init_slowroll()

!  call nest_sample_slowroll()

  call nest_init_aspic('lfi')
  call nest_sample_aspic()




end program bayaspic
