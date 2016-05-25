program bayslowroll
  use samplprec
  use nestwrap, only : nest_init_slowroll, nest_sample_slowroll
  
  call nest_init_slowroll()
  call nest_sample_slowroll()

end program bayslowroll
