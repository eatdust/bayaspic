module real_precision
  use rbfprec, only : fp
!  integer, parameter :: r8=kind(1._8)
  integer, parameter :: r8 = fp
end module real_precision

module shepprec
  use rbfprec, only : fp,ip
end module shepprec
