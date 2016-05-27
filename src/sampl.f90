module sampl
  implicit none

  public

  integer, parameter :: imn = kind(4)
  integer, parameter :: fmn = kind(1._8)

  integer, parameter :: lenmn = 100

  integer(imn), save :: samplNdim, samplNpars
  real(fmn), save, dimension(:), allocatable :: samplPmin, samplPmax

  character(len=lenmn), parameter :: rootDir = 'chains/'
  character(len=lenmn), parameter :: rootPrefix = 'bayesinf_'
  character(len=lenmn), save :: rootName = rootPrefix

contains

  subroutine init_samplparams(ndim)
    implicit none
    integer, intent(in) :: ndim

    samplNdim = ndim

    allocate(samplPmin(samplNdim))
    allocate(samplPmax(samplNdim))

  end subroutine init_samplparams



  subroutine free_samplparams()
    implicit none
   
    if (allocated(samplPmin)) deallocate(samplPmin)
    if (allocated(samplPmax)) deallocate(samplPmax)
    
  end subroutine free_samplparams
  


  function uncubize_samplparams(ndim,cube)
    implicit none
    integer(imn), intent(in) :: ndim
    real(fmn), dimension(ndim) :: cube
    real(fmn), dimension(ndim) :: uncubize_samplparams
    integer :: i

    if (.not.allocated(samplPmin).or..not.allocated(samplPmax)) then
       stop 'uncubize_samplparams: prior not allocated!'
    endif

    uncubize_samplparams = samplPmin + (samplPmax-samplPmin)*cube

  end function uncubize_samplparams



  subroutine sampl_print()
    implicit none

    if (allocated(samplPmin).and.allocated(samplPmax)) then
       write(*,*)
       write(*,*)'--------------PIORS ON SAMPLED PARAMS----------------'
       write(*,*)'MIN = ',samplPmin
       write(*,*)'MAX = ',samplPmax
       write(*,*)'-----------------------------------------------------'
       write(*,*)
    endif

  end subroutine sampl_print

  
end module sampl
