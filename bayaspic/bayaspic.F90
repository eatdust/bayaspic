program bayaspic
  use nestprec
  use nestwrap, only : nest_init_aspic, nest_sample_aspic, nest_free_aspic
  use scheduler, only : initialize_scheduler, free_scheduler
  use scheduler, only : start_scheduling, irq_scheduler, stop_scheduling
  use scheduler, only : restore_scheduler, scheduler_save_queue
  implicit none


#ifdef MPISCHED
  include "mpif.h"
  integer :: mpiCode
#endif  

  integer, save :: mpiRank = 0
  integer, save :: mpiSize = 1

  integer, parameter :: lmod = 7
  integer :: imodel, nmodels
  character(len=lmod), dimension(:), allocatable :: ModelNames
  character(len=lmod) :: name

  logical, parameter :: display = .true.

!restore all processes from previous run
  logical, parameter :: cpRestart = .false.

!save queues each time one element is done
  logical, parameter :: cpSave = .false.


!if zero, assume same mpiSize between runs.
  integer, save :: mpiPrevSize = 1


!  call initialize_manymodels()
  call initialize_onemodel('pli')
!  call initialize_filemodels('list_models.dat')


#ifdef MPISCHED
  call MPI_INIT(mpiCode)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiRank,mpiCode)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSize,mpiCode)
#endif


  if (mpiPrevSize.eq.0) mpiPrevSize = mpiSize

  if (cpRestart) then
     call restore_scheduler(mpiPrevSize)
  else
     call initialize_scheduler(nmodels)
  endif


  do

     call start_scheduling(imodel)

!     call irq_scheduler()
 
     name = modelNames(imodel)

     if (display) then
        write(*,*)
        write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(*,*)'MODEL: ',name,'        RANK= ',mpiRank
        write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(*,*)
     end if
        
     call nest_init_aspic(trim(name))
     call nest_sample_aspic()
     call nest_free_aspic()

     if (cpSave) then
        call scheduler_save_queue(mpiRank)
     endif

     if (stop_scheduling()) exit

  enddo

  call free_scheduler()
  call free_listmodels()

#ifdef MPISCHED   
  call MPI_BARRIER(MPI_COMM_WORLD,mpiCode)
  call MPI_FINALIZE(mpiCode)
#endif      

  if (mpiRank.ne.0) then
     write(*,*)'process stopped: mpiRank= ',mpiRank
     stop
  endif
  

contains



  subroutine initialize_onemodel(mname)
    implicit none
    character(len=*), intent(in) :: mname
    nmodels = 1

    allocate(ModelNames(0:0))

    ModelNames(0) = mname

  end subroutine initialize_onemodel



  subroutine initialize_manymodels()
    implicit none
    
    nmodels = 3

    allocate(ModelNames(0:nmodels-1))

    ModelNames(0) = 'hi'
    Modelnames(1) = 'lfi'
    Modelnames(2) = 'sfi'
        
  end subroutine initialize_manymodels



  subroutine initialize_filemodels(filename)
    implicit none
    character(len=*), intent(in) :: filename

    character(len=lmod) :: mname
    integer, parameter :: nunit = 100
    integer :: i,counter,ioerr

    counter = 0
    open(unit=nunit,file=filename,status='old')
    do 
       read(nunit,*,iostat=ioerr) mname
       if (ioerr.ne.0) exit
       counter = counter + 1
    enddo
    close(nunit)

    nmodels = counter
    allocate(ModelNames(0:nmodels-1))

    open(unit=nunit,file=filename,status='old')
    do i=0,nmodels-1
       read(nunit,*) ModelNames(i)
    enddo
    close(nunit)

    if (display) then
       write(*,*)
       write(*,*)'models read from file: ',trim(filename)
       write(*,*)'found nmodels=         ',nmodels
       write(*,*)'models list:           '
       do i=0,nmodels-1
          write(*,*)'            ',trim(ModelNames(i))
       enddo
       write(*,*)
    endif

  end subroutine initialize_filemodels


  subroutine free_listmodels()
    implicit none

    if (allocated(ModelNames)) deallocate(ModelNames)

  end subroutine free_listmodels


end program bayaspic
