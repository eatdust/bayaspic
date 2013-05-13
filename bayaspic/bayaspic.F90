program bayaspic
  use nestprec
  use nestwrap, only : nest_init_aspic, nest_sample_aspic
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

  integer, parameter :: lmod = 6
  integer :: imodel, nmodels
  character(len=lmod), dimension(:), allocatable :: ModelNames
  character(len=lmod) :: name

!restore all processes from previous run
  logical, parameter :: cpRestart = .false.

!save queues each time one element is done
  logical, parameter :: cpSave = .false.


!if zero, assume same mpiSize between runs.
  integer, save :: mpiPrevSize = 1


  call initialize_listmodels()

  

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

  print *,'testtt',nmodels

  do

     call start_scheduling(imodel)

!     call irq_scheduler()
 
     name = modelNames(imodel)
     print *,'test ', name
     call nest_init_aspic(trim(name))
     call nest_sample_aspic()

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


  subroutine initialize_listmodels()
    implicit none
    
    nmodels = 2

    allocate(ModelNames(0:nmodels-1))

    ModelNames(0) = 'sfi'
    Modelnames(1) = 'lfi'


  end subroutine initialize_listmodels


  subroutine free_listmodels()
    implicit none

    if (allocated(ModelNames)) deallocate(ModelNames)

  end subroutine free_listmodels


end program bayaspic
