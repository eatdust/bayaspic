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



program bayaspic
  use sampl
  use samplaspic, only : nest_init_aspic, nest_sample_aspic, nest_free_aspic
  use samplaspic, only : chord_init_aspic, chord_sample_aspic, chord_free_aspic
#ifdef MPISCHED
  use mpi
#endif
  use scheduler, only : initialize_scheduler, free_scheduler,scheduled_size
  use scheduler, only : start_scheduling, irq_scheduler, stop_scheduling
  use scheduler, only : restore_scheduler, scheduler_save_queue
  implicit none


  character(len=*), parameter :: sampler = 'multinest'
!  character(len=*), parameter :: sampler = 'polychord'

#ifdef MPISCHED
  integer :: mpiCode
#endif

  integer, save :: mpiRank = 0
  integer, save :: mpiSize = 1

  integer, parameter :: lmod = 16
  integer :: imodel, nmodels
  character(len=lmod), dimension(:), allocatable :: ModelNames
  character(len=lmod) :: name

  integer :: timein, timenow
  logical, parameter :: display = .true.

!restore all processes from previous run
  logical, parameter :: cpRestart = .false.

!save queues each time one element is done
  logical, parameter :: cpSave = .true.


!if zero, assume same mpiSize between runs.
  integer, save :: mpiPrevSize = 0



!  call initialize_onemodel('ccsi1')
!  call initialize_onemodel('kklti stg')
  call initialize_manymodels()

!  call initialize_filemodels('list_ootest.dat')
!  call initialize_filemodels('list_rreh_dpmodels.dat')



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

     call irq_scheduler()

     name = modelNames(imodel)

     if (display) then
        write(*,*)
        write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(*,*)'MODEL: ',name,'RANK= ',mpiRank,'   QSIZE= ',scheduled_size()
        write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(*,*)
     end if

     select case (sampler)

     case ('multinest')

        call nest_init_aspic(trim(name))
        call nest_sample_aspic()
        call nest_free_aspic()


     case ('polychord')

        call chord_init_aspic(trim(name))
        call chord_sample_aspic()
        call chord_free_aspic()

     case ('timer')

        timein = time()
        timenow = timein
        do while ((timenow-timein).lt.(mpiRank + 1))
           timenow = time()
        end do

     case default
        stop 'bayaspic: sampler not found!'

     end select



     if (cpSave) then
        call scheduler_save_queue(mpiRank)
     endif

     if (stop_scheduling()) exit

  enddo

  call free_scheduler()
  call free_listmodels()

#ifdef MPISCHED
  write(*,*)'process on barrier: mpiRank= ',mpiRank
  call MPI_BARRIER(MPI_COMM_WORLD,mpiCode)
  call MPI_FINALIZE(mpiCode)
#endif

  if (mpiRank.eq.0) then
     write(*,*)'all models done!'
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

    nmodels = 6

    allocate(ModelNames(0:nmodels-1))

    Modelnames(0) = 'saii1 n'
    Modelnames(1) = 'saii1 p'
    Modelnames(2) = 'saii1 f'
    Modelnames(3) = 'saii2 n'
    Modelnames(4) = 'saii2 p'
    Modelnames(5) = 'saii2 f'
!    Modelnames(6) = 'lpi2 6'
!    Modelnames(7) = 'lpi3 2'
!    Modelnames(8) = 'lpi3 4'
!    Modelnames(9) = 'lpi3 6'

  end subroutine initialize_manymodels



  subroutine initialize_filemodels(filename)
    implicit none
    character(len=*), intent(in) :: filename

    character(len=lmod) :: mname
    character(len=2) :: cmod
    character(len=len(cmod)+3) :: fmod

    integer, parameter :: nunit = 100
    integer :: i,counter,ioerr

    write(cmod,'(I2)')lmod
    fmod = '(A'//trim(cmod)//')'

    counter = 0
    open(unit=nunit,file=filename,status='old')
    do
       read(nunit,fmod,iostat=ioerr) mname
       if (ioerr.ne.0) exit
       counter = counter + 1
    enddo
    close(nunit)

    nmodels = counter
    allocate(ModelNames(0:nmodels-1))

    open(unit=nunit,file=filename,status='old')
    do i=0,nmodels-1
       read(nunit,fmod) ModelNames(i)
    enddo
    close(nunit)

    if ((display)) then
       write(*,*)
       write(*,*)'models read from file: ',trim(filename)
       write(*,*)'found nmodels=         ',nmodels
       write(*,*)'models list:           '
!       do i=0,nmodels-1
!          write(*,*)'            ',trim(ModelNames(i))
!       enddo
       write(*,*)
    endif

  end subroutine initialize_filemodels


  subroutine free_listmodels()
    implicit none

    if (allocated(ModelNames)) deallocate(ModelNames)

  end subroutine free_listmodels


end program bayaspic
