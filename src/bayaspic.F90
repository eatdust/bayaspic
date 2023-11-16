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



program bayaspic
  use sampl
  use samplaspic, only : nest_init_aspic, nest_sample_aspic, nest_free_aspic
  use samplaspic, only : chord_init_aspic, chord_sample_aspic, chord_free_aspic
#ifdef MPISCHED
  use mpi
#endif
  use scheduler, only : initialize_scheduler, free_scheduler,scheduled_size
  use scheduler, only : start_scheduling, irq_scheduler, stop_scheduling
  use scheduler, only : restore_scheduler, scheduler_save_queue, check_saved_queue
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


  call initialize_onemodel('si mc')
!  call initialize_manymodels()

!  call initialize_filemodels('list_rreh_dpmodels.dat')
!  call initialize_filemodels('list_rreh_qpmodels.dat')
  


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

  if ((cpSave).and.(.not.cpRestart)) then
     if (check_saved_queue(mpiRank)) stop 'previous saved files present!'
  endif

  do

     call start_scheduling(imodel)

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

     call irq_scheduler()

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

    nmodels = 20

    allocate(ModelNames(0:nmodels-1))

    Modelnames(0) = 'saiii1 n'
    Modelnames(1) = 'rclfi1 pm'
    Modelnames(2) = 'rclfi1 mm'
    Modelnames(3) = 'rclfi1 mp'
    Modelnames(4) = 'rclfi2 pm'
    Modelnames(5) = 'rclfi2 mm'
    Modelnames(6) = 'rclfi2 mp'
    Modelnames(7) = 'rclfi3 pp'
    Modelnames(8) = 'rclfi3 pm'
    Modelnames(9) = 'rclfi3 mp'
    Modelnames(10) = 'rcipi 1m'
    Modelnames(11) = 'rcipi 1tune2p'
    Modelnames(12) = 'rcipi 1tune4p'
    Modelnames(13) = 'rcipi 1tune4m'
    Modelnames(14) = 'rcipi 2'
    Modelnames(15) = 'rcipi 22'
    Modelnames(16) = 'rcipi 23'
    Modelnames(17) = 'rcipi 24'
    Modelnames(18) = 'rcipi 25'
    Modelnames(19) = 'rcipi 26'
    

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
