!   This file is part of bayaspic
!
!   Copyright (C) 2021 C. Ringeval
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
!   along with Foobar.  If not, see <https://www.gnu.org/licenses/>.


module chordparams
  use sampl
  implicit none
    
  public

!max no. of live points !20000
  integer(imn), parameter :: chordNlive = 20000

!evidence tolerance factor !1e-4
  real(fmn), parameter :: chordZTol = 1e-4

!dimension
  integer(imn), save :: chordNdim = 0

!number of derived parameters:
  integer(imn), save :: chordNderived = 0

!total number of parameters
  integer(imn), save :: chordNpars = 0

  real(fmn), parameter :: chordBoost = 5._fmn

!points with loglike < chordlogZero will be ignored (not disfavoured)
  real(fmn), parameter :: chordLogZero = -1d30

!Whether to calculate weighted posteriors
  logical, parameter :: chordWeightedPost = .true.

!Whether to calculate equally weighted posteriors
  logical, parameter :: chordEqualWeightedPost = .false.

!Whether to calculate equally weighted cluster posteriors
  logical, parameter :: chordClusterPost = .false.


!if clustering should be performed
  logical, save :: chordCluster = .false.

!file output name
  character(len=lenmn), save :: chordDir = 'chains'
  character(len=lenmn), save :: chordName = 'bayesinf_'

!need update on sampling progress?
  integer(imn), parameter :: chordFeedBack = 3

!whether to resume/save from a previous run
  logical, parameter :: chordSave = .true.
  logical, parameter :: chordResume = .false.

!whether to write output files
  logical, parameter :: chordOutLive = .true.
  logical, parameter :: chordOutDead = .false.
  logical, parameter :: chordOutStats = .true.
  logical, parameter :: chordOutParamNames = .false. 


contains


  subroutine chord_print()
    implicit none

    write(*,*)
    write(*,*)'-----------------------------------------------------'
    write(*,*)'Initializing polychord with:'
    write(*,*)'chordNdim     =         ',chordNdim
    write(*,*)'chordNpars    =         ',chordNpars
    write(*,*)'chordNderived =         ',chordNderived
    write(*,*)'chordNrepeats =         ',5*chordNdim
    write(*,*)'chordNlive    =         ',chordNlive
    write(*,*)'chordZtol     =         ',chordZtol
    write(*,*)'chordFeedBack =         ',chordFeedBack
    write(*,*)'chordBoost    =         ',chordBoost
    write(*,*)'chordLogZero  =         ',chordLogZero
    write(*,*)'chordSave     =         ',chordSave
    write(*,*)'chordResume   =         ',chordResume
    write(*,*)'chordCluster  =         ',chordCluster
    write(*,*)'chordWPost    =         ',chordWeightedPost
    write(*,*)'chordEWPost   =         ',chordEqualWeightedPost
    write(*,*)'chordCPost    =         ',chordClusterPost
    write(*,*)'chordOutLive  =         ',chordOutLive
    write(*,*)'chordOutDead  =         ',chordOutDead
    write(*,*)'chordOutStats =         ',chordOutStats
    write(*,*)'chordOutPNames=         ',chordOutParamNames
    write(*,*)'chordDir      =         ',trim(chordDir)
    write(*,*)'chordName     =         ',trim(chordName)
    write(*,*)'-----------------------------------------------------'
    write(*,*)

  end subroutine chord_print


  subroutine chord_settings(set)
    use settings_module, only : program_settings
    use settings_module, only : initialise_settings
    implicit none

    type(program_settings), intent(out) :: set

    set%ndims = chordNdim
    set%nderived = chordNderived
    set%nlive = chordNlive
    set%num_repeats = chordNdim*5
    set%do_clustering = chordCluster
    set%precision_criterion = chordZtol
    set%logzero = chordLogZero
    set%base_dir = trim(chordDir)
    set%file_root = trim(chordName)
    set%write_resume = chordSave
    set%read_resume = chordResume
    set%write_live = chordOutLive
    set%write_dead = chordOutDead
    set%write_stats=  chordOutStats
    set%equals = chordEqualWeightedPost
    set%posteriors = chordWeightedPost
    set%cluster_posteriors = chordClusterPost
    set%feedback = chordFeedBack
    set%boost_posterior = chordBoost
    

  end subroutine chord_settings


end module chordparams
