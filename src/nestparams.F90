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



module nestparams
  use sampl
  implicit none

  public

! Parameters for Nested Sampler

!whether to do Importance Nested Sampling
  logical, parameter :: nestINS = .false.

!whether to do multimodal sampling (false)
  logical, parameter :: nestMmodal = .false.

!max no. of live points !20000
  integer(imn), parameter :: nestNlive = 80000

!sample with constant efficiency (false). Should be true if hard prior
!regions are large and the chains remain trapped inside
!(evidence accuracy ensured in INS mode only)
  logical, parameter :: nestCteEff = .true.

!evidence tolerance factor (<0.5)
  real(fmn), parameter :: nestTol = 0.05_fmn

!sampling efficiency (enlargement factor reduction parameter, around 5%)
  real(fmn), parameter :: nestSampEff = 0.05_fmn

!dimension
  integer(imn), save :: nestNdim = 0

!total number of parameters nestNdim + extra
  integer(imn), save :: nestNpars = 0

!no. of parameters on which clustering should be performed (read below)
  integer(imn), save :: nestCdim = 0

!max modes expected, for memory allocation
  integer(imn), parameter :: nestMaxModes = 5

!after how many iterations feedback is required & the output files
!should be updated note: posterior files are updated & dumper routine
!is called after every updInt*10 iterations
  integer(imn), parameter :: nestUpdInt = 2000

!null evidence (set it to very high negative no. if null evidence is unknown)
  real(fmn), parameter :: nestNullZ = -1.0d99

!file output name
  character(len=lenmn), parameter :: nestRootPrefix = 'chains/bayesinf_'
  character(len=lenmn), save :: nestRootName = nestRootPrefix

!seed for nested sampler, -ve means take it from sys clock
  integer(imn), parameter ::  nestSeed = -1

!In order to sample from parameters with periodic boundary conditions
!(or wraparound parameters), set pWrap[i] where i is the index of the
!parameter to be wraparound, to a non-zero value. If pWrap[i] = 0,
!then the ith parameter is not wraparound.
  integer(imn), dimension(:), allocatable, save :: nestPWrap
  
!need update on sampling progress?
  logical, parameter :: nestFeedBack = .false.

!whether to resume from a previous run (within one model, bayaspic may
!still be resumed with this option set to false)
  logical, parameter :: nestResume = .false.

!whether to write output files
  logical, parameter ::  nestOutfile = .true.

!initialize MPI routines?, relevant only for -DMPI
  logical, parameter :: nestInitMPI = .true.

!points with loglike < nestlogZero will be ignored (not disfavoured)
  real(fmn), parameter :: nestLogZero = -1.0d50

 !max no. of iterations, a non-positive value means
 !infinity. MultiNest will terminate if either it has done max no. of
 !iterations or convergence criterion (defined through nest_tol) has
 !been satisfied
  integer(imn), parameter :: nestMaxIter = 1000000


contains


  subroutine nest_print()
    implicit none

    write(*,*)
    write(*,*)'-----------------------------------------------------'
    write(*,*)'Initializing multinest with:'
    write(*,*)'nestNdim     =         ',nestNdim
    write(*,*)'nestNpars    =         ',nestNpars
    write(*,*)'nestINS      =         ',nestINS
    write(*,*)'nestMmodal   =         ',nestMmodal
    write(*,*)'nestNlive    =         ',nestNlive
    write(*,*)'nestCteEff   =         ',nestCteEff
    write(*,*)'nestTol      =         ',nestTol
    write(*,*)'nestFeedBack =         ',nestFeedBack
    write(*,*)'nestResume   =         ',nestResume
    write(*,*)'nestRootName =         ',trim(nestRootName)
    write(*,*)'-----------------------------------------------------'
    write(*,*)

  end subroutine nest_print



  subroutine nest_dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr &
       , maxLogLike, logZ, INSlogZ, logZerr, context)
#ifdef MPISCHED
    use mpi
#endif    
    implicit none
    ! number of samples in posterior array
    integer :: nSamples				
    ! number of live points
    integer :: nlive					
    ! number of parameters saved (physical plus derived)
    integer :: nPar					
    ! array containing the last set of live points
    real(fmn), pointer :: physLive(:,:)	
    ! array with the posterior distribution
    real(fmn), pointer :: posterior(:,:)	
    ! array with mean, sigmas, maxlike & MAP parameters
    real(fmn), pointer :: paramConstr(:)	
    ! max loglikelihood value
    real(fmn) :: maxLogLike			
    ! log evidence
    real(fmn) :: logZ, INSlogZ			
    ! error on log evidence
    real(fmn) :: logZerr			
    ! not required by MultiNest, any additional information user wants to pass
    integer(imn) :: context

    integer :: mpiRank, mpiCode, mpiSize
        
#ifdef MPISCHED
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiRank,mpiCode)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSize,mpiCode)
#else
  mpiRank = 0
  mpiSize = 1
#endif
    
    write(*,*)
    write(*,*)'**************************************************************************'
    write(*,*)'nest_dumper: (rank= size=)',mpiRank,mpiSize
    write(*,*)'nestRoot: ',trim(nestRootName)
    write(*,*)'nSamples= ',nSamples
    write(*,*)'logZ= logZerr=        ',logZ, logZerr
    write(*,*)'maxLogLike= INSlogZ = ',maxLogLike, INSlogZ
    write(*,*)'**************************************************************************'
    write(*,*)


  end subroutine nest_dumper


end module nestparams
