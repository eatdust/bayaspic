module nestparams
  use sampl
  implicit none

  public

! Parameters for Nested Sampler

!whether to do Importance Nested Sampling
  logical, parameter :: nestINS = .false.

!whether to do multimodal sampling
  logical, parameter :: nestMmodal = .false.

!max no. of live points !20000
  integer(imn), parameter :: nestNlive = 30000

!sample with constant efficiency
  logical, parameter :: nestCteEff = .false.

!evidence tolerance factor !1e-4
  real(fmn), parameter :: nestZTol = 1e-4

!sampling efficiency (enlargement factor reduction parameter)
  real(fmn), parameter :: nestSampEff = 0.5

!dimension
  integer(imn), save :: nestNdim = 0

!total number of parameters nestNdim + extra
  integer(imn), save :: nestNpars = 0

!no. of parameters on which clustering should be performed (read below)
  integer(imn), save :: nestCdim = 0

!max modes expected, for memory allocation
  integer(imn), parameter :: nestMaxModes = 10

!after how many iterations feedback is required & the output files
!should be updated note: posterior files are updated & dumper routine
!is called after every updInt*10 iterations
  integer(imn), parameter :: nestUpdInt = 1000

!null evidence (set it to very high negative no. if null evidence is unknown)
  real(fmn), parameter :: nestNullZ = -1d99

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

!whether to resume from a previous run
  logical, parameter :: nestResume = .true.

!whether to write output files
  logical, parameter ::  nestOutfile = .true.

!initialize MPI routines?, relevant only for -DMPINEST
  logical, parameter :: nestInitMPI = .true.

!points with loglike < nestlogZero will be ignored (not disfavoured)
  real(fmn), parameter :: nestLogZero = -1d99


 !max no. of iterations, a non-positive value means
 !infinity. MultiNest will terminate if either it has done max no. of
 !iterations or convergence criterion (defined through nest_tol) has
 !been satisfied
  integer(imn), parameter :: nestMaxIter = 0


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
    write(*,*)'nestZtol     =         ',nestZtol
    write(*,*)'nestFeedBack =         ',nestFeedBack
    write(*,*)'nestResume   =         ',nestResume
    write(*,*)'nestRootName =         ',trim(nestRootName)
    write(*,*)'-----------------------------------------------------'
    write(*,*)

  end subroutine nest_print



  subroutine nest_dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr &
       , maxLogLike, logZ, INSlogZ, logZerr, context)

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

    write(*,*)
    write(*,*)'*****************************************************'
    write(*,*)'nest_dumper: '
    write(*,*)'nestRoot: ',trim(nestRootName)
    write(*,*)'nSamples= ',nSamples
    write(*,*)'logZ= logZerr=        ',logZ, logZerr
    write(*,*)'maxLogLike= INSlogZ = ',maxLogLike, INSlogZ
    write(*,*)'*****************************************************'
    write(*,*)


  end subroutine nest_dumper


end module nestparams
