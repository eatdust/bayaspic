module nestparams
  use nestprec
  implicit none

! Parameters for Nested Sampler

!whether to do multimodal sampling
  logical, parameter :: nestMmodal = .false.

!max no. of live points
  integer(imn), parameter :: nestNlive = 20000

!sample with constant efficiency
  logical, parameter :: nestCteEff = .false.

!evidence tolerance factor
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
  real(fmn), parameter :: nestNullZ = -1d90

!file output name
  character(len=lenmn), save :: nestRootName = 'chains/bayesinf-'

!seed for nested sampler, -ve means take it from sys clock
  integer(imn), parameter ::  nestSeed = -1

!In order to sample from parameters with periodic boundary conditions
!(or wraparound parameters), set pWrap[i] where i is the index of the
!parameter to be wraparound, to a non-zero value. If pWrap[i] = 0,
!then the ith parameter is not wraparound.
  integer(imn), dimension(:), allocatable, save :: nestPWrap
  
!need update on sampling progress?
  logical, parameter :: nestFeedBack = .true.

!whether to resume from a previous run
  logical, parameter :: nestResume = .false.

!whether to write output files
  logical, parameter ::  nestOutfile = .true.

!initialize MPI routines?, relevant only for -DMPI
  logical, parameter :: nestInitMPI = .true.

!points with loglike < logZero will be ignored 
  real(fmn), parameter :: nestLogZero = -huge(1d0)

 !max no. of iterations, a non-positive value means
 !infinity. MultiNest will terminate if either it has done max no. of
 !iterations or convergence criterion (defined through nest_tol) has
 !been satisfied
  integer(imn), parameter :: nestMaxIter = 0

end module nestparams
