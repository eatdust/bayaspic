module nestparams
  use nestprec
  implicit none

! Parameters for Nested Sampler

!wheter to do Importance Nested Sampling
  logical, parameter :: nestINS = .true.

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

!should match the zero cut in rbffit, or shepfit, that the minimal
!like achievable numerically
  real(fmn), save :: fitLogZero = 0._fmn

!points with loglike < nestlogZero will be ignored (not disfavoured)
  real(fmn), parameter :: nestLogZero = -1d99


 !max no. of iterations, a non-positive value means
 !infinity. MultiNest will terminate if either it has done max no. of
 !iterations or convergence criterion (defined through nest_tol) has
 !been satisfied
  integer(imn), parameter :: nestMaxIter = 0

end module nestparams
