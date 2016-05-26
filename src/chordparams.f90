module polyparams
  use sampl
  implicit none


!max no. of live points !20000
  integer(imn), parameter :: polyNlive = 30000

!evidence tolerance factor !1e-4
  real(fmn), parameter :: polyZTol = 1e-4

!dimension
  integer(imn), save :: polyNdim = 0

!num_repeats
  integer(imn), save :: polyNrepeats = 5*polyNdim

!total number of parameters nestNdim + extra
  integer(imn), save :: polyNpars = 0

!Whether to calculate weighted posteriors
  logical, parameter :: polyWeightedPost = .false.
!Whether to calculate equally weighted posteriors
  logical, parameter :: polyEqualWeightedPost = .true.
!if clustering should be performed
  logical, save :: polyCluster = .false.

!after how many iterations feedback is required & the output files
!should be updated note: posterior files are updated & dumper routine
!is called after every updInt*10 iterations
  integer(imn), parameter :: polyUpdInt = 1000

!file output name
  character(lenlenmn), parameter :: polyBaseDir = 'chains'
  character(len=lenmn), parameter :: polyRootName = 'bayesinf_'

!need update on sampling progress?
  integer(imn), parameter :: polyFeedBack = 2

!whether to resume/save from a previous run
  logical, parameter :: polySave = .true.
  logical, parameter :: polyResume = .false.
!How often to update the resume file
  integer(imn), parameter :: polyUpdSave = -1

!whether to write output files
  logical, parameter :: polyOutLive = .true.
  logical, parameter :: polyOutDead = .false.
  logical, parameter :: polyOutStats = .true.
  logical, parameter :: polyParamNames = .true. 


!should match the zero cut in rbffit, or shepfit, that the minimal
!like achievable numerically
  real(fmn), save :: fitLogZero = 0._fmn


end module polyparams
