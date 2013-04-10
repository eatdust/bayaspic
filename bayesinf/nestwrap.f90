module nestwrapper
  use nestprec
  implicit none

  private
 
  character(len=*), parameter :: fileweights = 'rbfdata/weights.dat'
  character(len=*), parameter :: filecentres = 'rbfdata/centres.dat'
  character(len=*), parameter :: filebounds = 'rbfdata/bounds.dat'

  public nest_init, nest_sample


contains

  subroutine nest_init()
    use rbflike, only : initialize_rbf_like
    use rbflike, only : get_rbf_ndim
    use nestparams, only : nestNdim, nestNpars, nestCdim
    use nestparams, only : nestPWrap
    implicit none
   

    call initialize_rbf_like(fileweights, filecentres, filebounds)

    nestNdim = get_rbf_ndim()
    nestNpars = nestNdim
    nestCdim = nestNdim

    allocate(nestPwrap(nestNdim))
    nestPwrap = 0

    call nest_print()

  end subroutine nest_init


  subroutine nest_print()
    use nestparams
    implicit none

    write(*,*)'-----------------------------------------------'
    write(*,*)'Initializing multinest with:'
    write(*,*)'nestNdim     =         ',nestNdim
    write(*,*)'nestNpars    =         ',nestNpars
    write(*,*)'nestMmodal   =         ',nestMmodal
    write(*,*)'nestNlive    =         ',nestNlive
    write(*,*)'nestCteEff   =         ',nestCteEff
    write(*,*)'nestZtol     =         ',nestZtol
    write(*,*)'nestResume   =         ',nestResume
    write(*,*)'nestRootName =         ',nestResume
    write(*,*)'-----------------------------------------------'

  end subroutine nest_print



  subroutine nest_sample()
    use nested, only : nestRun
    use nestparams
    implicit none

    integer(imn) :: nclusters			
    integer(imn) :: context
    integer(imn) :: maxNode 			

    call nestRun(nestMmodal,nestCteEff,nestNlive,nestZtol,nestSampEff,nestNdim,nestNpars, &
         nestCdim,nestMaxModes,nestUpdInt,nestNullZ,nestRootName,nestSeed,nestPwrap, &
         nestFeedBack,nestResume,nestOutfile,nestInitMPI,nestLogZero,nestMaxIter &
         ,rbf_multinest_loglike,nest_dumper,context)

  end subroutine nest_sample


  

  subroutine rbf_multinest_loglike(cube,nestdim,nestpars,lnew,context)
    use rbfprec, only : fp
    use rbflike, only : uncubize_params, check_rbf
    use rbflike, only : rbflike_eval
    implicit none   
    integer(imn) :: nestdim, nestpars
    real(fmn), dimension(nestpars) :: cube
    real(fmn) :: lnew
    integer(imn) :: context
    real(fp), dimension(nestdim) :: rbfcube

    if (.not.check_rbf()) stop 'rbf_multinest_loglike: not initialized!'
    
    rbfcube(1:nestdim) = cube(1:nestdim)

    lnew = rbflike_eval(rbfcube)
    
    cube(1:nestdim) = real(uncubize_params(nestdim,rbfcube),fmn)
    

  end subroutine rbf_multinest_loglike





  subroutine nest_dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr &
       , maxLogLike, logZ, logZerr, context)

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
    real(fmn) :: logZ				
! error on log evidence
    real(fmn) :: logZerr			
! not required by MultiNest, any additional information user wants to pass
    integer(imn) :: context

    
    write(*,*)
    write(*,*)'***********************************'
    write(*,*)'nest_dumper: '
    write(*,*)'nSamples= logZ= logZerr= ',nSamples, logZ, logZerr
    write(*,*)'maxLogLike= ',maxLogLike
    write(*,*)'***********************************'
    write(*,*)
   

  end subroutine nest_dumper


end module nestwrapper
