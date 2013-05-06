module nestwrap
  use nestprec
  implicit none

  private
 
  character(len=*), parameter :: fileweights = 'rbfdata/weights.dat'
  character(len=*), parameter :: filecentres = 'rbfdata/centres.dat'
  character(len=*), parameter :: filebounds = 'rbfdata/bounds.dat'

  integer(imn), save :: rbfNdim = 0
  integer(imn), parameter :: pstarpos = 1
  integer(imn), parameter :: eps1pos = 2

  real(fmn), save, dimension(:), allocatable :: nestPmin, nestPmax


  public nest_init_slowroll, nest_sample_slowroll
  public nest_init_aspic, nest_sample_aspic


contains


   subroutine nest_init_slowroll()
    use rbflike, only : initialize_rbf_like
    use rbflike, only : get_rbf_ndim
    use nestparams, only : nestNdim, nestNpars, nestCdim
    use nestparams, only : nestPWrap
    implicit none
       
    call initialize_rbf_like(fileweights, filecentres, filebounds)

    rbfNdim = get_rbf_ndim()
    nestNdim = rbfNdim
    nestNpars = nestNdim
    nestCdim = nestNdim
        
    allocate(nestPwrap(nestNdim))
    nestPwrap = 0

    call nest_print()

  end subroutine nest_init_slowroll



  subroutine nest_sample_slowroll()
    use nested, only : nestRun
    use nestparams
    implicit none

    integer(imn) :: nclusters			
    integer(imn) :: context
    integer(imn) :: maxNode 			

    call nestRun(nestMmodal,nestCteEff,nestNlive,nestZtol,nestSampEff,nestNdim,nestNpars, &
         nestCdim,nestMaxModes,nestUpdInt,nestNullZ,nestRootName,nestSeed,nestPwrap, &
         nestFeedBack,nestResume,nestOutfile,nestInitMPI,nestLogZero,nestMaxIter &
         ,rbf_multinest_slowroll_loglike,nest_dumper,context)

  end subroutine nest_sample_slowroll

  
  subroutine rbf_multinest_slowroll_loglike(cube,nestdim,nestpars,lnew,context)
    use rbfprec, only : fp
    use rbflike, only : uncubize_rbfparams, check_rbf
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
    
    cube(1:nestdim) = real(uncubize_rbfparams(nestdim,rbfcube),fmn)
    
  end subroutine rbf_multinest_slowroll_loglike



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
    write(*,*)'nestRootName =         ',nestRootName
    write(*,*)
    write(*,*)'rbfNdim      =         ',rbfNdim
    write(*,*)'-----------------------------------------------'

    if (allocated(nestPmin).and.allocated(nestPmax)) then

       write(*,*)
       write(*,*)'-------------Uniform Prior bounds--------------'
       write(*,*)'MIN = ',nestPmin
       write(*,*)'MAX = ',nestPmax
       write(*,*)'-----------------------------------------------'

    endif

  end subroutine nest_print


  subroutine nest_init_aspic(modelname)
    use aspicwrap, only : set_aspic_model, get_aspic_ntot
    use aspicwrap, only : get_aspic_allprior
    use rbflike, only : initialize_rbf_like
    use rbflike, only : get_rbf_ndim, get_rbf_xpmin, get_rbf_xpmax
    use nestparams, only : nestNdim, nestNpars, nestCdim
    use nestparams, only : nestPWrap, nestRootName
    implicit none    
    character(len=*), intent(in) :: modelname

    call set_aspic_model(modelname)

    nestNdim = get_aspic_ntot()
    nestNpars = nestNdim
    nestCdim = nestNdim
    nestRootName = nestRootName//modelname
    
    allocate(nestPmin(nestNdim))
    allocate(nestPmax(nestNdim))

    call get_aspic_allprior(nestPmin, nestPmax)
        
    call initialize_rbf_like(fileweights, filecentres, filebounds)

!cut prior of P* by the one encoded in the likelihood
    nestPmin(pstarpos) = max(get_rbf_xpmin(pstarpos),nestPmin(pstarpos))
    nestPmax(pstarpos) = min(get_rbf_xpmax(pstarpos),nestPmax(pstarpos))

    rbfNdim = get_rbf_ndim()

    if (nestNdim.ne.rbfNdim) then
       write(*,*)'nest_inif: '
       write(*,*)'nestdim= rbfdim= ',nestNdim, rbfNdim    
    endif    
       
    allocate(nestPwrap(nestNdim))
    nestPwrap = 0

    call nest_print()

  end subroutine nest_init_aspic

  
 
  subroutine nest_sample_aspic()
    use nested, only : nestRun
    use nestparams
    implicit none

    integer(imn) :: nclusters			
    integer(imn) :: context
    integer(imn) :: maxNode 			

    call nestRun(nestMmodal,nestCteEff,nestNlive,nestZtol,nestSampEff,nestNdim,nestNpars, &
         nestCdim,nestMaxModes,nestUpdInt,nestNullZ,nestRootName,nestSeed,nestPwrap, &
         nestFeedBack,nestResume,nestOutfile,nestInitMPI,nestLogZero,nestMaxIter &
         ,rbf_multinest_aspic_loglike,nest_dumper,context)

  end subroutine nest_sample_aspic



  subroutine rbf_multinest_aspic_loglike(cube,nestdim,nestpars,lnew,context)
    use rbfprec, only : fp
    use rbflike, only : cubize_rbfparams, uncubize_rbfparams, check_rbf
    use rbflike, only : rbflike_eval, cutmin_rbfparams
    use aspicwrap, only : get_aspic_slowroll, get_aspic_derived
    implicit none   
    integer(imn) :: nestdim, nestpars
    real(fmn), dimension(nestpars) :: cube
    real(fmn) :: lnew
    integer(imn) :: context,i
    real(fmn), dimension(nestdim) :: mnpars
    real(fp), dimension(rbfNdim) :: rbfcube, rbfpars, rbfcuts


    if (.not.check_rbf()) stop 'rbf_multinest_loglike: not initialized!'

!get the aspic parameters we are sampling on
    mnpars = uncubize_nestparams(nestdim,cube)

!use aspic to get the slowroll parameters
    rbfpars = get_aspic_slowroll(rbfNdim,mnpars)

!if eps1<eps1min, the likelihood is flat
    rbfcuts = cutmin_rbfparams(rbfNdim,eps1pos,rbfpars)

!go into cubic space for the rbf
    rbfcube = cubize_rbfparams(rbfNdim,rbfcuts)

!get the likelihood
    lnew = rbflike_eval(rbfcube)

!dump the physical params into cube
    cube(1:nestdim) = mnpars(1:nestdim)

!+ derived parameters
    do i=1,nestpars-nestdim
       cube(nestdim+i) = get_aspic_derived(i)
    enddo

  end subroutine rbf_multinest_aspic_loglike


  function uncubize_nestparams(nestdim,nestcube)
    implicit none
    integer(imn), intent(in) :: nestdim
    real(fmn), dimension(nestdim) :: nestcube
    real(fmn), dimension(nestdim) :: uncubize_nestparams
    integer :: i

    if (.not.allocated(nestPmin).or..not.allocated(nestPmax)) then
       stop 'uncubize_nestparams: prior not allocated!'
    endif

    uncubize_nestparams = nestPmin + (nestPmax-nestPmin)*nestcube
    
  end function uncubize_nestparams

 

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


end module nestwrap
