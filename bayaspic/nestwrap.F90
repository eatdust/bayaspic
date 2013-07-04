module nestwrap
  use nestprec
  implicit none

  private
 
  character(len=*), parameter :: fileweights = 'rbfdata/weights.dat'
  character(len=*), parameter :: filecentres = 'rbfdata/centres.dat'
  character(len=*), parameter :: filerbfbounds = 'rbfdata/bounds.dat'

  character(len=*), parameter :: fileshep = 'shepdata/shepdata.dat'
  character(len=*), parameter :: filepost = 'shepdata/postcubed.dat'
  character(len=*), parameter :: fileshepbounds = 'shepdata/bounds.dat'

  integer(imn), parameter :: numDerivedParams = 10

  integer(imn), save :: fitNdim = 0
  integer(imn), parameter :: pstarpos = 1
  integer(imn), parameter :: eps1pos = 2

  real(fmn), save, dimension(:), allocatable :: nestPmin, nestPmax

  logical, parameter :: display = .false.

!the name of the fastlike method we call

! radial basis functions
!  character(len=*), parameter :: fastLikeName = 'rbf'

! inverse shepard method
  character(len=*), parameter :: fastLikeName = 'shep'


  public nest_free_slowroll
  public nest_init_slowroll, nest_sample_slowroll
#ifdef ASPIC
  public nest_init_aspic, nest_sample_aspic, nest_free_aspic
#endif

contains


   subroutine nest_init_slowroll()
    use rbflike, only : initialize_rbf_like
    use rbflike, only : get_rbf_ndim, get_rbf_fmin
    use sheplike, only : initialize_shep_like
    use sheplike, only : get_shep_ndim, get_shep_fmin
    use nestparams, only : nestNdim, nestNpars, nestCdim
    use nestparams, only : nestPWrap, nestRootName, nestRootPrefix
    use nestparams, only : fitLogZero
    implicit none
    
    select case(fastLikeName)
    
    case ('rbf')

       call initialize_rbf_like(fileweights, filecentres, filerbfbounds)
       fitNdim = get_rbf_ndim()
       fitLogZero = get_rbf_fmin()

    case ('shep')

       call initialize_shep_like(fileshep, filepost, fileshepbounds)
       fitNdim = get_shep_ndim()
       fitLogZero = get_shep_fmin()

    case default

       stop 'nest_init_slowroll: fast like not found!'

    end select
       
    nestNdim = fitNdim
    nestNpars = nestNdim
    nestCdim = nestNdim
    if (fitNdim.eq.3) then
       nestRootName = trim(nestRootPrefix)//'sr2'
    elseif (fitNdim.eq.4) then
       nestRootName = trim(nestRootPrefix)//'sr3'
    endif

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

    select case (fastLikeName)

    case ('rbf')

       call nestRun(nestMmodal,nestCteEff,nestNlive,nestZtol,nestSampEff,nestNdim,nestNpars, &
            nestCdim,nestMaxModes,nestUpdInt,nestNullZ,nestRootName,nestSeed,nestPwrap, &
            nestFeedBack,nestResume,nestOutfile,nestInitMPI,nestLogZero,nestMaxIter &
            ,rbf_multinest_slowroll_loglike,nest_dumper,context)

    case ('shep')

       call nestRun(nestMmodal,nestCteEff,nestNlive,nestZtol,nestSampEff,nestNdim,nestNpars, &
            nestCdim,nestMaxModes,nestUpdInt,nestNullZ,nestRootName,nestSeed,nestPwrap, &
            nestFeedBack,nestResume,nestOutfile,nestInitMPI,nestLogZero,nestMaxIter &
            ,shep_multinest_slowroll_loglike,nest_dumper,context)

    case default

       stop 'nest_sample_slowroll: fast like not found!'

    end select


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



  subroutine shep_multinest_slowroll_loglike(cube,nestdim,nestpars,lnew,context)
    use shepprec, only : fp
    use sheplike, only : uncubize_shepparams, check_shep
    use sheplike, only : sheplike_eval
    implicit none   
    integer(imn) :: nestdim, nestpars
    real(fmn), dimension(nestpars) :: cube
    real(fmn) :: lnew
    integer(imn) :: context
    real(fp), dimension(nestdim) :: shepcube

    if (.not.check_shep()) stop 'shep_multinest_loglike: not initialized!'
    
    shepcube(1:nestdim) = cube(1:nestdim)

    lnew = sheplike_eval(shepcube)
    
    cube(1:nestdim) = real(uncubize_shepparams(nestdim,shepcube),fmn)
    
  end subroutine shep_multinest_slowroll_loglike


  subroutine nest_print()
    use nestparams
    implicit none

    write(*,*)'-----------------------------------------------------'
    write(*,*)'Initializing multinest with:'
    write(*,*)'nestNdim     =         ',nestNdim
    write(*,*)'nestNpars    =         ',nestNpars
    write(*,*)'nestMmodal   =         ',nestMmodal
    write(*,*)'nestNlive    =         ',nestNlive
    write(*,*)'nestCteEff   =         ',nestCteEff
    write(*,*)'nestZtol     =         ',nestZtol
    write(*,*)'nestFeedBack =         ',nestFeedBack
    write(*,*)'nestResume   =         ',nestResume
    write(*,*)'nestRootName =         ',trim(nestRootName)
    write(*,*)
    write(*,*)'fast like is :         ',fastLikeName
    write(*,*)'lnZeroMin    =         ',fitLogZero
    write(*,*)'fitNdim      =         ',fitNdim
    write(*,*)'-----------------------------------------------------'

    if (allocated(nestPmin).and.allocated(nestPmax)) then

       write(*,*)
       write(*,*)'--------------PIORS ON SAMPLED PARAMS----------------'
       write(*,*)'MIN = ',nestPmin
       write(*,*)'MAX = ',nestPmax
       write(*,*)'-----------------------------------------------------'

    endif

  end subroutine nest_print


  subroutine nest_dump_priors(extname)
    use nestparams, only : nestNdim, nestRootname
    implicit none
    character(len=*), intent(in) :: extname
    integer, parameter :: nunit = 412
    
    if ((.not.allocated(nestPmin)).or.(.not.allocated(nestPmax))) then
       stop 'nest_dump_priors: prior not allocated!'
    endif

    open(unit=nunit,file=trim(nestRootName)//'.ini', status='unknown')
    write(nunit,*)'limits[lnA]=',nestPmin(1),nestPmax(1)
    write(nunit,*)'limits[lnRreh]=',nestPmin(2),nestPmax(2)

    select case (nestNdim-2)
    case (0)
    case (1)
       write(nunit,*)'limits[c1]=',nestPmin(3),nestPmax(3)
    case(2)
       write(nunit,*)'limits[c1]=',nestPmin(3),nestPmax(3)
       write(nunit,*)'limits[c2]=',nestPmin(4),nestPmax(4)
    case(3)
       write(nunit,*)'limits[c1]=',nestPmin(3),nestPmax(3)
       write(nunit,*)'limits[c2]=',nestPmin(4),nestPmax(4)
       write(nunit,*)'limits[c3]=',nestPmin(5),nestPmax(5)
    case default
       stop 'nest_dump_priors: case not implemented!'
    end select

    close(nunit)


  end subroutine nest_dump_priors



  subroutine nest_free_slowroll()
    use nestparams, only : nestPwrap    
    implicit none
   
    if (allocated(nestPwrap)) deallocate(nestPwrap)
    if (allocated(nestpmin)) deallocate(nestpmin)
    if (allocated(nestpmax)) deallocate(nestpmax)

  end subroutine nest_free_slowroll

  

#ifdef ASPIC
  subroutine nest_init_aspic(modelname)
    use wraspic, only : set_model, get_ntot
    use wraspic, only : get_allprior
    use rbflike, only : initialize_rbf_like, check_rbf, get_rbf_fmin
    use rbflike, only : get_rbf_ndim, get_rbf_xpmin, get_rbf_xpmax
    use sheplike, only : initialize_shep_like, check_shep, get_shep_fmin
    use sheplike, only : get_shep_ndim, get_shep_xpmin, get_shep_xpmax
    use nestparams, only : nestNdim, nestNpars, nestCdim
    use nestparams, only : nestPWrap, nestRootName, nestRootPrefix
    use nestparams, only : fitLogZero
    implicit none    
    integer :: cpos, lenmod
    character(len=*), intent(in) :: modelname
    character(len=len(modelname)) :: name, subname

    character, parameter :: separator=' '

    lenmod = len(modelname)

    cpos = scan(modelname,separator)

    if (cpos.eq.0) then
       name = trim(adjustl(modelname))
       nestRootName = trim(nestRootPrefix)//name
       call set_model(name)

    else
       name = trim(adjustl(modelname(1:cpos-1)))
       subname = trim(adjustl(modelname(cpos:lenmod)))
       nestRootName = trim(nestRootPrefix)//trim(name)//subname
       call set_model(name,subname)

    endif

   
    nestNdim = get_ntot()
    nestNpars = nestNdim + numDerivedParams
    nestCdim = nestNdim

    allocate(nestPmin(nestNdim))
    allocate(nestPmax(nestNdim))

    call get_allprior(nestPmin, nestPmax)

    select case (fastLikeName)

    case ('rbf')

       if (.not.check_rbf()) then
          call initialize_rbf_like(fileweights, filecentres, filerbfbounds)
       endif

!cut prior of P* by the one encoded in the likelihood
       nestPmin(pstarpos) = max(get_rbf_xpmin(pstarpos),nestPmin(pstarpos))
       nestPmax(pstarpos) = min(get_rbf_xpmax(pstarpos),nestPmax(pstarpos))
    
       fitNdim = get_rbf_ndim()
       fitLogZero = get_rbf_fmin()

    case ('shep')

       if (.not.check_shep()) then
          call initialize_shep_like(fileshep, filepost, fileshepbounds)
       endif

!cut prior of P* by the one encoded in the likelihood
       nestPmin(pstarpos) = max(get_shep_xpmin(pstarpos),nestPmin(pstarpos))
       nestPmax(pstarpos) = min(get_shep_xpmax(pstarpos),nestPmax(pstarpos))
    
       fitNdim = get_shep_ndim()
       fitLogZero = get_shep_fmin()

    case default

       stop 'nest_init_aspic: fast like not found!'

    end select


    if (nestNdim.ne.fitNdim) then
       write(*,*)'nest_init_aspic: '
       write(*,*)'nestdim= fastlikedim= ',nestNdim, fitNdim    
    endif    
       
    allocate(nestPwrap(nestNdim))
    nestPwrap = 0

    call nest_dump_priors(trim(name)//subname)

    call nest_print()
    
  end subroutine nest_init_aspic

  
  subroutine nest_free_aspic()
    use nestparams, only : nestPwrap
    use wraspic, only : free_model
    implicit none

    call free_model()
    
    if (allocated(nestPwrap)) deallocate(nestPwrap)
    if (allocated(nestpmin)) deallocate(nestpmin)
    if (allocated(nestpmax)) deallocate(nestpmax)

  end subroutine nest_free_aspic


 
  subroutine nest_sample_aspic()
    use nested, only : nestRun
    use nestparams
    implicit none

    integer(imn) :: nclusters			
    integer(imn) :: context
    integer(imn) :: maxNode 			

    select case (fastLikeName)

    case ('rbf')

       call nestRun(nestMmodal,nestCteEff,nestNlive,nestZtol,nestSampEff,nestNdim,nestNpars, &
            nestCdim,nestMaxModes,nestUpdInt,nestNullZ,nestRootName,nestSeed,nestPwrap, &
            nestFeedBack,nestResume,nestOutfile,nestInitMPI,nestLogZero,nestMaxIter &
            ,rbf_multinest_aspic_loglike,nest_dumper,context)

    case ('shep')

       call nestRun(nestMmodal,nestCteEff,nestNlive,nestZtol,nestSampEff,nestNdim,nestNpars, &
            nestCdim,nestMaxModes,nestUpdInt,nestNullZ,nestRootName,nestSeed,nestPwrap, &
            nestFeedBack,nestResume,nestOutfile,nestInitMPI,nestLogZero,nestMaxIter &
            ,shep_multinest_aspic_loglike,nest_dumper,context)

    case default

       stop 'nest_sample_aspic: fast like not found!'

    end select

  end subroutine nest_sample_aspic



  subroutine rbf_multinest_aspic_loglike(cube,nestdim,nestpars,lnew,context)
    use rbfprec, only : fp
    use rbflike, only : cubize_rbfparams, uncubize_rbfparams, check_rbf
    use rbflike, only : rbflike_eval, cutmin_rbfparams
    use wraspic, only : get_slowroll, get_derived
    use wraspic, only : test_aspic_hardprior, test_reheating_hardprior
    use nestparams, only : fitLogZero, nestLogZero
    implicit none   
    integer(imn) :: nestdim, nestpars
    real(fmn), dimension(nestpars) :: cube
    real(fmn) :: lnew
    integer(imn) :: context,i
    real(fmn), dimension(nestdim) :: mnpars
    real(fp), dimension(fitNdim) :: rbfcube, rbfpars, rbfcuts


    if (.not.check_rbf()) stop 'rbf_multinest_loglike: not initialized!'

!get the physical parameters we are sampling on
    mnpars = uncubize_nestparams(nestdim,cube)

!check for any hard prior in aspic model parameters
    if (test_aspic_hardprior(mnpars)) then

       lnew = nestLogZero

    else
       
!use aspic to get the slowroll parameters
       rbfpars = get_slowroll(fitNdim,mnpars)

!if eps1<eps1min, the likelihood is flat
       rbfcuts = cutmin_rbfparams(fitNdim,eps1pos,rbfpars)

!go into cubic space for the rbf likelihood
       rbfcube = cubize_rbfparams(fitNdim,rbfcuts)
       
       if (any(rbfcube.gt.1._fp).or.any(rbfcube.lt.0._fp)) then
!if outside rbffits box, the real likelihood is so small that we
!cannot calculate it numerically, but we can define a junk one smaller
!than fitLogZero
          lnew = fitLogZero * 4._fp*sum((rbfcube(:)-0.5_fp)**2)
       elseif (test_reheating_hardprior(mnpars)) then
! reheating hardprior, ignoring those points
          lnew = nestLogZero
       else
          lnew = rbflike_eval(rbfcube)
       endif

    endif

!dump the physical params into cube for file output
    cube(1:nestdim) = mnpars(1:nestdim)

!+ derived parameters we want to dump too
    do i=1,nestpars-nestdim
       cube(nestdim+i) = get_derived(i)
    enddo

    if (display) then
       write(*,*)
       write(*,*)'rbf_multinest_aspic_loglike: '
       write(*,*)'lnA= log(eps1)= eps_i=       ',rbfpars
       write(*,*)'ln(like) =                   ',lnew
       write(*,*)
    end if

  end subroutine rbf_multinest_aspic_loglike


  subroutine shep_multinest_aspic_loglike(cube,nestdim,nestpars,lnew,context)
    use shepprec, only : fp
    use sheplike, only : cubize_shepparams, uncubize_shepparams, check_shep
    use sheplike, only : sheplike_eval, cutmin_shepparams
    use wraspic, only : get_slowroll, get_derived
    use wraspic, only : test_aspic_hardprior, test_reheating_hardprior
    use nestparams, only : nestLogZero,fitLogZero
    implicit none   
    integer(imn) :: nestdim, nestpars
    real(fmn), dimension(nestpars) :: cube
    real(fmn) :: lnew
    integer(imn) :: context,i
    real(fmn), dimension(nestdim) :: mnpars
    real(fp), dimension(fitNdim) :: shepcube, sheppars, shepcuts


    if (.not.check_shep()) stop 'shep_multinest_loglike: not initialized!'

!get the physical parameters we are sampling on
    mnpars = uncubize_nestparams(nestdim,cube)

!check for any hard prior in aspic model parameters
    if (test_aspic_hardprior(mnpars)) then

       lnew = nestLogZero

    else
       
!use aspic to get the slowroll parameters
       sheppars = get_slowroll(fitNdim,mnpars)

!if eps1<eps1min, the likelihood is flat
       shepcuts = cutmin_shepparams(fitNdim,eps1pos,sheppars)

!go into cubic space for the shep likelihood
       shepcube = cubize_shepparams(fitNdim,shepcuts)
       
       if (any(shepcube.gt.1._fp).or.any(shepcube.lt.0._fp)) then
!if outside shepfits box, the real likelihood is so small that we
!cannot calculate it numerically, but we can define a junk one smaller
!than fitLogZero
          lnew = fitLogZero * 4._fp*sum((shepcube(:)-0.5_fp)**2)
       elseif (test_reheating_hardprior(mnpars)) then
! reheating hardprior, those points are ignored
          lnew = nestLogZero
       else
          lnew = sheplike_eval(shepcube)
       endif

    endif

!dump the physical params into cube for file output
    cube(1:nestdim) = mnpars(1:nestdim)

!+ derived parameters we want to dump too
    do i=1,nestpars-nestdim
       cube(nestdim+i) = get_derived(i)
    enddo

    if (display) then
       write(*,*)
       write(*,*)'shep_multinest_aspic_loglike: '
       write(*,*)'lnA= log(eps1)= eps_i=       ',sheppars
       write(*,*)'ln(like) =                   ',lnew
       write(*,*)
    end if

  end subroutine shep_multinest_aspic_loglike


#endif

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
    write(*,*)'*****************************************************'
    write(*,*)'nest_dumper: '
    write(*,*)'nSamples= logZ= logZerr= ',nSamples, logZ, logZerr
    write(*,*)'maxLogLike= ',maxLogLike
    write(*,*)'*****************************************************'
    write(*,*)
   

  end subroutine nest_dumper


end module nestwrap
