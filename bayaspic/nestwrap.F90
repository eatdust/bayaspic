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

  integer(imn), parameter :: numDerivedParams = 11

  integer(imn), save :: fitNdim = 0
  integer(imn), parameter :: pstarpos = 1
  integer(imn), parameter :: eps1pos = 2

  real(fmn), save, dimension(:), allocatable :: nestPmin, nestPmax

  logical, parameter :: display = .false.

!the name of the fastlike method we call

!  character(len=*), parameter :: fastLikeName = 'null'

! radial basis functions
!  character(len=*), parameter :: fastLikeName = 'rbf'

! inverse shepard method
  character(len=*), parameter :: fastLikeName = 'shep'

!nothing


  public nest_free_slowroll
  public nest_init_slowroll, nest_sample_slowroll
#if defined ASPIC || ASPICQ
  public nest_init_aspic, nest_sample_aspic, nest_free_aspic
#endif

contains


   subroutine nest_init_slowroll()
    use rbflike, only : initialize_rbf_like
    use rbflike, only : get_rbf_ndim, get_rbf_fmin
    use rbflike, only : get_rbf_xpmin, get_rbf_xpmax
    use sheplike, only : initialize_shep_like
    use sheplike, only : get_shep_ndim, get_shep_fmin
    use sheplike, only : get_shep_xpmin, get_shep_xpmax
    use nestparams, only : nestNdim, nestNpars, nestCdim
    use nestparams, only : nestPWrap, nestRootName, nestRootPrefix
    use nestparams, only : fitLogZero
    implicit none

    integer, parameter :: srlen = 3
    character(len=srlen) :: name

    integer :: i

    select case(fastLikeName)
    
    case ('rbf')

       call initialize_rbf_like(fileweights, filecentres, filerbfbounds)
       fitNdim = get_rbf_ndim()
       fitLogZero = get_rbf_fmin()
       nestNdim = fitNdim
       
       allocate(nestPmin(nestNdim))
       allocate(nestPmax(nestNdim))
       
       do i=1,fitNdim
          nestPmin(i) = get_rbf_xpmin(i)
          nestPmax(i) = get_rbf_xpmax(i)
       enddo

    case ('shep')

       call initialize_shep_like(fileshep, filepost, fileshepbounds)
       fitNdim = get_shep_ndim()
       fitLogZero = get_shep_fmin()
       nestNdim = fitNdim

       allocate(nestPmin(nestNdim))
       allocate(nestPmax(nestNdim))

       do i=1,fitNdim
          nestPmin(i) = get_shep_xpmin(i)
          nestPmax(i) = get_shep_xpmax(i)
       enddo

    case default

       stop 'nest_init_slowroll: fast like not found!'

    end select
       
    nestNpars = nestNdim
    nestCdim = nestNdim

    select case (fitNdim)
    case (3)
       name = 'sr2'
    case (4)
       name = 'sr3'
    case default
       stop 'nest_init_slowroll: sr model not known'
    end select

    

    allocate(nestPwrap(nestNdim))
    nestPwrap = 0

    nestRootName = trim(nestRootPrefix)//name
    call nest_dump_slowroll_priors(name)

    call nest_print()

  end subroutine nest_init_slowroll


  subroutine nest_dump_slowroll_priors(extname)
    use nestparams, only : nestNdim, nestRootname
    implicit none
    character(len=*), intent(in) :: extname
    integer, parameter :: nunit = 414, nname = 415, nrange=416
    
    if ((.not.allocated(nestPmin)).or.(.not.allocated(nestPmax))) then
       stop 'nest_dump_slowroll_priors: prior not allocated!'
    endif

    open(unit=nunit,file=trim(nestRootName)//'.ini', status='unknown')
    open(unit=nname,file=trim(nestRootName)//'.paramnames', status='unknown')
    open(unit=nrange,file=trim(nestRootName)//'.ranges', status='unknown')

    write(nunit,*)'limits[lnA]=',nestPmin(1),nestPmax(1)    
    write(nname,*)'lnA         \ln(10^{10} P_*)'
    write(nrange,*)'lnA                    ',nestPmin(1),nestPmax(1)

    select case (extname)

    case ('sr2')

       write(nunit,*)'limits[sr1]=',nestPmin(2),nestPmax(2)
       write(nname,*)'sr1         \log(\epsilon_1)'
       write(nunit,*)'sr1                    ',nestPmin(2),nestPmax(2)
       write(nunit,*)'limits[sr2]=',nestPmin(3),nestPmax(3)
       write(nname,*)'sr2         \epsilon_2'
       write(nunit,*)'sr2                    ',nestPmin(3),nestPmax(3)

    case ('sr3')
       write(nunit,*)'limits[sr1]=',nestPmin(2),nestPmax(2)
       write(nname,*)'sr1         \log(\epsilon_1)'
       write(nunit,*)'limits[sr2]=',nestPmin(3),nestPmax(3)
       write(nname,*)'sr2         \epsilon_2'
       write(nunit,*)'limits[sr3]=',nestPmin(4),nestPmax(4)
       write(nname,*)'sr3         \epsilon_3'
       write(nunit,*)'sr3                    ',nestPmin(4),nestPmax(4)

    case default

       stop 'nest_dump_slowrollpriors: internal error!'

    end select

    close(nunit)
    close(nname)
    close(nrange)

  end subroutine nest_dump_slowroll_priors


  subroutine nest_sample_slowroll()
    use nested, only : nestRun
    use nestparams
    implicit none

    integer(imn) :: nclusters			
    integer(imn) :: context
    integer(imn) :: maxNode 			

    select case (fastLikeName)

    case ('rbf')

       call nestRun(nestINS,nestMmodal,nestCteEff,nestNlive,nestZtol,nestSampEff,nestNdim,nestNpars, &
            nestCdim,nestMaxModes,nestUpdInt,nestNullZ,nestRootName,nestSeed,nestPwrap, &
            nestFeedBack,nestResume,nestOutfile,nestInitMPI,nestLogZero,nestMaxIter &
            ,rbf_multinest_slowroll_loglike,nest_dumper,context)

    case ('shep')

       call nestRun(nestINS,nestMmodal,nestCteEff,nestNlive,nestZtol,nestSampEff,nestNdim,nestNpars, &
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
    write(*,*)'nestINS      =         ',nestINS
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



  subroutine nest_free_slowroll()
    use nestparams, only : nestPwrap    
    implicit none
   
    if (allocated(nestPwrap)) deallocate(nestPwrap)
    if (allocated(nestpmin)) deallocate(nestpmin)
    if (allocated(nestpmax)) deallocate(nestpmax)

  end subroutine nest_free_slowroll

  

#if defined ASPIC || ASPICQ
  subroutine nest_init_aspic(modelname)
    use wraspic, only : set_model, get_ntot
    use wraspic, only : get_allprior
    use rbflike, only : initialize_rbf_like, check_rbf, get_rbf_fmin
    use rbflike, only : get_rbf_ndim, get_rbf_xpmin, get_rbf_xpmax
    use sheplike, only : initialize_shep_like, check_shep, get_shep_fmin
    use sheplike, only : get_shep_ndim, get_shep_xpmin, get_shep_xpmax
    use nestparams, only : nestNdim, nestNpars, nestCdim
    use nestparams, only : nestPWrap, nestRootName, nestRootPrefix
    use nestparams, only : fitLogZero, nestLogZero
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


    case ('null')

       fitNdim = 0
       fitLogZero = nestLogZero

    case default

       stop 'nest_init_aspic: fast like not found!'

    end select


    if (nestNdim.ne.fitNdim) then
       write(*,*)'nest_init_aspic: '
       write(*,*)'nestdim= fastlikedim= ',nestNdim, fitNdim    
    endif    
       
    allocate(nestPwrap(nestNdim))
    nestPwrap = 0

    call nest_dump_aspic_priors(trim(name)//subname)

    call nest_print()
    
  end subroutine nest_init_aspic

  

  subroutine nest_dump_aspic_priors(extname)
    use nestparams, only : nestNdim, nestRootname
    use wraspic, only : get_nextra, ReheatModel
    implicit none
    character(len=*), intent(in) :: extname
    integer, parameter :: nunit = 412, nname = 413, nrange = 414
    integer :: nextra
    
    if ((.not.allocated(nestPmin)).or.(.not.allocated(nestPmax))) then
       stop 'nest_dump_aspic_priors: prior not allocated!'
    endif

    nextra = get_nextra()

    open(unit=nunit,file=trim(nestRootName)//'.ini', status='unknown')
    open(unit=nname,file=trim(nestRootName)//'.paramnames', status='unknown')
    open(unit=nrange,file=trim(nestRootName)//'.ranges', status='unknown')

    write(nunit,*)'limits[lnA]=',nestPmin(1),nestPmax(1)
    write(nname,*)'lnA            \ln(10^{10} P_*)'
    write(nrange,*)'lnA                    ',nestPmin(1),nestPmax(1)

    if (nextra.eq.2) then

       select case (ReheatModel)
       case ('Rreh')
          write(nunit,*)'limits[lnRreh]=',nestPmin(2),nestPmax(2)
          write(nname,*)'lnRreh         \ln(R)'
          write(nrange,*)'lnRreh                 ',nestPmin(2),nestPmax(2)
       case ('Rrad')
          write(nunit,*)'limits[lnRrad]=',nestPmin(2),nestPmax(2)
          write(nname,*)'lnRrad         \ln(R_{\rm rad})'
          write(nrange,*)'lnRrad                 ',nestPmin(2),nestPmax(2)
       case default
          stop 'nest_dump_aspic_priors: internal error!'
       end select
          
    elseif (nextra.eq.3) then

       if (Reheatmodel.ne.'Rhow') stop 'ReheatModel is not Rhow!'

       write(nunit,*)'limits[lnRhoReh]=',nestPmin(2),nestPmax(2)
       write(nunit,*)'limits[wreh]=',nestPmin(3),nestPmax(3)
       write(nname,*)'lnRhoReh       \ln(\rho_{\rm reh})'
       write(nname,*)'wreh           \bar{w}_{\rm reh}'
       write(nrange,*)'lnRhoReh               ',nestPmin(2),nestPmax(2)
       write(nrange,*)'wreh                   ',nestPmin(3),nestPmax(3)

    else
       stop 'nest_dump_aspic_priors: nextra not found!'
    endif


    select case (nestNdim-nextra)
    case (0)
    case (1)
       write(nunit,*)'limits[c1]=',nestPmin(nextra+1),nestPmax(nextra+1)
       write(nname,*)'c1             c_1'
       write(nrange,*)'c1                     ',nestPmin(nextra+1),nestPmax(nextra+1)
    case(2)
       write(nunit,*)'limits[c1]=',nestPmin(nextra+1),nestPmax(nextra+1)
       write(nunit,*)'limits[c2]=',nestPmin(nextra+2),nestPmax(nextra+2)
       write(nname,*)'c1             c_1'
       write(nname,*)'c2             c_2'
       write(nrange,*)'c1                     ',nestPmin(nextra+1),nestPmax(nextra+1)
       write(nrange,*)'c2                     ',nestPmin(nextra+2),nestPmax(nextra+2)
    case (3)
       write(nunit,*)'limits[c1]=',nestPmin(nextra+1),nestPmax(nextra+1)
       write(nunit,*)'limits[c2]=',nestPmin(nextra+2),nestPmax(nextra+2)
       write(nunit,*)'limits[c3]=',nestPmin(nextra+3),nestPmax(nextra+3)
       write(nname,*)'c1             c_1'
       write(nname,*)'c2             c_2'
       write(nname,*)'c3             c_3'
       write(nrange,*)'c1                     ',nestPmin(nextra+1),nestPmax(nextra+1)
       write(nrange,*)'c2                     ',nestPmin(nextra+2),nestPmax(nextra+2)
       write(nrange,*)'c3                     ',nestPmin(nextra+3),nestPmax(nextra+3)
    case (4)
       write(nunit,*)'limits[c1]=',nestPmin(nextra+1),nestPmax(nextra+1)
       write(nunit,*)'limits[c2]=',nestPmin(nextra+2),nestPmax(nextra+2)
       write(nunit,*)'limits[c3]=',nestPmin(nextra+3),nestPmax(nextra+3)
       write(nunit,*)'limits[c4]=',nestPmin(nextra+4),nestPmax(nextra+4)
       write(nname,*)'c1             c_1'
       write(nname,*)'c2             c_2'
       write(nname,*)'c3             c_3'
       write(nname,*)'c4             c_4'
       write(nrange,*)'c1                     ',nestPmin(nextra+1),nestPmax(nextra+1)
       write(nrange,*)'c2                     ',nestPmin(nextra+2),nestPmax(nextra+2)
       write(nrange,*)'c3                     ',nestPmin(nextra+3),nestPmax(nextra+3)
       write(nrange,*)'c4                     ',nestPmin(nextra+4),nestPmax(nextra+4)
    case default
       stop 'nest_dump_aspic_priors: case not implemented!'
    end select

    close(nunit)
    close(nname)
    close(nrange)

  end subroutine nest_dump_aspic_priors



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

       call nestRun(nestINS,nestMmodal,nestCteEff,nestNlive,nestZtol,nestSampEff,nestNdim,nestNpars, &
            nestCdim,nestMaxModes,nestUpdInt,nestNullZ,nestRootName,nestSeed,nestPwrap, &
            nestFeedBack,nestResume,nestOutfile,nestInitMPI,nestLogZero,nestMaxIter &
            ,rbf_multinest_aspic_loglike,nest_dumper,context)

    case ('shep')

       call nestRun(nestINS,nestMmodal,nestCteEff,nestNlive,nestZtol,nestSampEff,nestNdim,nestNpars, &
            nestCdim,nestMaxModes,nestUpdInt,nestNullZ,nestRootName,nestSeed,nestPwrap, &
            nestFeedBack,nestResume,nestOutfile,nestInitMPI,nestLogZero,nestMaxIter &
            ,shep_multinest_aspic_loglike,nest_dumper,context)

    case ('null')

       call nestRun(nestINS,nestMmodal,nestCteEff,nestNlive,nestZtol,nestSampEff,nestNdim,nestNpars, &
            nestCdim,nestMaxModes,nestUpdInt,nestNullZ,nestRootName,nestSeed,nestPwrap, &
            nestFeedBack,nestResume,nestOutfile,nestInitMPI,nestLogZero,nestMaxIter &
            ,null_multinest_aspic_loglike,nest_dumper,context)

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

! reheating hardprior, ignoring those points
       if (test_reheating_hardprior(mnpars)) then
          
          lnew = nestLogZero

       else

!if eps1<eps1min, the likelihood is flat
          rbfcuts = cutmin_rbfparams(fitNdim,eps1pos,rbfpars)

!go into cubic space for the rbf likelihood
          rbfcube = cubize_rbfparams(fitNdim,rbfcuts)
       
          if (any(rbfcube.gt.1._fp).or.any(rbfcube.lt.0._fp)) then
!if outside rbffits box, the real likelihood is so small that we
!cannot calculate it numerically, but we can define a junk one smaller
!than fitLogZero
             lnew = fitLogZero * 4._fp*sum((rbfcube(:)-0.5_fp)**2)

          else

             lnew = rbflike_eval(rbfcube)

          endif

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

! reheating hardprior, those points are ignored
       if (test_reheating_hardprior(mnpars)) then

          lnew = nestLogZero

       else

!if eps1<eps1min, the likelihood is flat
          shepcuts = cutmin_shepparams(fitNdim,eps1pos,sheppars)

!go into cubic space for the shep likelihood
          shepcube = cubize_shepparams(fitNdim,shepcuts)
       
          if (any(shepcube.gt.1._fp).or.any(shepcube.lt.0._fp)) then
!if outside shepfits box, the real likelihood is so small that we
!cannot calculate it numerically, but we can define a junk one smaller
!than fitLogZero
             lnew = fitLogZero * 4._fp*sum((shepcube(:)-0.5_fp)**2)
       
          else

             lnew = sheplike_eval(shepcube)

          endif

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



  subroutine null_multinest_aspic_loglike(cube,nestdim,nestpars,lnew,context)    
    use nestparams, only : nestLogZero
    use wraspic, only : get_derived, test_aspic_hardprior, test_reheating_hardprior
    use wraspic, only : get_slowroll
    implicit none   
    integer(imn) :: nestdim, nestpars
    real(fmn), dimension(nestpars) :: cube
    real(fmn) :: lnew
    integer(imn) :: context,i
    real(fmn), dimension(nestdim) :: mnpars    
    real(fmn), dimension(fitNdim) :: nullpars

!get the physical parameters we are sampling on
    mnpars = uncubize_nestparams(nestdim,cube)

!check for any hard prior in aspic model parameters
    if (test_aspic_hardprior(mnpars)) then

       lnew = nestLogZero

    else

!to get derived parameters
       nullpars = get_slowroll(fitNdim,mnpars)

       if (test_reheating_hardprior(mnpars)) then
! reheating hardprior, those points are ignored
          lnew = nestLogZero
       else
          lnew = 1._fmn
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
       write(*,*)'null_aspic_loglike: '       
       write(*,*)'ln(like) =                   ',lnew
       write(*,*)
    end if

  end subroutine null_multinest_aspic_loglike


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
       , maxLogLike, logZ, INSlogZ, logZerr, context)
    use nestparams, only : nestRootName, nestINS
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


end module nestwrap
