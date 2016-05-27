module samplaspic
  use sampl, only : imn, fmn, lenmn
  use sampl, only : samplNdim, samplNpars
  use sampl, only : rootDir, rootName, rootPrefix
  use sampl, only : init_samplparams, free_samplparams
  use sampl, only : uncubize_samplparams
  use sampl, only : samplPmin, samplPmax
  implicit none

  private

  logical, parameter :: display = .true.
  logical, parameter :: debug = .false.


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

  real(fmn), save :: fitLogZero = 0._fmn


!the name of the fastlike method we call

!  character(len=*), parameter :: fastLikeName = 'null'

! radial basis functions
!  character(len=*), parameter :: fastLikeName = 'rbf'

! inverse shepard method
  character(len=*), parameter :: fastLikeName = 'shep'

  public nest_init_aspic, nest_sample_aspic, nest_free_aspic
  public chord_init_aspic, chord_sample_aspic, chord_free_aspic

contains


   
  subroutine init_aspic(modelname)
    use wraspic, only : set_model, get_ntot
    use wraspic, only : get_allprior
    use rbflike, only : initialize_rbf_like, check_rbf, get_rbf_fmin
    use rbflike, only : get_rbf_ndim, get_rbf_xpmin, get_rbf_xpmax
    use sheplike, only : initialize_shep_like, check_shep, get_shep_fmin
    use sheplike, only : get_shep_ndim, get_shep_xpmin, get_shep_xpmax
    implicit none    
    integer :: cpos, lenmod
    character(len=*), intent(in) :: modelname
    character(len=len(modelname)) :: name, subname

    character, parameter :: separator=' '

    lenmod = len(modelname)

    cpos = scan(modelname,separator)

    if (cpos.eq.0) then
       name = trim(adjustl(modelname))
       rootName = trim(rootPrefix)//name
       call set_model(name)

    else
       name = trim(adjustl(modelname(1:cpos-1)))
       subname = trim(adjustl(modelname(cpos:lenmod)))
       rootName = trim(rootPrefix)//trim(name)//subname
       call set_model(name,subname)

    endif

   
    samplNdim = get_ntot()
    samplNpars = samplNdim + numDerivedParams

    allocate(samplPmin(samplNdim))
    allocate(samplPmax(samplNdim))

    call get_allprior(samplPmin, samplPmax)

    select case (fastLikeName)

    case ('rbf')

       if (.not.check_rbf()) then
          call initialize_rbf_like(fileweights, filecentres, filerbfbounds)
       endif

!cut prior of P* by the one encoded in the likelihood
       samplPmin(pstarpos) = max(get_rbf_xpmin(pstarpos),samplPmin(pstarpos))
       samplPmax(pstarpos) = min(get_rbf_xpmax(pstarpos),samplPmax(pstarpos))
    
       fitNdim = get_rbf_ndim()
       fitLogZero = get_rbf_fmin()

    case ('shep')

       if (.not.check_shep()) then
          call initialize_shep_like(fileshep, filepost, fileshepbounds)
       endif

!cut prior of P* by the one encoded in the likelihood
       samplPmin(pstarpos) = max(get_shep_xpmin(pstarpos),samplPmin(pstarpos))
       samplPmax(pstarpos) = min(get_shep_xpmax(pstarpos),samplPmax(pstarpos))
    
       fitNdim = get_shep_ndim()
       fitLogZero = get_shep_fmin()


    case ('null')

       fitNdim = 0
       fitLogZero = -1._fmn

    case default

       stop 'nest_init_aspic: fast like not found!'

    end select


    if (samplNdim.ne.fitNdim) then
       write(*,*)'init_aspic: '
       write(*,*)'samplingdim= fastlikedim= ',samplNdim, fitNdim    
    endif    
       
    call dump_aspic_priors(trim(name)//subname)

    
    if (display) then
       write(*,*)
       write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(*,*)'fast like is :         ',fastLikeName
       write(*,*)'lnZeroMin    =         ',fitLogZero
       write(*,*)'fitNdim      =         ',fitNdim
       write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(*,*)
    end if
    
  end subroutine init_aspic

  

  subroutine dump_aspic_priors(extname)
    use wraspic, only : get_nextra, ReheatModel
    implicit none
    character(len=*), intent(in) :: extname
    character(len=lenmn) :: fileName
    integer, parameter :: nunit = 412, nname = 413, nrange = 414
    integer :: nextra
    
    if ((.not.allocated(samplPmin)).or.(.not.allocated(samplPmax))) then
       stop 'dump_aspic_priors: prior not allocated!'
    endif

    nextra = get_nextra()

    fileName=trim(rootDir)//trim(rootName)

    open(unit=nunit,file=trim(fileName)//'.ini', status='unknown')
    open(unit=nname,file=trim(fileName)//'.paramnames', status='unknown')
    open(unit=nrange,file=trim(fileName)//'.ranges', status='unknown')

    write(nunit,*)'limits[lnA]=',samplPmin(1),samplPmax(1)
    write(nname,*)'lnA            \ln(10^{10} P_*)'
    write(nrange,*)'lnA                    ',samplPmin(1),samplPmax(1)

    if (nextra.eq.2) then

       select case (ReheatModel)
       case ('Rreh')
          write(nunit,*)'limits[lnRreh]=',samplPmin(2),samplPmax(2)
          write(nname,*)'lnRreh         \ln(R)'
          write(nrange,*)'lnRreh                 ',samplPmin(2),samplPmax(2)
       case ('Rrad')
          write(nunit,*)'limits[lnRrad]=',samplPmin(2),samplPmax(2)
          write(nname,*)'lnRrad         \ln(R_{\rm rad})'
          write(nrange,*)'lnRrad                 ',samplPmin(2),samplPmax(2)
       case default
          stop 'dump_aspic_priors: internal error!'
       end select
          
    elseif (nextra.eq.3) then

       if (Reheatmodel.ne.'Rhow') stop 'ReheatModel is not Rhow!'

       write(nunit,*)'limits[lnRhoReh]=',samplPmin(2),samplPmax(2)
       write(nunit,*)'limits[wreh]=',samplPmin(3),samplPmax(3)
       write(nname,*)'lnRhoReh       \ln(\rho_{\rm reh})'
       write(nname,*)'wreh           \bar{w}_{\rm reh}'
       write(nrange,*)'lnRhoReh               ',samplPmin(2),samplPmax(2)
       write(nrange,*)'wreh                   ',samplPmin(3),samplPmax(3)

    else
       stop 'dump_aspic_priors: nextra not found!'
    endif


    select case (samplNdim-nextra)
    case (0)
    case (1)
       write(nunit,*)'limits[c1]=',samplPmin(nextra+1),samplPmax(nextra+1)
       write(nname,*)'c1             c_1'
       write(nrange,*)'c1                     ',samplPmin(nextra+1),samplPmax(nextra+1)
    case(2)
       write(nunit,*)'limits[c1]=',samplPmin(nextra+1),samplPmax(nextra+1)
       write(nunit,*)'limits[c2]=',samplPmin(nextra+2),samplPmax(nextra+2)
       write(nname,*)'c1             c_1'
       write(nname,*)'c2             c_2'
       write(nrange,*)'c1                     ',samplPmin(nextra+1),samplPmax(nextra+1)
       write(nrange,*)'c2                     ',samplPmin(nextra+2),samplPmax(nextra+2)
    case (3)
       write(nunit,*)'limits[c1]=',samplPmin(nextra+1),samplPmax(nextra+1)
       write(nunit,*)'limits[c2]=',samplPmin(nextra+2),samplPmax(nextra+2)
       write(nunit,*)'limits[c3]=',samplPmin(nextra+3),samplPmax(nextra+3)
       write(nname,*)'c1             c_1'
       write(nname,*)'c2             c_2'
       write(nname,*)'c3             c_3'
       write(nrange,*)'c1                     ',samplPmin(nextra+1),samplPmax(nextra+1)
       write(nrange,*)'c2                     ',samplPmin(nextra+2),samplPmax(nextra+2)
       write(nrange,*)'c3                     ',samplPmin(nextra+3),samplPmax(nextra+3)
    case (4)
       write(nunit,*)'limits[c1]=',samplPmin(nextra+1),samplPmax(nextra+1)
       write(nunit,*)'limits[c2]=',samplPmin(nextra+2),samplPmax(nextra+2)
       write(nunit,*)'limits[c3]=',samplPmin(nextra+3),samplPmax(nextra+3)
       write(nunit,*)'limits[c4]=',samplPmin(nextra+4),samplPmax(nextra+4)
       write(nname,*)'c1             c_1'
       write(nname,*)'c2             c_2'
       write(nname,*)'c3             c_3'
       write(nname,*)'c4             c_4'
       write(nrange,*)'c1                     ',samplPmin(nextra+1),samplPmax(nextra+1)
       write(nrange,*)'c2                     ',samplPmin(nextra+2),samplPmax(nextra+2)
       write(nrange,*)'c3                     ',samplPmin(nextra+3),samplPmax(nextra+3)
       write(nrange,*)'c4                     ',samplPmin(nextra+4),samplPmax(nextra+4)
    case default
       stop 'dump_aspic_priors: case not implemented!'
    end select

    close(nunit)
    close(nname)
    close(nrange)

  end subroutine dump_aspic_priors



  subroutine nest_init_aspic(modelname)
    use nestparams, only : nestNdim, nestNpars, nestCdim
    use nestparams, only : nestPWrap, nestRootName, nest_print
    use sampl, only : sampl_print
    implicit none
    character(len=*), intent(in) :: modelname
    
    call init_aspic(modelname)

    nestNdim = samplNdim
    nestNpars = samplNpars
    nestCdim = samplNdim

    nestRootName = trim(rootDir)//trim(rootName)

    allocate(nestPwrap(nestNdim))
    nestPwrap = 0

    call nest_print()

    call sampl_print()

  end subroutine nest_init_aspic



  
 
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
    use nestparams, only : nestLogZero
    implicit none   
    integer(imn) :: nestdim, nestpars
    real(fmn), dimension(nestpars) :: cube
    real(fmn) :: lnew
    integer(imn) :: context,i
    real(fmn), dimension(nestdim) :: mnpars
    real(fp), dimension(fitNdim) :: rbfcube, rbfpars, rbfcuts


    if (.not.check_rbf()) stop 'rbf_multinest_loglike: not initialized!'

!get the physical parameters we are sampling on
    mnpars = uncubize_samplparams(nestdim,cube)

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

    if (debug) then
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
    use nestparams, only : nestLogZero
    implicit none   
    integer(imn) :: nestdim, nestpars
    real(fmn), dimension(nestpars) :: cube
    real(fmn) :: lnew
    integer(imn) :: context,i
    real(fmn), dimension(nestdim) :: mnpars
    real(fp), dimension(fitNdim) :: shepcube, sheppars, shepcuts


    if (.not.check_shep()) stop 'shep_multinest_loglike: not initialized!'

!get the physical parameters we are sampling on
    mnpars = uncubize_samplparams(nestdim,cube)

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

    if (debug) then
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
    mnpars = uncubize_samplparams(nestdim,cube)

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

    if (debug) then
       write(*,*)
       write(*,*)'null_aspic_loglike: '       
       write(*,*)'ln(like) =                   ',lnew
       write(*,*)
    end if

  end subroutine null_multinest_aspic_loglike

 


  subroutine nest_free_aspic()
    use nestparams, only : nestPwrap
    use wraspic, only : free_model
    implicit none

    call free_model()
    
    if (allocated(nestPwrap)) deallocate(nestPwrap)

    call free_samplparams()

  end subroutine nest_free_aspic




  subroutine chord_init_aspic(modelname)    
    use chordparams, only : chordNdim, chordNpars, chordNderived
    use chordparams, only : chordName, chordDir, chord_print
    use sampl, only : sampl_print
    implicit none
    character(len=*), intent(in) :: modelname

    call init_aspic(modelname)

    chordNdim = samplNdim
    chordNderived = numDerivedParams
    chordNpars = samplNpars

    chordDir=rootDir
    chordName = trim(rootName)
   
    call chord_print()

    call sampl_print()

  end subroutine chord_init_aspic



  subroutine chord_sample_aspic()
    use interfaces_module, only : run_polychord
    use settings_module, only : program_settings
    use chordparams
    implicit none
    type(program_settings) :: runset

    call chord_settings(runset)

    select case (fastLikeName)

    case ('rbf')

       call run_polychord(rbf_polychord_aspic_loglike,polychord_prior,runset)

    case ('shep')

       call run_polychord(shep_polychord_aspic_loglike,polychord_prior,runset)


    case ('null')
       
       call run_polychord(null_polychord_aspic_loglike,polychord_prior,runset)

    case default

       stop 'chord_sample_aspic: fast like not found!'

    end select

  end subroutine chord_sample_aspic


  function polychord_prior(cube)
    implicit none
    real(fmn), dimension(:), intent(in) :: cube
    real(fmn), dimension(size(cube,1)) :: polychord_prior

    integer(imn) :: ndim
    ndim = size(cube,1)

    polychord_prior = uncubize_samplparams(ndim,cube)

  end function polychord_prior



  function rbf_polychord_aspic_loglike(theta,phi) result(lnew)
    use chordparams, only : chordLogZero, chordNderived, chordNdim
    use rbfprec, only : fp
    use rbflike, only : cubize_rbfparams, check_rbf
    use rbflike, only : rbflike_eval, cutmin_rbfparams
    use wraspic, only : get_slowroll, get_derived
    use wraspic, only : test_aspic_hardprior, test_reheating_hardprior
    implicit none   
    real(fmn) :: lnew
    real(fmn), dimension(:), intent(in) :: theta
    real(fmn), dimension(:), intent(out) :: phi

    real(fp), dimension(fitNdim) :: rbfcube, rbfcuts, rbfpars
    integer(imn) :: ndim, i

    ndim = size(theta,1)

    if (.not.check_rbf()) stop 'rbf_polychord_aspic_loglike: not initialized!'
    if (ndim.ne.chordNdim) stop 'rbf_polychord_aspic_loglike: dim mismatch!'

!check for any hard prior in aspic model parameters
    if (test_aspic_hardprior(theta)) then

       lnew = chordLogZero

    else
       
!use aspic to get the slowroll parameters
       rbfpars = get_slowroll(fitNdim,theta)

! reheating hardprior, ignoring those points
       if (test_reheating_hardprior(theta)) then
          
          lnew = chordLogZero

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

!+ derived parameters we want to dump too
    do i=1,chordNderived
       phi(i) = get_derived(i)
    enddo

    if (debug) then
       write(*,*)
       write(*,*)'rbf_polychord_aspic_loglike: '
       write(*,*)'lnA= log(eps1)= eps_i=       ',rbfpars
       write(*,*)'ln(like) =                   ',lnew
       write(*,*)
    end if

  end function rbf_polychord_aspic_loglike
  



  function shep_polychord_aspic_loglike(theta,phi) result(lnew)
    use chordparams, only : chordLogZero, chordNderived, chordNdim
    use shepprec, only : fp
    use sheplike, only : cubize_shepparams, check_shep
    use sheplike, only : sheplike_eval, cutmin_shepparams
    use wraspic, only : get_slowroll, get_derived
    use wraspic, only : test_aspic_hardprior, test_reheating_hardprior
    implicit none   
    real(fmn) :: lnew
    real(fmn), dimension(:), intent(in) :: theta
    real(fmn), dimension(:), intent(out) :: phi

    integer(imn) :: ndim, i
    real(fp), dimension(fitNdim) :: shepcube, sheppars, shepcuts

    ndim = size(theta,1)
    if (.not.check_shep()) stop 'shep_polychord_aspic_loglike: not initialized!'
    if (ndim.ne.chordNdim) stop 'shep_polychord_aspic_loglike: dim mismatch!'

!check for any hard prior in aspic model parameters
    if (test_aspic_hardprior(theta)) then

       lnew = chordLogZero

    else
       
!use aspic to get the slowroll parameters
       sheppars = get_slowroll(fitNdim,theta)

! reheating hardprior, those points are ignored
       if (test_reheating_hardprior(theta)) then

          lnew = chordLogZero

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


!+ derived parameters we want to dump too
    do i=1,chordNderived
       phi(i) = get_derived(i)
    enddo


    if (debug) then
       write(*,*)
       write(*,*)'shep_polychord_aspic_loglike: '
       write(*,*)'lnA= log(eps1)= eps_i=       ',sheppars
       write(*,*)'ln(like) =                   ',lnew
       write(*,*)
    end if

  end function shep_polychord_aspic_loglike



  

  function null_polychord_aspic_loglike(theta,phi)    
    use chordparams, only : chordLogZero, chordNderived, chordNdim
    use wraspic, only : get_derived, test_aspic_hardprior, test_reheating_hardprior
    use wraspic, only : get_slowroll
    implicit none   
    real(fmn) :: null_polychord_aspic_loglike
    real(fmn), dimension(:), intent(in) :: theta
    real(fmn), dimension(:), intent(out) :: phi

    real(fmn) :: lnew
    integer(imn) :: ndim, i
    real(fmn), dimension(size(theta,1)) :: nullpars
   
    ndim = size(theta,1)
    if (ndim.ne.chordNdim) stop 'null_polychord_aspic_loglike: dim mismatch!'

!check for any hard prior in aspic model parameters
    if (test_aspic_hardprior(theta)) then

       lnew = chordLogZero

    else

!to get derived parameters
       nullpars = get_slowroll(fitNdim,theta)

       if (test_reheating_hardprior(theta)) then
! reheating hardprior, those points are ignored
          lnew = chordLogZero
       else
          lnew = 1._fmn
       endif
       
    endif


!+ derived parameters we want to dump too
    do i=1,chordNderived
       phi(i) = get_derived(i)
    enddo

    if (debug) then
       write(*,*)
       write(*,*)'null_polychord_aspic_loglike: '
       write(*,*)'lnA= log(eps1)= eps_i=       ',nullpars
       write(*,*)'ln(like) =                   ',lnew
       write(*,*)
    end if


  end function  null_polychord_aspic_loglike




  subroutine chord_free_aspic()
    use wraspic, only : free_model
    implicit none

    call free_model()       
    call free_samplparams()

  end subroutine chord_free_aspic


end module samplaspic
