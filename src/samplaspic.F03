!   This file is part of bayaspic
!
!   Copyright (C) 2021 C. Ringeval
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
!   along with Foobar.  If not, see <https://www.gnu.org/licenses/>.


module samplaspic
  use sampl, only : imn, fmn, lenmn
  use sampl, only : samplNdim, samplNpars
  use sampl, only : rootDir, rootName, rootPrefix
  use sampl, only : init_samplparams, free_samplparams
  use sampl, only : uncubize_samplparams
  use sampl, only : samplPmin, samplPmax
  use sampl, only : mn_loglike, pc_loglike, pc_prior
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


  character(len=*), parameter :: filefnn = 'fnndata/fnndata.dat'
  character(len=*), parameter :: filefnnbounds = 'fnndata/bounds.dat'
  
  integer(imn), save :: numDerivedParams = 0

  integer(imn), save :: fitNdim = 0

  integer(imn), parameter :: pstarpos = 1
  integer(imn), parameter :: eps1pos = 2
  integer(imn), parameter :: eps3pos = 4
  
  
  real(fmn), save :: fitLogZero = 0._fmn

  integer, parameter :: lenlabel = 16
  character(len=lenlabel) :: label

  integer, parameter :: countmax = 2000000
  logical, parameter :: cuteps3 = .false.
  
  
!the name of the fastlike method we call

#ifdef LIKENULL
  character(len=*), parameter :: fastLikeName = 'null'

#elif LIKERBF
! radial basis functions
  character(len=*), parameter :: fastLikeName = 'rbf'
#elif LIKEFANN
!feedforward neural network
  character(len=*), parameter :: fastLikeName = 'fnn'
#else
! inverse shepard method
  character(len=*), parameter :: fastLikeName = 'shep'
#endif
    
  
  public nest_init_aspic, nest_sample_aspic, nest_free_aspic
  public chord_init_aspic, chord_sample_aspic, chord_free_aspic

contains


   
  subroutine init_aspic(modelname)
    use wraspic, only : set_model, get_ntot, get_nderived
    use wraspic, only : allocate_and_set_allprior
    use rbflike, only : initialize_rbf_like, check_rbf, get_rbf_fmin
    use rbflike, only : get_rbf_ndim, get_rbf_xpmin, get_rbf_xpmax
    use sheplike, only : initialize_shep_like, check_shep, get_shep_fmin
    use sheplike, only : get_shep_ndim, get_shep_xpmin, get_shep_xpmax
#ifndef NOFANN
    use fnnlike, only : initialize_fnn_like, check_fnn, get_fnn_fmin
    use fnnlike, only : get_fnn_ndim, get_fnn_xpmin, get_fnn_xpmax
#endif
    
    implicit none    
    integer :: cpos, lenmod
    character(len=*), intent(in) :: modelname
    character(len=len(modelname)) :: name, subname

    character, parameter :: separator=' '

    label = trim(modelname)
    
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

    call allocate_and_set_allprior(samplPmin,samplPmax)

!sanity checks    
    if (.not.(allocated(samplPmin)).or.(.not.(allocated(samplPmax)))) then
       stop 'init_aspic: samplPmin/max allocation failed!'
    endif

    samplNdim = get_ntot()
    
    if ((size(samplPmin,1).ne.samplNdim).or.(size(samplPmax,1).ne.samplNdim)) then
       stop 'init_aspic: samplPmin/max size not matching ntot!'
    endif

    numDerivedParams = get_nderived()
    samplNpars = samplNdim + numDerivedParams


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

#ifndef NOFANN
    case ('fnn')

       if (.not.check_fnn()) then
          call initialize_fnn_like(filefnn, filefnnbounds)
       endif

!cut prior of P* by the one encoded in the likelihood
       samplPmin(pstarpos) = max(get_fnn_xpmin(pstarpos),samplPmin(pstarpos))
       samplPmax(pstarpos) = min(get_fnn_xpmax(pstarpos),samplPmax(pstarpos))

       fitNdim = get_fnn_ndim()
       fitLogZero = get_fnn_fmin()
#endif

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
    use wraspic, only : get_nextra, ReheatModel, get_derived_name
    implicit none
    character(len=*), intent(in) :: extname
    character(len=lenmn) :: fileName
    integer, parameter :: nunit = 412, nname = 413, nrange = 414
    integer :: i, nextra, nderived
    
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
          write(nname,*)'lnRreh         \ln(R_{\mathrm{reh}})'
          write(nrange,*)'lnRreh                 ',samplPmin(2),samplPmax(2)
       case ('Rrad')
          write(nunit,*)'limits[lnRrad]=',samplPmin(2),samplPmax(2)
          write(nname,*)'lnRrad         \ln(R_{\mathrm{rad}})'
          write(nrange,*)'lnRrad                 ',samplPmin(2),samplPmax(2)
       case default
          stop 'dump_aspic_priors: internal error!'
       end select
          
    elseif (nextra.eq.3) then

       if (Reheatmodel.ne.'Rhow') stop 'ReheatModel is not Rhow!'

       write(nunit,*)'limits[lnRhoReh]=',samplPmin(2),samplPmax(2)
       write(nunit,*)'limits[wreh]=',samplPmin(3),samplPmax(3)
       write(nname,*)'lnRhoReh       \ln(\rho_{\mathrm{reh}}/M_{\mathrm{Pl}}^4)'
       write(nname,*)'wreh           \bar{w}_{\mathrm{reh}}'
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

    nderived = samplNpars - samplNdim
    
    do i=1,nderived
       write(nname,*) get_derived_name(i)
    enddo

    
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
    
    procedure(mn_loglike), pointer :: ptrnest_aspic_loglike => null()

    integer(imn) :: nclusters			
    integer(imn) :: context
    integer(imn) :: maxNode 			

    select case (fastLikeName)

    case ('rbf')
       ptrnest_aspic_loglike => rbf_multinest_aspic_loglike

    case ('shep')
       ptrnest_aspic_loglike => shep_multinest_aspic_loglike

#ifndef NOFANN       
    case ('fnn')
       ptrnest_aspic_loglike => fnn_multinest_aspic_loglike
#endif
       
    case ('null')
       ptrnest_aspic_loglike => null_multinest_aspic_loglike

    case default

       stop 'nest_sample_aspic: fast like not found!'

    end select

    call nestRun(nestINS,nestMmodal,nestCteEff,nestNlive,nestTol,nestSampEff,nestNdim,nestNpars, &
         nestCdim,nestMaxModes,nestUpdInt,nestNullZ,nestRootName,nestSeed,nestPwrap, &
         nestFeedBack,nestResume,nestOutfile,nestInitMPI,nestLogZero,nestMaxIter &
         ,ptrnest_aspic_loglike,nest_dumper,context)

    if (associated(ptrnest_aspic_loglike)) ptrnest_aspic_loglike => null()


  end subroutine nest_sample_aspic



  subroutine rbf_multinest_aspic_loglike(cube,nestdim,nestpars,lnew,context)
    use rbfprec, only : fp
    use rbflike, only : cubize_rbfparams, uncubize_rbfparams, check_rbf
    use rbflike, only : rbflike_eval, cutmin_rbfparams, cutmax_rbfparams
    use wraspic, only : get_hubbleflow, get_derived
    use wraspic, only : test_aspic_hardprior, test_reheating_hardprior
    use nestparams, only : nestLogZero
    implicit none   
    integer(imn) :: nestdim, nestpars
    real(fmn), dimension(nestpars) :: cube
    real(fmn) :: lnew
    integer(imn) :: context,i
    real(fmn), dimension(nestdim) :: mnpars
    real(fp), dimension(fitNdim) :: rbfcube, rbfpars, rbfcuts, rbfclamp

    integer, save :: counthard = 0
    logical :: disfavour
    
    if (.not.check_rbf()) stop 'rbf_multinest_loglike: not initialized!'
        
!get the physical parameters we are sampling on
    mnpars(1:nestdim) = uncubize_samplparams(nestdim,cube(1:nestdim))

!check for any hard prior in aspic model parameters
    if (test_aspic_hardprior(mnpars,disfavour)) then

       if (disfavour) then
          lnew = nestLogZero / 10._fmn
       else
          lnew = nestLogZero 
       endif

       counthard = counthard + 1
       
    else
       
!use aspic to get the slowroll parameters
       rbfpars = get_hubbleflow(fitNdim,mnpars)

! reheating hardprior, ignoring those points
       if (test_reheating_hardprior(mnpars)) then
          
          lnew = nestLogZero
          counthard = counthard + 1
          
       else

          counthard = 0
          
!if eps1<eps1min, the likelihood is flat
          rbfcuts = cutmin_rbfparams(fitNdim,eps1pos,rbfpars)

!In case the posterior of eps3 is not going to zero on the boundaries          
          if (cuteps3) then
             rbfcuts = cutmin_rbfparams(fitNdim,eps3pos,rbfcuts)
             rbfcuts = cutmax_rbfparams(fitNdim,eps3pos,rbfcuts)
          endif
             
!go into cubic space for the rbf likelihood
          rbfcube = cubize_rbfparams(fitNdim,rbfcuts)
       
          if (any(rbfcube.gt.1._fp).or.any(rbfcube.lt.0._fp)) then
!if outside rbffits box, the real likelihood is so small that we
!cannot calculate it numerically, but we can define a junk one smaller

             rbfclamp = max(0._fp,min(rbfcube,1._fp))
             
             lnew = rbflike_eval(rbfclamp) - sum((rbfcube-rbfclamp)**2)
             
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

    if (isnan(lnew)) stop 'rbf_multinest_aspic_loglike: like is NaN!'

    if (counthard.gt.countmax) then
       write(*,*)'rbf_multinest_aspic_loglike: trapped into hard prior regions'
       write(*,*)'model is: ',label
       write(*,*)'mnparams= ',mnpars(1:nestdim)
       stop
    endif
    
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
    use sheplike, only : sheplike_eval, cutmin_shepparams, cutmax_shepparams
    use wraspic, only : get_hubbleflow, get_derived
    use wraspic, only : test_aspic_hardprior, test_reheating_hardprior
    use nestparams, only : nestLogZero
    implicit none   
    integer(imn) :: nestdim, nestpars
    real(fmn), dimension(nestpars) :: cube
    real(fmn) :: lnew
    integer(imn) :: context,i
    real(fmn), dimension(nestdim) :: mnpars
    real(fp), dimension(fitNdim) :: shepcube, sheppars, shepcuts, shepclamp

    integer, save :: counthard = 0
    logical :: disfavour
    
    if (.not.check_shep()) stop 'shep_multinest_loglike: not initialized!'       
    
!get the physical parameters we are sampling on
    mnpars(1:nestdim) = uncubize_samplparams(nestdim,cube(1:nestdim),label)
    
!check for any hard prior in aspic model parameters
    if (test_aspic_hardprior(mnpars, disfavour)) then

       if (disfavour) then
          lnew = nestLogZero / 10._fmn
       else
          lnew = nestLogZero
       endif
       
       counthard = counthard + 1
       
    else
       
!use aspic to get the slowroll parameters
       sheppars = get_hubbleflow(fitNdim,mnpars)

! reheating hardprior, those points are ignored
       if (test_reheating_hardprior(mnpars)) then

          lnew = nestLogZero
          counthard = counthard + 1
       else

          counthard = 0
          
!if eps1<eps1min, the likelihood is flat
          shepcuts = cutmin_shepparams(fitNdim,eps1pos,sheppars)

!In case the posterior of eps3 is not going to zero on the boundaries          
          if (cuteps3) then
             shepcuts = cutmin_shepparams(fitNdim,eps3pos,shepcuts)
             shepcuts = cutmax_shepparams(fitNdim,eps3pos,shepcuts)
          endif
          
!go into cubic space for the shep likelihood
          shepcube = cubize_shepparams(fitNdim,shepcuts)
       
          if (any(shepcube.gt.1._fp).or.any(shepcube.lt.0._fp)) then
!if outside shepfits box, the real likelihood is so small that we
!cannot calculate it numerically, but we can define a junk one smaller

             shepclamp = max(0._fp,min(shepcube,1._fp))
             
             lnew = sheplike_eval(shepclamp) -  sum((shepcube-shepclamp)**2)
             
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

!sanity check
    if (isnan(lnew)) stop 'shep_multinest_aspic_loglike: like is NaN!'

    if (counthard.gt.countmax) then
       write(*,*)'shep_multinest_aspic_loglike: trapped into hard prior regions'
       write(*,*)'model is: ',label
       write(*,*)'mnparams= ',mnpars(1:nestdim)
       stop
    endif
    
    if (debug) then
       write(*,*)
       write(*,*)'shep_multinest_aspic_loglike: '
       write(*,*)'lnA= log(eps1)= eps_i=       ',sheppars
       write(*,*)'ln(like) =                   ',lnew
       write(*,*)
    end if

  end subroutine shep_multinest_aspic_loglike


#ifndef NOFANN
  subroutine fnn_multinest_aspic_loglike(cube,nestdim,nestpars,lnew,context)
    use fnnprec, only : fp
    use fnnlike, only : cubize_fnnparams, uncubize_fnnparams, check_fnn
    use fnnlike, only : fnnlike_eval, cutmin_fnnparams, cutmax_fnnparams
    use wraspic, only : get_hubbleflow, get_derived
    use wraspic, only : test_aspic_hardprior, test_reheating_hardprior
    use nestparams, only : nestLogZero
    implicit none   
    integer(imn) :: nestdim, nestpars
    real(fmn), dimension(nestpars) :: cube
    real(fmn) :: lnew
    integer(imn) :: context,i
    real(fmn), dimension(nestdim) :: mnpars
    real(fp), dimension(fitNdim) :: fnncube, fnnpars, fnncuts, fnnclamp

    integer, save :: counthard = 0
    logical :: disfavour

    if (.not.check_fnn()) stop 'fnn_multinest_loglike: not initialized!'

!get the physical parameters we are sampling on
    mnpars(1:nestdim) = uncubize_samplparams(nestdim,cube(1:nestdim))

!check for any hard prior in aspic model parameters
    if (test_aspic_hardprior(mnpars,disfavour)) then

       if (disfavour) then
          lnew = nestLogZero / 10._fmn
       else
          lnew = nestLogZero
       endif
      
       counthard = counthard + 1
       
    else
       
!use aspic to get the slowroll parameters
       fnnpars = get_hubbleflow(fitNdim,mnpars)

! reheating hardprior, those points are ignored
       if (test_reheating_hardprior(mnpars)) then

          lnew = nestLogZero
          counthard = counthard + 1
          
       else

          counthard=0
          
!if eps1<eps1min, the likelihood is flat
          fnncuts = cutmin_fnnparams(fitNdim,eps1pos,fnnpars)

!In case the posterior of eps3 is not going to zero on the boundaries          
          if (cuteps3) then
             fnncuts = cutmin_fnnparams(fitNdim,eps3pos,fnncuts)
             fnncuts = cutmax_fnnparams(fitNdim,eps3pos,fnncuts)
          endif          
          
!go into cubic space for the fnn likelihood
          fnncube = cubize_fnnparams(fitNdim,fnncuts)
       
          if (any(fnncube.gt.1._fp).or.any(fnncube.lt.0._fp)) then
!if outside fnnfits box, the real likelihood is so small that we
!cannot calculate it numerically, but we can define a junk one smaller
             
             fnnclamp = max(0._fp,min(fnncube,1._fp))
             
             lnew = fnnlike_eval(fnnclamp) - sum((fnncube-fnnclamp)**2)

             
          else

             lnew = fnnlike_eval(fnncube)

          endif

       endif

    endif

!dump the physical params into cube for file output
    cube(1:nestdim) = mnpars(1:nestdim)

!+ derived parameters we want to dump too
    do i=1,nestpars-nestdim
       cube(nestdim+i) = get_derived(i)
    enddo

!sanity check
    if (isnan(lnew)) stop 'fnn_multinest_aspic_loglike: like is NaN!'

    if (counthard.gt.countmax) then
       write(*,*)'fnn_multinest_aspic_loglike: trapped into hard prior regions'
       write(*,*)'model is: ',label
       write(*,*)'mnparams= ',mnpars(1:nestdim)
       stop
    endif
    
    if (debug) then
       write(*,*)
       write(*,*)'fnn_multinest_aspic_loglike: '
       write(*,*)'lnA= log(eps1)= eps_i=       ',fnnpars
       write(*,*)'ln(like) =                   ',lnew
       write(*,*)
    end if

  end subroutine fnn_multinest_aspic_loglike
#endif



  subroutine null_multinest_aspic_loglike(cube,nestdim,nestpars,lnew,context)    
    use nestparams, only : nestLogZero
    use wraspic, only : get_derived, test_aspic_hardprior, test_reheating_hardprior
    use wraspic, only : get_hubbleflow
    implicit none   
    integer(imn) :: nestdim, nestpars
    real(fmn), dimension(nestpars) :: cube
    real(fmn) :: lnew
    integer(imn) :: context,i
    real(fmn), dimension(nestdim) :: mnpars    
    real(fmn), dimension(fitNdim) :: nullpars

    logical :: disfavour
    
!get the physical parameters we are sampling on
    mnpars(1:nestdim) = uncubize_samplparams(nestdim,cube(1:nestdim))

!check for any hard prior in aspic model parameters
    if (test_aspic_hardprior(mnpars,disfavour)) then

       if (disfavour) then
          lnew = nestLogZero / 10._fmn
       else
          lnew = nestLogZero
       endif

    else

!to get derived parameters
       nullpars = get_hubbleflow(fitNdim,mnpars)

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

    procedure(pc_loglike), pointer :: ptrchord_aspic_loglike => null()
    

    type(program_settings) :: runset

    call chord_settings(runset)

    select case (fastLikeName)

    case ('rbf')
       ptrchord_aspic_loglike => rbf_polychord_aspic_loglike

    case ('shep')
       ptrchord_aspic_loglike => shep_polychord_aspic_loglike

#ifndef NOFANN       
    case ('fnn')
       ptrchord_aspic_loglike => fnn_polychord_aspic_loglike
#endif
       
    case ('null')
       ptrchord_aspic_loglike => null_polychord_aspic_loglike

    case default

       stop 'chord_sample_aspic: fast like not found!'

    end select

    call run_polychord(ptrchord_aspic_loglike,polychord_prior,runset)

    if (associated(ptrchord_aspic_loglike)) ptrchord_aspic_loglike => null()

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
    use wraspic, only : get_hubbleflow, get_derived
    use wraspic, only : test_aspic_hardprior, test_reheating_hardprior
    implicit none   
    real(fmn) :: lnew
    real(fmn), dimension(:), intent(in) :: theta
    real(fmn), dimension(:), intent(out) :: phi

    real(fp), dimension(fitNdim) :: rbfcube, rbfcuts, rbfpars, rbfclamp
    integer(imn) :: ndim, i

    logical :: disfavour
    
    ndim = size(theta,1)

    if (.not.check_rbf()) stop 'rbf_polychord_aspic_loglike: not initialized!'
    if (ndim.ne.chordNdim) stop 'rbf_polychord_aspic_loglike: dim mismatch!'

!check for any hard prior in aspic model parameters
    if (test_aspic_hardprior(theta,disfavour)) then

       if (disfavour) then
          lnew = chordLogZero / 10._fmn
       else
          lnew = chordLogZero
       endif       

    else
       
!use aspic to get the slowroll parameters
       rbfpars = get_hubbleflow(fitNdim,theta)

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

             rbfclamp = max(0._fp,min(rbfcube,1._fp))
             
             lnew = rbflike_eval(rbfclamp) - sum((rbfcube-rbfclamp)**2)
             
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
    use wraspic, only : get_hubbleflow, get_derived
    use wraspic, only : test_aspic_hardprior, test_reheating_hardprior
    implicit none   
    real(fmn) :: lnew
    real(fmn), dimension(:), intent(in) :: theta
    real(fmn), dimension(:), intent(out) :: phi

    integer(imn) :: ndim, i
    real(fp), dimension(fitNdim) :: shepcube, sheppars, shepcuts, shepclamp

    logical :: disfavour
    
    ndim = size(theta,1)
    if (.not.check_shep()) stop 'shep_polychord_aspic_loglike: not initialized!'
    if (ndim.ne.chordNdim) stop 'shep_polychord_aspic_loglike: dim mismatch!'

!check for any hard prior in aspic model parameters
    if (test_aspic_hardprior(theta,disfavour)) then

       if (disfavour) then
          lnew = chordLogZero / 10._fmn
       else
          lnew = chordLogZero
       endif              

    else
       
!use aspic to get the slowroll parameters
       sheppars = get_hubbleflow(fitNdim,theta)

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

             shepclamp = max(0._fp,min(shepcube,1._fp))             
             lnew = sheplike_eval(shepclamp) - sum((shepcube-shepclamp)**2)
             
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



#ifndef NOFANN
  function fnn_polychord_aspic_loglike(theta,phi) result(lnew)
    use chordparams, only : chordLogZero, chordNderived, chordNdim
    use fnnprec, only : fp
    use fnnlike, only : cubize_fnnparams, check_fnn
    use fnnlike, only : fnnlike_eval, cutmin_fnnparams
    use wraspic, only : get_hubbleflow, get_derived
    use wraspic, only : test_aspic_hardprior, test_reheating_hardprior
    implicit none   
    real(fmn) :: lnew
    real(fmn), dimension(:), intent(in) :: theta
    real(fmn), dimension(:), intent(out) :: phi

    integer(imn) :: ndim, i
    real(fp), dimension(fitNdim) :: fnncube, fnnpars, fnncuts, fnnclamp

    logical :: disfavour
    
    ndim = size(theta,1)
    if (.not.check_fnn()) stop 'fnn_polychord_aspic_loglike: not initialized!'
    if (ndim.ne.chordNdim) stop 'fnn_polychord_aspic_loglike: dim mismatch!'

!check for any hard prior in aspic model parameters
    if (test_aspic_hardprior(theta,disfavour)) then

       if (disfavour) then
          lnew = chordLogZero / 10._fmn
       else
          lnew = chordLogZero
       endif              
       
    else
       
!use aspic to get the slowroll parameters
       fnnpars = get_hubbleflow(fitNdim,theta)

! reheating hardprior, those points are ignored
       if (test_reheating_hardprior(theta)) then

          lnew = chordLogZero

       else

!if eps1<eps1min, the likelihood is flat
          fnncuts = cutmin_fnnparams(fitNdim,eps1pos,fnnpars)

!go into cubic space for the fnn likelihood
          fnncube = cubize_fnnparams(fitNdim,fnncuts)
       
          if (any(fnncube.gt.1._fp).or.any(fnncube.lt.0._fp)) then
!if outside fnnfits box, the real likelihood is so small that we
!cannot calculate it numerically, but we can define a junk one smaller
!than fitLogZero
             fnnclamp = max(0._fp,min(fnncube,1._fp))
             lnew = fnnlike_eval(fnnclamp) - sum((fnncube-fnnclamp)**2)

          else

             lnew = fnnlike_eval(fnncube)

          endif

       endif

    endif


!+ derived parameters we want to dump too
    do i=1,chordNderived
       phi(i) = get_derived(i)
    enddo


    if (debug) then
       write(*,*)
       write(*,*)'fnn_polychord_aspic_loglike: '
       write(*,*)'lnA= log(eps1)= eps_i=       ',fnnpars
       write(*,*)'ln(like) =                   ',lnew
       write(*,*)
    end if

  end function fnn_polychord_aspic_loglike
#endif

  

  function null_polychord_aspic_loglike(theta,phi)    
    use chordparams, only : chordLogZero, chordNderived, chordNdim
    use wraspic, only : get_derived, test_aspic_hardprior, test_reheating_hardprior
    use wraspic, only : get_hubbleflow
    implicit none   
    real(fmn) :: null_polychord_aspic_loglike
    real(fmn), dimension(:), intent(in) :: theta
    real(fmn), dimension(:), intent(out) :: phi

    real(fmn) :: lnew
    integer(imn) :: ndim, i
    real(fmn), dimension(size(theta,1)) :: nullpars

    logical :: disfavour
    
    ndim = size(theta,1)
    if (ndim.ne.chordNdim) stop 'null_polychord_aspic_loglike: dim mismatch!'

!check for any hard prior in aspic model parameters
    if (test_aspic_hardprior(theta,disfavour)) then

       if (disfavour) then
          lnew = chordLogZero / 10._fmn
       else
          lnew = chordLogZero
       endif                            

    else

!to get derived parameters
       nullpars = get_hubbleflow(fitNdim,theta)

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
