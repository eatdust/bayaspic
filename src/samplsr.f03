!sampling the slow-roll parameter space only
module samplsr
  use sampl, only : imn, fmn, lenmn
  use sampl, only : samplNdim, rootName, rootPrefix, rootDir
  use sampl, only : init_samplparams, free_samplparams
  use sampl, only : samplPmin, samplPmax
  use sampl, only : mn_loglike, pc_loglike, pc_prior
  implicit none

  private
 
  character(len=*), parameter :: fileweights = 'rbfdata/weights.dat'
  character(len=*), parameter :: filecentres = 'rbfdata/centres.dat'
  character(len=*), parameter :: filerbfbounds = 'rbfdata/bounds.dat'

  character(len=*), parameter :: fileshep = 'shepdata/shepdata.dat'
  character(len=*), parameter :: filepost = 'shepdata/postcubed.dat'
  character(len=*), parameter :: fileshepbounds = 'shepdata/bounds.dat'

  character(len=*), parameter :: filefnn = 'fnndata/fnndata.dat'
  character(len=*), parameter :: filefnnbounds = 'fnndata/bounds.dat'

  logical, parameter :: display = .true.
  
  integer, parameter :: lensr = 3
  character(len=lensr), save :: srname = 'ooo'

  integer(imn), save :: fitNdim = 0
  real(fmn), save :: fitLogZero = 0._fmn



!the name of the fastlike method we call

!  character(len=*), parameter :: fastLikeName = 'null'

! radial basis functions
!  character(len=*), parameter :: fastLikeName = 'rbf'

! inverse shepard method
!  character(len=*), parameter :: fastLikeName = 'shep'

!feedforward neural network
  character(len=*), parameter :: fastLikeName = 'fnn'


  public nest_init_slowroll,  nest_free_slowroll
  public nest_sample_slowroll

  public chord_init_slowroll, chord_free_slowroll
  public chord_sample_slowroll

contains

 
   subroutine init_slowroll()
    use rbflike, only : initialize_rbf_like
    use rbflike, only : get_rbf_ndim, get_rbf_fmin
    use rbflike, only : get_rbf_xpmin, get_rbf_xpmax
    use sheplike, only : initialize_shep_like
    use sheplike, only : get_shep_ndim, get_shep_fmin
    use sheplike, only : get_shep_xpmin, get_shep_xpmax
    use fnnlike, only : initialize_fnn_like
    use fnnlike, only : get_fnn_ndim, get_fnn_fmin
    use fnnlike, only : get_fnn_xpmin, get_fnn_xpmax
    implicit none

    integer :: i

    select case(fastLikeName)
    
    case ('rbf')

       call initialize_rbf_like(fileweights, filecentres, filerbfbounds)
       fitNdim = get_rbf_ndim()
       fitLogZero = get_rbf_fmin()

       call init_samplparams(fitNdim)
       
       do i=1,fitNdim
          samplPmin(i) = get_rbf_xpmin(i)
          samplPmax(i) = get_rbf_xpmax(i)
       enddo

    case ('shep')

       call initialize_shep_like(fileshep, filepost, fileshepbounds)
       fitNdim = get_shep_ndim()
       fitLogZero = get_shep_fmin()
      
       call init_samplparams(fitNdim)

       do i=1,fitNdim
          samplPmin(i) = get_shep_xpmin(i)
          samplPmax(i) = get_shep_xpmax(i)
       enddo

    case ('fnn')

       call initialize_fnn_like(filefnn,filefnnbounds)
       fitNdim = get_fnn_ndim()
       fitLogZero = get_fnn_fmin()
      
       print *,'fitNdim',fitndim

       call init_samplparams(fitNdim)

       do i=1,fitNdim
          samplPmin(i) = get_fnn_xpmin(i)
          samplPmax(i) = get_fnn_xpmax(i)
       enddo

    case default

       stop 'sampl_init_slowroll: fast like not found!'

    end select
         

    select case (fitNdim)
    case (3)
       srname = 'sr2'
    case (4)
       srname = 'sr3'
    case default
       stop 'sampl_init_slowroll: sr model not known'
    end select

    rootName = trim(rootPrefix)//srname

    call dump_slowroll_priors(srname)

    if (display) then
       write(*,*)
       write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(*,*)'fast like is :         ',fastLikeName
       write(*,*)'lnZeroMin    =         ',fitLogZero
       write(*,*)'fitNdim      =         ',fitNdim
       write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(*,*)
    end if

  end subroutine init_slowroll




  subroutine dump_slowroll_priors(extname)
    implicit none
    character(len=*), intent(in) :: extname
    character(len=lenmn) :: fileName
    integer, parameter :: nunit = 414, nname = 415, nrange=416
    
    if ((.not.allocated(samplPmin)).or.(.not.allocated(samplPmax))) then
       stop 'dump_slowroll_priors: prior not allocated!'
    endif
   
    fileName=trim(rootDir)//trim(rootPrefix)//extname

    open(unit=nunit,file=trim(fileName)//'.ini', status='unknown')
    open(unit=nname,file=trim(fileName)//'.paramnames', status='unknown')
    open(unit=nrange,file=trim(fileName)//'.ranges', status='unknown')

    write(nunit,*)'limits[lnA]=',samplPmin(1),samplPmax(1)    
    write(nname,*)'lnA         \ln(10^{10} P_*)'
    write(nrange,*)'lnA                    ',samplPmin(1),samplPmax(1)

    select case (extname)

    case ('sr2')

       write(nunit,*)'limits[sr1]=',samplPmin(2),samplPmax(2)
       write(nname,*)'sr1         \log(\epsilon_1)'
       write(nrange,*)'sr1                    ',samplPmin(2),samplPmax(2)
       write(nunit,*)'limits[sr2]=',samplPmin(3),samplPmax(3)
       write(nname,*)'sr2         \epsilon_2'
       write(nrange,*)'sr2                    ',samplPmin(3),samplPmax(3)

    case ('sr3')
       write(nunit,*)'limits[sr1]=',samplPmin(2),samplPmax(2)
       write(nname,*)'sr1         \log(\epsilon_1)'
       write(nrange,*)'sr1                    ',samplPmin(2),samplPmax(2)
       write(nunit,*)'limits[sr2]=',samplPmin(3),samplPmax(3)
       write(nname,*)'sr2         \epsilon_2'
       write(nrange,*)'sr2                    ',samplPmin(3),samplPmax(3)
       write(nunit,*)'limits[sr3]=',samplPmin(4),samplPmax(4)
       write(nname,*)'sr3         \epsilon_3'
       write(nrange,*)'sr3                     ',samplPmin(4),samplPmax(4)

    case default

       stop 'sampl_dump_slowroll_priors: internal error!'

    end select

    close(nunit)
    close(nname)
    close(nrange)

  end subroutine dump_slowroll_priors



  

  subroutine nest_init_slowroll()    
    use nestparams, only : nestNdim, nestNpars, nestCdim
    use nestparams, only : nestPWrap, nestRootName, nest_print
    use sampl, only : sampl_print
    implicit none

    call init_slowroll()

    nestNdim = samplNdim
    nestNpars = samplNdim
    nestCdim = samplNdim

    nestRootName = trim(rootDir)//trim(rootName)

    allocate(nestPwrap(nestNdim))
    nestPwrap = 0

    call nest_print()

    call sampl_print()

  end subroutine nest_init_slowroll


  subroutine nest_sample_slowroll()
    use nested, only : nestRun
    use nestparams
    implicit none

    procedure(mn_loglike), pointer :: ptrnest_slowroll_loglike => null()

    integer(imn) :: nclusters			
    integer(imn) :: context
    integer(imn) :: maxNode 			

    select case (fastLikeName)

    case ('rbf')
       ptrnest_slowroll_loglike => rbf_multinest_slowroll_loglike
       
    case ('shep')
       ptrnest_slowroll_loglike => shep_multinest_slowroll_loglike

    case ('fnn')
       ptrnest_slowroll_loglike => fnn_multinest_slowroll_loglike
          
    case default

       stop 'nest_sample_slowroll: fast like not found!'

    end select

    call nestRun(nestINS,nestMmodal,nestCteEff,nestNlive,nestZtol,nestSampEff,nestNdim,nestNpars, &
         nestCdim,nestMaxModes,nestUpdInt,nestNullZ,nestRootName,nestSeed,nestPwrap, &
         nestFeedBack,nestResume,nestOutfile,nestInitMPI,nestLogZero,nestMaxIter &
         ,ptrnest_slowroll_loglike,nest_dumper,context)

    if (associated(ptrnest_slowroll_loglike)) ptrnest_slowroll_loglike => null()

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



  subroutine fnn_multinest_slowroll_loglike(cube,nestdim,nestpars,lnew,context)
    use fnnprec, only : fp
    use fnnlike, only : uncubize_fnnparams, check_fnn
    use fnnlike, only : fnnlike_eval
    implicit none   
    integer(imn) :: nestdim, nestpars
    real(fmn), dimension(nestpars) :: cube
    real(fmn) :: lnew
    integer(imn) :: context
    real(fp), dimension(nestdim) :: fnncube

    if (.not.check_fnn()) stop 'fnn_multinest_loglike: not initialized!'
    
    fnncube(1:nestdim) = cube(1:nestdim)

    lnew = fnnlike_eval(fnncube)
    
    cube(1:nestdim) = real(uncubize_fnnparams(nestdim,fnncube),fmn)
    
  end subroutine fnn_multinest_slowroll_loglike




  subroutine nest_free_slowroll()
    use nestparams, only : nestPwrap    
    implicit none
   
    if (allocated(nestPwrap)) deallocate(nestPwrap)

    call free_samplparams()

  end subroutine nest_free_slowroll

  

  
  subroutine chord_init_slowroll()    
    use chordparams, only : chordNdim, chordNpars, chordNderived
    use chordparams, only : chordName, chordDir, chord_print
    use sampl, only : sampl_print
    implicit none

    call init_slowroll()

    chordNdim = samplNdim
    chordNderived = 0
    chordNpars = samplNdim

    chordDir=rootDir
    chordName = trim(rootName)
   
    call chord_print()

    call sampl_print()

  end subroutine chord_init_slowroll



  subroutine chord_sample_slowroll()
    use interfaces_module, only : run_polychord
    use settings_module, only : program_settings
    use chordparams
    implicit none

    procedure(pc_loglike), pointer :: ptrchord_slowroll_loglike => null()
    procedure(pc_prior), pointer :: ptrchord_prior => null()

    type(program_settings) :: runset

    call chord_settings(runset)

    select case (fastLikeName)

    case ('rbf')
       ptrchord_slowroll_loglike => rbf_polychord_slowroll_loglike
       ptrchord_prior => rbf_polychord_prior

    case ('shep')
       ptrchord_slowroll_loglike => shep_polychord_slowroll_loglike
       ptrchord_prior => shep_polychord_prior

    case ('fnn')
       ptrchord_slowroll_loglike => fnn_polychord_slowroll_loglike
       ptrchord_prior => fnn_polychord_prior

    case default
       stop 'chord_sample_slowroll: fast like not found!'

    end select

    call run_polychord(ptrchord_slowroll_loglike,ptrchord_prior,runset)

    if (associated(ptrchord_slowroll_loglike)) ptrchord_slowroll_loglike => null()
    if (associated(ptrchord_prior)) ptrchord_prior => null()

  end subroutine chord_sample_slowroll



  function rbf_polychord_prior(cube)
    use rbfprec, only : fp
    use rbflike, only :  uncubize_rbfparams
    implicit none
    real(fmn), dimension(:), intent(in) :: cube
    real(fmn), dimension(size(cube,1)) :: rbf_polychord_prior
    integer(imn) :: ndim
    ndim = size(cube,1)

    rbf_polychord_prior = real(uncubize_rbfparams(ndim,real(cube(1:ndim),fp)),fmn)

  end function rbf_polychord_prior



  function rbf_polychord_slowroll_loglike(theta,phi)
    use rbfprec, only : fp
    use rbflike, only : cubize_rbfparams, check_rbf
    use rbflike, only : rbflike_eval
    implicit none   
    real(fmn) :: rbf_polychord_slowroll_loglike
    real(fmn), dimension(:), intent(in) :: theta
    real(fmn), dimension(:), intent(out) :: phi

    integer(imn) :: ndim
    real(fp), dimension(size(theta,1)) :: rbfcube

    ndim = size(theta,1)

    if (.not.check_rbf()) stop 'rbf_polychord_slowroll_loglike: not initialized!'
    
    rbfcube(1:ndim) = cubize_rbfparams(ndim,real(theta(1:ndim),fp))

    rbf_polychord_slowroll_loglike = rbflike_eval(rbfcube)
        
  end function rbf_polychord_slowroll_loglike




  function shep_polychord_prior(cube)
    use shepprec, only : fp
    use sheplike, only :  uncubize_shepparams
    implicit none
    real(fmn), dimension(:), intent(in) :: cube
    real(fmn), dimension(size(cube,1)) :: shep_polychord_prior
    integer(imn) :: ndim
    ndim = size(cube,1)

    shep_polychord_prior = real(uncubize_shepparams(ndim,real(cube(1:ndim),fp)),fmn)

  end function shep_polychord_prior




  function shep_polychord_slowroll_loglike(theta,phi)
    use shepprec, only : fp
    use sheplike, only : cubize_shepparams, check_shep
    use sheplike, only : sheplike_eval
    implicit none   
    real(fmn) :: shep_polychord_slowroll_loglike
    real(fmn), dimension(:), intent(in) :: theta
    real(fmn), dimension(:), intent(out) :: phi

    integer(imn) :: ndim
    real(fp), dimension(size(theta,1)) :: shepcube

    ndim = size(theta,1)

    if (.not.check_shep()) stop 'shep_polychord_slowroll_loglike: not initialized!'
    
    shepcube(1:ndim) = cubize_shepparams(ndim,real(theta(1:ndim),fp))

    shep_polychord_slowroll_loglike = sheplike_eval(shepcube)

  end function shep_polychord_slowroll_loglike



  function fnn_polychord_prior(cube)
    use fnnprec, only : fp
    use fnnlike, only :  uncubize_fnnparams
    implicit none
    real(fmn), dimension(:), intent(in) :: cube
    real(fmn), dimension(size(cube,1)) :: fnn_polychord_prior
    integer(imn) :: ndim
    ndim = size(cube,1)

    fnn_polychord_prior = real(uncubize_fnnparams(ndim,real(cube(1:ndim),fp)),fmn)

  end function fnn_polychord_prior




  function fnn_polychord_slowroll_loglike(theta,phi)
    use fnnprec, only : fp
    use fnnlike, only : cubize_fnnparams, check_fnn
    use fnnlike, only : fnnlike_eval
    implicit none   
    real(fmn) :: fnn_polychord_slowroll_loglike
    real(fmn), dimension(:), intent(in) :: theta
    real(fmn), dimension(:), intent(out) :: phi

    integer(imn) :: ndim
    real(fp), dimension(size(theta,1)) :: fnncube

    ndim = size(theta,1)

    if (.not.check_fnn()) stop 'fnn_polychord_slowroll_loglike: not initialized!'

    fnncube(1:ndim) = cubize_fnnparams(ndim,real(theta(1:ndim),fp))

    fnn_polychord_slowroll_loglike = fnnlike_eval(fnncube)

  end function fnn_polychord_slowroll_loglike




  subroutine chord_free_slowroll()
    implicit none
   
    call free_samplparams()

  end subroutine chord_free_slowroll



end module samplsr
