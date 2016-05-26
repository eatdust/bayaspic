!sampling the slow-roll parameter space only
module samplsr
  use sampl, only : imn, fmn
  use sampl, only : samplNdim, rootName, rootPrefix
  use sampl, only : init_samplparams, free_samplparams
  use sampl, only : samplPmin, samplPmax
  implicit none

  private
 
  character(len=*), parameter :: fileweights = 'rbfdata/weights.dat'
  character(len=*), parameter :: filecentres = 'rbfdata/centres.dat'
  character(len=*), parameter :: filerbfbounds = 'rbfdata/bounds.dat'

  character(len=*), parameter :: fileshep = 'shepdata/shepdata.dat'
  character(len=*), parameter :: filepost = 'shepdata/postcubed.dat'
  character(len=*), parameter :: fileshepbounds = 'shepdata/bounds.dat'

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
  character(len=*), parameter :: fastLikeName = 'shep'



  public nest_init_slowroll,  nest_free_slowroll
  public nest_sample_slowroll

  

contains

 
   subroutine init_slowroll()
    use rbflike, only : initialize_rbf_like
    use rbflike, only : get_rbf_ndim, get_rbf_fmin
    use rbflike, only : get_rbf_xpmin, get_rbf_xpmax
    use sheplike, only : initialize_shep_like
    use sheplike, only : get_shep_ndim, get_shep_fmin
    use sheplike, only : get_shep_xpmin, get_shep_xpmax
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
       write(*,*)'fast like is :         ',fastLikeName
       write(*,*)'lnZeroMin    =         ',fitLogZero
       write(*,*)'fitNdim      =         ',fitNdim
       write(*,*)
    end if

  end subroutine init_slowroll




  subroutine dump_slowroll_priors(extname)
    implicit none
    character(len=*), intent(in) :: extname
    integer, parameter :: nunit = 414, nname = 415, nrange=416
    
    if ((.not.allocated(samplPmin)).or.(.not.allocated(samplPmax))) then
       stop 'dump_slowroll_priors: prior not allocated!'
    endif
   

    open(unit=nunit,file=trim(rootName)//'.ini', status='unknown')
    open(unit=nname,file=trim(rootName)//'.paramnames', status='unknown')
    open(unit=nrange,file=trim(rootName)//'.ranges', status='unknown')

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



  function rbf_slowroll_loglike(cube)
    use rbfprec, only : fp
    use rbflike, only : check_rbf
    use rbflike, only : rbflike_eval
    implicit none   
    real(fmn) :: rbf_slowroll_loglike
    real(fmn), dimension(:) :: cube
    real(fp), dimension(size(cube,1)) :: rbfcube

    integer :: ndim
    ndim = size(cube,1)

    if (.not.check_rbf()) stop 'rbf_loglike: not initialized!'
    
    rbfcube(1:ndim) = cube(1:ndim)

    rbf_slowroll_loglike = rbflike_eval(rbfcube)

  end function rbf_slowroll_loglike



  function shep_slowroll_loglike(cube)
    use shepprec, only : fp
    use sheplike, only : check_shep
    use sheplike, only : sheplike_eval
    implicit none   
    real(fmn) :: shep_slowroll_loglike
    real(fmn), dimension(:) :: cube
    real(fp), dimension(size(cube,1)) :: shepcube
    integer :: ndim

    if (.not.check_shep()) stop 'shep_loglike: not initialized!'

    ndim = size(cube,1)

    shepcube(1:ndim) = cube(1:ndim)

    shep_slowroll_loglike = sheplike_eval(shepcube)

  end function shep_slowroll_loglike

  


  subroutine nest_init_slowroll()    
    use nestparams, only : nestNdim, nestNpars, nestCdim
    use nestparams, only : nestPWrap, nestRootName, nest_print
    use sampl, only : sampl_print
    implicit none

    call init_slowroll()

    nestNdim = samplNdim
    nestNpars = samplNdim
    nestCdim = samplNdim

    nestRootName = rootName

    allocate(nestPwrap(nestNdim))
    nestPwrap = 0

    call nest_print()

    call sampl_print()

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
    use rbflike, only : uncubize_rbfparams
    implicit none   
    integer(imn) :: nestdim, nestpars
    real(fmn), dimension(nestpars) :: cube
    real(fmn) :: lnew
    integer(imn) :: context
    
    lnew = rbf_slowroll_loglike(cube(1:nestdim))
    
    cube(1:nestdim) = real(uncubize_rbfparams(nestdim,real(cube(1:nestdim),fp)),fmn)
    
  end subroutine rbf_multinest_slowroll_loglike



  subroutine shep_multinest_slowroll_loglike(cube,nestdim,nestpars,lnew,context)
    use shepprec, only : fp
    use sheplike, only : uncubize_shepparams
    implicit none   
    integer(imn) :: nestdim, nestpars
    real(fmn), dimension(nestpars) :: cube
    real(fmn) :: lnew
    integer(imn) :: context
    
    lnew = shep_slowroll_loglike(cube(1:nestdim))
    
    cube(1:nestdim) = real(uncubize_shepparams(nestdim,real(cube(1:nestdim),fp)),fmn)
    
  end subroutine shep_multinest_slowroll_loglike

  

  subroutine nest_free_slowroll()
    use nestparams, only : nestPwrap    
    implicit none
   
    if (allocated(nestPwrap)) deallocate(nestPwrap)

    call free_samplparams()

  end subroutine nest_free_slowroll

  

  
  subroutine chord_init_slowroll()    
    use chordparams, only : chordNdim, chordNpars, chordCdim
    use chordparams, only : chordPWrap, chordRootName, chord_print
    use sampl, only : sampl_print
    implicit none

    call init_slowroll()

    chordNdim = samplNdim
    chordNpars = samplNdim
    chordCdim = samplNdim

    chordRootName = rootName

    allocate(chordPwrap(chordNdim))
    chordPwrap = 0

    call chord_print()

    call sampl_print()

  end subroutine chord_init_slowroll


  subroutine chord_sample_slowroll()
    use chorded, only : chordRun
    use chordparams
    implicit none

    integer(imn) :: nclusters			
    integer(imn) :: context
    integer(imn) :: maxNode 			

    select case (fastLikeName)

    case ('rbf')

       call chordRun(chordINS,chordMmodal,chordCteEff,chordNlive,chordZtol,chordSampEff,chordNdim,chordNpars, &
            chordCdim,chordMaxModes,chordUpdInt,chordNullZ,chordRootName,chordSeed,chordPwrap, &
            chordFeedBack,chordResume,chordOutfile,chordInitMPI,chordLogZero,chordMaxIter &
            ,rbf_multichord_slowroll_loglike,chord_dumper,context)

    case ('shep')

       call chordRun(chordINS,chordMmodal,chordCteEff,chordNlive,chordZtol,chordSampEff,chordNdim,chordNpars, &
            chordCdim,chordMaxModes,chordUpdInt,chordNullZ,chordRootName,chordSeed,chordPwrap, &
            chordFeedBack,chordResume,chordOutfile,chordInitMPI,chordLogZero,chordMaxIter &
            ,shep_multichord_slowroll_loglike,chord_dumper,context)

    case default

       stop 'chord_sample_slowroll: fast like not found!'

    end select


  end subroutine chord_sample_slowroll

  

  subroutine rbf_multichord_slowroll_loglike(cube,chorddim,chordpars,lnew,context)
    use rbfprec, only : fp
    use rbflike, only : uncubize_rbfparams
    implicit none   
    integer(imn) :: chorddim, chordpars
    real(fmn), dimension(chordpars) :: cube
    real(fmn) :: lnew
    integer(imn) :: context
    
    lnew = rbf_slowroll_loglike(cube(1:chorddim))
    
    cube(1:chorddim) = real(uncubize_rbfparams(chorddim,real(cube(1:chorddim),fp)),fmn)
    
  end subroutine rbf_multichord_slowroll_loglike



  subroutine shep_multichord_slowroll_loglike(cube,chorddim,chordpars,lnew,context)
    use shepprec, only : fp
    use sheplike, only : uncubize_shepparams
    implicit none   
    integer(imn) :: chorddim, chordpars
    real(fmn), dimension(chordpars) :: cube
    real(fmn) :: lnew
    integer(imn) :: context
    
    lnew = shep_slowroll_loglike(cube(1:chorddim))
    
    cube(1:chorddim) = real(uncubize_shepparams(chorddim,real(cube(1:chorddim),fp)),fmn)
    
  end subroutine shep_multichord_slowroll_loglike

  

  subroutine chord_free_slowroll()
    use chordparams, only : chordPwrap    
    implicit none
   
    if (allocated(chordPwrap)) deallocate(chordPwrap)

    call free_samplparams()

  end subroutine chord_free_slowroll



end module samplsr
