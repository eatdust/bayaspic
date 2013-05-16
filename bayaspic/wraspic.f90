module wraspic
  use nestprec, only : fmn, imn
  use aspicvars, only : kp, infaspic, lname
  implicit none

  private


  integer, parameter :: nextra = 2
  integer, parameter :: ireh= 2
  integer, parameter :: ilnA = 1
  logical, parameter :: useRrad = .false.

  type(infaspic), save :: AspicModel

  integer, parameter :: naspmax = 3
  integer, parameter :: nepsmax = 3

  public set_model, check_model, get_allprior
  public get_slowroll, get_ntot, get_derived
  public test_aspic_hardprior, test_reheating_hardprior

  logical, parameter :: display = .true.

contains



  subroutine set_model(mname)
    use aspicmodels, only : initialize_aspic_ptrs
    use aspicmodels, only : get_aspic_numparams

    implicit none
    character(len=*), intent(in) :: mname

    call initialize_aspic_ptrs(trim(adjustl(mname)))

    AspicModel%name = trim(adjustl(mname))
    AspicModel%nasp = get_aspic_numparams()

    allocate(AspicModel%params(AspicModel%nasp))
    allocate(AspicModel%cmaps(AspicModel%nasp))

  end subroutine set_model



  function check_model()
    implicit none
    logical :: check_model

    check_model = associated(AspicModel%params)

  end function check_model
  

  function get_ntot()
    implicit none
    integer :: get_ntot

    if (.not.check_model()) stop 'get_ntot: aspic not set!'

    get_ntot = AspicModel%nasp + nextra
    
  end function get_ntot



  subroutine get_prior_lnA(lnAmin,lnAmax)
    implicit none
    real(fmn) , intent(out) :: lnAmin,lnAmax

    lnAmin = 2.7
    lnAmax = 4.
!    lnAmin = 3.000
!    lnAmax = 3.175


  end subroutine get_prior_lnA



  subroutine get_prior_lnRrad(lnRradmin,lnRradmax)
    use cosmopar, only : lnRhoNuc
    implicit none
    real(fmn) , intent(out) :: lnRradmin,lnRradmax

    lnRradMin = 0.25 * lnRhoNuc
    lnRradMax = -1._kp/12._kp * lnRhoNuc

  end subroutine get_prior_lnRrad



  subroutine get_prior_lnRreh(lnRmin,lnRmax)
    use cosmopar, only : lnRhoNuc
    implicit none
    real(fmn) , intent(out) :: lnRmin,lnRmax

!    lnRMin = 0.25_kp * lnRhoNuc
!which is ~
    lnRMin = -46

!    lnRMax = -1._kp/12._kp * lnRhoNuc
!which is ~
    lnRMax = 15

  end subroutine get_prior_lnRreh


  
  function get_derived(i)
    implicit none
    integer, intent(in) :: i
    real(fmn) :: get_derived

    select case(i)

    case (1)
       get_derived = AspicModel%lnM
   
    case(2)
       if (useRrad) then
          get_derived = AspicModel%lnRreh
       else
          get_derived = AspicModel%lnRrad
       endif

    case (3)       
       get_derived = AspicModel%lnRhoEnd

    case (4)
       get_derived = AspicModel%bfold

    case (5)
       get_derived = AspicModel%logeps

    case (6)
       get_derived = AspicModel%eps2

    case (7)
       get_derived = AspicModel%eps3

    case (8)
       get_derived = AspicModel%ns

    case (9)
       get_derived = AspicModel%logr

    case (10)
       get_derived = AspicModel%alpha

    

    case default
       stop 'get_derived: incorrect parameters number!'

    end select


  end function get_derived



  subroutine get_allprior(pmin,pmax)
    use aspicpriors, only : get_aspic_priors
    implicit none
    real(fmn), dimension(:), intent(out) :: pmin, pmax

    integer :: nsize, nasp, ntot
    real(kp), dimension(naspmax) :: aspmin,aspmax

    nsize = size(pmin,1)
    nasp = AspicModel%nasp
    ntot = get_ntot()

    if (nsize.ne.size(pmax,1) &
         .or.(nsize.ne.(nextra+size(AspicModel%params,1)))) then
       stop 'get_allprior: size mismatch!'
    endif

!ln[10^10 P*]
    call get_prior_lnA(pmin(1),pmax(1))

!lnRrad or lnR
    if (useRrad) then
       call get_prior_lnRrad(pmin(2),pmax(2))
    else
       call get_prior_lnRreh(pmin(2),pmax(2))
    endif

!model dependant

    call get_aspic_priors(AspicModel%name,aspmin,aspmax,AspicModel%cmaps)
    
    pmin(nextra+1:ntot) = aspmin(1:nasp)
    pmax(nextra+1:ntot) = aspmax(1:nasp)

  end subroutine get_allprior



  function map_power_amplitude(lnA)
    implicit none
    real(kp) :: map_power_amplitude
    real(kp), intent(in) :: lnA
    real(kp) :: Pstar

!    lnA = ln[10^10 P*]
    
    Pstar = exp(lnA)*1d-10
    map_power_amplitude = Pstar

  end function map_power_amplitude




  function map_aspic_params(nasp,inparams,mapnames) result(outparams)
    implicit none
    integer, intent(in) :: nasp
    real(fmn), intent(in), dimension(nasp) :: inparams
    character(len=*), dimension(nasp), intent(in) :: mapnames
    real(kp), dimension(nasp) :: outparams

    integer :: i

    do i=1,nasp
    
       select case (mapnames(i))

       case ('flat')

          outparams(i) = inparams(i)

       case ('log')

          outparams(i) = 10._kp**(inparams(i))

       case ('ln')

          outparams(i) = exp(inparams(i))

       case ('mlog')
          
          outparams(i) = -10._kp**(inparams(i))

       case ('mln')

          outparams(i) = -exp(inparams(i))

       case default
          
          stop 'map_aspic_params: not such functions!'

       end select

    end do
   
  end function map_aspic_params




  function get_slowroll(nstar,mnParams)
    use srreheat, only : ln_potential_normalization
    use srreheat, only : ln_rho_endinf
    use srreheat, only : get_lnrreh_rrad, get_lnrrad_rreh
    use srflow, only : scalar_spectral_index
    use srflow, only : tensor_to_scalar_ratio
    use srflow, only : scalar_running

    use aspicmodels

    implicit none
    integer, intent(in) :: nstar
    real(fmn), dimension(nstar) :: get_slowroll
    real(fmn), dimension(:), intent(in) :: mnParams

    real(kp), dimension(nepsmax) :: epsStar
    real(kp), dimension(naspmax) :: asparams
    character(len=lname) :: aspname
    character(len=lname), dimension(naspmax) :: mapnames

    real(kp) :: bfoldstar, lnA
    real(kp) :: Pstar, lnRrad, lnRreh, lnM
    real(kp) :: xstar, xend
    real(kp) :: epsOneEnd
    real(kp) :: Vstar, lnRhoEnd, Vend

    integer :: nasp, ntot, neps
    integer :: i

    ntot = get_ntot()
    nasp = AspicModel%nasp
    neps = nstar - 1


    if (size(mnParams,1).ne.ntot) then
       stop 'get_slowroll: size mismatch!'
    endif
    
    lnA = mnParams(ilnA)
    Pstar = map_power_amplitude(lnA)

    if (useRrad) then
       lnRrad = mnParams(ireh)
    else
       lnRreh = mnParams(ireh)
    endif

!let's get everything from libaspic
    aspname = trim(AspicModel%name)
    forall (i=1:nasp)
       mapnames(i) = trim(AspicModel%cmaps(i))
    end forall
!    asparams(1:nasp) = mnParams(nextra+1:ntot)


    asparams(1:nasp) = map_aspic_params(nasp,mnparams(nextra+1:ntot) &
         ,mapnames(1:nasp))

    if (useRrad) then
       xstar = aspic_x_rrad(aspname,asparams,lnRrad,Pstar,bfoldstar)
    else
       xstar = aspic_x_rreh(aspname,asparams,lnRreh,bfoldstar)
    endif

    epsStar(1) = aspic_epsilon_one(aspname,xstar,asparams)
    epsStar(2) = aspic_epsilon_two(aspname,xstar,asparams)
    epsStar(3) = aspic_epsilon_three(aspname,xstar,asparams)


    Vstar = aspic_norm_potential(aspname,xstar,asparams)       

!aspicmodels returns the last params if xend is itself a param
    xend = aspic_x_endinf(aspname,asparams)

    epsOneEnd = aspic_epsilon_one(aspname,xend,asparams)
    Vend = aspic_norm_potential(aspname,xend,asparams)
                         
    lnM = ln_potential_normalization(Pstar,epsStar(1),Vstar)
    lnRhoEnd = ln_rho_endinf(Pstar,epsStar(1) &
         ,epsOneEnd,Vend/Vstar)
    
    if (useRrad) then
       lnRreh = get_lnrreh_rrad(lnRrad,lnRhoEnd)
    else
       lnRrad = get_lnrrad_rreh(lnRreh,lnRhoEnd)
    endif

!update AspicModel shared variables
    AspicModel%Pstar = Pstar
    AspicModel%lnRrad = lnRrad
    AspicModel%params(1:nasp) = mnParams(nextra+1:ntot)
    AspicModel%lnM = lnM
    AspicModel%lnRreh = lnRreh
    AspicModel%logeps = log10(epsStar(1))
    AspicModel%eps2 = epsStar(2)
    AspicModel%eps3 = epsStar(3)
    AspicModel%lnRhoEnd = lnRhoEnd
    AspicModel%bfold = bfoldstar
    AspicModel%ns = scalar_spectral_index(epsStar(1:3))
    AspicModel%logr = log10(tensor_to_scalar_ratio(epsStar(1:2)))
    AspicModel%alpha = scalar_running(epsStar(1:3))

    if (display) then
       write(*,*)
       write(*,*)'get_slowroll:',neps
       call print_aspicmodel(AspicModel)
    end if
      

!output the slow-roll params for the likelihood
    get_slowroll(1) = lnA
    get_slowroll(2) = log10(epsStar(1))
    get_slowroll(3:nstar) = epsStar(2:neps)

  end function get_slowroll

  
  subroutine print_aspicmodel(model)
    implicit none
    type(infaspic), intent(in) :: model

    write(*,*)'AspicModel Params: '
    write(*,*)'Pstar= asparams=   ', Model%Pstar,Model%params
    write(*,*)'lnRrad= lnRreh=    ', Model%lnRrad, Model%lnRreh
    write(*,*)'lnM= lnRhoEnd=     ',Model%lnM,Model%lnRhoEnd
    write(*,*)'N*-Nend=           ',Model%bfold
    write(*,*)'log(eps1)= eps23=  ',Model%logeps, Model%eps2, Model%eps3
    write(*,*)'ns= log(r)= alpha= ',Model%ns,Model%logr,Model%alpha
    write(*,*)

  end subroutine print_aspicmodel



  function test_aspic_hardprior(mnParams)
    use aspicpriors, only : check_aspic_hardprior
    implicit none    
    logical :: test_aspic_hardprior
    real(fmn), dimension(:), intent(in) :: mnParams

    real(kp), dimension(naspmax) :: asparams
    character(len=lname) :: aspname
    integer :: nasp, ntot

    ntot = get_ntot()
    nasp = AspicModel%nasp

    if (size(mnParams,1).ne.ntot) then
       stop 'test_aspic_hardprior: size mismatch!'
    endif

    !aspic params (with xend)
    asparams(1:nasp) = mnParams(nextra+1:ntot)
    aspname = trim(AspicModel%name)

    test_aspic_hardprior = check_aspic_hardprior(aspname,asparams)

  end function test_aspic_hardprior



  function test_reheating_hardprior(mnParams)   
    use cosmopar, only : lnRhoNuc
    implicit none    
    logical :: test_reheating_hardprior
    real(fmn), dimension(:), intent(in) :: mnParams
    real(fmn) :: lnRrad, lnRreh

    if (useRrad) then
       lnRrad = mnParams(ireh)       
       test_reheating_hardprior = check_lnrrad_hardprior(lnRrad)
    else
       lnRreh = mnParams(ireh)       
       test_reheating_hardprior = check_lnrreh_hardprior(lnRreh)
    endif
    
  end function test_reheating_hardprior



  function check_lnrreh_hardprior(lnRreh) result(reject)
    use cosmopar, only : lnRhoNuc
    implicit none
    logical :: reject
    real(fmn), intent(in) :: lnRreh
    real(fmn) :: lnRrehMax, lnRrehMin
    real(fmn) :: lnRhoEnd

    if (AspicModel%lnRreh.ne.lnRreh) then
       stop 'test_lnrreh_hardprior: must be called after get_slowroll!'
    endif

    lnRhoEnd = AspicModel%lnRhoEnd

    lnRrehMax = - 1._kp/12._kp * lnRhoNuc + 1._kp/3._kp * lnRhoEnd

    lnRrehMin = 1._kp/4._kp * lnRhoNuc

    reject = (lnRreh.lt.lnRrehMin).or.(lnRreh.gt.lnRrehMax) &
         .or.(lnRhoEnd.lt.lnRhoNuc)

    if (display) then
       if (reject) then
          write(*,*)
          write(*,*)'check_lnrreh_hardprior:'
          write(*,*)'lnRhoEnd= lnRhoNuc= ',lnrhoEnd,lnRhoNuc
          write(*,*)'lnRreh= (max= min=)',lnRreh,lnRrehMax,lnRrehMin
          write(*,*)
       end if
    end if

  end function check_lnrreh_hardprior



  function check_lnrrad_hardprior(lnRrad) result(reject)
    use cosmopar, only : lnRhoNuc
    implicit none
    logical :: reject
    real(fmn), intent(in) :: lnRrad
    real(fmn) :: lnRhoEnd, lnRradMax, lnRradMin

    if (AspicModel%lnRrad.ne.lnRrad) then
       stop 'test_lnrrad_hardprior: must be called after get_slowroll!'
    endif

    lnRhoEnd = AspicModel%lnRhoEnd


    lnRradMax = - 1._kp/12._kp * (lnRhoNuc - lnRhoEnd)

    lnRradMin = 0.25_kp * ( lnRhoNuc - lnRhoEnd)

    reject = (lnRrad.lt.lnRradMin).or.(lnRrad.gt.lnRradMax) &
         .or.(lnRhoEnd.lt.lnRhoNuc)

    if (display) then
       if (reject) then
          write(*,*)
          write(*,*)'check_lnrrad_hardprior:'
          write(*,*)'lnRhoEnd= lnRhoNuc= ',lnRhoEnd,lnRhoNuc
          write(*,*)'lnRrad= (max= min=)',lnRrad,lnRradMax,lnRradMin
          write(*,*)
       end if
    end if

  end function check_lnrrad_hardprior


end module wraspic
