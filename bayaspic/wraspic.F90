module wraspic
  use nestprec, only : fmn, imn
  use aspicvars, only : kp, infaspic, lname
  implicit none

  private

  integer, parameter :: nextra = 2
  logical, parameter :: useRrad = .true.


  type(infaspic), save :: AspicModel

  integer, parameter :: naspmax = 3
  integer, parameter :: nepsmax = 3

  public set_model, check_model, get_allprior
  public get_slowroll, get_ntot, get_derived


contains

  
  subroutine set_model(mname)
    implicit none
    character(len=*), intent(in) :: mname

    AspicModel%name = trim(adjustl(mname))
    
    AspicModel%nasp = 1

    allocate(AspicModel%params(AspicModel%nasp))

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
    real(kp) , intent(out) :: lnAmin,lnAmax

    lnAmin = 2.7
    lnAmax = 4

  end subroutine get_prior_lnA



  subroutine get_prior_rrad(lnRradmin,lnRradmax)
    implicit none
    real(kp) , intent(out) :: lnRradmin,lnRradmax

    lnRradMin = -46
    lnRradMax = 10

  end subroutine get_prior_rrad



  subroutine get_prior_rreh(lnRmin,lnRmax)
    implicit none
    real(kp) , intent(out) :: lnRmin,lnRmax

    lnRMin = -46
    lnRMax = 10

  end subroutine get_prior_rreh


  
  function get_derived(i)
    implicit none
    integer, intent(in) :: i
    real(kp) :: get_derived

    select case(i)

    case (1)
       get_derived = AspicModel%lnM

    case (2)       
       get_derived = AspicModel%lnRhoEnd

    case(3)
       if (useRrad) then
          get_derived = AspicModel%lnRreh
       else
          get_derived = AspicModel%lnRrad
       endif
       
    case (4)
       get_derived = AspicModel%ns

    case (5)
       get_derived = AspicModel%r

    case (6)
       get_derived = AspicModel%alpha

    case (7)
       get_derived = AspicModel%bfold

    case default
       stop 'get_derived: incorrect parameters number!'

    end select


  end function get_derived


#ifdef ASPIC

  subroutine get_allprior(pmin,pmax)
    use aspicpriors, only : get_aspic_priors
    implicit none
    real(kp), dimension(:), intent(out) :: pmin, pmax

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
       call get_prior_rrad(pmin(2),pmax(2))
    else
       call get_prior_rreh(pmin(2),pmax(2))
    endif

!model dependant

    call get_aspic_priors(AspicModel%name,aspmin,aspmax)
    
    pmin(nextra+1:ntot) = aspmin(1:nasp)
    pmax(nextra+1:ntot) = aspmax(1:nasp)

  end subroutine get_allprior



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
    real(kp), dimension(nstar) :: get_slowroll
    real(fmn), dimension(:), intent(in) :: mnParams

    real(kp), dimension(nepsmax) :: epsStar
    real(kp), dimension(naspmax) :: asparams
    character(len=lname) :: aspname

    real(kp) :: bfoldstar
    real(kp) :: Pstar, lnRrad, lnRreh, lnM
    real(kp) :: xstar, xend
    real(kp) :: epsOneEnd
    real(kp) :: Vstar, lnRhoEnd, Vend

    integer :: nasp, ntot, neps


    ntot = get_ntot()
    nasp = AspicModel%nasp
    neps = nstar - 1

    if (size(mnParams,1).ne.ntot) then
       stop 'get_slowroll: size mismatch!'
    endif
    
    Pstar = mnParams(1)

    if (useRrad) then
       lnRrad = mnParams(2)
    else
       lnRreh = mnParams(2)
    endif

!let's get everything from libaspic
    asparams(1:nasp) = mnParams(nextra+1:ntot)
    aspname = trim(AspicModel%name)

    if (useRrad) then
       xstar = aspic_x_rrad(aspname,asparams,lnRrad,Pstar,bfoldstar)
    else
       xstar = aspic_x_rreh(aspname,asparams,lnRreh,bfoldstar)
    endif

    epsStar(1) = aspic_epsilon_one(aspname,xstar,asparams)
    epsStar(2) = aspic_epsilon_two(aspname,xstar,asparams)

    if (neps.eq.3) then
       epsStar(3) = aspic_epsilon_three(aspname,xstar,asparams)
    endif

    Vstar = aspic_norm_potential(aspname,xstar,asparams)       

    xend = aspic_x_endinf(aspname,asparams)
    epsOneEnd = aspic_epsilon_one(aspname,xend,asparams)
    Vstar = aspic_norm_potential(aspname,xend,asparams)
                         
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
    AspicModel%lnRhoEnd = lnRhoEnd
    AspicModel%lnM = lnM
    AspicModel%lnRreh = lnRreh
    AspicModel%bfold = bfoldstar
    AspicModel%ns = scalar_spectral_index(epsStar(1:neps))
    AspicModel%r = tensor_to_scalar_ratio(epsStar(1:2))
    AspicModel%alpha = scalar_running(epsStar)

!output the slow-roll params for the likelihood
    get_slowroll(1) = Pstar
    get_slowroll(2) = log10(epsStar(1))
    get_slowroll(3:nstar) = epsStar(2:neps)

  end function  get_slowroll

#else

  function get_slowroll(nstar,mnParams)
    use aspicvars, only : neps
    implicit none
    integer, intent(in) :: nstar
    real(kp), dimension(nstar) :: get_slowroll
    real(fmn), dimension(:), intent(in) :: mnParams

    write(*,*)'ASPIC NOT DEFINED...'
    get_slowroll = mnParams

  end function get_slowroll

  subroutine get_allprior(pmin,pmax)    
    implicit none
    real(kp), dimension(:), intent(out) :: pmin, pmax

    integer :: nsize, nasp
    
    nsize = size(pmin,1)
    nasp = AspicModel%nasp

    if (nsize.ne.size(pmax,1) &
         .or.(nsize.ne.(nextra+size(AspicModel%params,1)))) then
       stop 'get_allprior: size mismatch!'
    endif

!ln[10^10 P*]
    call get_prior_lnA(pmin(1),pmax(1))

!lnRrad or lnR
    if (useRrad) then
       call get_prior_rrad(pmin(2),pmax(2))
    else
       call get_prior_rreh(pmin(2),pmax(2))
    endif

!model dependant    
    pmin(nextra+1:ntot) = 0.1
    pmax(nextra+1:ntot) = 1

  end subroutine get_allprior

#endif

end module wraspic
