module aspicwrap
  use nestprec, only : fmn, imn
  use aspicvars, only : kp, infaspic
  implicit none

  private

  integer, parameter :: nextra = 2

  type(infaspic), save :: AspicModel

  public set_aspic_model, check_aspic, get_aspic_allprior
  public get_aspic_slowroll, get_aspic_ntot, get_aspic_derived


contains

  
  subroutine set_aspic_model(mname)
    implicit none
    character(len=*), intent(in) :: mname

    AspicModel%name = trim(adjustl(mname))
    
    AspicModel%nasp = 1

    allocate(AspicModel%params(AspicModel%nasp))

  end subroutine set_aspic_model



  function check_aspic()
    implicit none
    logical :: check_aspic

    check_aspic = associated(AspicModel%params)

  end function check_aspic
  

  function get_aspic_ntot()
    implicit none
    integer :: get_aspic_ntot

    if (.not.check_aspic()) stop 'get_aspic_ntot: aspic not set!'

    get_aspic_ntot = AspicModel%nasp + nextra
    
  end function get_aspic_ntot



  subroutine get_aspic_allprior(pmin,pmax)
    implicit none
    real(kp), dimension(:), intent(out) :: pmin, pmax

    integer :: nsize

    nsize = size(pmin,1)
    if (nsize.ne.size(pmax,1).or.(nsize.ne.(nextra+size(AspicModel%params,1)))) then
       stop 'get_aspic_allprior: size mismatch!'
    endif

!ln[10^10 P*]
    pmin(1) = 2.7
    pmax(1) = 4
!lnRrad
    pmin(2) = -45
    pmax(2) = 10
!model dependant

    select case (AspicModel%name)

    case ('lfi')
       pmin(3) = 0.1
       pmax(3) = 5

    case default
       stop 'get_aspic_allprior: model not found!'

    end select


  end subroutine get_aspic_allprior


#ifdef ASPIC
  function get_aspic_slowroll(nstar,mnParams)
    use aspicvars, only : neps
    include 'aspic.h'
    implicit none
    integer, intent(in) :: nstar
    real(kp), dimension(nstar) :: get_aspic_slowroll
    real(fmn), dimension(:), intent(in) :: mnParams

    integer :: nasp, ntot
    real(kp) :: bfoldstar
    real(kp) :: Pstar, lnRrad, lnR, lnM
    real(kp) :: xstar, xend

    real(kp), dimension(neps) :: epsStar
    real(kp) :: epsOneEnd

    real(kp) :: Vstar, lnRhoEnd, Vend

    ntot = get_aspic_ntot()
    nasp = AspicModel%nasp

    if (size(mnParams,1).ne.ntot) then
       stop 'get_aspic_slowroll: size mismatch!'
    endif
    
    Pstar = mnParams(1)
    lnRrad = mnParams(2)


    select case (AspicModel%name)

    case ('lfi')
       
       xstar = lfi_x_rrad(mnParams(1),lnRrad,Pstar,bfoldstar)
       epsStar(1) = lfi_epsilon_one(xstar,mnParams(1))
       epsStar(2) = lfi_epsilon_two(xstar,mnParams(1))
       epsStar(3) = lfi_epsilon_three(xstar,mnParams(1))
       Vstar = lfi_norm_potential(xstar,mnParams(1))
       lnM = ln_potential_normalization(Pstar,epsStar(1),Vstar)
       xend = lfi_x_endinf(mnParams(1))
       epsOneEnd = lfi_epsilon_one(xend,mnParams(1))
       Vstar = lfi_norm_potential(xend,mnParams(1))
              
       lnRhoEnd = ln_rho_endinf(Pstar,epsStar(1) &
            ,epsOneEnd,Vend/Vstar)
      
    case default
       stop 'get_aspic_slowroll: model not found!'
    end select
    
!update AspicModel shared variables    
    AspicModel%Pstar = Pstar
    AspicModel%lnRrad = lnRrad
    AspicModel%params(1:nasp) = mnParams(nextra+1:ntot)
    AspicModel%lnRhoEnd = lnRhoEnd
    AspicModel%lnM = lnM
    AspicModel%lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
    AspicModel%bfold = bfoldstar
    AspicModel%eps(:) = epsStar(:)
    AspicModel%ns = scalar_spectral_index(epsStar)
    AspicModel%r = tensor_to_scalar_ratio(epsStar(1:2))
    AspicModel%alpha = scalar_running(epsStar)

!out the slow-roll params
    get_aspic_slowroll(1) = mnParams(1)
    get_aspic_slowroll(2) = log10(epsStar(1))
    get_aspic_slowroll(3:4) = epsStar(2:3)


  end function  get_aspic_slowroll
#else

  function get_aspic_slowroll(nstar,mnParams)
    use aspicvars, only : neps
    implicit none
    integer, intent(in) :: nstar
    real(kp), dimension(nstar) :: get_aspic_slowroll
    real(fmn), dimension(:), intent(in) :: mnParams

    write(*,*)'ASPIC NOT DEFINED...'
    get_aspic_slowroll = mnParams

  end function get_aspic_slowroll

#endif

  function get_aspic_derived(i)
    implicit none
    integer, intent(in) :: i
    real(kp) :: get_aspic_derived

    select case(i)

    case (1)
       get_aspic_derived = AspicModel%lnM

    case (2)       
       get_aspic_derived = AspicModel%lnRhoEnd

    case(3)
       get_aspic_derived = AspicModel%lnRrad
       
    case (4)
       get_aspic_derived = AspicModel%ns

    case (5)
       get_aspic_derived = AspicModel%r

    case (6)
       get_aspic_derived = AspicModel%alpha

    case (7)
       get_aspic_derived = AspicModel%bfold

    case default
       stop 'get_aspic_derived: incorrect parameters number!'

    end select


  end function get_aspic_derived


end module aspicwrap
