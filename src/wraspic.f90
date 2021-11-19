module wraspic
  use sampl, only : fmn, imn
  use aspicvars, only : kp, infaspic, lname
  implicit none

  private


  integer, parameter :: ilnA = 1
  integer, parameter :: ireh= 2
  integer, parameter :: iw=3

  integer, save :: nextra = 0


  character(len=*), parameter :: ReheatModel = 'Rreh'
!  character(len=*), parameter :: ReheatModel = 'Rrad'
!  character(len=*), parameter :: ReheatModel = 'Rhow'



  type(infaspic), save :: AspicModel

  integer, parameter :: nparmax = 4
  integer, parameter :: nepsmax = 3

  public ReheatModel
  public set_model, check_model, free_model, allocate_and_set_allprior
  public get_hubbleflow, get_ntot, get_nextra, get_derived, get_derived_name
  public test_aspic_hardprior, test_reheating_hardprior

  logical, parameter :: display = .false.

contains



  subroutine set_model(mname,msubname)
    use aspicmodels, only : initialize_aspic_ptrs
    use aspicmodels, only : get_aspic_numparams

    implicit none
    character(len=*), intent(in) :: mname
    character(len=*), intent(in), optional :: msubname

    call initialize_aspic_ptrs(trim(adjustl(mname)))

    AspicModel%name = trim(adjustl(mname))
    if (present(msubname)) then
       AspicModel%extname = trim(adjustl(mname))//trim(adjustl(msubname))
    else
       AspicModel%extname = AspicModel%name
    endif

    AspicModel%nasp = get_aspic_numparams()
    AspicModel%nhid = 0
    
    if (AspicModel%nasp.gt.nparmax) stop 'set_model: nparmax too small!'

  end subroutine set_model



  function check_model()
    implicit none
    logical :: check_model

    check_model = associated(AspicModel%params)

  end function check_model
  

  subroutine free_model()
    use aspicmodels, only : free_aspic_ptrs
    implicit none

    if (.not.check_model()) stop 'free_model: nothing allocated!'
       
    call free_aspic_ptrs()

    if (associated(AspicModel%params)) deallocate(AspicModel%params)
    AspicModel%params => null()
    if (associated(AspicModel%cmaps)) deallocate(AspicModel%cmaps)
    AspicModel%cmaps => null()

  end subroutine free_model



  function get_ntot()
    implicit none
    integer :: get_ntot

    if (.not.check_model()) stop 'get_ntot: aspic not set!'    

    get_ntot = AspicModel%nasp + AspicModel%nhid + nextra
    
  end function get_ntot


  function get_nhid()
    implicit none
    integer :: get_nhid
    if (.not.check_model()) stop 'get_nhid: aspic not set!'    

    get_nhid = AspicModel%nhid

  end function get_nhid

  

  function get_nextra()
    implicit none
    integer :: get_nextra

    if (.not.check_model()) stop 'get_nextra: aspic not set!'    

    get_nextra = nextra

  end function get_nextra




  subroutine get_prior_lnA(lnAmin,lnAmax)
    implicit none
    real(fmn) , intent(out) :: lnAmin,lnAmax

    lnAmin = 2.4
    lnAmax = 3.6
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



  subroutine get_prior_lnRhoReh(lnRhoRehMin,lnRhoRehMax)
    use cosmopar, only : lnRhoNuc
    implicit none
    real(fmn), intent(out) :: lnRhoRehMin, lnRhoRehMax
    
    lnRhoRehMin = lnRhoNuc
    lnRhoRehMax = 0._kp
  end subroutine get_prior_lnRhoReh



  subroutine get_prior_wreh(wmin,wmax)
    use cosmopar, only : lnRhoNuc
    implicit none
    real(fmn), intent(out) :: wmin, wmax

    wmin = -1._kp/3._kp
    wmax = 1._kp

  end subroutine get_prior_wreh


  
  function get_derived(i)
    implicit none
    integer, intent(in) :: i
    real(fmn) :: get_derived

    select case(i)

    case (1)
       get_derived = AspicModel%lnM
            
    case (2)
       get_derived = AspicModel%bfold

    case (3)
       get_derived = AspicModel%logeps

    case (4)
       get_derived = AspicModel%eps2

    case (5)
       get_derived = AspicModel%eps3

    case (6)
       get_derived = AspicModel%ns

    case (7)
       get_derived = AspicModel%logr

    case (8)
       get_derived = AspicModel%alpha   

    case (9)       
       get_derived = AspicModel%lnRhoEnd
       
    case (10)
       get_derived = AspicModel%lnRreh
       
    case (11)
       get_derived = AspicModel%lnRrad
       
    case default
       stop 'get_derived: incorrect parameters number!'

    end select


  end function get_derived


!should be consistent with the get_derived, not check made.  
  function get_derived_name(i)
    implicit none
    integer, intent(in) :: i
    integer, parameter :: lenderived = 40
    character(len=lenderived) :: get_derived_name

    select case(i)

    case (1)
       get_derived_name = 'lnM*           \ln(M)                  '
          
    case (2)
       get_derived_name = 'bfold*         N_{\mathrm{end}}-N_*    '

    case (3)
       get_derived_name = 'logeps*        \log(\epsilon_1)        '

    case (4)
       get_derived_name = 'eps2*          \epsilon_2              '

    case (5)
       get_derived_name = 'eps3*          \epsilon_3              '

    case (6)
       get_derived_name = 'ns*            n_{\mathrm{S}}          '

    case (7)
       get_derived_name = 'logr*          \log(r_\epsilon)        '

    case (8)
       get_derived_name = 'alpha*         \alpha_{\mathrm{S}}     '

    case (9)       
       get_derived_name = 'lnRhoEnd*      \ln(\rho_{\mathrm{end}})'

    case (10)
       get_derived_name = 'lnRreh*        \ln(R_{\mathrm{reh}})   '
       
    case (11)
       get_derived_name = 'lnRrad*        \ln(R_{\mathrm{rad}})   '
       
    case default
       stop 'get_derived_name: incorrect parameters number!'

    end select


  end function get_derived_name

  

  subroutine allocate_and_set_allprior(pmin,pmax)
    use aspicpriors, only : get_aspic_priors, get_aspic_numpriors
    implicit none
    real(fmn), dimension(:), allocatable, intent(out) :: pmin, pmax

    real(kp), dimension(nparmax) :: aspmin,aspmax
    character(len=lname), dimension(nparmax) :: ascmaps

    integer :: npar, ntot, i
    
!model dependant

    call get_aspic_priors(AspicModel%extname,aspmin,aspmax,ascmaps)

    npar = get_aspic_numpriors()

    if (npar.ge.AspicModel%nasp) then
       AspicModel%nhid = npar - AspicModel%nasp
    else
       stop 'allocate_and_set_allprior: nprior < naspic!'
    endif
    
    
    if (npar.gt.nparmax) stop 'allocate_and_set_allprior: npar > nparmax!'
    
    allocate(AspicModel%cmaps(npar))
    allocate(AspicModel%params(npar))

    do i=1,npar
       AspicModel%cmaps(i) = ascmaps(i)
    enddo
       

    
!CMB amplitude
    nextra = 1

    select case (ReheatModel)

    case ('Rrad','Rreh')
       nextra = nextra + 1

    case ('Rhow')
       nextra = nextra + 2

    case default
       stop 'allocate_and_set_allprior: not such a reheating modelisation!'

    end select
       
    ntot = npar + nextra
    allocate(pmin(ntot),pmax(ntot))
    

!ln[10^10 P*]
    call get_prior_lnA(pmin(ilnA),pmax(ilnA))
    
!lnRrad or lnR
    select case (ReheatModel)

    case ('Rrad')

       call get_prior_lnRrad(pmin(ireh),pmax(ireh))

    case ('Rreh')

       call get_prior_lnRreh(pmin(ireh),pmax(ireh))
       
    case ('Rhow')

       call get_prior_lnRhoReh(pmin(ireh),pmax(ireh))
       call get_prior_wreh(pmin(iw),pmax(iw))
       
    case default

       stop 'allocate_and_set_allprior: not such a reheating modelisation!'

    end select


!model parameter are in the last bits    
    pmin(nextra+1:ntot) = aspmin(1:npar)
    pmax(nextra+1:ntot) = aspmax(1:npar)

  end subroutine allocate_and_set_allprior

  
  
  function map_power_amplitude(lnA)
    implicit none
    real(kp) :: map_power_amplitude
    real(kp), intent(in) :: lnA
    real(kp) :: Pstar

!    lnA = ln[10^10 P*]
    
    Pstar = exp(lnA)*1d-10
    map_power_amplitude = Pstar

  end function map_power_amplitude




  function map_aspic_params(npar,inparams,mapnames) result(outparams)
    use aspicpriors, only : redefine_aspic_params
    implicit none
    integer, intent(in) :: npar
    real(fmn), intent(in), dimension(npar) :: inparams
    character(len=*), dimension(npar), intent(in) :: mapnames
    real(kp), dimension(npar) :: outparams, scalparams

    integer :: i

    do i=1,npar
    
       select case (trim(mapnames(i)))

       case ('flat')

          scalparams(i) = inparams(i)

       case ('log')

          scalparams(i) = 10._kp**(inparams(i))

       case ('ln')

          scalparams(i) = exp(inparams(i))

       case ('mlog')
          
          scalparams(i) = -10._kp**(inparams(i))

       case ('mln')

          scalparams(i) = -exp(inparams(i))

       case ('inv')

          scalparams(i) = 1._kp/inparams(i)

       case ('invsqrt')
          stop 'this one is weird!'
          scalparams(i) = 1._kp/sqrt(inparams(i))

       
       case default

          write(*,*)'mapname= ',mapnames(i)
          stop 'map_aspic_params: not such priors!'

       end select

    end do

!some models have redefined parameters, because this relationships are
!done on the fundamental aspic parameters, this is called after the
!sampled parameter scaling
    outparams =  redefine_aspic_params(AspicModel%extname,scalparams)
   
   
  end function map_aspic_params



  function get_hubbleflow(nstar,mnParams)
    use srreheat, only : potential_normalization
    use srreheat, only : ln_rho_endinf
    use srreheat, only : get_lnrreh_rrad, get_lnrrad_rreh
    use srreheat, only : get_lnrreh_rhow, get_lnrrad_rhow
    use srflow, only : slowroll_to_hubble
    use srflow, only : scalar_spectral_index
    use srflow, only : tensor_to_scalar_ratio
    use srflow, only : scalar_running

    use aspicmodels

    implicit none
    integer, intent(in) :: nstar
    real(fmn), dimension(nstar) :: get_hubbleflow
    real(fmn), dimension(:), intent(in) :: mnParams

    real(kp), dimension(nepsmax) :: epsVStar, epsHStar
    real(kp), dimension(nparmax) :: asparams
    character(len=lname) :: aspname
    character(len=lname), dimension(nparmax) :: mapnames

    real(kp) :: bfoldstar, lnA
    real(kp) :: Pstar, lnRrad, lnRreh, lnRhoReh, w, lnM
    real(kp) :: xstar, xend
    real(kp) :: epsOneEnd
    real(kp) :: Vstar, lnRhoEnd, Vend

    integer :: npar, ntot, neps
    integer :: i

    ntot = get_ntot()
    npar = AspicModel%nasp + AspicModel%nhid
    neps = nstar - 1


    if (size(mnParams,1).ne.ntot) then
       stop 'get_hubbleflow: size mismatch!'
    endif
    
    lnA = mnParams(ilnA)
    Pstar = map_power_amplitude(lnA)
   

!let's get everything from libaspic
    aspname = trim(AspicModel%name)
    do i=1,npar
       mapnames(i) = AspicModel%cmaps(i)
    enddo
!    asparams(1:npar) = mnParams(nextra+1:ntot)

    asparams(1:npar) = map_aspic_params(npar,mnparams(nextra+1:ntot) &
         ,mapnames(1:npar))


    select case (ReheatModel)

    case ('Rrad')

       lnRrad = mnParams(ireh)
       xend = aspic_x_endinf(aspname,asparams)
       xstar = aspic_x_rrad(aspname,asparams,xend,lnRrad,Pstar,bfoldstar)

    case ('Rreh')

       lnRreh = mnParams(ireh)
       xend = aspic_x_endinf(aspname,asparams)
       xstar = aspic_x_rreh(aspname,asparams,xend,lnRreh,bfoldstar)

    case ('Rhow')
       lnRhoReh = mnParams(ireh)
       w = mnParams(iw)
       xend = aspic_x_endinf(aspname,asparams)
       xstar = aspic_x_rhow(aspname,asparams,xend,w,lnRhoReh,Pstar,bfoldstar)

    case default
       stop 'get_hubbleflow: not such a reheating modelisation!'

    end select

    epsVStar(1) = aspic_epsilon_one(aspname,xstar,asparams)
    epsVStar(2) = aspic_epsilon_two(aspname,xstar,asparams)
    epsVStar(3) = aspic_epsilon_three(aspname,xstar,asparams)

    epsHstar = slowroll_to_hubble(epsVstar)
    
    Vstar = aspic_norm_potential(aspname,xstar,asparams)       

    epsOneEnd = aspic_epsilon_one(aspname,xend,asparams)
    Vend = aspic_norm_potential(aspname,xend,asparams)
                         
    lnM = log(potential_normalization(Pstar,epsHStar(1),Vstar))
    lnRhoEnd = ln_rho_endinf(Pstar,epsHStar(1) &
         ,epsOneEnd,Vend/Vstar)

    select case (ReheatModel)
    case ('Rrad')
       lnRreh = get_lnrreh_rrad(lnRrad,lnRhoEnd)
    case ('Rreh')
       lnRrad = get_lnrrad_rreh(lnRreh,lnRhoEnd)
    case ('Rhow')
       lnRreh = get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd)
       lnRrad = get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd)
    case default
       stop 'get_hubbleflow: not such a reheating modelisation!'
    end select

!update AspicModel shared variables
    AspicModel%Pstar = Pstar
    AspicModel%lnRrad = lnRrad
!better displaying the aspic params rather than the mnparams
    AspicModel%params(1:npar) = asparams(1:npar)
    AspicModel%lnM = lnM
    AspicModel%lnRreh = lnRreh
    AspicModel%logeps = log10(epsVStar(1))
    AspicModel%eps2 = epsVStar(2)
    AspicModel%eps3 = epsVStar(3)
    AspicModel%lnRhoEnd = lnRhoEnd
    AspicModel%bfold = bfoldstar
    AspicModel%ns = scalar_spectral_index(epsVStar(1:3))
    AspicModel%logr = log10(tensor_to_scalar_ratio(epsVStar(1:2)))
    AspicModel%alpha = scalar_running(epsVStar(1:3))

    if (display) then
       write(*,*)
       write(*,*)'get_hubbleflow:',neps
       call print_aspicmodel(AspicModel)
    end if
      

!output the slow-roll params for the likelihood
    get_hubbleflow(1) = lnA
    get_hubbleflow(2) = log10(epsHStar(1))
    get_hubbleflow(3:nstar) = epsHStar(2:neps)

  end function get_hubbleflow

  
  subroutine print_aspicmodel(model)
    implicit none
    type(infaspic), intent(in) :: model

    write(*,*)'AspicModel Params:   '
    write(*,*)'Pstar= asparams=     ',Model%Pstar,Model%params
    write(*,*)'lnRrad= lnRreh=      ',Model%lnRrad, Model%lnRreh
    write(*,*)'lnM= lnRhoEnd=       ',Model%lnM,Model%lnRhoEnd
    write(*,*)'N*-Nend=             ',Model%bfold
    write(*,*)'log(epsV1)= epsV23=  ',Model%logeps, Model%eps2, Model%eps3
    write(*,*)'ns= log(r)= alpha=   ',Model%ns,Model%logr,Model%alpha
    write(*,*)

  end subroutine print_aspicmodel



  function test_aspic_hardprior(mnParams)
    use aspicpriors, only : check_aspic_hardprior
    implicit none    
    logical :: test_aspic_hardprior
    real(fmn), dimension(:), intent(in) :: mnParams

    real(kp), dimension(nparmax) :: asparams
    character(len=lname), dimension(nparmax) :: mapnames
    character(len=lname) :: extname
    integer :: npar, ntot, i

    ntot = get_ntot()
    npar = AspicModel%nasp + AspicModel%nhid

    if (size(mnParams,1).ne.ntot) then
       stop 'test_aspic_hardprior: size mismatch!'
    endif

    extname = trim(AspicModel%extname)

    do i=1,npar
       mapnames(i) = AspicModel%cmaps(i)
    end do

    asparams(1:npar) = map_aspic_params(npar,mnparams(nextra+1:ntot) &
         ,mapnames(1:npar))

    test_aspic_hardprior = check_aspic_hardprior(extname,asparams(1:npar))

  end function test_aspic_hardprior



  function test_reheating_hardprior(mnParams)   
    use cosmopar, only : lnRhoNuc
    implicit none    
    logical :: test_reheating_hardprior
    real(fmn), dimension(:), intent(in) :: mnParams
    real(fmn) :: lnRrad, lnRreh, lnRhoReh,w

    select case (ReheatModel)
    case ('Rrad')
       lnRrad = mnParams(ireh)       
       test_reheating_hardprior = check_lnrrad_hardprior(lnRrad)
    case ('Rreh')
       lnRreh = mnParams(ireh)       
       test_reheating_hardprior = check_lnrreh_hardprior(lnRreh)
    case ('Rhow')
       lnRhoReh = mnParams(ireh)
       w = mnParams(iw)
       test_reheating_hardprior = check_rhow_hardprior(lnRhoReh,w)
    case default
       stop 'test_reheating_hardprior: not such a reheating modelisation!'
    end select

    
  end function test_reheating_hardprior



  function check_lnrreh_hardprior(lnRreh) result(reject)
    use cosmopar, only : lnRhoNuc
    implicit none
    logical :: reject
    real(fmn), intent(in) :: lnRreh
    real(fmn) :: lnRrehMax, lnRrehMin
    real(fmn) :: lnRhoEnd

    if (AspicModel%lnRreh.ne.lnRreh) then
       stop 'test_lnrreh_hardprior: must be called after get_hubbleflow!'
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
       stop 'test_lnrrad_hardprior: must be called after get_hubbleflow!'
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



  function check_rhow_hardprior(lnRhoReh,w) result(reject)
    use cosmopar, only : lnRhoNuc
    implicit none
    logical :: reject
    real(fmn), intent(in) :: lnRhoReh,w
    real(fmn) :: lnRhoEnd
    
    lnRhoEnd = AspicModel%lnRhoEnd
    
    reject = (lnRhoReh.lt.lnRhoNuc).or.(lnRhoReh.gt.lnRhoEnd) &
         .or.(lnRhoEnd.lt.lnRhoNuc)

    if (display) then
       if (reject) then
          write(*,*)
          write(*,*)'check_rhow_hardprior:'
          write(*,*)'lnRhoEnd= lnRhoNuc= ',lnRhoEnd,lnRhoNuc
          write(*,*)'wreh= lnRhoReh= ',w, lnRhoReh
          write(*,*)
       end if
    end if

  end function check_rhow_hardprior

end module wraspic
