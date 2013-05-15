module aspicpriors
  use aspicvars, only : kp, lname
  implicit none

  private

  public get_aspic_priors, check_aspic_hardprior

contains

  subroutine get_aspic_priors(name,pmin,pmax,mapname)
    use aspicmodels, only : get_aspic_numparams
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(out) :: pmin, pmax
    character(len=lname), dimension(:), intent(out) :: mapname

    integer :: naspic

    naspic = get_aspic_numparams()

    select case (name)

#define ZEROPRIOR(fooi) \
    case ('fooi')

!reminder: cpp is only used in "traditional mode" for fortran
#define ONEPRIOR(fooi,foomin,foomax,foomap) \
    case ('fooi'); \
       pmin(1) = foomin ; \
       pmax(1) = foomax ; \
       mapname(1) = 'foomap' ; \
       if (naspic.ne.1) stop 'get_aspic_priors: wrong npriors'

#define TWOPRIORS(fooi,foo1min,foo1max,foo1map,foo2min,foo2max,foo2map) \
    case ('fooi'); \
       pmin(1) = foo1min ; \
       pmax(1) = foo1max ; \
       mapname(1) = 'foo1map' ; \
       pmin(2) = foo2min ; \
       pmax(2) = foo2max ; \
       mapname(2) = 'foo2map' ; \
       if (naspic.ne.2) stop 'get_aspic_priors: wrong npriors'

#define THREEPRIORS(fooi,foo1min,foo1max,foo1map,foo2min,foo2max,foo2map \
       ,foo3min,foo3max,foo3map) \
    case ('fooi'); \
       pmin(1) = foo1min ; \
       pmax(1) = foo1max ; \
       mapname(1) = 'foo1map' ; \
       pmin(2) = foo2min ; \
       pmax(2) = foo2max ; \
       mapname(2) = 'foo2map' ; \
       pmin(3) = foo3min ; \
       pmax(3) = foo3max ; \
       mapname(3) = 'foo3map' ; \
       if (naspic.ne.3) stop 'get_aspic_priors: wrong npriors'

#include "aspicpriors.pp"
#undef ZEROPRIOR
#undef ONEPRIOR
#undef TWOPRIORS
#undef THREEPRIORS

    case default
       stop 'get_aspic_priors: model not found!'

    end select

  end subroutine get_aspic_priors


!case by case model if the prior space is not cubic
  function check_aspic_hardprior(name,params) result(reject)

    include 'aspicpriors.h'

    implicit none
    logical :: check_aspic_hardprior
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params
    logical :: reject

    real(kp), parameter :: efoldMax = 120._kp
    real(kp) :: q

    reject = .false.

    select case (name)

    case ('ii')
       reject = params(2).lt.ii_xendmin(efoldMax,params(1))

    case ('bsusybi')
       reject = params(2).gt.bsusybi_xendmax(efoldMax,params(1))

    case ('csi')
       reject = params(2).gt.csi_xendmax(efoldMax,params(1))

    case ('cnci')
       reject = params(2).lt.cnci_xendmin(efoldMax,params(1))

    case ('imi')
       reject = params(2).lt.imi_xendmin(efoldMax,params(1))

    case ('kmiii')
       reject = (params(1).lt.kmiii_alphamin(params(2))) &
            .or.(params(1).gt.kmiii_alphamax(params(2)))

    case ('sbi')
       reject = params(1).lt.sbi_alphamin(params(2))

    case ('ssbi1')
       reject = params(1).lt.ssbi1_alphamin(params(2))

    case ('ssbi3')
       reject = params(1).lt.ssbi3_alphamin(params(2))

    case ('ssbi5')
       reject = params(1).gt.ssbi5_alphamax(params(2))

    case ('ssbi6')
       reject = params(1).gt.ssbi6_alphamax(params(2))
       
    case ('vhi')
       reject = (params(3).gt.vhi_xendmax(efoldMax,params(1),params(2))) &
            .or.(params(3).lt.vhi_xendmin(params(1),params(2)))
       
    case ('dsi')
       q = 8._kp
       reject = (params(2).gt.dsi_mumax(efoldMax,params(1),q)) &
            .or.(params(3).gt.dsi_xendmax(params(1),params(2),q)) &
            .or.(params(3).lt.dsi_xendmin(efoldMax,params(1),params(2)))

    case ('cndi')
       reject = params(3).gt.cndi_xendmax(efoldMax,params(1),params(2))

    end select

  end function check_aspic_hardprior


end module aspicpriors
