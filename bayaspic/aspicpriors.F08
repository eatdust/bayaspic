module aspicpriors
  use aspicvars, only : kp, pi, lname
  implicit none

  private

  public get_aspic_priors, check_aspic_hardprior
  public redefine_aspic_params

!hard prior: all model are assumed to support, at lest, this number of efolds
  real(kp), parameter :: efoldNum = 120._kp

!dsi correction power
  real(kp), parameter :: qDSI = 8._kp

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
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params
    logical :: reject
    real(kp) :: eps2max

    reject = .false.

    select case (name)

    case ('ii', 'iif', 'iilambda')
       reject = params(2).lt.ii_xendmin(efoldNum,params(1))

    case ('lmi1o', 'lmi1p')
       eps2max=0.2_kp
       reject = lmi_epstwo_potmax(params(1),params(2)).gt.eps2max      

    case ('lmi2o','lmi2p')
!Such that exp(-beta xend^gamma)>10^(-100) to make it numerically
!tractable
       reject = params(3).gt.(100._kp/params(2))**(1._kp/params(1))
       
    case ('bsusybif', 'bsusybil')
       reject = params(2).gt.bsusybi_xendmax(efoldNum,params(1))

    case ('gmssmiB', 'gmssmiD', 'gmssmiF', 'gmssmiH')
       reject = abs(1._kp-params(1)).gt.params(2)**4*pi**2/(900._kp*60._kp**2)

    case ('gripiB', 'gripiD', 'gripiF', 'gripiH')
       reject = abs(1._kp-params(1)).gt.params(2)**4*pi**2/(576._kp*60._kp**2)

    case ('tiB', 'tiD', 'tiF')
       reject = abs(1._kp-2._kp*params(1)).gt.0.1_kp*params(2)**2

    case ('psniB')
       reject = params(1).gt.0.1_kp*params(2)**2

    case ('csi')
       reject = params(2).gt.csi_xendmax(efoldNum,params(1))

    case ('cnci')
       reject = params(2).lt.cnci_xendmin(efoldNum,params(1))

    case ('sbi')
       reject = params(1).lt.sbi_alphamin(params(2))*(1._kp+epsilon(1._kp))

    case ('imi')
       reject = params(2).lt.imi_xendmin(efoldNum,params(1))

    case ('kmiii')
       reject = (params(1).lt.kmiii_alphamin(params(2))) &
            .or.(params(1).gt.kmiii_alphamax(params(2)))

    case ('ssbi1')
       reject = params(1).lt.ssbi1_alphamin(params(2))

    case ('ssbi3')
       reject = params(1).lt.ssbi3_alphamin(params(2))

    case ('ssbi5')
       reject = params(1).gt.ssbi5_alphamax(params(2))

    case ('ssbi6')
       reject = params(1).gt.ssbi6_alphamax(params(2))
              
    case ('rmi1')
!       reject = params(3).gt.rmi1_numacc_xendmax(70._kp,params(1),params(2))
       reject = params(3).gt.rmi1_numacc_xendmax(efoldNum,params(1),params(2))

    case ('rmi2')
!       reject = params(3).lt.rmi2_numacc_xendmin(70._kp,params(1),params(2))
       reject = params(3).lt.rmi2_numacc_xendmin(efoldNum,params(1),params(2))

    case ('rmi3')
       reject = params(3).gt.rmi3_numacc_xendmax(params(1),params(2))

    case ('rmi4')
       reject = params(3).lt.rmi4_numacc_xendmin(params(1),params(2))

    case ('dsio', 'dsi2')
       reject = params(2) .gt. dsi_mumax(efoldNum,params(1),qDsi)

    end select

  end function check_aspic_hardprior



!redefinition of aspic params (called after change of priors, see wraspic)
  function redefine_aspic_params(name,nasp,newparams) result(asparams)
    include 'aspicpriors.h'
    implicit none   
    character(len=*), intent(in) :: name
    integer, intent(in) :: nasp
    real(kp), dimension(nasp), intent(in) :: newparams
    real(kp), dimension(nasp) :: asparams

    asparams = newparams

    select case (name)

    case ('lmi1o')
       asparams(1) = 2._kp/(newparams(1)+1._kp)
       asparams(2) = 2._kp*((newparams(1)+1._kp)/(2._kp* &
            sqrt(2._kp*newparams(1)*newparams(2))))** &
            (2._kp/(newparams(1)+1._kp))

    case ('lmi2o')
       asparams(1) = 2._kp/(newparams(1)+1._kp)
       asparams(2) = 2._kp*((newparams(1)+1._kp)/(2._kp* &
            sqrt(2._kp*newparams(1)*newparams(2))))** &
            (2._kp/(newparams(1)+1._kp))
       asparams(3) = newparams(3) * lmi2_xendmin(efoldNum,asparams(1),asparams(2))

    case ('lmi2p')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2)
       asparams(3) = newparams(3) * lmi2_xendmin(efoldNum,asparams(1),asparams(2))

    case ('twi')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2)*newparams(1)
 
    case ('csi')
       asparams(1) = newparams(1)
!       asparams(2) = newparams(2) * csi_xendmax(59._kp,asparams(1))
       asparams(2) = newparams(2) * csi_xendmax(efoldNum,asparams(1))

    case ('cnci')
       asparams(1) = newparams(1)
!       asparams(2) = newparams(2) * cnci_xendmin(58._kp,asparams(1))
       asparams(2) = newparams(2) * cnci_xendmin(efoldNum,asparams(1))

    case ('imi')
       asparams(1) = newparams(1)
!       asparams(2) = newparams(2) * imi_xendmin(65._kp,asparams(1))
       asparams(2) = newparams(2) * imi_xendmin(efoldNum,asparams(1))

    case ('gmssmiA', 'gmssmiB', 'gmssmiC', 'gmssmiD')
       asparams(1) = 1._kp + newparams(1)
       asparams(2) = newparams(2)

    case ('gmssmiE', 'gmssmiF', 'gmssmiG', 'gmssmiH')
       asparams(1) = 1._kp - newparams(1)
       asparams(2) = newparams(2)

    case ('gripiA', 'gripiB', 'gripiC', 'gripiD')
       asparams(1) = 1._kp + newparams(1)
       asparams(2) = newparams(2)

    case ('gripiE', 'gripiF', 'gripiG', 'gripiH')
       asparams(1) = 1._kp - newparams(1)
       asparams(2) = newparams(2)

    case ('iif')
       asparams(1) = 4._kp * (newparams(1) - 1._kp)
       asparams(2) = newparams(2)

    case ('iilambda')
       asparams(1) = 2._kp/(newparams(1) - 1._kp)
       asparams(2) = newparams(2)

    case ('cndi','cndi+')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2)
       asparams(3) = newparams(3) &
            * cndi_xendmax(efoldNum,asparams(1),asparams(2))

    case ('psnift')
       asparams(1) = newparams(1)*newparams(2)**2
       asparams(2) = newparams(2)

    case ('vhi')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2)
       asparams(3) = vhi_xendmin(asparams(1),asparams(2)) &
            + newparams(3) * (vhi_xendmax(efoldNum,asparams(1),asparams(2)) &
            - vhi_xendmin(asparams(1),asparams(2)))

    case ('dsi')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2) * dsi_mumax(efoldNum,asparams(1),qDSI)
       asparams(3) = dsi_xendmin(efoldNum,asparams(1),asparams(2)) &
            + newparams(3) * ( dsi_xendmax(asparams(1),asparams(2),qDSI) &
            - dsi_xendmin(efoldNum,asparams(1),asparams(2)) )

    case ('dsio', 'dsi2')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2)
       asparams(3) = dsi_xendmin(efoldNum,asparams(1),asparams(2)) &
            + newparams(3) * ( dsi_xendmax(asparams(1),asparams(2),qDSI) &
            - dsi_xendmin(efoldNum,asparams(1),asparams(2)) )
            
    end select

  end function redefine_aspic_params


end module aspicpriors
