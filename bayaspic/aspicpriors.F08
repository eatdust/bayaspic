module aspicpriors
  use aspicvars, only : kp, pi, lname
  implicit none

  private

  public get_aspic_priors, check_aspic_hardprior
  public redefine_aspic_params

!hard prior: all models are assumed to support, at least, this number
!of efolds which allows for all possible reheating history
  real(kp), parameter :: efoldNum = 120._kp

!but some models are so fine-tuned that they cannot even support 120
!efolds of inflation and would be ruled out otherwise. In that case,
!we ask for the minimum possible efold number
  real(kp), parameter :: efoldMini = 60._kp

!hard prior: all models inflating in a region with eps2>eps2minMax are not counted
  real(kp), parameter :: epsilonMax = 0.2_kp

!dsi correction power
  real(kp), parameter :: qDSI = 8._kp


contains

  subroutine get_aspic_priors(name,pmin,pmax,mapname)
    use aspicmodels, only : get_aspic_numparams
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(out) :: pmin, pmax
    character(len=lname), dimension(:), intent(inout) :: mapname

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

#define FOURPRIORS(fooi,foo1min,foo1max,foo1map,foo2min,foo2max,foo2map \
       ,foo3min,foo3max,foo3map,foo4min,foo4max,foo4map) \
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
       pmin(4) = foo4min ; \
       pmax(4) = foo4max ; \
       mapname(4) = 'foo4map' ; \
       if (naspic.ne.4) stop 'get_aspic_priors: wrong npriors'


#include "aspicpriors.pp"
#undef ZEROPRIOR
#undef ONEPRIOR
#undef TWOPRIORS
#undef THREEPRIORS
#undef FOURPRIORS

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

    reject = .false.

    select case (name)

    case ('ii','iif','iilambda')
       reject = params(2).lt.ii_xendmin(efoldNum,params(1))

    case ('lmi1o', 'lmi1p')
       reject = lmi_epstwo_potmax(params(1),params(2)).gt.epsilonMax      

    case ('lmi2o','lmi2p')
!Such that exp(-beta xend^gamma) is computable; and epsilon1>numacc;
!and epsilon1< 0.2 during inflation
       reject = (params(3).gt.lmi_numacc_x_big(params(1),params(2))) &
            .or. (params(2).lt.lmi2_numacc_betamin(params(1))) &
            .or. (lmi2_epsilon_one(params(3),params(1),params(2)) &
            .gt.epsilonMax)

    case('twiA1','twiA2','twiB1','twiB2')
       reject = (params(2).gt.1._kp)
       
    case ('bsusybif', 'bsusybil')
       reject = params(2).gt.bsusybi_xendmax(efoldNum,params(1))

    case ('gmssmi')
       reject = (abs(gmssmi_epstwomin(params(1),params(2))).gt.epsilonMax) &
            .or. (params(1).lt.gmssmi_alphamin(efoldNum,params(2)) &
            .or. (abs(params(1)-1._kp).lt.epsilon(1._kp)))

    case ('gmssmiopA','gmssmiopB','gmssmiomA','gmssmiomB')
       reject = params(1).lt.gmssmi_alphamin(efoldMini,params(2))

    case ('gripi')
       reject = (abs(gripi_epstwomin(params(1),params(2))).gt.epsilonMax) &
            .or. (params(1).lt.gripi_alphamin(efoldNum,params(2)) &
            .or. (abs(params(1)-1._kp).lt.epsilon(1._kp)))

    case ('gripiopA','gripiopB','gripiomA','gripiomB')
       reject = params(1).lt.gripi_alphamin(efoldNum,params(2))

    case ('lin')
       reject = params(1).lt.li_alphamin(efoldNum)

    case ('li')
       reject = params(1) .gt. -0.1_kp .and. params(1) .lt. 0.003_kp 

    case ('tiB1', 'tiD1', 'tiF1')
       reject = abs(1._kp-2._kp*params(1)).gt.0.1_kp*params(2)**2

    case ('tiB2')
       reject = abs(1._kp-2._kp*params(1)).gt.0.2_kp*params(2)**2

    case ('tiB3', 'tiD2', 'tiF2')
       reject = abs(1._kp-2._kp*params(1)).gt.1._kp*params(2)**2

    case ('tiD3', 'tiF3')
       reject = abs(1._kp-2._kp*params(1)).gt.10._kp*params(2)**2

    case ('psnioA', 'psnioB', 'psnioC')
       reject = (params(1).gt.0.5_kp*epsilonMax*params(2)**2)

    case ('sbi')
       reject = params(1).lt.(sbi_alphamin(params(2)) +epsilon(1._kp))

    case ('ssbi1', 'ssbi1f')
!otherwise inflation does not stop
      reject = ssbi1_epsilon_one(ssbi1_x_epsonemax(params(1),params(2)), & 
               params(1),params(2)).lt.1._kp

    case ('ssbi3', 'ssbi3f')
      reject = (ssbi3_epsilon_one(ssbi3_x_epsonemax(params(1),params(2)), & 
               params(1),params(2)).lt.1._kp) .or. &
               (ssbi3_epstwo_potmax(params(1),params(2)) .gt. epsilonMax)

    case ('ssbi4', 'ssbi4f')
      reject = ssbi4_epstwo_potmax(params(1),params(2)) .gt. epsilonMax

    case ('ssbi5','ssbi5f')
       reject = params(1).gt.ssbi5_alphamax(params(2))

    case ('ssbi6','ssbi6f')
       reject = params(1).gt.ssbi6_alphamax(params(2))
              
    case ('rmi1', 'rmi1l')
       reject = params(3).gt.rmi1_numacc_xendmax(efoldNum,params(1),params(2))

    case ('rmi2', 'rmi2l')
       reject = params(3).lt.rmi2_numacc_xendmin(efoldNum,params(1),params(2))

    case ('rmi3', 'rmi3l')
       reject = (params(3).gt.rmi3_numacc_xendmax(params(1),params(2))) &
            .or. (params(3).lt.rmi3_numacc_xendmin(efoldNum,params(1),params(2)))

    case ('rmi4', 'rmi4l')
       reject = params(3).lt.rmi4_numacc_xendmin(params(1),params(2))

    case ('dsio', 'dsi2')
       reject = params(2) .gt. dsi_mumax(efoldNum,params(1),qDsi)

    case ('bistg')
!enforce the hard prior bound xini<xuv (asparams(4) is xuv)
       reject = params(4).lt.bi_x_trajectory(-efoldNum,params(3),params(1),params(2))

    case ('kkltistg')
!enforce the hard prior bound xini<xuv (asparams(4) is xuv)
       reject = params(4).lt.kklti_x_trajectory(-efoldNum,params(3),params(1),params(2))

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

    case ('sbialphamin')
       asparams(1) = sbi_alphamin(newparams(2))
       asparams(2) = newparams(2)

    case ('rpi2')
       asparams(2) = newparams(1)
       asparams(2) = newparams(2)*rpi_x_potmax(asparams(1))

    case ('kmiii')
       asparams(1) = newparams(1)*newparams(2)**(5._kp/3._kp)
       asparams(2) = newparams(2)**(2._kp/3._kp)

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

    case ('twiA1','twiA2')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2)*newparams(1)

    case ('twiB1','twiB2')
       asparams(1) = 10._kp**(-5._kp)/sqrt(newparams(1))
       asparams(2) = newparams(2)*asparams(1)
 
    case ('csi')
       asparams(1) = newparams(1)
       asparams(2) = csi_x_epsonenumacc(asparams(1)) &
            + newparams(2) * (csi_xendmax(efoldNum,asparams(1)) &
            - csi_x_epsonenumacc(asparams(1)))
  
    case ('cnci')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2) * cnci_xendmin(efoldNum,asparams(1))

    case ('imi', 'imi1', 'imi2', 'imi3', 'imi4', 'imi5', 'imi6')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2) * imi_xendmin(efoldNum,asparams(1))

    case ('gmssmiopA','gmssmiopB')
       asparams(1) = 1._kp + newparams(1)
       asparams(2) = newparams(2)    

    case ('gmssmiomA','gmssmiomB')
       asparams(1) = 1._kp - newparams(1)
       asparams(2) = newparams(2)

    case ('gripiopA','gripiopB')
       asparams(1) = 1._kp + newparams(1)
       asparams(2) = newparams(2)      

    case ('gripiomA','gripiomB')
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

    case ('psnift1', 'psnift2', 'psnift3', 'psnift4')
       asparams(1) = newparams(1)*newparams(2)**2
       asparams(2) = newparams(2)

    case ('rmi4', 'rmi4l')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2)
       asparams(3) = 1._kp + newparams(3) &
            * (rmi4_xendmax(efoldNum,asparams(1),asparams(2)) - 1._kp)

    case ('vhi','vhi12','vhi1','vhi2','vhi3','vhi4','vhil1')
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

    case ('dsio','dsi2')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2)
       asparams(3) = dsi_xendmin(efoldNum,asparams(1),asparams(2)) &
            + newparams(3) * ( dsi_xendmax(asparams(1),asparams(2),qDSI) &
            - dsi_xendmin(efoldNum,asparams(1),asparams(2)) )
            
    case ('bi','bis','bi1s','bi2s','bi3s','bi4s','bi5s','bi6s')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2)
!overload xend with xeps1 for the pheno bi models
       asparams(3) = bi_x_epsoneunity(asparams(1),asparams(2))
!unused
       asparams(4) = 0._kp

    case ('kklti','kkltis','kklti1s','kklti2s','kklti3s','kklti4s','kklti5s','kklti6s')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2)
!overload xend with xeps1 for the pheno kklti models
       asparams(3) = kklti_x_epsoneunity(asparams(1),asparams(2))
!unused
       asparams(4) = 0._kp
        
    case ('bistg')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2)
!overload xend with max(xeps1,xstg)
       asparams(3) = max(bi_x_epsoneunity(asparams(1),asparams(2)),newparams(3))
!multinest param is phiuv, we map it to xuv = phiuv/mu
       asparams(4) = newparams(4)/newparams(2)

    case ('kkltistg')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2)
!overload xend with max(xeps1,xstg)
       asparams(3) = max(kklti_x_epsoneunity(asparams(1),asparams(2)),newparams(3))
!multinest param is phiuv, we map it to xuv = phiuv/mu
       asparams(4) = newparams(4)/newparams(2)

    end select

  end function redefine_aspic_params


end module aspicpriors
