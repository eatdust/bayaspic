!   This file is part of bayaspic
!
!   Copyright (C) 2013-2023 C. Ringeval
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
!   along with bayaspic.  If not, see <https://www.gnu.org/licenses/>.


module aspicpriors
  use aspicvars, only : kp, pi, lname
  use aspicvars, only : AspicModel

  include 'aspicpriors.h'
  
  implicit none


  
  private

  public get_aspic_priors, check_aspic_hardprior
  public get_aspic_numpriors
  public redefine_aspic_params

  integer, save :: numAspicPriors = 0


!hard prior: all models are assumed to support, at least, this number
!of efolds which allows for all possible reheating history
  real(kp), parameter :: efoldNum = 120._kp

!but some models are so fine-tuned that they cannot even support 120
!efolds of inflation and would be ruled out otherwise. In that case,
!we ask for the minimum possible efold number
  real(kp), parameter :: efoldMini = 60._kp

!conversely, in some models the potential runs over so many orders of magnitude that a larger number of efolds need to be realised
  real(kp), parameter :: efoldMaxi = 180._kp

!hard prior: all the parameter space inflating in a domain with
!eps2>eps2minMax is ignored
  real(kp), parameter :: epsilonMax = 0.2_kp

!dsi correction power
  real(kp), parameter :: qDSI = 8._kp
  
contains

  subroutine get_aspic_priors(name,pmin,pmax,mapname)
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(out) :: pmin, pmax
    character(len=lname), dimension(:), intent(inout) :: mapname


    numAspicPriors = 0

    select case (name)

#define ZEROPRIOR(fooi) \
    case ('fooi')

!reminder: cpp is only used in "traditional mode" for fortran
#define ONEPRIOR(fooi,foomin,foomax,foomap) \
    case ('fooi'); \
       pmin(1) = foomin ; \
       pmax(1) = foomax ; \
       mapname(1) = 'foomap' ; \
       numAspicPriors = 1

#define TWOPRIORS(fooi,foo1min,foo1max,foo1map,foo2min,foo2max,foo2map) \
    case ('fooi'); \
       pmin(1) = foo1min ; \
       pmax(1) = foo1max ; \
       mapname(1) = 'foo1map' ; \
       pmin(2) = foo2min ; \
       pmax(2) = foo2max ; \
       mapname(2) = 'foo2map' ; \
       numAspicPriors = 2

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
       numAspicPriors = 3

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
       numAspicPriors = 4


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


  function get_aspic_numpriors()
    implicit none
    integer :: get_aspic_numpriors

    get_aspic_numpriors = numAspicPriors

  end function get_aspic_numpriors




!case by case model if the prior space is not cubic
  function check_aspic_hardprior(name,params) result(answer)
    implicit none
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params
    logical :: ignore, disfavour
    logical, dimension(2) :: answer

!ignore is for hard prior, these regions are simply ignored
    ignore = .false.

!disfavour is for ruling out, these regions will be associated to
!logZero (this is changing the evidence)
    disfavour = .false.

    select case (name)

!this is SI assuming a hard prior as a minimal coupling to gravity    
    case ('si mc')
       ignore = .false.
       disfavour = .false.
       call si_set_minimal_coupling(.true.)
       
    case ('rcipi1tune2p','rcipi1tune2m','rcipi1tune4p','rcipi1tune4m','rcipi1m')
       ignore = .not.rcipi_check_params(params(1),params(2),params(3))
       disfavour = .false.
       if (.not.ignore) then
          ignore = rcipi_efoldmax(params(1),params(2),params(3)).lt.efoldNum
          disfavour = ignore
       endif
       
    case ('rclfi1pm','rclfi1mm')
       ignore = .not.rclfi1_check_params(params(1),params(2),params(3))
       disfavour = .false.

    case ('rclfi1mp')
       ignore = .not.rclfi1_check_params(params(1),params(2),params(3))
       if (.not.ignore) ignore = params(2).gt.rclfi1_numacc_alphamax(params(1))
       disfavour = .false.

    case ('rclfi2pm','rclfi2mm')
       ignore = .not.rclfi2_check_params(params(1),params(2),params(3))
       disfavour = .false.

    case ('rclfi2mp')
       ignore = .not.rclfi2_check_params(params(1),params(2),params(3))
       if (.not.ignore) ignore = params(2).gt.rclfi2_numacc_alphamax(params(1))
       disfavour = .false.

    case ('rclfi3pp','rclfi3mp')
       ignore = .not.rclfi3_check_params(params(1),params(2),params(3))
       disfavour = .false.

    case ('rclfi3pm')
       ignore = .not.rclfi3_check_params(params(1),params(2),params(3))
       if (.not.ignore) ignore = (params(1).lt.rclfi3_numacc_pmin()) .or. &
            (params(2).lt.rclfi3_numacc_alphamin(params(1)))
       disfavour = .false.

    case ('rclfi4p','rclfi4m')
       ignore = .not.rclfi4_check_params(params(1),params(2),params(3))
       disfavour = .false.
       
    case ('nmlfi1s','nmlfi1l','nmlfi11','nmlfi12','nmlfi13')
       ignore = .not.nmlfi1_check_params(params(1),params(2))
       if (.not.ignore) ignore = nmlfi1_numacc_efoldmax(params(1),params(2)).lt.efoldNum
       if (.not.ignore) ignore = nmlfi_epsilon_one_infinity(params(1),params(2)).gt.epsilonMax
       disfavour = .false.

    case ('nmlfi3s','nmlfi3l','nmlfi31','nmlfi32','nmlfi33')
       ignore = .not.nmlfi3_check_params(params(1),params(2))
       if (.not.ignore) ignore = params(3).lt.nmlfi3_numacc_hbarendmin(efoldNum,params(1),params(2))
       disfavour = .false.
       
    case ('saiii1p','saiii1n')
       ignore = .not.saiii1_check_params(params(1),params(2),params(3))
       if (.not.ignore) then
          ignore = saiii1_numacc_efoldmax(params(1),params(2),params(3)).lt.efoldNum
       endif
       disfavour = .false.

    case ('saiii2p','saiii2n')

       ignore = .not.saiii2_check_params(params(1),params(2),params(3))

       if (.not.ignore) then
          ignore = saiii2_numacc_efoldmax(params(1),params(2),params(3)).lt.efoldNum
       endif
       disfavour = .false.
       
       
    case ('saiii3p')
       ignore = .not.saiii3_check_params(params(1),params(2),params(3))

!not too close to the model boundary (trajectory integration
!takes ages)
       ignore = ignore .or. params(1).lt.1.1_kp*saiii_alpha_two(params(2))
       disfavour = .false.

    case ('saiii3n')
       ignore = .not.saiii3_check_params(params(1),params(2),params(3))

!not too close to the model boundary (trajectory integration
!takes ages)
       ignore = ignore .or. params(1).gt.1.1_kp*saiii_alpha_one(params(2))
       if ((.not.ignore).and.(params(2).gt.beta1)) then
          ignore = params(1).lt.0.98_kp*saiii_alpha_three(params(2))
       endif
       disfavour = .false.

       
    case('saii1n', 'saii1p', 'saii1f')
       ignore = params(2).lt.saii1_numacc_mumin(efoldNum,params(1))
       disfavour = .false.

    case('saii2n', 'saii2p', 'saii2f')
       ignore = params(2).lt.saii1_numacc_mumin(efoldNum,params(1))
       disfavour = .false.

    case('sdi')
       ignore = params(2).lt.sdi_numacc_xendmin(efoldNum,params(1)) .or. &
            params(2).gt.sdi_numacc_xendmax(params(1))
       disfavour = .false.

    case('hni1f')
       ignore = params(1).le.hni1_alphamin(params(2))
       if (.not.ignore) ignore = hni1_numacc_efoldmax(params(1),params(2)).le.efoldNum
       disfavour = .false.

    case('hni2f','hni2l')
       ignore = params(3).ge.hni2_xendmax(params(1),params(2))
       if (.not.ignore) ignore = hni2_numacc_efoldmax(params(1),params(2)).le.efoldNum
       endif
       disfavour = .false.

    case('ncli1','ncli2','ncli3','ncli4', 'ncli')
       disfavour = .not.(ncli_check_params(efoldMini,epsilonMax, &
            params(1),params(2),params(3)))
       ignore = disfavour
       
    case('fi0')
       disfavour = .not.(fi_check_params(100._kp,epsilonMax,params(1),params(2)))
       ignore = disfavour

    case('fi1')
       disfavour = .not.(fi_check_params(80._kp,epsilonMax,params(1),params(2)))
       ignore = disfavour

    case('fi2')
       disfavour = .not.(fi_check_params(80._kp,epsilonMax,params(1),params(2)))
       ignore = disfavour

    case('sbkin','sbkip')
      disfavour = (sbki_efoldmax(params(1)).lt.efoldNum) &
              .or. (sbki_epsilon_one_min(params(1)).gt.efoldMini)
      ignore = disfavour

    case ('ii','iif','iilambda')
       disfavour = params(2).lt.ii_xendmin(efoldNum,params(1))
       ignore = disfavour

    case ('lmi1o', 'lmi1p')
       disfavour = lmi_epstwo_potmax(params(1),params(2)).gt.epsilonMax
       ignore = disfavour

    case ('lmi2o','lmi2p')
!Such that exp(-beta xend^gamma) is computable; and epsilon1>numacc;
!and epsilon1< 0.2 during inflation
       ignore = (params(3).gt.lmi_numacc_x_big(params(1),params(2))) &
            .or. (params(2).lt.lmi2_numacc_betamin(params(1)))

!eps1>0.2 all along inflation is ruled out
       if (lmi2_epsilon_one(params(3),params(1),params(2)) &
            .gt.epsilonMax) then
          ignore = .true.
          disfavour = .true.
       endif

    case('twiA1','twiA2','twiB1','twiB2')
       ignore = (params(2).gt.1._kp)
       disfavour = .false.
       
    case ('bsusybif', 'bsusybil')
       disfavour = params(2).gt.bsusybi_xendmax(efoldNum,params(1))
       ignore = disfavour
!this is roughly -95
       ignore = ignore .or. (params(2).lt.bsusybi_numacc_xendmin(1d-99,params(1)))

    case ('gmssmi')
       ignore = (params(1).lt.gmssmi_alphamin(efoldNum,params(2)) &
            .or. (abs(params(1)-1._kp).lt.epsilon(1._kp)))
       if (abs(gmssmi_epstwomin(params(1),params(2))).gt.epsilonMax) then
          ignore = .true.
          disfavour = .true.
       endif

    case ('gmssmiopA','gmssmiopB','gmssmiomA','gmssmiomB')
       disfavour = params(1).lt.gmssmi_alphamin(efoldMini,params(2))
       ignore = disfavour
       
    case ('gripi')
       ignore = (params(1).lt.gripi_alphamin(efoldNum,params(2)) &
            .or. (abs(params(1)-1._kp).lt.epsilon(1._kp)))
       if (abs(gmssmi_epstwomin(params(1),params(2))).gt.epsilonMax) then
          ignore = .true.
          disfavour = .true.
       endif

    case ('gripiopA','gripiopB','gripiomA','gripiomB')
       disfavour = params(1).lt.gripi_alphamin(efoldNum,params(2))
       ignore = disfavour
       
    case ('lin')
       disfavour = params(1).lt.li_alphamin(efoldNum)
       ignore = disfavour
       
    case ('li')
       ignore = params(1) .gt. -0.1_kp .and. params(1) .lt. 0.003_kp
       disfavour = .false.
       
    case ('tiB1', 'tiD1', 'tiF1')
       ignore = abs(1._kp-2._kp*params(1)).gt.0.1_kp*params(2)**2
       disfavour = .false.
       
    case ('tiB2')
       ignore = abs(1._kp-2._kp*params(1)).gt.0.2_kp*params(2)**2
       disfavour = .false.
       
    case ('tiB3', 'tiD2', 'tiF2')
       ignore = abs(1._kp-2._kp*params(1)).gt.1._kp*params(2)**2
       disfavour = .false.
       
    case ('tiD3', 'tiF3')
       ignore = abs(1._kp-2._kp*params(1)).gt.10._kp*params(2)**2
       disfavour = .false.
       
    case ('psnioA', 'psnioB', 'psnioC')
       ignore = (params(1).gt.0.5_kp*epsilonMax*params(2)**2)
       disfavour = .false.
       
    case ('sbi')
       ignore = params(1).lt.(sbi_alphamin(params(2)) +epsilon(1._kp))
       disfavour = .false.
       
    case ('ssbi1', 'ssbi1f')
!otherwise inflation does not stop
      ignore = ssbi1_epsilon_one(ssbi1_x_epsonemax(params(1),params(2)), &
               params(1),params(2)).lt.1._kp
      disfavour = .false.
      
    case ('ssbi3', 'ssbi3f')
      ignore = (ssbi3_epsilon_one(ssbi3_x_epsonemax(params(1),params(2)), &
           params(1),params(2)).lt.1._kp)

      if (ssbi3_epstwo_potmax(params(1),params(2)) .gt. epsilonMax) then
         ignore = .true.
         disfavour = .true.
      endif

    case ('ssbi4', 'ssbi4f')
      disfavour = ssbi4_epstwo_potmax(params(1),params(2)) .gt. epsilonMax
      ignore = disfavour

    case ('ssbi5','ssbi5f')
       ignore = params(1).gt.ssbi5_alphamax(params(2))
       disfavour = .false.
       
    case ('ssbi6','ssbi6f')
       ignore = params(1).gt.ssbi6_alphamax(params(2))
       disfavour = .false.
       
    case ('rmi1', 'rmi1l')
       ignore = params(3).gt.rmi1_numacc_xendmax(efoldNum,params(1),params(2))
       disfavour = .false.
       
    case ('rmi2', 'rmi2l')
       ignore = params(3).lt.rmi2_numacc_xendmin(efoldNum,params(1),params(2))
       disfavour = .false.
       
    case ('rmi3', 'rmi3l')
       ignore = (params(3).gt.rmi3_numacc_xendmax(params(1),params(2))) &
            .or. (params(3).lt.rmi3_numacc_xendmin(efoldNum,params(1),params(2)))
       disfavour = .false.
       
    case ('rmi4', 'rmi4l')
       ignore = params(3).lt.rmi4_numacc_xendmin(params(1),params(2))
       disfavour = .false.
       
    case ('rpi2')
       ignore = params(2).lt.rpi2_numacc_xendmin(efoldNum,params(1))
       disfavour = .false.
       
!equivalent to an hard prior on mumin to get efoldnum
    case ('vhi12')
       ignore = vhi_xendmax(efoldNum,params(1),params(2)) &
            .eq.vhi_xendmin(params(1),params(2))
       disfavour = .false.
       
    case ('dsio', 'dsi2')
       ignore = params(2) .gt. dsi_mumax(efoldNum,params(1),qDsi)
       disfavour = .false.
       
    case ('bistg')
!enforce the hard prior bound xini<xuv (asparams(4) is xuv)
       ignore = params(4).lt.bi_x_trajectory(-efoldNum,params(3),params(1),params(2))
       disfavour = .false.
       
    case ('kkltistg')
!enforce the hard prior bound xini<xuv (asparams(4) is xuv)
       ignore = params(4).lt.kklti_x_trajectory(-efoldNum,params(3),params(1),params(2))
       disfavour = .false.
       
    case('nfi1')
       ignore = params(1).lt.nfi1_numacc_amin(params(2)) .or. &
            (params(2).lt.2._kp .and. &
            params(1).gt.nfi1_amax(efoldNum,params(2))) .or. &
            params(1).gt.nfi1_numacc_amax(efoldNum,params(2))
       disfavour = params(1).lt.nfi1_numacc_amin(params(2)) .or. &
            params(2).lt.2._kp .and. &
            params(1).gt.nfi1_amax(efoldNum,params(2))

    case ('nfi2')
       ignore = (params(1).lt.nfi2_numacc_amin(efoldMaxi,params(2))) .or. &
            nfi2_efold_max(params(1),params(2)) .lt. efoldMaxi
       disfavour = (params(1).lt.nfi2_numacc_amin(efoldMaxi,params(2))) .or. &
            nfi2_efold_max(params(1),params(2)) .lt. efoldMaxi
       if (nfi2_epsilon_one(params(3),params(1),params(2)).gt.epsilonMax) then
          ignore = .true.
          disfavour = .true.
       endif

    case ('nfi3n','nfi3p')
       ignore = abs(params(1)).gt.0.1_kp*nfi3_numacc_absamax(params(2))
       disfavour = .false.

    case ('nfi4n', 'nfi4p')
       ignore = nfi4_xendmin(efoldMaxi,params(1),params(2)) &
            .ge.nfi4_numacc_xendmax(params(1),params(2))
       disfavour = .false.

    case ('vfmi', 'vfmis')
!Such that exp{-3 beta exp[-x/sqrt(3 beta)]} is computable, the
!potential should not become too large.
       ignore = params(2).ge.vfmi_numacc_betamax(efoldNum,params(1))
       disfavour = .false.

    case ('hbio','hbip')
       if (hbi_epsonemin(params(1),params(2)).gt.epsilonMax) then
          ignore = .true.
          disfavour = .true.
       endif

    case ('shi')
       if (shi_epstwomin(params(1),params(2)).gt.epsilonMax) then
          ignore = .true.
          disfavour = .true.
       endif

    case ('ccsi1')
       ignore = ccsi1_numacc_efoldmax(params(1)).lt.efoldNum
       disfavour = .false.
       
    end select

!Currently, even though some models have many regions
!trivially ruled out, we just ignore them and give these model a
!chance... (see wraspic & samplaspic modules)  
!let's never disfavour
    disfavour = .false.

    answer = (/ignore,disfavour/)

  end function check_aspic_hardprior



!redefinition of aspic params (called after change of priors, see wraspic)
  function redefine_aspic_params(name,newparams) result(asparams)
    implicit none
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: newparams
    real(kp), dimension(size(newparams,1)) :: asparams

    asparams = newparams

    select case (name)

    case ('sbialphamin')
       asparams(1) = sbi_alphamin(newparams(2))
       asparams(2) = newparams(2)

    case ('rpi2')
       asparams(1) = newparams(1)
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
       asparams(2) = csi_numacc_x_epsonesmall(asparams(1),1d-5) &
            + newparams(2) * ( csi_xendmax(efoldNum,asparams(1)) &
            - csi_numacc_x_epsonesmall(asparams(1),1d-5) )

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


    case ('kklti','kkltis','kklti1s','kklti2s','kklti3s','kklti4s','kklti5s','kklti6s')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2)
!overload xend with xeps1 for the pheno kklti models
       asparams(3) = kklti_x_epsoneunity(asparams(1),asparams(2))

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

    case ('nfi2')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2)
       asparams(3) = nfi2_numacc_xendmin(asparams(1),asparams(2)) &
            + newparams(3)*(nfi2_numacc_xendmax(efoldMaxi,asparams(1),asparams(2)) &
            - nfi2_numacc_xendmin(asparams(1),asparams(2)))

    case ('nfi4p','nfi4n')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2)
       asparams(3) = nfi4_xendmin(efoldMaxi,asparams(1),asparams(2)) &
            +newparams(3)*(nfi4_numacc_xendmax(asparams(1),asparams(2)) &
            -nfi4_xendmin(efoldMaxi,asparams(1),asparams(2)))

    case ('hbio')
!n and mu as a function of nphi and nf
       asparams(1) = 2._kp*(newparams(1)+newparams(2))/newparams(2)
       asparams(2) = 2._kp*sqrt(newparams(1)+newparams(2))/newparams(2)

!we sample on alpha/ccsi3_alphamin(efoldNum) where alphamin < 0
    case ('ccsi3')
       asparams(1) = newparams(1)*ccsi3_alphamin(efoldNum)

!first param is alpha, second is xend/xendmin with xendmin the closest possible to xVmax
    case ('ccsi2')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2) * ccsi2_numacc_xendmin(efoldNum,asparams(1))
       
    case ('rclfi1pm')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2) * rclfi_alpha_one(asparams(1))
!for test       AspicModel%derived(1) = rclfi1_numacc_mumin(efoldNum,asparams(1),asparams(2))
       asparams(3) = newparams(3) * rclfi1_numacc_mumin(efoldNum,asparams(1),asparams(2))

    case ('rclfi1mm')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2)
       asparams(3) = newparams(3) * rclfi1_numacc_mumin(efoldNum,asparams(1),asparams(2))

    case ('rclfi1mp')
       asparams(1) = newparams(1) * rclfi1_numacc_pmax()
       asparams(2) = newparams(2) * rclfi_alpha_one(asparams(1))
       asparams(3) = newparams(3) * rclfi1_numacc_mumin(efoldNum,asparams(1),asparams(2))       

    case ('rclfi2pm')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2) * rclfi_alpha_zero(asparams(1))
       asparams(3) = newparams(3) * rclfi2_numacc_mumin(efoldNum,asparams(1),asparams(2))

    case ('rclfi2mm')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2)
       asparams(3) = newparams(3) * rclfi2_numacc_mumin(efoldNum,asparams(1),asparams(2))
       
    case ('rclfi2mp')
       asparams(1) = newparams(1) * rclfi2_numacc_pmax()
       asparams(2) = newparams(2) * rclfi_alpha_zero(asparams(1))
       asparams(3) = newparams(3) * rclfi2_numacc_mumin(efoldNum,asparams(1),asparams(2))

    case ('rclfi3pm','rclfi3mp')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2) * rclfi_alpha_zero(asparams(1))
       asparams(3) = newparams(3)

    case ('rclfi4p','rclfi4m')
       asparams(1) = newparams(1)
       asparams(2) = newparams(2) * rclfi_alpha_one(asparams(1))
       asparams(3) = newparams(3)

    case ('rcipi1tune2p','rcipi1tune4p')
       asparams(1) = newparams(1)
       asparams(3) = newparams(3)
       asparams(2) = newparams(2) + rcipi_alpha_zero(asparams(1),asparams(3))

    case ('rcipi1tune2m','rcipi1tune4m')
       asparams(1) = newparams(1)
       asparams(3) = newparams(3)
       asparams(2) = newparams(2) - rcipi_alpha_zero(asparams(1),asparams(3))
       
    case ('rcipi1m')
       asparams(1) = newparams(1)
       asparams(3) = newparams(3)
       asparams(2) = -2._kp*sqrt(asparams(3)) + newparams(2) &
            * ( -rcipi_alpha_zero(asparams(1),asparams(3)) + 2._kp*sqrt(asparams(3)) )
       
    case ('rcipi22','rcipi23','rcipi24','rcipi25','rcipi26')
       asparams(1) = newparams(1)
       asparams(3) = newparams(3)
       asparams(2) = newparams(2) * rcipi_alpha_zero(asparams(1),asparams(3))

    case ('rcipi2')
       asparams(1) = newparams(1)
       asparams(3) = newparams(3) * sqrt(asparams(1))
       asparams(2) = newparams(2) * rcipi_alpha_zero(asparams(1),asparams(3))
       
    end select

  end function redefine_aspic_params


end module aspicpriors
