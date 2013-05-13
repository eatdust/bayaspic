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
#undef ONEPRIOR
#undef TWOPRIORS
#undef THREEPRIORS

    case default
       stop 'get_aspic_priors: model not found!'

    end select

  end subroutine get_aspic_priors


!case by case model if the prior space is not cubic
  function check_aspic_hardprior(name,params)

    include 'aspicpriors.h'

    implicit none
    logical :: check_aspic_hardprior
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(in) :: params
    logical :: reject

    real(kp), parameter :: efoldMax = 120._kp

    reject = .false.

    select case (name)
       
    case ('vhi')
       reject = (params(3).gt.vhi_xendmax(efoldMax,params(1),params(2))) &
            .or.(params(3).lt.vhi_xendmin(params(1),params(2)))
       
    end select

    check_aspic_hardprior = reject

  end function check_aspic_hardprior


end module aspicpriors
