module aspicpriors
  use aspicvars, only : kp
  implicit none

  private

  public get_aspic_priors

contains

  subroutine get_aspic_priors(name,pmin,pmax)
    character(len=*), intent(in) :: name
    real(kp), dimension(:), intent(out) :: pmin, pmax

    select case (name)

!reminder: cpp is only used in "traditional mode" for fortran
#define ONEPRIOR(fooi,foomin,foomax) \
    case ('fooi'); \
       pmin(1) = foomin ; \
       pmax(1) = foomax

#define TWOPRIOR(fooi,foo1min,foo1max,foo2min,foo2max) \
    case ('fooi'); \
       pmin(1) = foo1min ; \
       pmax(1) = foo1max ; \
       pmin(2) = foo2min ; \
       pmax(2) = foo2max

#include "aspicpriors.pp"
#undef ONEPRIOR
#undef TWOPRIOR

    case default
       stop 'get_aspic_priors: model not found!'

    end select

  end subroutine get_aspic_priors



end module aspicpriors
