module rbflike
  use rbfprec
  use rbfnd, only : rbf_polyharmonic_two
  implicit none


  private
   
 
  integer(ip), save :: ndim, nctrs
  real(fp), save :: scale
  real(fp), save, dimension(:,:), pointer :: xctrs => null()
  real(fp), save, dimension(:), pointer :: weights => null()
  real(fp), save, dimension(:), pointer :: xpmin => null(), xpmax=>null()

  logical, parameter :: display = .false.

  public initialize_rbf_like
  public rbflike_eval, uncubize_rbfparams, cubize_rbfparams
  public posterior_boundaries, cubize_rbfparamspace
  public check_rbf, get_rbf_ndim, get_rbf_nctrs
  public cutmin_rbfparams, cutmax_rbfparams
  public get_rbf_xpmin, get_rbf_xpmax

contains


  function check_rbf()
    implicit none
    logical :: check_rbf
    
    check_rbf = associated(xctrs) .and. associated(weights) &
         .and. associated(xpmin) .and. associated(xpmax)

  end function check_rbf


  function get_rbf_ndim()
    implicit none
    integer :: get_rbf_ndim

    if (.not.check_rbf()) then
       get_rbf_ndim = 0
       return
    endif

    get_rbf_ndim = ndim

  end function get_rbf_ndim


  function get_rbf_nctrs()
    implicit none
    integer :: get_rbf_nctrs

    if (.not.check_rbf()) then
       get_rbf_nctrs = 0
       return
    endif

    get_rbf_nctrs = nctrs

  end function get_rbf_nctrs


  function get_rbf_xpmin(i)
    implicit none
    integer, intent(in) :: i
    real(fp) :: get_rbf_xpmin

    if (.not.associated(xpmin)) then
       stop 'get_rbf_xpmin: not allocated!'
    elseif ((i.gt.ndim).or.(i.le.0)) then
       stop 'get_rbf_xpmin: i out of range!'
    endif

    get_rbf_xpmin = xpmin(i)

  end function get_rbf_xpmin


  function get_rbf_xpmax(i)
    implicit none
    integer, intent(in) :: i
    real(fp) :: get_rbf_xpmax

    if (.not.associated(xpmax)) then
       stop 'get_rbf_xpmax: not allocated!'
    elseif ((i.gt.ndim).or.(i.le.0)) then
       stop 'get_rbf_xpmax: i out of range!'
    endif

    get_rbf_xpmax = xpmax(i)

  end function get_rbf_xpmax



  subroutine initialize_rbf_like(fileweights, filecentres, filebounds)
    use iorbf, only : load_weights, load_centres, read_boundaries
    implicit none   
    character(len=*), intent(in) :: fileweights, filecentres
    character(len=*), intent(in), optional :: filebounds
    
    call load_weights(fileweights,scale,weights)
    call load_centres(filecentres,xctrs)
    if (present(filebounds)) call read_boundaries(filebounds,xpmin,xpmax)
    if (size(weights,1).ne.size(xctrs,2)) stop 'weights/centres mismatch!'

    ndim = size(xctrs,1)
    nctrs = size(weights,1)

  end subroutine initialize_rbf_like



  function rbflike_eval(x)
    use rbfnd, only : rbf_svd_eval
    implicit none
    real(fp) :: rbflike_eval
    real(fp), dimension(:), intent(in) :: x
    if (size(x,1).ne.ndim) stop 'rbflike_eval: x dim screwed!'


    rbflike_eval = rbf_svd_eval(ndim,nctrs,scale,rbf_polyharmonic_two,xctrs,weights,x)

  end function rbflike_eval




  subroutine posterior_boundaries(xdata,xmin,xmax)
    implicit none
    real(fp), dimension(:,:), intent(in) :: xdata
    real(fp), dimension(:), intent(out) :: xmin,xmax

    integer(ip) :: ndim, i
    
    ndim = size(xdata,1)
    
    if ((size(xmin,1).ne.ndim).or.(size(xmax,1).ne.ndim)) then
       stop 'posterior_boundaries: array mismatch!'
    endif

    do i=1,ndim
       xmin(i) = minval(xdata(i,:))
       xmax(i) = maxval(xdata(i,:))
    enddo

  end subroutine posterior_boundaries


  function cutmin_rbfparams(ndim,icut,uncubed)
    implicit none
    integer(ip), intent(in) :: ndim,icut
    real(fp), dimension(ndim) :: cutmin_rbfparams
    real(fp), dimension(ndim), intent(in) :: uncubed

    integer(ip) :: i

    if (icut.gt.ndim) stop 'cutmin_rbfparams: icut > ndim!'

    cutmin_rbfparams = uncubed

    if (uncubed(icut).lt.xpmin(icut)) then
       cutmin_rbfparams(icut) = xpmin(icut)
       if (display) then
          write(*,*)'cutmin_rbfparams: ',icut,uncubed(icut),xpmin(icut)
       endif
    endif    

  end function cutmin_rbfparams



  function cutmax_rbfparams(ndim,icut,uncubed)
    implicit none
    integer(ip), intent(in) :: ndim,icut
    real(fp), dimension(ndim) :: cutmax_rbfparams
    real(fp), dimension(ndim), intent(in) :: uncubed

    integer(ip) :: i

    if (icut.gt.ndim) stop 'cutmax_rbfparams: icut > ndim!'

    cutmax_rbfparams = uncubed

    if (uncubed(icut).gt.xpmax(icut)) then
       cutmax_rbfparams(icut) = xpmax(icut)
       if (display) then
          write(*,*)'cutmax_rbfparams: ',icut,uncubed(icut),xpmax(icut)
       endif
    endif    

  end function cutmax_rbfparams




  function cubize_rbfparams(ndim,uncubed)
    implicit none
    integer(ip), intent(in) :: ndim
    real(fp), dimension(ndim) :: cubize_rbfparams
    real(fp), dimension(ndim), intent(in) :: uncubed

    if (ndim.ne.size(uncubed,1)) then
       stop 'uncubize_rbfparams: sizes do not match!'
    endif

    cubize_rbfparams = (uncubed - xpmin)/(xpmax-xpmin)

  end function cubize_rbfparams


  function uncubize_rbfparams(ndim,cubed)
    implicit none
    integer(ip), intent(in) :: ndim
    real(fp), dimension(ndim) :: uncubize_rbfparams
    real(fp), dimension(ndim), intent(in) :: cubed

    if (ndim.ne.size(cubed,1)) then
       stop 'uncubize_rbfparams: sizes do not match!'
    endif

    uncubize_rbfparams = xpmin + (xpmax-xpmin)*cubed

  end function uncubize_rbfparams


 
  subroutine cubize_rbfparamspace(xdata,xcubes)
    implicit none
    real(fp), dimension(:,:), intent(in) :: xdata
    real(fp), dimension(:,:), intent(out) :: xcubes
    
    integer(ip) :: ndim, ndata, i
    real(fp) :: xmin,xmax

    ndim = size(xdata,1)
    ndata = size(xdata,2)

    if ((size(xcubes,1).ne.ndim).or.(size(xcubes,2).ne.ndata)) then
       stop 'cube_rbfparamspace: mismatch arrays!'
    endif

    do i=1,ndim
       xmin = minval(xdata(i,:))
       xmax = maxval(xdata(i,:))
       xcubes(i,:) = (xdata(i,:) - xmin)/(xmax-xmin)
    enddo


  end subroutine cubize_rbfparamspace



end module rbflike
