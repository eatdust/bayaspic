module rbflike
  use rbfprec
  use rbfnd, only : rbf_polyharmonic_two
  implicit none

  integer, parameter :: imn = kind(4)
  integer, parameter :: fmn = kind(1._8)

  private
   
 
  integer(ip), save :: ndim, nctrs
  real(fp), save :: scale
  real(fp), save, dimension(:,:), pointer :: xctrs => null()
  real(fp), save, dimension(:), pointer :: weights => null()

  real(fp), save, dimension(:), pointer :: xpmin => null(), xpmax=>null()


  public initialize_rbf_like, rbf_multinest_loglike
  public rbflike_eval
  public posterior_boundaries, cubize_paramspace


contains


  function check_rbf()
    implicit none
    logical :: check_rbf
    
    check_rbf = associated(xctrs) .and. associated(weights) &
         .and. associated(xpmin) .and. associated(xpmax)

  end function check_rbf

  
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


  subroutine rbf_multinest_loglike(cube,ndimmn,npars,lnew)
    use rbfnd, only : rbf_svd_eval
    implicit none   
    integer(imn), intent(in) :: ndimmn, npars
    real(fmn), dimension(ndim+npars), intent(in) :: cube
    real(fmn) :: lnew

    integer(ip) :: ndimrbf
    real(fp), dimension(ndim) :: xcuberbf

    if (.not.check_rbf()) stop 'rbf_multinest_loglike: not initialized!'
    if (ndimmn.ne.ndim) stop 'rbf_multinest_loglike: dim mismatch!'

    ndimrbf = ndimmn

    xcuberbf(1:ndimrbf) = cube(1:ndimmn)
    

    lnew = rbf_svd_eval(ndimrbf,nctrs,scale,rbf_polyharmonic_two,xctrs,weights,xcuberbf)

  end subroutine rbf_multinest_loglike



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


  function cubize_params(ndim,pmin,pmax,uncubed)
    implicit none
    integer(ip), intent(in) :: ndim
    real(fp), dimension(ndim) :: cubize_params
    real(fp), dimension(ndim), intent(in) :: pmin,pmax,uncubed

    cubize_params = (uncubed - pmin)/(pmax-pmin)

  end function cubize_params


  function uncubize_params(ndim,pmin,pmax,cubed)
    implicit none
    integer(ip), intent(in) :: ndim
    real(fp), dimension(ndim) :: uncubize_params
    real(fp), dimension(ndim), intent(in) :: pmin,pmax,cubed

    uncubize_params = pmin + (pmax-pmin)*cubed

  end function uncubize_params


 
  subroutine cubize_paramspace(xdata,xcubes)
    implicit none
    real(fp), dimension(:,:), intent(in) :: xdata
    real(fp), dimension(:,:), intent(out) :: xcubes
    
    integer(ip) :: ndim, ndata, i
    real(fp) :: xmin,xmax

    ndim = size(xdata,1)
    ndata = size(xdata,2)

    if ((size(xcubes,1).ne.ndim).or.(size(xcubes,2).ne.ndata)) then
       stop 'cube_paramspace: mismatch arrays!'
    endif

    do i=1,ndim
       xmin = minval(xdata(i,:))
       xmax = maxval(xdata(i,:))
       xcubes(i,:) = (xdata(i,:) - xmin)/(xmax-xmin)
    enddo


  end subroutine cubize_paramspace



end module rbflike
