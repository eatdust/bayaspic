module iond
  use rbfprec
  implicit none

  private

  integer(ip), parameter :: nrecmax = 1000000
  integer, parameter :: ndimmax = 5

  public read_binned_posterior, posterior_boundaries
  public save_posterior, load_posterior, cubize_paramspace
  public save_boundaries, read_boundaries

  public save_weights, load_weights
  public save_centres, load_centres
  public save_shepdata, load_shepdata

contains


  subroutine read_binned_posterior(filename,posterior,params)
    implicit none
    character(len=*), intent(in) :: filename   

    logical, parameter :: logpost = .true.
    integer, parameter :: nzeroskip = 0

    real(fp), save :: eps = exp(-10._fp)
    real(fp), save :: minNonZero = huge(1._fp)


    integer, parameter :: nunit = 210

    real(fp), dimension(:), intent(inout), pointer :: posterior
    real(fp), dimension(:,:), intent(inout), pointer :: params

    real(fp), save :: fpbuffer
    character(LEN=200) :: cbuffer
    
    real(fp), dimension(ndimmax), save :: statbuffer

    integer(ip) :: i,j,ioerr,nrec,ndim,nnzero,count

    if (associated(posterior).or.associated(params)) then
       stop 'read_binned_posterior: data already loaded!'
    endif

    write(*,*)'read_binned_posterior:'         
    write(*,*)'opening ',filename
   
!counting columns
    ndim = count_column(filename,'E') - 1
    
    
!counting number records and number of non-zero posterior values
    open(unit=nunit,file=filename,status='old')
    nnzero=0
    nrec=0
    count = 0
    do i=1,nrecmax
       read(nunit,*,iostat=ioerr) statbuffer(1:ndim+1)
       if (ioerr.ne.0) exit
       nrec = nrec + 1

       if (statbuffer(1).eq.0._fp) then
          count = count + 1
          if (count.le.nzeroskip) cycle
          count = 0
       else
          minNonZero = min(minNonZero,statbuffer(1))
       endif

       nnzero = nnzero + 1
    enddo
    close(nunit)

    write(*,*)'number of params:',ndim
    write(*,*)'number or records:',nrec
    write(*,*)'number of bins kept:',nnzero
    write(*,*)'ln(minNonZero)=    ',log(minNonZero)

    eps = exp(0._fp+real(int(log(minNonZero)),fp))


!reading non-zero records
    allocate(posterior(nnzero))
    allocate(params(ndim,nnzero))
         
    open(unit=nunit,file=filename,status='old')
    j=0
    count = 0
    do i=1,nrec
       read(nunit,*,iostat=ioerr) statbuffer(1:ndim+1)

       if (statbuffer(1).eq.0._fp) then
          count = count + 1
          if (count.le.nzeroskip) cycle
          count = 0
       end if

       j=j+1
       if (logpost) then
          posterior(j) = log(statbuffer(1)+eps)
       else
          posterior(j) = statbuffer(1)
       endif
       params(1:ndim,j) = statbuffer(2:ndim+1)
!       print *,'test',posterior(j),params(1:ndim,j)
       if (ioerr.ne.0) stop 'counting screwed!'
    enddo
    close(nunit)

  end subroutine read_binned_posterior

  

  subroutine save_posterior(filename,f,x)
    implicit none
    character(len=*), intent(in) :: filename   
    real(fp), dimension(:), intent(in) :: f
    real(fp), dimension(:,:), intent(in) :: x

    integer, parameter :: nunit = 211

    integer :: i,j, ndim, ndata

    ndim = size(x,1)
    ndata = size(x,2)
    
    if (size(f,1).ne.ndata) stop 'save_posterior: size mismatch!'

    open(unit=nunit,file=filename,status='unknown')
    
    write(nunit,*) ndim, ndata
    do i=1,ndata
       write(nunit,*) f(i),(x(j,i),j=1,ndim)
    enddo
    
    close(nunit)
  end subroutine save_posterior


  subroutine load_posterior(filename,f,x)
    implicit none
    character(len=*), intent(in) :: filename   
    real(fp), dimension(:), pointer :: f
    real(fp), dimension(:,:), pointer :: x

    integer, parameter :: nunit = 211

    integer :: i,j, ndim, ndata


    if (associated(f).or.associated(x)) then
       stop 'load_posterior: data already loaded!'
    endif

    open(unit=nunit,file=filename,status='old')
    
    read(nunit,*) ndim, ndata

    allocate(f(ndata),x(ndim,ndata))

    do i=1,ndata
       read(nunit,*) f(i),(x(j,i),j=1,ndim)
    enddo
    
    close(nunit)

  end subroutine load_posterior



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



  subroutine cubize_paramspace(xdata,xcubes)
    implicit none
    real(fp), dimension(:,:), intent(in) :: xdata
    real(fp), dimension(:,:), intent(out) :: xcubes
    
    integer(ip) :: ndim, ndata, i
    real(fp), dimension(:), allocatable :: xmin,xmax

    ndim = size(xdata,1)
    ndata = size(xdata,2)

    if ((size(xcubes,1).ne.ndim).or.(size(xcubes,2).ne.ndata)) then
       stop 'cube_rbfparamspace: mismatch arrays!'
    endif

    allocate(xmin(ndim), xmax(ndim))

    call posterior_boundaries(xdata,xmin,xmax)

    do i=1,ndim       
       xcubes(i,:) = (xdata(i,:) - xmin(i))/(xmax(i)-xmin(i))
    enddo

    deallocate(xmin,xmax)

  end subroutine cubize_paramspace



  function count_column(filename,delimiter)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: delimiter

    logical, parameter :: display = .false.

    integer :: count_column

    integer :: stat, j, num_delim
    character :: single_byte, CR, LF, column_delimiter

    integer, parameter :: nunit = 200

    LF = char(10) ! Line Feed
    CR = char(13) ! Carriage Return

    column_delimiter = delimiter    
    
    open (nunit, file=filename, form='unformatted', &
         access='direct', status='old', recl = 1, iostat=stat)
    if (stat .ne. 0) stop 'read_hearder: missing file!'
    
    ! process header line of the file
    j = 0
    num_delim = 0
    single_byte='o'

    do while ((single_byte .ne. CR) .and. (single_byte .ne. LF))
       j = j + 1
       read(nunit, rec = j) single_byte
       if (single_byte .eq. column_delimiter) then
          num_delim = num_delim + 1
       end if
       !write (*,'(I3,5x,I3,5x,A)') j, ichar(single_byte), single_byte
    end do
    close(nunit)

    if (display) then
       write (*,*)'delimiter ',delimiter,' found ',num_delim,'times'
    endif

    count_column = num_delim    
  end function count_column



  subroutine save_weights(filename,scale,weight)
    implicit none
    character(len=*), intent(in) :: filename
    real(fp), intent(in) :: scale
    real(fp), dimension(:), intent(in) :: weight
    
    integer :: nunit
    integer(ip) :: i, nw
   
    nw = size(weight,1)

    nunit = 310

    open(nunit,file=filename,status='unknown')
    write(nunit,*) nw, scale
    do i=1,nw
       write(nunit,*) weight(i)
    enddo
    close(nunit)
  end subroutine save_weights


  



  subroutine load_weights(filename,scale,weight)
    implicit none
    character(len=*), intent(in) :: filename
    real(fp), intent(out) :: scale
    real(fp), dimension(:), pointer, intent(inout) :: weight

    integer :: nunit
    integer(ip) :: i, nw

    if (associated(weight)) then
       stop 'load_weight: weight already loaded!'
    endif
    
    nunit = 311

    open(nunit,file=filename,status='old')
    read(nunit,*) nw, scale

    write(*,*)'load_weights: nw= scale= ',nw, scale

    allocate(weight(nw))    
    do i=1,nw
       read(nunit,*) weight(i)
    enddo
    close(nunit)

  end subroutine load_weights



  subroutine save_centres(filename,xctrs)
    implicit none
    character(len=*), intent(in) :: filename
    real(fp), dimension(:,:), intent(in) :: xctrs

    integer :: nunit
    integer(ip) :: i, j, nw, ndim

    ndim = size(xctrs,1)
    nw = size(xctrs,2)

    nunit = 320

    open(nunit,file=filename,status='unknown')
    write(nunit,*) ndim, nw
    do j=1,nw
       write(nunit,*) (xctrs(i,j),i=1,ndim)
    enddo
    close(nunit)

  end subroutine save_centres



   subroutine load_centres(filename,xctrs)
    implicit none
    character(len=*), intent(in) :: filename    
    real(fp), dimension(:,:), pointer, intent(inout) :: xctrs

    integer :: nunit
    integer(ip) :: i, j,nw ,ndim

    if (associated(xctrs)) then
       stop 'load_centres: centres already loaded!'
    endif
    
    nunit = 341

    open(nunit,file=filename,status='old')
    read(nunit,*) ndim, nw
      
    write(*,*)'load_centres: ndim= nw= ',ndim,nw

    allocate(xctrs(ndim,nw))    
    do j=1,nw
       read(nunit,*) (xctrs(i,j),i=1,ndim)
    enddo
    close(nunit)

  end subroutine load_centres



  subroutine save_boundaries(filename,xmin,xmax)
    implicit none
    character(len=*), intent(in) :: filename
    real(fp), dimension(:), intent(in) :: xmin,xmax

    integer :: nunit
    integer(ip) :: j, ndim

    ndim = size(xmin,1)
    if (ndim.ne.size(xmax,1)) stop 'save_boundaries: xmin/xmax dim!'

    nunit = 410

    open(nunit,file=filename,status='unknown')
    write(nunit,*) ndim
    write(nunit,*)(xmin(j),j=1,ndim)
    write(nunit,*)(xmax(j),j=1,ndim)
    close(nunit)

  end subroutine save_boundaries


  subroutine read_boundaries(filename,xmin,xmax)
    implicit none
    character(len=*), intent(in) :: filename
    real(fp), dimension(:), pointer, intent(inout) :: xmin,xmax

    integer :: nunit
    integer(ip) :: j, ndim

    if (associated(xmin).or.associated(xmax)) stop 'read_boundaries: already done!'

    nunit = 411

    open(nunit,file=filename,status='old')
    read(nunit,*) ndim

    allocate(xmin(ndim),xmax(ndim))

    read(nunit,*)(xmin(j),j=1,ndim)
    read(nunit,*)(xmax(j),j=1,ndim)

    close(nunit)

  end subroutine read_boundaries



  subroutine save_shepdata(filename,rmax)
    use qshepmdata, only : A, IW
    use qshepmdata, only : DX, XMIN, RSQ, WS
    use qshepmdata, only : LCELL, LNEXT

    implicit none
    character(len=*), intent(in) :: filename
    real(fp), intent(in) :: rmax
    integer :: nunit

    integer :: i,j
    integer :: n, ntt, m, nlcell

    if (.not.allocated(A)) stop 'save_shepdata: not found!'

    n = size(A,1)
    ntt = size(A,2)
    m = size(xmin)
    nlcell = size(LCELL,1)

    nunit = 511
    open(nunit,file=filename,status='unknown')
    write(nunit,*) n,ntt,m,nlcell
    write(nunit,*) rmax

    do i=1,m
       write(nunit,*) DX(i),XMIN(i)
       write(nunit,*) (IW(i,j),j=1,5)
    enddo

    do i=1,n
       write(nunit,*) RSQ(i),LNEXT(i)
       write(nunit,*) (A(i,j),j=1,ntt)
    enddo

    do i=1,ntt*ntt
       write(nunit,*) WS(i)
    enddo

    do i=1,nlcell
       write(nunit,*) LCELL(i)
    enddo

    close(nunit)

  end subroutine save_shepdata

 

  subroutine load_shepdata(filename,rmax)
    use qshepmdata, only : A, IW
    use qshepmdata, only : DX, XMIN, RSQ, WS
    use qshepmdata, only : LCELL, LNEXT

    implicit none
    character(len=*), intent(in) :: filename
    real(fp), intent(out) :: rmax
    integer :: nunit

    integer :: i,j
    integer :: n, ntt, m, nlcell

    if (allocated(A)) stop 'load_shepdata: already loaded'

    n = size(A,1)
    ntt = size(A,2)
    m = size(xmin)
    nlcell = size(LCELL,1)

    nunit = 511
    open(nunit,file=filename,status='old')
    read(nunit,*) n,ntt,m,nlcell
    read(nunit,*) rmax

    allocate(A(n,ntt))
    allocate(DX(m),XMIN(m),IW(m,5))
    allocate(RSQ(N),LNEXT(N))
    allocate(LCELL(nlcell))
    allocate(WS(ntt*ntt))
    

    do i=1,m
       read(nunit,*) DX(i),XMIN(i)
       read(nunit,*) (IW(i,j),j=1,5)
    enddo

    do i=1,n
       read(nunit,*) RSQ(i),LNEXT(i)
       read(nunit,*) (A(i,j),j=1,ntt)
    enddo

    do i=1,ntt*ntt
       read(nunit,*) WS(i)
    enddo

    do i=1,nlcell
       read(nunit,*) LCELL(i)
    enddo

    close(nunit)

  end subroutine load_shepdata


end module iond
