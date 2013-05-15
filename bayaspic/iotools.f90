module iotools

  implicit none
  
  integer, parameter :: lenIoRank = 4
  integer, parameter :: lenIoMax = 64
  integer, parameter :: sp=kind(1._4)
  integer, parameter :: dp=kind(1._8)
  
 
  interface livewrite
     module procedure sp_livewrite, dp_livewrite
  end interface

  interface allwrite
     module procedure sp_allwrite, dp_allwrite
  end interface

  private

  public lenIoRank, lenIoMax
  public replace_char,step2char,run2char
  public delete_file,livewrite,allwrite

contains

  

  subroutine replace_char(strg,charold,charnew)
    implicit none
    character(len=*), intent(inout) :: strg
    character, intent(in) :: charold, charnew

    integer :: position
    position = 1

    do 
       position = scan(strg,charold)
       if (position.ne.0) then
          strg(position:position) = charnew
       else
          exit
       endif
    enddo
        
  end subroutine replace_char



  subroutine step2char(igstep,clen,cgstep)
    implicit none
    integer, intent(in) :: igstep
    integer, intent(in) :: clen
    character(len=clen), intent(out) :: cgstep

    character(len=lenIoMax) :: numToStrg
    character(len=clen) :: strg

    if (clen.gt.lenIoMax) then
       stop 'step2char: clen > lenIoMax!'
    endif

    write(numToStrg,*) igstep
    strg = trim(adjustl(numToStrg)) 
    strg = adjustr(strg)
    call replace_char(strg,' ','0')

    cgstep(1:clen) = strg(1:clen)

  end subroutine step2char



  subroutine run2char(irun,clen,crun)
    implicit none
    integer, intent(in) :: irun
    integer, intent(in) :: clen
    character(len=clen), intent(out) :: crun

    character(len=lenIoMax) :: numToStrg
    character(len=clen) :: strg

    if (clen.gt.lenIoMax) then
       stop 'step2char: clen > lenIoMax!'
    endif

    write(numToStrg,*) irun
    strg = trim(adjustl(numToStrg)) 
!    strg = adjustr(strg)
!    call replace_char(strg,' ','0')

    crun(1:clen) = strg(1:clen)

  end subroutine run2char



   subroutine delete_file(name)
    implicit none
    character(len=*) :: name
    logical :: isthere

    inquire(file=name,exist=isthere)

    if (isthere) then
       open(unit=10,file=name)
       close(unit=10,status='delete')
    endif

  end subroutine delete_file


  subroutine sp_livewrite(name,x,a,b,c,d,e,f,g)
    implicit none
    character(len=*) :: name    
    real(sp) :: x,a
    real(sp), optional :: b,c,d,e,f,g
      
    open(10,file=name,position='append',status='unknown')
    
    if (.not.present(b)) then
       write(10,100) x,a         
    elseif (.not.present(c)) then           
       write(10,100) x,a,b         
    elseif (.not.present(d)) then             
       write(10,100) x,a,b,c                    
    elseif (.not.present(e)) then         
       write(10,100) x,a,b,c,d                    
    elseif (.not.present(f)) then
       write(10,100) x,a,b,c,d,e            
    elseif (.not.present(g)) then
       write(10,100) x,a,b,c,d,e,f            
    else
       write(10,100) x,a,b,c,d,e,f,g            
    endif
    
    close(10)

100 format(8(ES25.16E3))

  end subroutine sp_livewrite


  subroutine dp_livewrite(name,x,a,b,c,d,e,f,g)
    implicit none
    character(len=*) :: name    
    real(dp) :: x,a
    real(dp), optional :: b,c,d,e,f,g
    
    open(10,file=name,position='append',status='unknown')
    
    if (.not.present(b)) then
       write(10,100) x,a         
    elseif (.not.present(c)) then           
       write(10,100) x,a,b         
    elseif (.not.present(d)) then             
       write(10,100) x,a,b,c                    
    elseif (.not.present(e)) then         
       write(10,100) x,a,b,c,d                    
    elseif (.not.present(f)) then
       write(10,100) x,a,b,c,d,e            
    elseif (.not.present(g)) then
       write(10,100) x,a,b,c,d,e,f            
    else
       write(10,100) x,a,b,c,d,e,f,g            
    endif
    
    close(10)
    
100 format(8(ES25.16E3))      
    
  end subroutine dp_livewrite


  subroutine sp_allwrite(name,x,a,b,c,d,e,f,g)
    implicit none
      character(*) :: name
      integer :: j,npts
      real(sp) :: x(:),a(:)
      real(sp), optional :: b(:),c(:),d(:),e(:),f(:),g(:)

      npts=ubound(x,1)
      
      if (ubound(a,1).ne.npts) then
         write(*,*)'WARNING: vectors length differ'
      endif

!      write(*,*)'__write: save in ',name
      open(10,file=name,status='unknown')
      
      if (.not.present(b)) then
         do j=1,npts      
            write(10,100) x(j),a(j)
         enddo
      elseif (.not.present(c)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j)
         enddo
      elseif (.not.present(d)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j)            
         enddo
      elseif (.not.present(e)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j)            
         enddo
      elseif (.not.present(f)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j),e(j)            
         enddo
      elseif (.not.present(g)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j),e(j),f(j)       
         enddo
      else
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j),e(j),f(j),g(j)            
         enddo
      endif
      
      close(10)

100   format(8(ES24.16))      

    end subroutine sp_allwrite


    subroutine dp_allwrite(name,x,a,b,c,d,e,f,g)
      implicit none
      character(*) :: name
      integer :: j,npts
      real(dp) :: x(:),a(:)
      real(dp), optional :: b(:),c(:),d(:),e(:),f(:),g(:)

      npts=ubound(x,1)
      
      if (ubound(a,1).ne.npts) then
         write(*,*)'WARNING: vectors length differ'
      endif

!      write(*,*)'__write: save in ',name
      open(10,file=name,status='unknown')
      
      if (.not.present(b)) then
         do j=1,npts      
            write(10,100) x(j),a(j)
         enddo
      elseif (.not.present(c)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j)
         enddo
      elseif (.not.present(d)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j)            
         enddo
      elseif (.not.present(e)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j)            
         enddo
      elseif (.not.present(f)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j),e(j)            
         enddo
      elseif (.not.present(g)) then
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j),e(j),f(j)       
         enddo
      else
         do j=1,npts      
            write(10,100) x(j),a(j),b(j),c(j),d(j),e(j),f(j),g(j)            
         enddo
      endif
      
      close(10)

100   format(8(ES24.16))      

    end subroutine dp_allwrite




end module iotools

