module inoutfile

implicit none

private

integer, parameter :: sp = kind(1.0_4)
integer, parameter :: dp = kind(1.0_8)


interface livewrite
   module procedure slivewrite, dlivewrite
end interface

interface allwrite
   module procedure sallwrite, dallwrite
end interface

interface vecallwrite
   module procedure svecallwrite, dvecallwrite
end interface

interface vecbinallwrite
   module procedure svecbinallwrite, dvecbinallwrite
end interface

interface allread
   module procedure slecture,dlecture
end interface

interface vecallread
   module procedure sveclecture,dveclecture
end interface

interface vecbinallread
   module procedure svecbinlecture,dvecbinlecture
end interface

interface lecture
   module procedure slecture,dlecture
end interface

interface veclecture
   module procedure sveclecture,dveclecture
end interface

interface vecbinlecture
   module procedure svecbinlecture,dvecbinlecture
end interface


interface livecriture
   module procedure slivewrite, dlivewrite
end interface
 
interface ecriture 
   module procedure sallwrite, dallwrite
end interface

public livewrite, allwrite, vecallwrite
public allread,vecallread, vecbinallread
public veclecture, livecriture, ecriture
public vecbinallwrite,vecbinlecture, lecture
public delete_file

contains

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



  subroutine slivewrite(name,x,a,b,c,d,e,f,g)
    implicit none
    character(*) :: name
    integer :: istat
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
    
  end subroutine slivewrite


  subroutine dlivewrite(name,x,a,b,c,d,e,f,g)
    implicit none
    character(*) :: name
    integer :: istat
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

  end subroutine dlivewrite



  subroutine sallwrite(name,x,a,b,c,d,e,f,g)
    implicit none
    character(*) :: name
    integer :: j,npts
    real(sp) :: x(:),a(:)
    real(sp), optional :: b(:),c(:),d(:),e(:),f(:),g(:)
    
    npts=ubound(x,1)
      
    if (ubound(a,1).ne.npts) then
       write(*,*)'WARNING: vectors length differ'
    endif
    
    write(*,*)'__write: save in ',name
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
    
100 format(8(ES25.16E3))      

  end subroutine sallwrite



  subroutine dallwrite(name,x,a,b,c,d,e,f,g)
    implicit none
    character(*) :: name
    integer :: j,npts
    real(dp) :: x(:),a(:)
    real(dp), optional :: b(:),c(:),d(:),e(:),f(:),g(:)

    npts=ubound(x,1)
      
    if (ubound(a,1).ne.npts) then
       write(*,*)'WARNING: vectors length differ'
    endif
    
    write(*,*)'__write: save in ',name
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
    
100 format(8(ES25.16E3))      
    
  end subroutine dallwrite




  subroutine svecallwrite(name,x,vec)
    implicit none
    character(*) :: name
    integer :: i,j,npts,ncol
    real(sp), dimension(:) :: x
    real(sp), dimension(:,:) :: vec
    
    npts=size(x,1)
    ncol = size(vec,1)
      
    if (size(vec,2).ne.npts) then
       write(*,*)'inoutfile: WARNING: vectors length differ'
       read(*,*)
    endif
    
    write(*,*)'__write: save in ',name
    open(10,file=name,status='unknown')
    
    
    do j=1,npts      
       write(10,100) x(j),(vec(i,j),i=1,ncol)
    enddo
    
    close(10)
    
100 format(50(ES25.16E3))      

  end subroutine svecallwrite



  subroutine dvecallwrite(name,x,vec)
    implicit none
    character(*) :: name
    integer :: i,j,npts,ncol
    real(dp), dimension(:) :: x
    real(dp), dimension(:,:) :: vec

    npts=size(x,1)
    ncol = size(vec,1)
      
    if (size(vec,2).ne.npts) then
       write(*,*)'WARNING: vectors length differ'
    endif
    
    write(*,*)'__write: save in ',name         

    open(10,file=name,status='unknown')
    do j=1,npts      
       write(10,100) x(j),(vec(i,j),i=1,ncol)
    enddo

    close(10)
    
100 format(50(ES25.16E3))      

  end subroutine dvecallwrite

  

 subroutine svecbinallwrite(name,x,vec)
    implicit none
    character(*) :: name
    integer :: i,j,npts,ncol
    real(sp), dimension(:) :: x
    real(sp), dimension(:,:) :: vec
    
    npts=size(x,1)
    ncol = size(vec,1)
      
    if (size(vec,2).ne.npts) then
       write(*,*)'inoutfile: WARNING: vectors length differ'
       read(*,*)
    endif
    
    write(*,*)'__write: save in ',name
    open(10,file=name,status='unknown',form='unformatted')
    
    
    do j=1,npts      
       write(10) x(j),(vec(i,j),i=1,ncol)
    enddo
    
    close(10)
    


  end subroutine svecbinallwrite



  subroutine dvecbinallwrite(name,x,vec)
    implicit none
    character(*) :: name
    integer :: i,j,npts,ncol
    real(dp), dimension(:) :: x
    real(dp), dimension(:,:) :: vec

    npts=size(x,1)
    ncol = size(vec,1)
      
    if (size(vec,2).ne.npts) then
       write(*,*)'WARNING: vectors length differ'
    endif
    
    write(*,*)'__write: save in ',name         

    open(10,file=name,status='unknown',form='unformatted')
    do j=1,npts      
       write(10) x(j),(vec(i,j),i=1,ncol)
    enddo

    close(10)
    


  end subroutine dvecbinallwrite



  subroutine slecture(name,npts,x,a,b,c,d,e,f,g)
    implicit none
    
    character (*) :: name
    integer :: ncoldata
    real(sp) :: x(:),a(:)
    real(sp), optional :: b(:),c(:),d(:),e(:),f(:),g(:)
    integer :: jcol,jlin,npts
    
    
    jlin=0
    
    open(11,file=name,status='old')
    
    do while (.true.)
       jlin=jlin+1
       
       if (.not.present(b)) then
          read(11,*,end=100) x(jlin),a(jlin)
       elseif (.not.present(c)) then
          read(11,*,end=100) x(jlin),a(jlin),b(jlin)
       elseif (.not.present(d)) then
          read(11,*,end=100) x(jlin),a(jlin),b(jlin),c(jlin)            
       elseif (.not.present(e)) then
          read(11,*,end=100) x(jlin),a(jlin),b(jlin),c(jlin),d(jlin)
       elseif (.not.present(f)) then
          read(11,*,end=100) x(jlin),a(jlin),b(jlin),c(jlin),d(jlin)&
               ,e(jlin)            
       elseif (.not.present(g)) then
          read(11,*,end=100) x(jlin),a(jlin),b(jlin),c(jlin),d(jlin)&
               ,e(jlin),f(jlin)            
       else
          read(11,*,end=100) x(jlin),a(jlin),b(jlin),c(jlin),d(jlin)&
               ,e(jlin),f(jlin),g(jlin)        
       endif
       
       if (jlin.gt.ubound(x,1)) then
          write(*,*)'__read: array length too small'
          stop
       endif
    enddo
    
100 close(11)
    
    npts=jlin-1
    write(*,*)'__read:',npts,'lines read'
    
  end subroutine slecture
  
  
  subroutine dlecture(name,npts,x,a,b,c,d,e,f,g)
    implicit none
    
    character (*) :: name
    integer :: ncoldata
    real(dp) :: x(:),a(:)
    real(dp), optional :: b(:),c(:),d(:),e(:),f(:),g(:)
    integer :: jcol,jlin,npts
    
    
    jlin=0
    
    open(11,file=name,status='old')
    
    do while (.true.)
       jlin=jlin+1
       
       if (.not.present(b)) then
          read(11,*,end=100) x(jlin),a(jlin)
       elseif (.not.present(c)) then
          read(11,*,end=100) x(jlin),a(jlin),b(jlin)
       elseif (.not.present(d)) then
          read(11,*,end=100) x(jlin),a(jlin),b(jlin),c(jlin)            
       elseif (.not.present(e)) then
          read(11,*,end=100) x(jlin),a(jlin),b(jlin),c(jlin),d(jlin)
       elseif (.not.present(f)) then
          read(11,*,end=100) x(jlin),a(jlin),b(jlin),c(jlin),d(jlin)&
               ,e(jlin)            
       elseif (.not.present(g)) then
          read(11,*,end=100) x(jlin),a(jlin),b(jlin),c(jlin),d(jlin)&
               ,e(jlin),f(jlin)            
       else
          read(11,*,end=100) x(jlin),a(jlin),b(jlin),c(jlin),d(jlin)&
               ,e(jlin),f(jlin),g(jlin)        
       endif
       
       if (jlin.gt.ubound(x,1)) then
          write(*,*)'__read: array length too small'
          stop
       endif
    enddo
    
100 close(11)
    
    npts=jlin-1
    write(*,*)'__read:',npts,'lines read'
    
  end subroutine dlecture


  subroutine sveclecture(name,npts,x,vec)
    implicit none
    
    character (*) :: name
    integer :: ncoldata
    real(sp), dimension(:) :: x
    real(sp), dimension(:,:) :: vec
    
    integer :: jcol,jlin,npts,ncol
    
    ncol = size(vec,1)
    
    jlin=0
    
    open(11,file=name,status='old')
    
    do while (.true.)
       jlin=jlin+1
       
       read(11,*,end=100) x(jlin),(vec(jcol,jlin),jcol=1,ncol)    

       if (jlin.gt.size(x,1)) then
          write(*,*)'__read: array length too small'
          stop
       endif
    enddo
    
100 close(11)
    
    npts=jlin-1
    write(*,*)'__read:',npts,'lines read'
    
  end subroutine sveclecture



  subroutine dveclecture(name,npts,x,vec)
    implicit none
    
    character (*) :: name
    integer :: ncoldata
    real(dp), dimension(:) :: x
    real(dp), dimension(:,:) :: vec
    
    integer :: jcol,jlin,npts,ncol
    
    ncol = size(vec,1)
    
    jlin=0
    
    open(11,file=name,status='old')
    
    do while (.true.)
       jlin=jlin+1
       
       read(11,*,end=100) x(jlin),(vec(jcol,jlin),jcol=1,ncol)    

       if (jlin.gt.size(x,1)) then
          write(*,*)'__read: array length too small'
          stop
       endif
    enddo
    
100 close(11)
    
    npts=jlin-1
    write(*,*)'__read:',npts,'lines read'
    
  end subroutine dveclecture


  subroutine svecbinlecture(name,npts,x,vec)
    implicit none
    
    character (*) :: name
    integer :: ncoldata
    real(sp), dimension(:) :: x
    real(sp), dimension(:,:) :: vec
    
    integer :: jcol,jlin,npts,ncol
    
    ncol = size(vec,1)
    
    jlin=0
    
    open(11,file=name,status='old',form='unformatted')
    
    do while (.true.)
       jlin=jlin+1
       
       read(11,end=100) x(jlin),(vec(jcol,jlin),jcol=1,ncol)    

       if (jlin.gt.size(x,1)) then
          write(*,*)'__read: array length too small'
          stop
       endif
    enddo
    
100 close(11)
    
    npts=jlin-1
    write(*,*)'__read:',npts,'lines read'
    
  end subroutine svecbinlecture



  subroutine dvecbinlecture(name,npts,x,vec)
    implicit none
    
    character (*) :: name
    integer :: ncoldata
    real(dp), dimension(:) :: x
    real(dp), dimension(:,:) :: vec
    
    integer :: jcol,jlin,npts,ncol
    
    ncol = size(vec,1)
    
    jlin=0
    
    open(11,file=name,status='old',form='unformatted')
    
    do while (.true.)
       jlin=jlin+1
       
       read(11,end=100) x(jlin),(vec(jcol,jlin),jcol=1,ncol)    

       if (jlin.gt.size(x,1)) then
          write(*,*)'__read: array length too small'
          stop
       endif
    enddo
    
100 close(11)
    
    npts=jlin-1
    write(*,*)'__read:',npts,'lines read'
    
  end subroutine dvecbinlecture



end module inoutfile


