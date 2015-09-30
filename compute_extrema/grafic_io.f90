module grafic_io
  ! Modifs PC le 08/08/08
  use ISO_C_BINDING
! END Modifs PC le 08/08/08
  
  use grafic_types
  implicit none
!  integer, parameter :: nblock = 10
  integer, parameter :: RECORD_FLAG=4 !bytes 
  
  type taille
     integer :: nx
     integer :: ny
     integer :: nz
     real(kind=sp) :: dx
     real(kind=sp) :: lx
     real(kind=sp) :: ly
     real(kind=sp) :: lz
  end type taille

  type cosmo
     real(kind=sp) :: astart
     real(kind=sp) :: omegam
     real(kind=sp) :: omegav
     real(kind=sp) :: h0
  end type cosmo

  interface grafic_write
     module procedure grafic_write_single, grafic_write_double
  end interface

  interface grafic_read
     module procedure grafic_read_single, grafic_read_double
  end interface

contains

  
  ! This routine write a fftw_mpi slice of a cube 
  ! This slice has interior dimension nz, exterior nx.
  ! Beware that interior dimension is padded: nz <- 2(nz/2+1)
  ! Padded region will not be written.
  subroutine grafic_write_double(buffer,local_nz,local_z_start,ny,nx,filename,padding_in, white_in)
    
    ! Arguments
    implicit none
    include 'mpif.h'
    real(dp), dimension(:), intent(in) :: buffer
    integer, intent(in) :: local_nz, local_z_start, ny, nx
    character(len=128), intent(in) :: filename
    !integer, intent(in) :: filenamelen
    logical, optional :: padding_in
    logical, optional :: white_in

    ! Local variables
    real(sp), allocatable, dimension(:) :: tampon
    integer :: record_size
    integer(RECORD_FLAG) :: taille_tampon
    integer :: n2x
    ! Modifs PC le 08/08/08
    ! integer(i8b) :: index, index2, offset, length, toto
    integer(i8b) :: index, index2, toto
    integer(i8b) :: offset
    integer(C_LONG) :: length

    ! END Modifs PC le 08/08/08
    integer :: fd
    
    integer :: idummy=0
    integer :: i,j,k,myid,ierr
    logical :: padding, white
    character(len=5)::nchar
    integer::iocheck

    padding=.false.
    white=.false.
    if (present(padding_in)) padding = padding_in
    if (present(white_in)) white = white_in
    ! 2*4 for leading and trailing record sizes 
    n2x = 2*(nx/2+1)
    if (padding) then
       allocate(tampon(ny*n2x))
       taille_tampon = ny*n2x*sp ! in bytes
    else
       allocate(tampon(ny*nx))
       taille_tampon = ny*nx*sp ! in bytes
    endif

#ifndef SLICES
    if (white) then
       offset = 4*4 + 2*RECORD_FLAG
    else
       offset = 3*4+8*sp + 2*RECORD_FLAG ! First record
    endif

    if (padding) then
       offset = offset + local_z_start*int(ny*n2x*sp + 2*RECORD_FLAG,8)
    else
       offset = offset + local_z_start*int(ny*nx*sp + 2*RECORD_FLAG,8)
    endif
#endif

#ifndef SLICES
    call f77_parallel_openwrite(trim(filename),len_trim(filename),fd)
#endif

    do k=1,local_nz
#ifdef SLICES     
       call title(k+local_z_start,nchar)
       open(111,file=trim(filename)//'_'//nchar,form='unformatted',status='unknown',iostat=iocheck)
       if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_WRITE_DOUBLE OPEN',iocheck
          stop
       endif
#endif
       
       index2=1
       ! This loop to (eventually) get rid of fftw padding zone
       if (padding) then
       do j=1,ny
          do i=1,n2x
             index = ((k-1)*ny+j-1)*n2x+i
             tampon(index2) = real(buffer(index),kind=sp)
             index2 = index2 + 1
          enddo
       enddo
       else
          do j=1,ny
             do i=1,nx
                index = ((k-1)*ny+j-1)*n2x+i
                tampon(index2) = real(buffer(index),kind=sp)
                index2 = index2 + 1
             enddo
          enddo
       endif

#ifndef SLICES       
       ! First write f... Fortran record length
       length=RECORD_FLAG
       call f77_parallel_write(fd, length, offset, taille_tampon)
       offset = offset+RECORD_FLAG
       if (padding) then
          length = ny*n2x*sp
       else
          length = ny*nx*sp
       endif
       call f77_parallel_write(fd, length, offset, tampon)
       offset = offset+length
       ! Again ..
       length=RECORD_FLAG
       call f77_parallel_write(fd, length, offset, taille_tampon)
       offset = offset+RECORD_FLAG
#endif

#ifdef SLICES
       write(111,iostat=iocheck)tampon
       if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_WRITE_DOUBLE WRITE',iocheck
          stop
       endif
#endif


#ifdef SLICES     
       close(111,iostat=iocheck)
       if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_WRITE_DOUBLE CLOSE',iocheck
          stop
       endif
#endif

    enddo
    deallocate(tampon)

#ifndef SLICES
    call f77_parallel_close(fd)
#endif
  end subroutine grafic_write_double
       
  subroutine grafic_write_single(buffer,local_nz,local_z_start,ny,nx,filename,padding_in, white_in)
    
    ! Arguments
    implicit none
    include 'mpif.h'
    real(sp), dimension(:), intent(in) :: buffer
    integer, intent(in) :: local_nz, local_z_start, ny, nx
    character(len=128), intent(in) :: filename
    !integer, intent(in) :: filenamelen
    logical, optional :: padding_in
    logical, optional :: white_in


    ! Local variables
    real(sp), allocatable, dimension(:) :: tampon
    integer :: record_size
    integer(RECORD_FLAG) :: taille_tampon
    integer :: n2x
    ! Modifs PC le 08/08/08
    !    integer(i8b) :: index, index2, offset, length, toto
    integer(i8b) :: index, index2, toto
    integer(i8b) :: offset
    integer(C_LONG) :: length
    ! END Modifs PC le 08/08/08
    integer :: fd
    integer :: idummy=0
    integer :: i,j,k,myid,ierr
    logical :: padding, white
    character(len=5)::nchar
    integer::iocheck

    padding=.false.
    white=.false.
    if (present(padding_in)) padding = padding_in
    if (present(white_in)) white = white_in
    ! 2*4 for leading and trailing record sizes 
    n2x = 2*(nx/2+1)
    if (padding) then
       allocate(tampon(ny*n2x))
       taille_tampon = ny*n2x*sp ! in bytes
    else
       allocate(tampon(ny*nx))
       taille_tampon = ny*nx*sp ! in bytes
    endif
#ifndef SLICES
    if (white) then
       offset = 4*4 + 2*RECORD_FLAG
    else
       offset = 3*4+8*sp + 2*RECORD_FLAG ! First record
    endif

    if (padding) then
       offset = offset + local_z_start*int(ny*n2x*sp + 2*RECORD_FLAG,8)
    else
       offset = offset + local_z_start*int(ny*nx*sp + 2*RECORD_FLAG,8)
    endif
#endif
#ifndef SLICES
    call f77_parallel_openwrite(trim(filename),len_trim(filename),fd)
#endif

    do k=1,local_nz
#ifdef SLICES     
       call title(k+local_z_start,nchar)
       open(111,file=trim(filename)//'_'//nchar,form='unformatted',status='unknown',iostat=iocheck)
       if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_WRITE_SINGLE OPEN',iocheck
          stop
       endif
#endif

       index2=1
       ! This loop to (eventually) get rid of fftw padding zone
       if (padding) then
       do j=1,ny
          do i=1,n2x
             index = ((k-1)*ny+j-1)*n2x+i
             tampon(index2) = real(buffer(index),kind=sp)
             index2 = index2 + 1
          enddo
       enddo
       else
          do j=1,ny
             do i=1,nx
                index = ((k-1)*ny+j-1)*n2x+i
                tampon(index2) = real(buffer(index),kind=sp)
                index2 = index2 + 1
             enddo
          enddo
       endif
#ifndef SLICES
       ! First write f... Fortran record length
       length=RECORD_FLAG
       call f77_parallel_write(fd, length, offset, taille_tampon)
       offset = offset+RECORD_FLAG
       if (padding) then
          length = ny*n2x*sp
       else
          length = ny*nx*sp
       endif
       call f77_parallel_write(fd, length, offset, tampon)
       offset = offset+length
       ! Again ..
       length=RECORD_FLAG
       call f77_parallel_write(fd, length, offset, taille_tampon)
       offset = offset+RECORD_FLAG
#endif

#ifdef SLICES
       write(111,iostat=iocheck)tampon
       if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_WRITE_SINGLE WRITE',iocheck
          stop
       endif
#endif
#ifdef SLICES     
       close(111,iostat=iocheck)
       if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_WRITE_SINGLE CLOSE',iocheck
          stop
       endif
#endif

    enddo
    deallocate(tampon)
#ifndef SLICES
    call f77_parallel_close(fd)
#endif
  end subroutine grafic_write_single

  subroutine grafic_read_double(buffer,local_nz,local_z_start,ny,nx,filename, padding_in, white_in)

    ! Arguments
    implicit none
    include 'mpif.h'
    integer, intent(in) :: local_nz, local_z_start, ny, nx
    real(dp), dimension(local_nz*ny*2*(nx/2+1)), intent(out) :: buffer
    character(len=128), intent(in) :: filename
    logical, optional :: padding_in ! Read padded zone or not
    logical, optional :: white_in ! Different header for white noise file
    ! Local variables
    real(sp), allocatable, dimension(:) :: tampon
    integer :: record_size
    integer(RECORD_FLAG) :: taille_tampon
    integer :: n2x
    
    ! Modifs PC le 08/08/08
    !    integer(i8b) :: index, index2, offset, length
    integer(i8b) :: index, index2
    integer(i8b) :: offset
    integer(C_LONG) :: length


    ! END Modifs PC le 08/08/08
    integer :: fd

    integer :: idummy=0
    integer :: i,j,k
    integer :: myid,ierr
    logical :: padding
    logical :: white
    character(len=5)::nchar

    integer::sizearray
    integer(i8b)::count
    integer::iocheck

    count=0
    sizearray=local_nz*ny*2*(nx/2+1)
    
    buffer=0.!!!!!RY 25/11/11 (to avoid undefined value)

    padding = .false.
    white = .false.
    if (present(padding_in)) padding = padding_in
    if (present(white_in)) white = white_in

    ! 2*4 for leading and trailing record sizes 
    n2x = 2*(nx/2+1)
    if (padding) then
       allocate(tampon(ny*n2x))
    else 
       allocate(tampon(ny*nx))
    endif
    taille_tampon = 0

#ifndef SLICES    
    if (white) then
       offset = 4*4 + 2*RECORD_FLAG
    else
       offset = 3*4+8*sp + 2*RECORD_FLAG ! First record
    endif
    if (padding) then
       offset = offset + local_z_start*int(ny*n2x*sp + 2*RECORD_FLAG,8)
    else 
       offset = offset + local_z_start*int(ny*nx*sp+2*RECORD_FLAG,8)
    endif
#endif
#ifndef SLICES    
    call f77_parallel_openread(trim(filename),len_trim(filename),fd)
#endif
    do k=1,local_nz
#ifdef SLICES     
       call title(k+local_z_start,nchar)
       open(111,file=trim(filename)//'_'//nchar,form='unformatted',status='old',iostat=iocheck)
       if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_READ_DOUBLE OPEN',iocheck
          stop
       endif
#endif

#ifndef SLICES         
       length = RECORD_FLAG 
       ! First read f... Fortran record length
       call f77_parallel_read(fd, length, offset, taille_tampon)
       call f77_parallel_read(fd, length, offset, tampon)
       if (padding) then
           if (taille_tampon /= ny*n2x*sp) then
             print*,'Record size ',taille_tampon,' is different from expected', &
                  & ny*n2x*sp
             stop
          endif
       else
           if (taille_tampon /= ny*nx*sp) then
             print*,'Record size ',taille_tampon,' is different from expected', &
                  & ny*nx*sp
             print*,trim(filename),len_trim(filename)
             print*,length,offset
             print*,i8b
             stop
          endif
       endif
       offset = offset+RECORD_FLAG
       if (padding) then
          length = ny*n2x*sp
       else 
          length = ny*nx*sp
       endif
     
       call f77_parallel_read(fd, length, offset, tampon)
       offset = offset+length
       ! Again ..
       length=RECORD_FLAG
       call f77_parallel_read(fd, length, offset, taille_tampon)
       offset = offset+RECORD_FLAG
#endif

#ifdef SLICES
       read(111,iostat=iocheck)tampon
       if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_READ_DOUBLE READ',iocheck
          stop
       endif
#endif
#ifdef SLICES     
       close(111,iostat=iocheck)
       if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_READ_DOUBLE CLOSE',iocheck
          stop
       endif
#endif


       ! This loop to (eventually) get rid of fftw padding zone
       index2=1
       if (padding) then
          do j=1,ny
             do i=1,n2x ! Read padding (Nyquist) : reading in k space ...
                index=((k-1)*ny+j-1)*n2x+i

                !!!!RY check array index
                if(index<1.or.index>sizearray) then!!!RY
                   write(*,*)'PB INDEX',index,sizearray!!!RY
                   stop!!!RY
                endif!!!RY

                buffer(index) = real(tampon(index2),kind=dp)
                index2 = index2 + 1
             enddo
          enddo
       else
          do j=1,ny
             do i=1,nx ! Do not read padding : reading in real space ...
                index = ((k-1)*ny+j-1)*n2x+i
                
                !!!!RY check array index
                if(index<1.or.index>sizearray) then!!!RY
                   write(*,*)'PB INDEX',index,sizearray!!!RY
                   stop!!!RY
                endif!!!RY
                
                if(buffer(index).ne.0.)write(*,*)'WARNING not zero in buffer before fill',i,j,k,index,index2, buffer(index),tampon(index2)
                buffer(index) = real(tampon(index2),kind=dp)
                if(buffer(index)==0.)write(*,*)'WARNING zero in buffer',i,j,k,index,index2, buffer(index),tampon(index2)
                index2 = index2 + 1
                count=count+1
             enddo
          enddo

          !additionnal check yann to be removed
          do j=1,ny
             do i=nx+1,n2x ! Do not read padding : reading in real space ...                                                                                                                                       
                index = ((k-1)*ny+j-1)*n2x+i

                if(buffer(index).ne.0.)write(*,*)'WARNING non zero in padding in grafic_io',i,j,k,index,buffer(index)
             enddo
          enddo
          !end additionnal check yann to be removed


       endif


       

    enddo
    if(count.ne.nx*ny*local_nz.and..not.padding)write(*,*)'bizarre count',count,nx*ny*local_nz,nx,ny,local_nz

    deallocate(tampon)
#ifndef SLICES  
    call f77_parallel_close(fd)
#endif
  end subroutine grafic_read_double

 

  subroutine grafic_read_single(buffer,local_nz,local_z_start,ny,nx,filename, padding_in, white_in)

    ! Arguments
    implicit none
    include 'mpif.h'
    integer, intent(in) :: local_nz, local_z_start, ny, nx
    real(sp), dimension(local_nz*ny*2*(nx/2+1)), intent(out) :: buffer
    character(len=128), intent(in) :: filename
    logical, optional :: padding_in ! Read padded zone or not
    logical, optional :: white_in ! Different header for white noise file
    ! Local variables
    real(sp), allocatable, dimension(:) :: tampon
    integer :: record_size
    integer(RECORD_FLAG) :: taille_tampon
    integer :: n2x
    
    ! Modifs PC le 08/08/08
    !    integer(i8b) :: index, index2, offset, length
    integer(i8b) :: index, index2
    integer(i8b) :: offset
    integer(C_LONG) :: length


    ! END Modifs PC le 08/08/08
    integer :: fd

    integer :: idummy=0
    integer :: i,j,k
    integer :: myid,ierr
    logical :: padding
    logical :: white
    character(len=5)::nchar

    integer::sizearray
    integer::iocheck

    sizearray=local_nz*ny*2*(nx/2+1)
        
    buffer=0.!!!!!RY 25/11/11 (to avoid undefined value)

    padding = .false.
    white = .false.
    if (present(padding_in)) padding = padding_in
    if (present(white_in)) white = white_in
    

    ! 2*4 for leading and trailing record sizes 
    n2x = 2*(nx/2+1)
    if (padding) then
       allocate(tampon(ny*n2x))
    else 
       allocate(tampon(ny*nx))
    endif
    taille_tampon = 0

#ifndef SLICES    
    if (white) then
       offset = 4*4 + 2*RECORD_FLAG
    else
       offset = 3*4+8*sp + 2*RECORD_FLAG ! First record
    endif
    if (padding) then
       offset = offset + local_z_start*int(ny*n2x*sp + 2*RECORD_FLAG,8)
    else 
       offset = offset + local_z_start*int(ny*nx*sp+2*RECORD_FLAG,8)
    endif
#endif
#ifndef SLICES    
    call f77_parallel_openread(trim(filename),len_trim(filename),fd)
#endif
    do k=1,local_nz
#ifdef SLICES     
       call title(k+local_z_start,nchar)
       open(111,file=trim(filename)//'_'//nchar,form='unformatted',status='old',iostat=iocheck)
       if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_READ_SINGLE OPEN',iocheck
          stop
       endif
#endif

#ifndef SLICES         
       length = RECORD_FLAG 
       ! First read f... Fortran record length
       call f77_parallel_read(fd, length, offset, taille_tampon)
       call f77_parallel_read(fd, length, offset, tampon)
       if (padding) then
           if (taille_tampon /= ny*n2x*sp) then
             print*,'Record size ',taille_tampon,' is different from expected', &
                  & ny*n2x*sp
             stop
          endif
       else
           if (taille_tampon /= ny*nx*sp) then
             print*,'Record size ',taille_tampon,' is different from expected', &
                  & ny*nx*sp
             print*,trim(filename),len_trim(filename)
             print*,length,offset
             print*,i8b
             stop
          endif
       endif
       offset = offset+RECORD_FLAG
       if (padding) then
          length = ny*n2x*sp
       else 
          length = ny*nx*sp
       endif
     
       call f77_parallel_read(fd, length, offset, tampon)
       offset = offset+length
       ! Again ..
       length=RECORD_FLAG
       call f77_parallel_read(fd, length, offset, taille_tampon)
       offset = offset+RECORD_FLAG
#endif

#ifdef SLICES
       read(111,iostat=iocheck)tampon
       if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_READ_SINGLE READ',iocheck
          stop
       endif
#endif
#ifdef SLICES     
       close(111,iostat=iocheck)
       if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_READ_SINGLE CLOSE',iocheck
          stop
       endif
#endif


       ! This loop to (eventually) get rid of fftw padding zone
       index2=1
       if (padding) then
          do j=1,ny
             do i=1,n2x ! Read padding (Nyquist) : reading in k space ...
                index=((k-1)*ny+j-1)*n2x+i

                !!!!RY check array index
                if(index<1.or.index>sizearray) then!!!RY
                   write(*,*)'PB INDEX',index,sizearray!!!RY
                   stop!!!RY
                endif!!!RY

                buffer(index) = real(tampon(index2),kind=dp)
                index2 = index2 + 1
             enddo
          enddo
       else
          do j=1,ny
             do i=1,nx ! Do not read padding : reading in real space ...
                index = ((k-1)*ny+j-1)*n2x+i
                
                !!!!RY check array index
                if(index<1.or.index>sizearray) then!!!RY
                   write(*,*)'PB INDEX',index,sizearray!!!RY
                   stop!!!RY
                endif!!!RY

                buffer(index) = real(tampon(index2),kind=dp)
                index2 = index2 + 1
             enddo
          enddo
       endif

       

    enddo
    deallocate(tampon)
#ifndef SLICES  
    call f77_parallel_close(fd)
#endif
  end subroutine grafic_read_single


  subroutine grafic_read_header(filename,head_taille,head_cosmo)
    
    type(taille), intent(out) :: head_taille
    type(cosmo), intent(out) :: head_cosmo
    character(len=128) :: filename
    character(len=150) :: file
    character(len=15) :: ext
    logical :: ok
    integer::iocheck
    
    ext='_header'
    
#ifndef SLICES
    file=filename
#else
    file=trim(filename)//trim(ext)
#endif   
    
    inquire(file=trim(file),exist=ok)
    if (.not. ok) then
       print*,'File '
       write(*,'(a150)')trim(file)
       print*,' does not exist, aborting.'
       stop
    endif

    open(2,file=trim(file),form='unformatted',status='old',iostat=iocheck)
    if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_READ_HEADER OPEN',iocheck
          stop
       endif
    read(2,iostat=iocheck) head_taille%nx, head_taille%ny, head_taille%nz, head_taille%dx, &
         head_taille%lx, head_taille%ly, head_taille%lz, head_cosmo%astart, &
         head_cosmo%omegam, head_cosmo%omegav, head_cosmo%h0
    if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_READ_HEADER READ',iocheck
          stop
       endif
    close(2,iostat=iocheck)
    if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_READ_HEADER CLOSE',iocheck
          stop
       endif

  end subroutine grafic_read_header

  subroutine grafic_read_header_white(filename,nx,ny,nz,iseed)
    
    integer, intent(out) :: nx,ny,nz,iseed
    character(len=128), intent(in) :: filename
     character(len=150) :: file
    logical :: ok
    integer::iocheck
#ifndef SLICES
    file=filename
#else
    file=trim(filename)//'_header'
#endif

    inquire(file=file,exist=ok)
    if (.not. ok) then
       print*,'File ',trim(file),' does not exist, aborting.'
       stop
    endif

    open(2,file=file,form='unformatted',status='old',iostat=iocheck)
    if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_READ_HEADER_WHITE OPEN',iocheck
          stop
       endif
    read(2,iostat=iocheck) nx, ny, nz, iseed
     if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_READ_HEADER_WHITE READ',iocheck
          stop
       endif
    close(2,iostat=iocheck)
     if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_READ_HEADER_WHITE CLOSE',iocheck
          stop
       endif


  end subroutine grafic_read_header_white
  
  subroutine grafic_write_header(filename,head_taille,head_cosmo)
    
    type(taille), intent(in) :: head_taille
    type(cosmo), intent(in) :: head_cosmo
    character(len=128), intent(in) :: filename
    character(len=150) :: file
    integer::iocheck

#ifndef SLICES
    file=filename
#else
    file=trim(filename)//'_header'
#endif

    open(2,file=file,form='unformatted',status='unknown',iostat=iocheck)
    if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_WRITE_HEADER OPEN',iocheck
          stop
       endif
    write(2,iostat=iocheck) head_taille%nx, head_taille%ny, head_taille%nz, head_taille%dx, &
         head_taille%lx, head_taille%ly, head_taille%lz, head_cosmo%astart, &
         head_cosmo%omegam, head_cosmo%omegav, head_cosmo%h0
    if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_WRITE_HEADER WRITE',iocheck
          stop
       endif
    close(2,iostat=iocheck)
    if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_WRITE_HEADER CLOSE',iocheck
          stop
       endif

  end subroutine grafic_write_header

  subroutine grafic_write_header_white(filename,nx,ny,nz,iseed)
    
    integer, intent(in) :: nx,ny,nz,iseed
    character(len=128), intent(in) :: filename
     character(len=150) :: file
     integer::iocheck
#ifndef SLICES
    file=filename
#else
    file=trim(filename)//'_header'
#endif

    open(2,file=file,form='unformatted',status='unknown',iostat=iocheck)
    if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_WRITE_HEADER_WHITE OPEN',iocheck
          stop
       endif
    write(2,iostat=iocheck) nx, ny, nz, iseed
    if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_WRITE_HEADER_WHITE WRITE',iocheck
          stop
       endif
    close(2,iostat=iocheck)
    if(iocheck.ne.0)then
          write(*,*)'ERREUR IOSTAT NON ZERO GRAFIC_WRITE_HEADER_WHITE CLOSE',iocheck
          stop
       endif

  end subroutine grafic_write_header_white

       
end module grafic_io


  
    subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character(LEN=5)::nchar

  character(LEN=1)::nchar1
  character(LEN=2)::nchar2
  character(LEN=3)::nchar3
  character(LEN=4)::nchar4
  character(LEN=5)::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif


end subroutine title


