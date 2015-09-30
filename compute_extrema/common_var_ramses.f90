Module modmem
contains
subroutine writemem(usedmem)
  real::usedmem
  integer::getpagesize

#ifdef NOSYSTEM
  !  call PXFSYSCONF(_SC_PAGESIZE,ipagesize,ierror)
  ipagesize=4096
#else
  !  ipagesize = getpagesize()
  ipagesize=4096
#endif
  usedmem=dble(usedmem)*dble(ipagesize)

  if(usedmem>1024.**3.)then
     write(*,999)usedmem/1024.**3.
  else if (usedmem>1024.**2.) then
     write(*,998)usedmem/1024.**2
  else if (usedmem>1024.) then
     write(*,997)usedmem/1024.
  endif

997 format(' Used memory:',F8.3,' kB')
998 format(' Used memory:',F8.3,' MB')
999 format(' Used memory:',F8.3,' GB')

end subroutine writemem



subroutine getmem(outmem)
  real::outmem
  character(len=300) :: dir, dir2,  cmd, file
  file='/proc/self/stat'
  open(unit=1,file=file,form='formatted')
  read(1,'(A300)')dir
  close(1)
  ind=300
  j=0
  do while (j<23)
     ind=index(dir,' ')
     dir2=dir(ind+1:300)
     j=j+1
     dir=dir2
  end do
  ind=index(dir,' ')
  dir2=dir(1:ind)
  read(dir2,'(I12)')nmem

  outmem=dble(nmem)

end subroutine getmem


end Module modmem


Module modconst
  use grafic_types

! Define LONGINT in makefile to analyze simulations with more than 2**31 particles.
#ifdef LONGINT
  Integer, parameter :: PRI = 8
  Integer, parameter :: MPI_PRI = Mpi_Integer8
#else
  Integer, parameter :: PRI = 4
  Integer, parameter :: MPI_PRI = Mpi_Integer
#endif

  ! Output Units 
  Integer, parameter :: Ulog=50, Ucub=51, Umas=52, Ustr=53, Uopa=54

  
 
  

End Module modconst

Module modvariable

  Use modconst
  Integer(kind=4)   :: mynpart
  Real(   kind=SP ), dimension(:,:), allocatable :: xx   ! particles positions - Always simple precision
End Module modvariable


! timing variables
Module modtiming

  Use modconst
  Real(kind=8) :: time0, timeInt
  Real(kind=8) :: tReadFile, tTailPart, tInitRead, tReadRA
  Real(kind=8) :: tObs,tOut,tFoF,tFoFloc, tFoFinit,tRaccord

End Module modtiming


! Input parameters
Module modparam

  Use modconst
  Character(len=90) :: root               ! base name for output
  Character(len=90) :: pathinput          ! path to the output groups
  Character(len=90) :: nameinfo           ! info file name
  Character(len=90) :: namepart           ! particles files base name
  Character(len=3)  :: code_index         ! input file type: RA2= RAMSES v2, RA3= RAMSES v3
  Integer(kind=4)   :: Mmin, Mmax         ! min and max structure mass
  Integer(kind=4)   :: grpsize            ! number of particles files in each group_N directory
  Real(kind=SP)     :: perco              ! percolation parameter for Friends of Friends halo detection
  Logical           :: outcube            ! should there be an output after the particles reading/sharing?
  Logical           :: star               ! is there information about stars in the RAMSES files?
  Logical           :: metal              ! and about metallicity?
  Logical           :: dofof              ! should the structures be detected?
  Logical           :: readfromcube       ! should the particles be read from cube files?
  Logical           :: dotimings          ! should there be timings?

End Module modparam

module pathramses
  use grafic_types
  integer(i4b) :: ramses_read_part=2
  character(LEN=5)::nchar
  character(len=300)::rep
  integer::iogroupsize=0
  character(len=filenamelen) :: inputfiles,infofile,outputfile
  integer(i4b) :: nfiles=0
end module pathramses



Module modmpicom

    Use modconst

    Character, dimension(:),allocatable :: header     ! buffer for pack/unpack
    Integer(kind=4) :: mpierr                         ! mpi error code
    Integer(kind=4) :: mpireqs1,mpireqs2,mpireqs3,mpireqs4,mpireqr1,mpireqr2,mpireqr3,mpireqr4 ! mpi non blocking communications request
    Integer(kind=4) :: procID, procNB                 ! process ID and number of processes
    Integer(kind=4) :: h_length, h_pos                ! buffer length and position in the buffer
    
    Contains

    ! Abort execution when a bug is detected
    Subroutine EmergencyStop(message,errcode)

        Character(len=*), Intent(in) :: message  ! Error message
        Integer, Intent(in) :: errcode           ! Error code: this code is not meaningful

        Write(*,1000) message,procID

        1000 Format('*** ',A,' on process ',I5.5,'. ***')

        Call Mpi_Abort(MPICube,errcode,mpierr)
        Stop

    End Subroutine EmergencyStop

End Module modmpicom


Module modio

  Use modconst

  Integer(kind=4)   :: reccnt, dircnt       ! local particle number
  Integer(kind=4)   :: ndim          ! number of dimensions
  Integer(kind=4)   :: lmin          ! minimum mesh refinement level
  Integer(kind=4)   :: nres          ! 1-D resolution: number of grid points in each dimension ; nres = 2 ** lmin
  Integer(kind=PRI) :: nptot,nptotout         !   particle number: nptot = nres ** 3
  Integer(kind=PRI) :: ngrid         ! total number of grid points: ngrid = nres ** 3
  Real   (kind=SP)  :: xmin, xmax, ymin, ymax, zmin, zmax  ! min and max (x,y,z) for each process
  Real   (kind=SP)  :: rescale, boxlen_size
  Real   (kind=SP)  :: aexpt,unit_t,unit_l,unit_d,Omm0

  integer(kind=4)  ::noutmin,noutmax
  Integer(kind=4)  :: nprocramses       ! process number  read in RAMSES info file

Contains



  ! Read particles files created by RAMSES
  ! RAMSES writes as many particles files as it used process to run the simulation.
  ! Each process will read a part of the files.
  ! Particles are spread among the different files and have to be sorted by position and distributed between the processes.
  ! Each process in the parallel FoF will analize a cube of particles, this means the particles have to be geographically
  ! distributed.
  Subroutine ramses_lecture(xminin,xmaxin,yminin,ymaxin,zminin,zmaxin,zmininbuf,zmaxinbuf,nxin,nyin,nzin,local_nz)
    Use modvariable
    Use modparam
    Use modmpicom
    Use modtiming
    use pathramses
    Implicit none

    ! Local variables
    Character(len=5)               :: ncharcpu
    Character(len=9)               :: tmpstr1, tmpstr2
    Character(len=200)             :: nomfich
    Character(len=13)              :: dumchar
    Character(len=11)              :: grpchar

    Integer(kind=4)                :: i, j, icpu, icputmp,idim   ! loop variables
    Integer(kind=4)                :: destCoord(3)       ! coords of the destination MPI process in MPI process cart
    Integer(kind=4)                :: nrecv              ! number of elements received in a Mpi_Recv
    Integer(kind=4)                :: recvpoint          ! address of the 1st element received in the local vector
    Integer(kind=4)                :: mynbfile            ! number of RAMSES part files read by local process
    Integer(kind=4)                :: firstp, lastp      ! id of 1st and last RAMSES part file read
    Integer(kind=4), allocatable   :: npartvloc(:), npartv(:)  ! temp and global table of particle numbers for each process
    Integer(kind=4)                :: n_i, n_j, n_k, nsd, ind
    Integer(kind=4)                :: nprocramses       ! process number  read in RAMSES info file
    Integer(kind=4)                :: ncpu2       ! process number  read in RAMSES part files
    Integer(kind=4)                :: ndim2       ! dimension       read in RAMSES part files
    Integer(kind=4)                :: npartloc    ! particle number read in RAMSES part files
    Integer(kind=4)                :: prov, dest  ! provenance and destination process number for p2p MPI communications
    Integer(kind=4)                :: mpistat(Mpi_Status_Size)   ! status of MPI communication
    Integer(kind=PRI)              :: tmplongint              ! temp integer8 variable
    Integer(kind=4)                :: errcode
    Integer(kind=4)                :: grpnb

    Real(kind=SP), allocatable     :: tmpsimple(:)            ! temporary variable for Ramses v2 output
    Real(kind=SP), allocatable     :: tmpsimplebis(:,:)       ! temporary variable for fofcube output
    Real(kind=DP), allocatable     :: tmpdouble(:)            ! temporary variable for Ramses v3 output
    Real(kind=SP), allocatable     :: tmpsendx(:,:)
    Real(kind=SP), allocatable     :: tmpx(:,:)
    Real(kind=SP)                  :: deltasd
    integer(kind=4)                ::nmod
    
    real(8)::xminin,xmaxin,yminin,ymaxin,zminin,zmaxin,zmininbuf,zmaxinbuf
    integer::nxin,nyin,nzin,local_nz,local_z_start_dest
    real(8),allocatable::zminintab(:),zmaxintab(:),zminmaxtmp(:)
    integer::iproc
    integer(8)::nptottmp,nptottmp2
    logical::periodic
    logical::verbose=.false.
    
    periodic=.true.!periodic boundary conditions

    !!!!!parameters to be removed (test)
    Call Mpi_Comm_size(Mpi_Comm_World, procNB, mpierr)
    Call Mpi_Comm_rank  (MPI_Comm_World, procID, mpierr)
    
    pathinput=trim(rep)
    nameinfo=infofile
    namepart=inputfiles
    
    if(ramses_read_part==1)then
       code_index='RA3'
    else if (ramses_read_part==2)then
       code_index='CUB'
    else if (ramses_read_part==3)then
       code_index='RA2'
    else
       print*,'Format not supported ramses_read_part=', ramses_read_part
       Call Mpi_Finalize(mpierr)
       Stop
    endif
    
    grpsize=iogroupsize
    dotimings=.false.

    !!!!!end parameters to be removed (test)


    ! buffer for Mpi_Pack
    ! 3*bit_size(nres)/8 : 3x Integer(kind=4) nprocramses, ndim, nres
    h_length = 3* bit_size(nres)/8
    Allocate (header(0:h_length-1))
    h_pos = 0

    ! Timer initialization
    time0 = Mpi_Wtime()

    ! group directory base name: particles files are dispatched in different 'group' directories
    grpchar = 'group_00001'

    ! Process 0 reads parameters in info file
    If(procID == 0) Then
       ! RAMSES v2 particles files
       If(code_index.eq.'RA2') Then
          Print *,'Reading Ramses v2 output...'
!          Write(Ulog,*) 'Reading Ramses v2 output...'
          ! RAMSES v3 particles files
       Else if(code_index.eq.'RA3') Then
          Print *,'Reading Ramses v3 output...'
!          Write(Ulog,*) 'Reading Ramses v3 output...'
       Else if(code_index.eq.'CUB') Then
          Print *,'Reading fof cube output...'
!          Write(Ulog,*) 'Reading fof cube output...'   
       End If

       ! info file
       ! Use of groups in case a large simulations
       if(grpsize>0) then
          nomfich = trim(pathinput)//'/'//trim(grpchar)//'/'//trim(nameinfo)
       else
          nomfich = trim(pathinput)//'/'//trim(nameinfo)
       endif
       Print *,'Reading RAMSES info file:',trim(nomfich)

       Open(11,file=nomfich, form='formatted', Status='Old', Iostat=errcode)
       If(errcode > 0) Then
          Call EmergencyStop('Error opening '//trim(nameinfo)//' file',errcode)
       End If
       Rewind 11

       ! Read number of particles files, nb of dimensions and level of refinement in the coarse grid (RAMSES is based on AMR)
       Read(11,'(A13,I11)') dumchar,nprocramses
       Read(11,'(A13,I11)') dumchar,ndim
       Read(11,'(A13,I11)') dumchar,lmin

       Close(11)
       
       !Number of files to read is different if fof cube format
       if(code_index=='CUB')nprocramses=nfiles
       
       ! nb of grid point in each dimension = 2^lmin
       nres = 2**lmin

       Write(*,*) 'Number of:'
       Write(*,'(A25,I6)') ' - files for each output:',nprocramses
       Write(*,'(A25,I6)') ' - dimensions:           ',ndim
       Write(*,'(A25,I6)') ' - grid points:          ',nres
!       Write(Ulog,*) 'nb_proc = ',nprocramses,'ndim = ',ndim,'nres = ',nres

       ! Pack these 3 parameters for broadcast
       Call Mpi_Pack(nprocramses, 1, Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack( ndim, 1, Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack( nres, 1, Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)

    End If

    ! Broadcast
    Call Mpi_Bcast(header,h_length,Mpi_Packed,0,Mpi_Comm_World,mpierr)

    ! Other processes unpack the 3 parameters
    If(procID /= 0) Then

       Call Mpi_Unpack(header, h_length, h_pos, nprocramses, 1, Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(header, h_length, h_pos,  ndim, 1, Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(header, h_length, h_pos,  nres, 1, Mpi_Integer, Mpi_Comm_World, mpierr)

    End If

    ! total nb of grid points
    ngrid = int(nres,kind=8)**3

    If(allocated(header)) Deallocate(header)

    If(procID==0) Print *,'Reading positions...'
!    If(procID==0) Write(Ulog,*) 'Reading positions...'

    ! total nb of particles is initialized to 0
    nptot = 0
    ! each process read mynbfile particles files
    mynbfile = nprocramses / procNB
    ! index of first and last particles file read by current process
    firstp  = procID * mynbfile + 1
    lastp   = (procID+1) * mynbfile

    ! Various Ramses proc number taken into account 
    nmod=mod(nprocramses,procNB)
    if (nmod.ne.0)then
       if(procID.le.nmod-1)then
          mynbfile=mynbfile+1
          firstp  = procID * mynbfile + 1
          lastp   = (procID+1) * mynbfile
       else
          firstp  = procID * mynbfile + 1+nmod
          lastp   = (procID+1) * mynbfile+nmod
       endif
    endif

    ! nb of particles in files read by current process
    mynpart = 0

    ! particles will  be distributed "geographicaly" among the process, i.e. at the end of this subroutine
    ! each process will keep only the particles located in a particular subdomain of the whole simulation domain
    ! allocate array containing the nb of particles that will remain on each process at the end of this subroutine (npartv)
    ! and array containing the nb of particles destined to each process but read by current process (npartvloc)
    Allocate(npartv(procNB))
    Allocate(npartvloc(procNB))

    npartv = 0
    npartvloc = 0

    ! nsd = nb of subdomains in each direction, i.e. there are nsd**3 subdomains
    nsd = int(procNB**(1./3.))
    ! deltasd = dimension of a subdomain, i.e. each subdomain is a cube of size deltasd
    deltasd = 1./nsd

    If(procID == 0) Then
       Write(*,*) 'Number of subdomains in each dimension:',nsd
       Write(*,*) 'Size of each subdomain:',deltasd
    End If

    ! definition of the current subdomain by min and max coordinates
    !xmin =  CubeCoord(1)      * deltasd
    !xmax = (CubeCoord(1) + 1) * deltasd
    !ymin =  CubeCoord(2)      * deltasd
    !ymax = (CubeCoord(2) + 1) * deltasd
    !zmin =  CubeCoord(3)      * deltasd
    !zmax = (CubeCoord(3) + 1) * deltasd

    xmin=xminin
    xmax=xmaxin
    ymin=yminin
    ymax=ymaxin
    zmin=zmininbuf
    zmax=zmaxinbuf
    
    allocate(zminintab(procNB))
    allocate(zmaxintab(procNB))
    zminintab=1e30
    zmaxintab=-1e30

    allocate(zminmaxtmp(procNB))
    zminmaxtmp=1e30
    zminmaxtmp(procID+1)=zmin
    Call Mpi_AllReduce(zminmaxtmp,zminintab,procNB,MPI_REAL8,Mpi_Min,Mpi_Comm_World,mpierr)

    zminmaxtmp=-1e30
    zminmaxtmp(procID+1)=zmax
    Call Mpi_AllReduce(zminmaxtmp,zmaxintab,procNB,MPI_REAL8,Mpi_Max,Mpi_Comm_World,mpierr)
    deallocate(zminmaxtmp)
   
    ! loop over the "file/cpu" index characterizing the input particles files to be read by current process
    Do icpu = firstp,lastp
       if(code_index=='CUB')then
          icputmp=icpu-1
       else
          icputmp=icpu
       endif
       Write(ncharcpu(1:5),'(I5.5)') icputmp
       ! the RAMSES files are written in group_XXXXX directories, where XXXXX is the group number
       ! the grpsize first files are in group number 1, the grpsize following are in group number 2, etc...
       !grpnb = (icpu-1)/grpsize + 1
       !Write(grpchar(7:11),'(I5.5)') grpnb
       !nomfich = trim(pathinput)//'/'//trim(grpchar)//'/'//trim(namepart)//trim(ncharcpu)

       ! Use of groups in case a large simulations
       if(grpsize >0)then
          grpnb = (icpu-1)/grpsize + 1
          Write(grpchar(7:11),'(I5.5)') grpnb
          nomfich = trim(pathinput)//'/'//trim(grpchar)//'/'//trim(namepart)//trim(ncharcpu)
       else
          nomfich = trim(pathinput)//'/'//trim(namepart)//trim(ncharcpu)
       endif

       ! read number of files, dimensions and particles
       Open(unit=1,file=nomfich,status='old',form='unformatted',Iostat=errcode)
       If(errcode > 0) Then
          Call EmergencyStop('Error opening '//trim(nomfich)//' file',errcode)
       End If
       if(code_index.ne.'CUB')then
          Read(1,Iostat=errcode) ncpu2
          If(errcode/=0) Print *,'Erreur de lecture: ',errcode
          Read(1,Iostat=errcode) ndim2
          If(errcode/=0) Print *,'Erreur de lecture: ',errcode
       else
          ncpu2=nfiles
          ndim2=3
       endif
       Read(1,Iostat=errcode) npartloc
       If(errcode/=0) Print *,'Erreur de lecture: ',errcode
       Close(1)

       ! check if nb of files and dimensions are consistent with what is written in info file
       If((ncpu2/=nprocramses).Or.(ndim2/=ndim)) Then
          Call EmergencyStop('Files'//trim(nomfich)// ' and '//trim(nameinfo)//' are not consistent for ncpu and/or ndim',2)
       End If

       ! nb of particles read by current process
       mynpart = mynpart + npartloc
    End Do

    ! we use a temporary long integer in case of simulations larger than 1024**3
    tmplongint = mynpart
    ! each mynpart is added -the resulst is nptot, total nb of particles - and broadcasted in a single collective communication
    Call Mpi_AllReduce(tmplongint,nptot,1,MPI_PRI,Mpi_Sum,Mpi_Comm_World,mpierr)

    If(procID == 0) Then
       Write(* ,*)'There are ',nptot,' DM particles'
!       Write(Ulog,*)'There are ',nptot,' DM particles'
    End If

    ! memory allocation for temporary arrays which will contain positions, velocities and id's
    ! of the particles read by current process
    Allocate(tmpx(3,mynpart))

    ! intialization
    tmpx=0.
    mynpart = 0

    ! Begin timings if requested
    If(dotimings) Then
       Call Mpi_Barrier(MPI_comm_world,mpierr)
       timeInt = Mpi_Wtime()
       tInitRead = timeInt - time0
    End If

    ! current process reads its particles files again, but this time it reads position, velocities and id's for each particle
    ! in the file
    Do icpu = firstp,lastp
       !   Write(ncharcpu(1:5),'(I5.5)') icpu
       !   grpnb = (icpu-1)/grpsize + 1
       !   Write(grpchar(7:11),'(I5.5)') grpnb
       !   nomfich = trim(pathinput)//'/'//trim(grpchar)//'/'//trim(namepart)//trim(ncharcpu)

       ! Use of groups in case a large simulations
       if(code_index=='CUB')then
          icputmp=icpu-1
       else
          icputmp=icpu
       endif
       Write(ncharcpu(1:5),'(I5.5)') icputmp
       if(grpsize >0) then
          grpnb = (icpu-1)/grpsize + 1
          Write(grpchar(7:11),'(I5.5)') grpnb
          nomfich = trim(pathinput)//'/'//trim(grpchar)//'/'//trim(namepart)//trim(ncharcpu)
       else
          nomfich = trim(pathinput)//'/'//trim(namepart)//trim(ncharcpu)
       endif

       Open(unit=1,file=nomfich,status='old',form='unformatted')
       if(code_index.ne.'CUB')then
          Read(1) ncpu2
          Read(1) ndim2
       endif
       Read(1) npartloc

       ! There are some differences between RAMSES v2 and RAMSES v3 files
       ramsesV2 : If(code_index.eq.'RA2') Then
          ! allocate an array of real simple precision to read position/velocities
          Allocate(tmpsimple(1:npartloc))
          ! Read positions
          Do idim = 1,ndim
             Read(1) tmpsimple
             ! put all positions in tmpx vector
             tmpx(idim,mynpart+1:mynpart+npartloc) = tmpsimple
          End Do
          Deallocate(tmpsimple)
       End If ramsesV2

       ! There are some differences between RAMSES v2 and RAMSES v3 files
       ramsesV3 : If(code_index.eq.'RA3') Then
          ! allocate an array of real double precision to read position/velocities
          ! these double precision values will be stored in simple precision variables
          Allocate(tmpdouble(1:npartloc))
          ! Read some dummy variables
          Read(1)
          Read(1)
          Read(1)
          Read(1)
          Read(1)

          ! Read positions
          Do idim = 1,ndim
             Read(1) tmpdouble
             ! put all positions in tmpx vector
             tmpx(idim,mynpart+1:mynpart+npartloc) = tmpdouble
          End Do
          Deallocate(tmpdouble)
       End If ramsesV3
       ! End of the differences between the two versions of RAMSES
       
       ! There are some differences between RAMSES v2 and RAMSES v3 files
       fofcube : If(code_index.eq.'CUB') Then
          ! allocate an array of real simple precision to read position/velocities
          Allocate(tmpsimplebis(1:ndim,1:npartloc))
          Read(1)
          Read(1)
          ! Read positions
          Read(1) tmpsimplebis
          ! put all positions in tmpx vector
          tmpx(1:ndim,mynpart+1:mynpart+npartloc) = tmpsimplebis
          Deallocate(tmpsimplebis)
       End If fofcube



       ! Read particle id

       ! loop over the particles read in the current particles file
       Do j = mynpart+1,mynpart+npartloc
          ! particle positions must be >= 0 and <1.0
          ! as the domain is periodic, every position=1.0 is set to 0.0
          Do idim = 1,ndim
             If(tmpx(idim,j)==1.0) tmpx(idim,j) = 0.0
          End Do
          ! in which subdomain is located this particle?
          !n_i = int(tmpx(1,j)/deltasd)
          ! n_j = int(tmpx(2,j)/deltasd)
          !n_k = int(tmpx(3,j)/deltasd)
          !ind = nsd**2 *n_i + nsd*n_j + n_k + 1
          ! this means subdomain number 'ind' has one more particle to take care of
          !npartvloc(ind) = npartvloc(ind)+1
          

          
!          n_i=int(real(tmpx(3,j),kind=dp)*real(nzin,kind=dp)/real(local_nz,kind=dp))
!          ind=n_i+1
!          npartvloc(ind) = npartvloc(ind)+1
          if(periodic) then
             do iproc=1,procNB !!!Might be slow (to be improved if too slow)
                If(tmpx(1,j)>=xmin .and. tmpx(1,j)<xmax .and. &
                     tmpx(2,j)>= ymin .and. tmpx(2,j) < ymax .and. &
                     ((tmpx(3,j)>= zminintab(iproc) .and. tmpx(3,j) < zmaxintab(iproc)).or.&
                     (tmpx(3,j)+1.>= zminintab(iproc) .and. tmpx(3,j)+1. < zmaxintab(iproc)).or.&
                     (tmpx(3,j)-1.>= zminintab(iproc) .and. tmpx(3,j)-1. < zmaxintab(iproc))))then
                   npartvloc(iproc) = npartvloc(iproc)+1
                endif
             end do
          else
             do iproc=1,procNB !!!Might be slow (to be improved if too slow)
                If(tmpx(1,j)>=xmin .and. tmpx(1,j)<xmax .and. &
                     tmpx(2,j)>= ymin .and. tmpx(2,j) < ymax .and. &
                     tmpx(3,j)>= zminintab(iproc) .and. tmpx(3,j) < zmaxintab(iproc))then
                   npartvloc(iproc) = npartvloc(iproc)+1
                endif
             end do
          endif
       End Do


       ! close particles file
       Close(1)
       ! current process has read mynpart particles
       mynpart = mynpart+npartloc
    End Do

    

    ! call writememfilestep(1, 0, 2, 0, ' After allocate tmp array ', 'memall.txt',974, 0.0, 0.0)

    ! add and broadcast npartvloc so that each process knows how many particles it and its friends will have to analize
    Call Mpi_AllReduce(npartvloc,npartv,procNB,Mpi_Integer,Mpi_Sum,Mpi_Comm_World,mpierr)

    
    ! timing if requested
    If(dotimings) Then
       tReadfile = Mpi_Wtime() - timeInt
       timeInt = Mpi_Wtime()
    End If

    ! ------------------------------------------------
    ! sharing out of the particles between processes
    ! ------------------------------------------------

    ! allocate final positions/velocities/id's arrays
    Allocate (xx(3,npartv(procID+1)))

    ! index of the position in the arrays where the data received from other processes must be stored
    recvpoint = 1

    ! loop on processes: each process will received from the other processes the positions/velocities/id's of the particles
    ! it has to analyze
    processus : Do i = 1,procNB - 1
       ! ID of the process current process has to send data to
       dest = mod(procID + i,procNB)
       ! ID of the process current process has to receive data from
       prov = mod(procID + procNB - i, procNB)
       ! Coordinates of dest process in the process grid
       
       ! current process sends and receives the size of the data it will send/receive
       Call Mpi_Isend(npartvloc(dest+1),1,Mpi_Integer,dest,procID,Mpi_comm_world,mpireqs1,mpierr)
       Call Mpi_Irecv(nrecv,1,Mpi_Integer,prov,prov,Mpi_comm_world,mpireqr1,mpierr)
       ! current process determines the boundaries of the region it has to send to dest
       !xmin =  destCoord(1)      * deltasd
       !xmax = (destCoord(1) + 1) * deltasd
       !ymin =  destCoord(2)      * deltasd
       !ymax = (destCoord(2) + 1) * deltasd
       !zmin =  destCoord(3)      * deltasd
       !zmax = (destCoord(3) + 1) * deltasd
       
       local_z_start_dest=dest*local_nz
       xmin = xminin
       xmax = xmaxin
       ymin = yminin
       ymax = ymaxin
       zmin = zminintab(dest+1)  
       zmax = zmaxintab(dest+1)


       ! there are npartvloc(dest+1) particles on the current process that it has to send to dest
       ! if this nb is  greater than 0 then current process proceeds
       If(npartvloc(dest+1)>0) Then
          ! allocation of temporary arrays
          Allocate(tmpsendx(3,npartvloc(dest+1)))

          ind = 1
          
          ! loop over local particles
          Do j=1,mynpart
             ! if the particle is located in the cube of dest, then it is added to the send buffers
             if(periodic) then
                If(tmpx(1,j)>=xmin .and. tmpx(1,j)<xmax .and. &
                     tmpx(2,j)>= ymin .and. tmpx(2,j) < ymax .and. &
                     ((tmpx(3,j)>= zmin .and. tmpx(3,j) < zmax).or.&
                     (tmpx(3,j)+1.>= zmin .and. tmpx(3,j)+1. < zmax).or.&
                     (tmpx(3,j)-1.>= zmin .and. tmpx(3,j)-1. < zmax))) Then
                   tmpsendx(:,ind) = tmpx(:,j)
                   ind=ind+1
                End If
             else
                If(tmpx(1,j)>=xmin .and. tmpx(1,j)<xmax .and. &
                     tmpx(2,j)>= ymin .and. tmpx(2,j) < ymax .and. &
                     tmpx(3,j)>= zmin .and. tmpx(3,j) < zmax) Then
                   tmpsendx(:,ind) = tmpx(:,j)
                   ind=ind+1
                End If
             endif
          End Do

          ! if current process has not found the same nb of particles in dest cube as it had found while
          ! reading from files, then there is a bug somewhere
          If(ind/=npartvloc(dest+1)+1) Then
             Print *,'Erreur dans la repartition des particules'
             Call Mpi_Finalize(mpierr)
             Stop
          End If

       End If

       ! wait for size of send/recv arrays communication to complete
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

       ! send particles to dest process if needed
       If(npartvloc(dest+1)/=0) Then
          Call Mpi_Isend(tmpsendx,3*npartvloc(dest+1),Mpi_Real,dest,procID,Mpi_comm_world,mpireqs1,mpierr)
       End If
       ! receive particles from recv process if needed
       If(nrecv/=0) Then
          Call Mpi_Irecv(xx(1,recvpoint),3*nrecv,Mpi_Real,prov,  prov,Mpi_comm_world,mpireqr1,mpierr)
       End If
       ! received particles are stored in final arrays at position recvpoint in the arrays.
       ! recvpoint is updated after each reception
       recvpoint=recvpoint+nrecv

       ! wait for send communication to complete and deallocate temporary send arrays
       If(npartvloc(dest+1)/=0) Then
          Call Mpi_Wait(mpireqs1,mpistat,mpierr)
          Deallocate(tmpsendx)
       End If
       ! wait for recv communication to complete
       If(nrecv/=0) Then
          Call Mpi_Wait(mpireqr1,mpistat,mpierr)
       End If
       ! end of loop over the processes
    End Do processus

    ! last step: current process copy the particles he read from files which are located in its subdomain to
    ! the final positions/velocities/id's arrays
   ! xmin =  CubeCoord(1)      * deltasd
   ! xmax = (CubeCoord(1) + 1) * deltasd
   ! ymin =  CubeCoord(2)      * deltasd
   ! ymax = (CubeCoord(2) + 1) * deltasd
   ! zmin =  CubeCoord(3)      * deltasd
   ! zmax = (CubeCoord(3) + 1) * deltasd
   
    xmin =xminin
    xmax =xmaxin 
    ymin =yminin 
    ymax =ymaxin
    zmin =zmininbuf 
    zmax =zmaxinbuf 



    ind = 0
    Do j=1,mynpart
       if(periodic)then
          If(tmpx(1,j)>=xmin .and. tmpx(1,j)<xmax .and. &
               tmpx(2,j)>= ymin .and. tmpx(2,j) < ymax .and. &
               ((tmpx(3,j)>= zmin .and. tmpx(3,j) < zmax).or.&
               (tmpx(3,j)+1.>= zmin .and. tmpx(3,j)+1. < zmax).or.&
               (tmpx(3,j)-1.>= zmin .and. tmpx(3,j)-1. < zmax) )) Then
             xx(:,recvpoint+ind) = tmpx(:,j)
             ind = ind+1
          End If
       else
          If(tmpx(1,j)>=xmin .and. tmpx(1,j)<xmax .and. &
               tmpx(2,j)>= ymin .and. tmpx(2,j) < ymax .and. &
               tmpx(3,j)>= zmin .and. tmpx(3,j) < zmax) Then
             xx(:,recvpoint+ind) = tmpx(:,j)
             ind = ind+1
          endif
       endif
    End Do

    ! if current process has not the same number of particles at the end of this subroutine than the number it found
    ! while reading, then there's a bug somewhere
    If(recvpoint+ind /= npartv(procID+1)+1) Then
       Write(tmpstr1,'(I9.9)') recvpoint+ind
       Write(tmpstr2,'(I9.9)') npartv(procID+1)+1
       Call EmergencyStop('Wrong particles number found after send/recv swaps:'//tmpstr1//' ; '//tmpstr2,2)
    End If

    ! call writememfilestep(1, 0, 2, 0, ' After allocate final arrays ', 'memall.txt',974, 0.0, 0.0)

    ! nb of particles for current process = number of particles in its subdomain
    mynpart = npartv(procID+1)



    ! deallocate temporary arrays
    Deallocate(tmpx)
    Deallocate(npartv, npartvloc)
    deallocate(zminintab,zmaxintab)

    ! timings if requested
    If(dotimings) Then
       tTailPart = Mpi_Wtime() - timeInt
       tReadRA = Mpi_Wtime() - time0
    End If


    call mpi_barrier(MPI_comm_world,mpierr)
    
    nptottmp=mynpart
    Call Mpi_AllReduce(nptottmp,nptottmp2,1,Mpi_Integer8,Mpi_Sum,Mpi_Comm_World,mpierr)
    
    if(procID==0)write(*,*)'TOT number of part',nptottmp2
    if(verbose)then
       write(*,*)procid,'mynpart  ',mynpart
       write(*,*)procid,'minmax(x)',minval(xx(1,:)),maxval(xx(1,:))
       write(*,*)procid,'minmax(y)',minval(xx(2,:)),maxval(xx(2,:))
       write(*,*)procid,'minmax(z)',minval(xx(3,:)),maxval(xx(3,:))
    endif
  End Subroutine ramses_lecture
end Module modio






module common_var_ramses
  use grafic_types
  use pathramses
  type info_type
     integer::ncpu
     integer::ndim
     integer::levelmin
     integer::levelmax
     integer::ngridmax
     integer::nstep_coarse
     real(8)::boxlen
     real(8)::time
     real(8)::aexp
     real(8)::omega_m
     real(8)::omega_l
     real(8)::omega_k
     real(8)::omega_b
     real(8)::unit_l
     real(8)::unit_d
     real(8)::unit_t
     real(8)::h0
     real(8),dimension(:),allocatable::bound_key
  end type info_type
  type(info_type)::info 

  integer,dimension(:),allocatable::cpu_list
  logical,dimension(:),allocatable::cpu_read
  integer::ncpu_read
#ifdef DOUB
  real(dp),dimension(:,:,:),allocatable::cubeout
#else
  real(sp),dimension(:,:,:),allocatable::cubeout
#endif
  !Version RAMSES
  integer::versionramses=3
  
contains


  subroutine cubestoread(ncubestot,xxmin,xxmax,yymin,yymax,zzmin,zzmax)
    !Give the list of cube files to read which intersects
    !the volume xxmin,xxmax,yymin,yymax,zzmin,zzmax

    implicit none
    real(8) xxmin,xxmax,yymin,yymax,zzmin,zzmax
    integer::nxcpu,nycpu,nzcpu,ncubestot
    integer::imintoread,imaxtoread,jmintoread,jmaxtoread,kmintoread,kmaxtoread
    integer::i,j,k,num

    nxcpu=int(real(ncubestot)**(1./3.))
    nycpu=int(real(ncubestot)**(1./3.))
    nzcpu=int(real(ncubestot)**(1./3.))
    imintoread=max(0,int(xxmin*nxcpu))
    imaxtoread=min(int(xxmax*nxcpu),nxcpu-1)
    jmintoread=max(0,int(yymin*nycpu))
    jmaxtoread=min(int(yymax*nycpu),nycpu-1)
    kmintoread=max(0,int(zzmin*nzcpu))
    kmaxtoread=min(int(zzmax*nzcpu),nzcpu-1)


    do i=imintoread,imaxtoread
       do j=jmintoread,jmaxtoread
          do k=kmintoread,kmaxtoread
             num=k+j*nycpu+i*nycpu*nxcpu
             if(.not. cpu_read(num+1))then
                ncpu_read=ncpu_read+1
                cpu_list(ncpu_read)=num
                cpu_read(num+1)=.true.
             endif
          end do
       end do
    end do

  end subroutine cubestoread



  subroutine filestoread(xxmin,xxmax,yymin,yymax,zzmin,zzmax)
    !CPU and their associated output files are ordered using Hilbert curve. 
    !This program determine the list of cpu files to read corresponding to 
    !the region defined by the user xr, yr , zr. It's better to have a total 
    !volume corresponding to a 512^3 at most (for memory purpose). 
    !From amr2cube.f90 Romain Teyssier 10/06 with slight modifications.
    !use common_var_ramses,only:info,cpu_list,cpu_read,ncpu_read,nchar,rep
    implicit none


    real(8) xxmin,xxmax,yymin,yymax,zzmin,zzmax

    integer:: ipos,ncpu,ndim,nx2,ny2,nz2,nlevelmax,ngridmax,nstep_coarse,twotondim

    character(LEN=128)::nomfich
    logical::ok
    real::boxlen,t,aexp,hexp,omega_m,omega_l,omega_k,omega_b,scale_l,scale_d,scale_t

    real(8),dimension(:),allocatable::bound_key

    integer::imin,imax,jmin,jmax,kmin,kmax
    integer:: impi,lmax=0,lmin,ilevel,bit_length,maxdom,ndom,i,j
    real(8) dmax,dx,dkey
    real(8),dimension(1)::order_min
    integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
    integer::nx_full,ny_full,nz_full
    real(8),dimension(1:8)::bounding_min,bounding_max


    ncpu=info%ncpu
    nlevelmax=info%levelmax
    allocate(bound_key(0:ncpu))
    bound_key=info%bound_key
    ndim=info%ndim

    !-----------------------
    ! Range --> files to read
    !-----------------------
    if(lmax==0)then
       lmax=info%levelmax
    endif
    lmin=-1
    dmax=max(xxmax-xxmin,yymax-yymin,zzmax-zzmin)
    do ilevel=1,lmax
       dx=0.5d0**ilevel
       if(dx.lt.dmax)exit
    end do
    if (lmin<0)lmin=ilevel
    bit_length=lmin-1
    maxdom=2**bit_length
    imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
    if(bit_length>0)then
       imin=int(xxmin*dble(maxdom))
       imax=imin+1
       jmin=int(yymin*dble(maxdom))
       jmax=jmin+1
       kmin=int(zzmin*dble(maxdom))
       kmax=kmin+1
    endif

    dkey=(dble(2**(nlevelmax+1)/dble(maxdom)))**ndim
    ndom=1
    if(bit_length>0)ndom=8
    idom(1)=imin; idom(2)=imax
    idom(3)=imin; idom(4)=imax
    idom(5)=imin; idom(6)=imax
    idom(7)=imin; idom(8)=imax
    jdom(1)=jmin; jdom(2)=jmin
    jdom(3)=jmax; jdom(4)=jmax
    jdom(5)=jmin; jdom(6)=jmin
    jdom(7)=jmax; jdom(8)=jmax
    kdom(1)=kmin; kdom(2)=kmin
    kdom(3)=kmin; kdom(4)=kmin
    kdom(5)=kmax; kdom(6)=kmax
    kdom(7)=kmax; kdom(8)=kmax

    do i=1,ndom
       if(bit_length>0)then
          call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
       else
          order_min=0.0d0
       endif
       bounding_min(i)=(order_min(1))*dkey
       bounding_max(i)=(order_min(1)+1.0D0)*dkey
    end do


    cpu_min=0; cpu_max=0
    do impi=1,ncpu
       do i=1,ndom
          if (   bound_key(impi-1).le.bounding_min(i).and.&
               & bound_key(impi  ).gt.bounding_min(i))then
             cpu_min(i)=impi
          endif
          if (   bound_key(impi-1).lt.bounding_max(i).and.&
               & bound_key(impi  ).ge.bounding_max(i))then
             cpu_max(i)=impi
          endif
       end do
    end do

    !Informations about files (=cpu) to read. 
    !ncpu_read is the number of cpu to read. cpu_list is the list of the cpu to read. cpu_read tells if a cpu has to be read
    !note that lmax (above) is important too since it is used when one makes the tree
    do i=1,ndom
       do j=cpu_min(i),cpu_max(i)
          if(.not. cpu_read(j))then
             ncpu_read=ncpu_read+1
             cpu_list(ncpu_read)=j
             cpu_read(j)=.true.
          endif
       enddo
    enddo

    deallocate(bound_key)
  end subroutine filestoread


  subroutine rd_info(myid)
    ! use common_var_ramses,only:info,rep,nchar

    implicit none
    integer::myid

    character(len=300)::name,file,filetmp

    integer::info_ncpu,info_ndim,info_levelmin,info_levelmax,info_ngridmax,info_nstep_coarse
    real(8)::info_boxlen,info_time,info_aexp,info_omega_m,info_omega_l,info_omega_k
    real(8)::info_omega_b,info_unit_l,info_unit_d,info_unit_t,info_h0
    real(8),dimension(:),allocatable::bound_key
    real(8)::tmp1,tmp2
    integer::impi,i

    if(ramses_read_part==3)then
       versionramses=2
    else
       versionramses=3
    endif
    
    if(myid==0)write(*,*)'Rd_info'
    if(iogroupsize>0) then
       filetmp=TRIM(rep)//'/group_00001'
    else
       filetmp=TRIM(rep)
    endif
    file=TRIM(filetmp)//'/info_'//trim(nchar)//'.txt'
    ! file=TRIM(rep)//'/info_'//trim(nchar)//'.txt'

    open(unit=1,file=file,status='old',form='formatted')
    read (1,'(a13,I11)')name,info_ncpu
    if(myid==0)write(*,'(a13,I11)')name,info_ncpu

    read (1,'(a13,I11)')name,info_ndim
    if(myid==0)write(*,'(a13,I11)')name,info_ndim

    read (1,'(a13,I11)')name,info_levelmin
    if(myid==0)write(*,'(a13,I11)')name,info_levelmin

    read (1,'(a13,I11)')name,info_levelmax
    if(myid==0)write(*,'(a13,I11)')name,info_levelmax

    read (1,'(a13,I11)')name,info_ngridmax
    if(myid==0)write(*,'(a13,I11)')name,info_ngridmax

    read (1,'(a13,I11)')name,info_nstep_coarse
    if(myid==0)write(*,'(a13,I11)')name,info_nstep_coarse

    read (1,*)

    read (1,'(a13,E23.15)')name,info_boxlen
    if(myid==0)write(*,'(a13,E23.15)')name,info_boxlen

    read (1,'(a13,E23.15)')name,info_time
    if(myid==0)write(*,'(a13,E23.15)')name,info_time

    read (1,'(a13,E23.15)')name,info_aexp
    if(myid==0)write(*,'(a13,E23.15)')name,info_aexp

    if (versionramses>=3)then
       read (1,'(a13,E23.15)')name,info_h0 !v3 only
       if(myid==0)write(*,'(a13,E23.15)')name,info_h0 !v3 only
    endif

    read (1,'(a13,E23.15)')name,info_omega_m
    if(myid==0)write(*,'(a13,E23.15)')name,info_omega_m

    read (1,'(a13,E23.15)')name,info_omega_l
    if(myid==0)write(*,'(a13,E23.15)')name,info_omega_l

    read (1,'(a13,E23.15)')name,info_omega_k
    if(myid==0)write(*,'(a13,E23.15)')name,info_omega_k

    read (1,'(a13,E23.15)')name,info_omega_b
    if(myid==0)write(*,'(a13,E23.15)')name,info_omega_b

    read (1,'(a13,E23.15)')name,info_unit_l
    if(myid==0)write(*,'(a13,E23.15)')name,info_unit_l

    read (1,'(a13,E23.15)')name,info_unit_d
    if(myid==0)write(*,'(a13,E23.15)')name,info_unit_d

    read (1,'(a13,E23.15)')name,info_unit_t
    if(myid==0)write(*,'(a13,E23.15)')name,info_unit_t

    read(1,*)
    read(1,*)
    read(1,*)

    allocate(bound_key(0:info_ncpu))
    do impi=1,info_ncpu
       read(1,'(I8,1X,E23.15,1X,E23.15)')i,tmp1,tmp2
       bound_key(impi-1)=tmp1
       bound_key(impi)=tmp2
    end do

    close(1)

    info%ncpu=info_ncpu
    info%ndim=info_ndim
    info%levelmin=info_levelmin
    info%levelmax=info_levelmax
    info%ngridmax=info_ngridmax
    info%nstep_coarse=info_nstep_coarse
    info%boxlen=info_boxlen
    info%time=info_time
    info%aexp=info_aexp
    info%omega_m=info_omega_m
    info%omega_l=info_omega_l
    info%omega_k=info_omega_k
    info%omega_b=info_omega_b
    info%unit_l=info_unit_l
    info%unit_d=info_unit_d
    info%unit_t=info_unit_t
    allocate(info%bound_key(0:info_ncpu))
    info%bound_key=bound_key
    deallocate(bound_key)
  end subroutine rd_info

  subroutine hilbert3d(x,y,z,order,bit_length,npoint)
    implicit none

    integer     ,INTENT(IN)                     ::bit_length,npoint
    integer     ,INTENT(IN) ,dimension(1:npoint)::x,y,z
    real(8),INTENT(OUT),dimension(1:npoint)::order

    logical,dimension(0:3*bit_length-1)::i_bit_mask
    logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask
    integer,dimension(0:7,0:1,0:11)::state_diagram
    integer::i,ip,cstate,nstate,b0,b1,b2,sdigit,hdigit

    if(bit_length>bit_size(bit_length))then
       write(*,*)'Maximum bit length=',bit_size(bit_length)
       write(*,*)'stop in hilbert3d'
       stop
    endif

    state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
         &   0, 1, 3, 2, 7, 6, 4, 5,&
         &   2, 6, 0, 7, 8, 8, 0, 7,&
         &   0, 7, 1, 6, 3, 4, 2, 5,&
         &   0, 9,10, 9, 1, 1,11,11,&
         &   0, 3, 7, 4, 1, 2, 6, 5,&
         &   6, 0, 6,11, 9, 0, 9, 8,&
         &   2, 3, 1, 0, 5, 4, 6, 7,&
         &  11,11, 0, 7, 5, 9, 0, 7,&
         &   4, 3, 5, 2, 7, 0, 6, 1,&
         &   4, 4, 8, 8, 0, 6,10, 6,&
         &   6, 5, 1, 2, 7, 4, 0, 3,&
         &   5, 7, 5, 3, 1, 1,11,11,&
         &   4, 7, 3, 0, 5, 6, 2, 1,&
         &   6, 1, 6,10, 9, 4, 9,10,&
         &   6, 7, 5, 4, 1, 0, 2, 3,&
         &  10, 3, 1, 1,10, 3, 5, 9,&
         &   2, 5, 3, 4, 1, 6, 0, 7,&
         &   4, 4, 8, 8, 2, 7, 2, 3,&
         &   2, 1, 5, 6, 3, 0, 4, 7,&
         &   7, 2,11, 2, 7, 5, 8, 5,&
         &   4, 5, 7, 6, 3, 2, 0, 1,&
         &  10, 3, 2, 6,10, 3, 4, 4,&
         &   6, 1, 7, 0, 5, 2, 4, 3 /), &
         & (/8 ,2, 12 /) )

    do ip=1,npoint

       ! convert to binary
       do i=0,bit_length-1
          x_bit_mask(i)=btest(x(ip),i)
          y_bit_mask(i)=btest(y(ip),i)
          z_bit_mask(i)=btest(z(ip),i)
       enddo

       ! interleave bits
       do i=0,bit_length-1
          i_bit_mask(3*i+2)=x_bit_mask(i)
          i_bit_mask(3*i+1)=y_bit_mask(i)
          i_bit_mask(3*i  )=z_bit_mask(i)
       end do

       ! build Hilbert ordering using state diagram
       cstate=0
       do i=bit_length-1,0,-1
          b2=0 ; if(i_bit_mask(3*i+2))b2=1
          b1=0 ; if(i_bit_mask(3*i+1))b1=1
          b0=0 ; if(i_bit_mask(3*i  ))b0=1
          sdigit=b2*4+b1*2+b0
          nstate=state_diagram(sdigit,0,cstate)
          hdigit=state_diagram(sdigit,1,cstate)
          i_bit_mask(3*i+2)=btest(hdigit,2)
          i_bit_mask(3*i+1)=btest(hdigit,1)
          i_bit_mask(3*i  )=btest(hdigit,0)
          cstate=nstate
       enddo

       ! save Hilbert key as double precision real
       order(ip)=0.
       do i=0,3*bit_length-1
          b0=0 ; if(i_bit_mask(i))b0=1
          order(ip)=order(ip)+dble(b0)*dble(2)**i
       end do

    end do

  end subroutine hilbert3d



  subroutine filestoreadrectangular(nregions,xxmin,xxmax,yymin,yymax,zzmin,zzmax)
    !Compute the number of cpu to read in a rectangular region of size nx=n,ny=n,nz=n/nregions
    !The idea is to split the rectangular regions in cube of length (xxmax-xxmin)/nregions=
    !(yymax-yymin)/nregions=zzmax-zzmin and then run filestoread wich works better for cubic
    !regions
    ! use common_var_ramses,only:info,cpu_list,cpu_read,ncpu_read
    implicit none
    real(8) xxmin,xxmax,yymin,yymax,zzmin,zzmax
    real(8) xmintmp,xmaxtmp,ymintmp,ymaxtmp,zmintmp,zmaxtmp
    real(8)dx,dy,dz
    integer::i,j,nregions

    dx=(xxmax-xxmin)/nregions
    dy=(yymax-yymin)/nregions
    dz=(zzmax-zzmin)!

    if ((dx.ne.dy).or.(dy.ne.dz)) then
       write(*,*)'Error size dx dy dz are different',dx,dy,dz
       stop
    endif

    zmintmp=zzmin
    zmaxtmp=zzmin+dz

    xmintmp=xxmin
    xmaxtmp=xxmin+dx
    do i=1,nregions
       ymintmp=yymin
       ymaxtmp=yymin+dy
       do j=1,nregions
          call filestoread(xmintmp,xmaxtmp,ymintmp,ymaxtmp,zmintmp,zmaxtmp)
          ymintmp= ymintmp+dy
          ymaxtmp= ymaxtmp+dy
       end do
       xmintmp=xmintmp+dx
       xmaxtmp=xmaxtmp+dx
    end do

  end subroutine filestoreadrectangular


  !=======================================================================
  subroutine title(n,nchar)
    !=======================================================================
    implicit none
    integer::n
    character*5::nchar

    character*1::nchar1
    character*2::nchar2
    character*3::nchar3
    character*4::nchar4
    character*5::nchar5

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

  !================================================================
  !================================================================
  !================================================================
  !================================================================


  subroutine cube2buffer(buffer,nx,ny,local_nz)
    !use common_var_ramses,only:cubeout
    implicit none
    integer::nx,ny,local_nz
#ifdef DOUB
    real(dp), dimension(local_nz*ny*2*(nx/2+1)), intent(out) :: buffer
#else
    real(sp), dimension(local_nz*ny*2*(nx/2+1)), intent(out) :: buffer
#endif

    integer::i,j,k,index,n2x


    n2x=2*(nx/2+1)
    do k=1,local_nz
       do j=1,ny
          do i=1,nx 
             index = ((k-1)*ny+j-1)*n2x+i
             buffer(index) = real(cubeout(i,j,k),kind=4)
          enddo
       enddo
    end do



  end subroutine cube2buffer



  subroutine filestoread_per(ncpu_mpi,xmin,xmax,ymin,ymax,zmin,zmax)
    !Not used because very general. We made a version for a decomposition along z...
    !Files to read taking into account periodic boundary conditions and perpendicular regions along slices along z
    implicit none
    real(8) xxmin,xxmax,yymin,yymax,zzmin,zzmax, xmin,xmax,ymin,ymax,zmin,zmax
    integer::j
    real(8),dimension(0:1)::xr,yr,zr,xrT,yrT,zrT
    real(8),dimension(0:2)::ctre
    integer::ncpu_mpi
    integer::ncpu
    logical::testx,testy,testz
    integer::lmin
    integer::nzcoarse

    !Allocate array to store which files (=cpus) to read
    ncpu_read=0
    ncpu=info%ncpu
    nzcoarse=2**info%levelmin
    if (allocated(cpu_read)) deallocate(cpu_read) 
    allocate(cpu_read(1:ncpu))
    cpu_read=.false.
    if (allocated(cpu_list)) deallocate(cpu_list)
    allocate(cpu_list(1:ncpu))
    cpu_list=-1


    !Main domain
    xxmin=max(0.d0,xmin)
    xxmax=min(1.d0,xmax)

    yymin=max(0.d0,ymin)
    yymax=min(1.d0,ymax)

    zzmin=max(0.d0,zmin)
    zzmax=min(1.d0,zmax)

    lmin=-1
    call filestoreadrectangular(ncpu_mpi,xxmin,xxmax,yymin,yymax,zzmin,zzmax)
    !assume lmin is determined by the main region (and not the one coming from periodic boundary conditions)


    !Periodic boundary conditions
    ctre=(/(xmin+xmax)/2.,(ymin+ymax)/2.,(zmin+zmax)/2./)
    xr=(/xmin,xmax/)
    yr=(/ymin,ymax/)
    zr=(/zmin,zmax/)

    do j=1,7 
       !translate xyzrange   
       call range(j,ctre,xr,yr,zr,xrT,yrT,zrT)

       !Test intersection with [0.,1.]
       testx=(xrT(1)>=0. .and. xrT(1)<=1.) .or. (xrT(0)<=1. .and. xrT(0)>=0.)
       testy=(yrT(1)>=0. .and. yrT(1)<=1.) .or. (yrT(0)<=1. .and. yrT(0)>=0.)
       testz=(zrT(1)>=0. .and. zrT(1)<=1.) .or. (zrT(0)<=1. .and. zrT(0)>=0.)

       if (testx .and. testy .and. testz) then          
          !Add cpu of these regions
          xxmin=max(0.d0,xrT(0))
          xxmax=min(1.d0,xrT(1))
          yymin=max(0.d0,yrT(0))
          yymax=min(1.d0,yrT(1))
          zzmin=max(0.d0,zrT(0))
          zzmax=min(1.d0,zrT(1))
          write(*,*)'zzminzzmaxinper',zzmin,zzmax
          call filestoreadrectangular(ncpu_mpi,xxmin,xxmax,yymin,yymax,zzmin,zzmax)
       endif
    end do

  end subroutine filestoread_per

  !Translate xyzrange to take into account periodic boundary condition
  subroutine range(j,ctre,xr,yr,zr,xrT,yrT,zrT)
    !j,ctre,xr,yr,zr->xrT,yrT,zrT
    implicit none
    integer::j
    real(8),dimension(0:1)::xr,yr,zr,xrT,yrT,zrT
    real(8),dimension(0:2)::ctre

    !translation selon x
    if (j == 1) then
       if (ctre(0) > 0.5) then 
          xrT=xr-1
       else
          xrT=xr+1
       endif
       yrT=yr
       zrT=zr
    endif

    !translation selon y
    if (j == 2) then
       if (ctre(1) > 0.5) then 
          yrT=yr-1
       else
          yrT=yr+1
       endif
       xrT=xr
       zrT=zr
    endif

    !translation selon z
    if (j == 3) then
       if (ctre(2) > 0.5) then 
          zrT=zr-1
       else
          zrT=zr+1
       endif
       yrT=yr
       xrT=xr
    endif


    !translation selon x et y
    if (j == 4) then
       if (ctre(0) > 0.5) then 
          xrT=xr-1
       else
          xrT=xr+1
       endif
       if (ctre(1) > 0.5) then 
          yrT=yr-1
       else
          yrT=yr+1
       endif
       zrT=zr
    endif

    !translation selon x et z
    if (j == 5) then
       if (ctre(0) > 0.5) then 
          xrT=xr-1
       else
          xrT=xr+1
       endif
       if (ctre(2) > 0.5) then 
          zrT=zr-1
       else
          zrT=zr+1
       endif
       yrT=yr
    endif

    !translation selon y et z
    if (j == 6) then
       if (ctre(2) > 0.5) then 
          zrT=zr-1
       else
          zrT=zr+1
       endif
       if (ctre(1) > 0.5) then 
          yrT=yr-1
       else
          yrT=yr+1
       endif
       xrT=xr
    endif

    !translation selon x et y et z
    if (j == 7) then
       if (ctre(0) > 0.5) then 
          xrT=xr-1
       else
          xrT=xr+1
       endif
       if (ctre(1) > 0.5) then 
          yrT=yr-1
       else
          yrT=yr+1
       endif
       if (ctre(2) > 0.5) then 
          zrT=zr-1
       else
          zrT=zr+1
       endif
    endif

  end subroutine range

end module common_var_ramses



