
!--------------------------------------------------------------------
!
! The aim of this program is to compute a normalized power spectrum
! from (a collection of) grafic files.
! The output of this program is supposed to be directly comparable
! to the "power.dat" file generated by the grafic code.
!  
!-------------------------------------------------------------------- 

program compute_extrema

    !  use grafic_types
    !  use grafic_io
    !  use transform
    use modpart2extrema    
    use common_var_ramses
    use modio
    use modmem
    implicit none
#ifdef ADD1US
#define  rfftw3d_f77_create_plan  rfftw3d_f77_create_plan_
#define  rfftwnd_f77_destroy_plan rfftwnd_f77_destroy_plan_
    #define  rfftwnd_f77_one_real_to_complex rfftwnd_f77_one_real_to_complex_
#define  rfftw3d_f77_mpi_create_plan  rfftw3d_f77_mpi_create_plan_
#define  rfftwnd_f77_mpi_destroy_plan rfftwnd_f77_mpi_destroy_plan_
    #define  rfftwnd_f77_mpi rfftwnd_f77_mpi_
#define  rfftwnd_f77_mpi_local_sizes rfftwnd_f77_mpi_local_sizes_
#endif
    
  ! Code identification
    character(len = *), parameter :: code = 'COMPUTE EXTREMA'
    character(len = *), parameter :: version = '1.0'

    ! Namelist (parameters) declaration
    integer(i4b), parameter :: nfilesmax = 20
    character(len = filenamelen) :: namelistfile
    logical(lgt) :: hanning = .false., shot = .false., count_out = .false.
    integer(i4b) :: cic = 0
    real :: h0
    real :: filterScale! such as filtering is on exp(-1/2*(x/FilterScale)^2)    
    integer :: nfact
    character(len = 128) :: filename

    namelist /parameters/ nfact, nfiles, inputfiles, infofile, outputfile, count_out, ramses_read_part, rep, nchar, h0, iogroupsize,filterScale

    ! Mpi stuff
    integer(i4b) :: myid, nproc, ierr

    ! Shift variable
    real(dp), dimension(3) :: shift

    ! Miscellaneous
    integer(i4b) :: narg
    integer(i4b) :: i, j, k, ifiles
    integer(i4b) :: iargc

    integer(i8b) :: plan
#ifdef DOUB
    real(dp), allocatable, dimension(:) :: buffer
#else
    real(sp), allocatable, dimension(:) :: buffer
#endif
    integer(i4b) :: nx, ny, nz, npow
    integer(i4b) :: local_nz, local_ny, total_local_size
    integer(i4b) :: local_z_start, local_y_start
    real(dp) :: fact
    character(len = filenamelen) :: outfile = ''
    character(len = 4) :: dummy = ''

    !For ramses read
    real(8) xxmin, xxmax, yymin, yymax, zzmin, zzmax
    real(8) dzbuffer, dzcoarse
    real(8) zminwithdzbuffer, zminwithdzcoarse, zmaxwithdzbuffer, zmaxwithdzcoarse, zmincube, zmaxcube
    integer :: ii, jj, ncpu
    real(8) :: nfactcoarse, nfactbuffer
    integer :: nregions
    integer :: nfilesreadtot
    logical :: periodic = .true.

    real(8) :: tfull1, tfull2, tfileselec1, tfileselec2, tfft1, tfft2, tpower1, tpower2, tpart2cube1, tpart2cube2
    real :: memfileselec2, memfft2, mempower2, mempart2cube2, real_mem, membuffer
    !------------------------------------------------------------------

    !Initialize mpi
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world, myid, ierr)
    call mpi_comm_size(mpi_comm_world, nproc, ierr)

    tfull1 = Mpi_Wtime()
    tfileselec1 = Mpi_Wtime()

    ! Announce program, and read namelist
    if (myid == 0) then
        print *, ' '
        print *, '         ****************************************************** '
        print *, '                           ' // code // ' ' // version
        print *, '                      Compute extrema from cube files        '
        print *, '                            Based on powergrid        '
        print *, '                   Using currently ', nproc, ' processors        '
        print *, '         ****************************************************** '
        print *, ' '
    endif

    narg = iargc()
    IF (narg .LT. 1)THEN
        write(*, *) 'You should type: powergrid parameters.nml'
        write(*, *) 'File parameters.nml should contain a parameter namelist'
        call mpi_finalize(ierr)
        stop
    END IF
    CALL getarg(1, namelistfile)
    namelistfile = TRIM(namelistfile)
    open(1, file = namelistfile)
    read(1, nml = parameters)
    close(1)

    ! Output namelist content
    if (myid == 0) then
        print*, ' '
        print*, ' COMPUTE EXTREMA has been called with the following parameters:'
        print*, 'nfact(ignored for mpgrafic files)=', nfact
        print*, 'cic=', cic, ' hanning=', hanning, ' count_out', count_out
        print*, 'h0=', h0
        print*, 'filterScale', filterScale        
        print*, 'ramses_read_part(not for mpgrafic)', ramses_read_part
        print*, 'rep', trim(rep)
        print*, 'Name of info file ', trim(infofile)
        print*, 'Name of file ', trim(inputfiles)
        print*, 'Snapshot number (not for mpgrafic)', nchar
        print*, 'Number of files (ignored for ramses) : ', nfiles
        print*, 'Number of files per group :', iogroupsize
        print*, 'Output file ', trim(outputfile)
    endif

    if (trim(inputfiles) == trim(outputfile))then
        print*, 'Name of inpufiles and outputfile are too similar, risk of erasing the inputfiles, please choose different name'
        print*, trim(inputfiles), trim(outputfile)
        call mpi_finalize(ierr)
        stop
    endif

    !Start processing file(s)

    if (myid == 0) then
        print*, ' '
        print*, ' ************************************************************** '
        print*, ' Processing file(s) in ', trim(rep)
    endif

    ! Get file header and content

    if (myid == 0) then
        write(*, *) 'READ RAMSES PARTICLES FILES MODE'
        write(*, *) 'for ramses version', versionramses
    endif
    call rd_info(myid)
    nx = 2**info % levelmin * nfact
    ny = 2**info % levelmin * nfact
    nz = 2**info % levelmin * nfact
    if (myid == 0.and.nfact .ne. 1) then
        write(*, *) 'NX,NY,NZ FFT=', nfact, 'times NX,NY,NZ RAMSES COARSE'
    endif

    npow = nx
    fact = (1.0_dp * nx) * ny * nz
    ! Get local buffer sizes, and prepare fft forward plan
    call rfftw3d_f77_mpi_create_plan(plan, mpi_comm_world, nx, ny&
    &, nz, fftw_real_to_complex, fftw_estimate)
    call rfftwnd_f77_mpi_local_sizes(plan, local_nz, local_z_start&
    &, local_ny, local_y_start, total_local_size)

#ifdef DOUB
    membuffer = total_local_size * 8./1024./1024.
#else        
    membuffer = total_local_size * 4./1024./1024.
#endif


    !Define regions of interest for x and y direction
    xxmin = 0.
    xxmax = 1.
    yymin = 0.
    yymax = 1.

    !Define z range where to compute density by this proc [zmincube,zmaxcube]
    !as well as additionnal regions for selecting particles for CIC [zminwithdzcoarse,zmaxwithdzcoarse]

    nfactcoarse = 1. !You can choose 1.,1./2,1./4... 1 is safe... Be careful to stay below or equal to one
    dzcoarse = 1.d0/real(nz, kind = dp)/nfactcoarse
    if (nfactcoarse > 1.)then
        write(*, *) 'Please adjust nfactcoarse in powergrid (should be >=1)', nfactcoarse
        call mpi_finalize(ierr)
        stop
    endif

    zmincube = real(local_z_start, kind = dp)/real(nz, kind = dp)
    zmaxcube = real(local_z_start + local_nz, kind = dp)/real(nz, kind = dp)
    zminwithdzcoarse = zmincube - dzcoarse
    zmaxwithdzcoarse = zmaxcube + dzcoarse

    tfileselec2 = Mpi_Wtime()

    call getmem(real_mem)
    call MPI_ALLREDUCE(real_mem, memfileselec2, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr)

    tpart2cube1 = Mpi_Wtime()
    !Compute density using CIC
    if (myid == 0)write(*, *) 'part2extrema: Read particles and compute density in slices using CIC'
  
    call part2extrema(myid, xxmin, xxmax, yymin, yymax, zmincube, zmaxcube, zminwithdzcoarse, zmaxwithdzcoarse, nx, ny, nz, local_nz, nproc, filterScale)
    call getmem(real_mem)
    call MPI_ALLREDUCE(real_mem,mempart2cube2 ,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
    
    tpart2cube2 = Mpi_Wtime()        
    tfull2 = Mpi_Wtime()
        
    if (myid == 0)then
        print*, ''
        print*, 'CPU TIME'
        print*, 'tfull(s)          =', tfull2 - tfull1
        print*, ''
        print*, 'tfileselec(s)     =', tfileselec2 - tfileselec1
        print*, 'tpart2extrema+read(s)=', tpart2cube2 - tpart2cube1
    endif

    if (myid == 0)then
        print*, ''
        print*, 'MEMORY USAGE'
        print*, 'End of fileselec'
        call writemem(memfileselec2)
        print*, 'End of part2extrema+read'
        call writemem(mempart2cube2)        
        print*, 'to be compared to buffer size (MB)', membuffer
    endif

    if (myid == 0) print*, 'Run completed'
    call mpi_finalize(ierr)


end program compute_extrema
