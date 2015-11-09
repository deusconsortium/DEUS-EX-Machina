Module modpart2extrema

    Use modconst
    use common_var_ramses

   type extremaList
        type(extremaList), pointer :: next
        real, dimension(6) :: extremum
    end type extremaList


contains


    subroutine part2extrema(myid, xmin, xmax, ymin, ymax, zmincube, zmaxcube, zminwithdzcoarse, zmaxwithdzcoarse, nx, ny, nz, local_nz, nproc, filterScale)
        !  use common_var_ramses,only:cpu_list,cpu_read,ncpu_read,nchar,rep,cubeout
        use modvariable
        use modio, only:ramses_lecture
        implicit none
        Integer(kind = 4) :: mpierr
        integer :: myid, nproc
        integer :: ncpu, ndim, npart, ngrid, n, i, j, k, l, icpu, ipos, igroup
        integer, allocatable :: all_minima(:), all_seuil0(:)
        integer total_minima
        integer :: nb_minima, nb_seuil0, nb_dens_pos, nb_dens_neg, nb_dens_nul
        integer :: ncpu2, npart2, ndim2, levelmin, levelmax, ilevel
        integer :: nbA, nbB
        integer :: nx, ny, nz, ix, iy, iz, ixp1, iyp1, izp1, ixm1, iym1, izm1, idim, jdim, kdim
        integer :: ixm2, iym2, izm2, ixp2, iyp2, izp2
        real(8) :: mtot, ddx, ddy, ddz, dex, dey, dez, t
        real(8) :: xmin, xmax, ymin, ymax, zmin, zmax, zmincube, zmaxcube, zminwithdzcoarse, zmaxwithdzcoarse
        real(8) :: zminwithdzcoarseold, zmaxwithdzcoarseold
        integer :: imin, imax, jmin, jmax, kmin, kmax, lmin
        real(8) :: xxmin, xxmax, yymin, yymax, zzmin, zzmax, dx, dy, dz, deltax
#ifdef DOUB    
        real(dp), dimension(:,:,:), allocatable :: cube
        !Smoothed cube
        real(dp), dimension(:,:,:), allocatable :: Smoothedcube
        real(dp), dimension(:), allocatable :: tempBufferCube
#else
        real(sp), dimension(:,:,:), allocatable :: cube
        real(sp), dimension(:,:,:), allocatable :: Smoothedcube
        real(sp), dimension(:), allocatable :: tempBufferCube
#endif

        !les parametres pour le filtrage gaussien
        !-----------------
        real(sp), dimension(:,:,:), allocatable :: FilterMatrix
        integer, dimension(:), allocatable :: shiftx
        integer, dimension(:), allocatable :: shifty
        integer, dimension(:), allocatable :: shiftz
        real :: filterScale! such as filtering is on exp(-1/2*(x/filterScale)^2)
        integer :: SmoothBuffer, totalBuffer! the number of buffered cells of +/- z direction used for smoothing
        !-----------------

        real(sp), dimension(:,:), allocatable :: x !v3 
        real(sp), dimension(:), allocatable :: tmpsimple
        real(dp), dimension(:), allocatable :: tmpdouble
        real(sp) :: mp !assumed constant here to save memory
        character(LEN = 1) :: proj = 'z'
        character(LEN = 5) :: ncharcpu, nchargroup
        character(LEN = 80) :: ordering
        character(LEN = 128) :: nomfich, outfich, filetmp
        logical :: ok, ok_part, ok_part2
        integer :: impi, ndom, bit_length, maxdom
        integer, dimension(1:8) :: idom, jdom, kdom, cpu_min, cpu_max
        real(8), dimension(1:8) :: bounding_min, bounding_max
        real(8) :: dkey, order_min, dmax
        real(8), dimension(:), allocatable :: bound_key
        integer :: local_nz, usable_local_nz
        integer :: npartin
        real(8), dimension(:), allocatable :: rho
        integer(8), dimension(:), allocatable :: idp, idptmp
        !seuil pour le calcul du minimum
        real :: xminproc, xmaxproc, yminproc, ymaxproc, zminproc, zmaxproc, seuil, seuil_moyen
        integer :: procid
        real(KIND = 4), dimension(:,:), allocatable :: tmp
        integer :: ierr
        integer :: ipart
        integer :: status(MPI_STATUS_SIZE)
        integer :: MPI_Type
        
        ! Density Histogram variables
        real(sp), dimension(:), allocatable :: density_histo
        integer :: nb_histo, rank_histo
        real(sp), dimension(:,:), allocatable :: all_density_histo
        
        type(extremaList), pointer :: list, curr
        real(sp), dimension(:,:), allocatable :: extremumToSave
        
        real, dimension(26):: neighbor_density
        real local_seuil
        
        
        
#ifdef DOUB    
      MPI_Type = MPI_Double
#else
      MPI_Type = MPI_Float
#endif
      
      ! to gather results in multiprocs
      allocate(all_minima(0:nproc-1))
      call title(myid, ncharcpu)
        
        !Define region
        if (myid == 0)write(*, *) 'Working cube =', nx, ny, local_nz
        idim = 1
        jdim = 2
        kdim = 3
        xxmin = xmin; xxmax = xmax
        yymin = ymin; yymax = ymax
        zzmin = zmincube; zzmax = zmaxcube !!!

!        if (myid == 0) write(*, *) 'input infos =', myid, xmin, xmax, ymin, ymax, zmincube, zmaxcube, zminwithdzcoarse, zmaxwithdzcoarse, nx, ny, nz, local_nz,nproc

        dx = 1.d0/real(nx, kind = 8)
        dy = 1.d0/real(ny, kind = 8)
        dz = 1.d0/real(nz, kind = 8)

        !!! Assume constant mass (to save memory)
        mp = 1.d0  !1.d0/2.d0**(3. * info % levelmin)

        !----------Smoothing part------------------

        !SmoothBuffer ; number of neighbours cells to take into account for the smoothing (in each direction)
        SmoothBuffer = ceiling(3.0 * filterScale)

        !CUSTOM to manage Smooth + minimum search
        totalBuffer = SmoothBuffer + 1
        usable_local_nz = local_nz
        local_nz = local_nz + totalBuffer * 2        
        zminwithdzcoarseold = zminwithdzcoarse
        zzmin = zmincube - float(totalBuffer)/nz
        zminwithdzcoarse = zmincube - (float(totalBuffer) + 1.0)/nz
        zmaxwithdzcoarseold = zmaxwithdzcoarse
        zzmax = zmaxcube + float(totalBuffer)/nz
        zmaxwithdzcoarse = zmaxcube + (float(totalBuffer) + 1.0)/nz
        
        !The region to compute the cube is different from the region to select particule!!!!
        zmin = zminwithdzcoarse
        zmax = zmaxwithdzcoarse

!        if (myid == 0) write(*, *) 'input infos custom =', myid, xmin, xmax, ymin, ymax, zmincube, zmaxcube, zminwithdzcoarse, zmaxwithdzcoarse, nx, ny, nz,local_nz,nproc

        !Allocate the smoothing matrix
        allocate(FilterMatrix(1:2 * SmoothBuffer + 1, 1:2 * SmoothBuffer + 1, 1:2 * SmoothBuffer + 1))
        !Smoothedcube includes the extra +/-1 buffer for extremum computation
        allocate(Smoothedcube(0:nx, 0:ny, 0:local_nz - 2 * SmoothBuffer))
        Smoothedcube = 0.
        !Allocate the shift_xyz tabs 
        allocate(shiftx(1:1 + 2 * SmoothBuffer))
        allocate(shifty(1:1 + 2 * SmoothBuffer))
        allocate(shiftz(1:1 + 2 * SmoothBuffer))

        seuil = 0.
        do i = 1, 2 * SmoothBuffer + 1
            do j = 1, 2 * SmoothBuffer + 1
                do k = 1, 2 * SmoothBuffer + 1
                    !Gaussian filtering exp(-(r-r0)^2/(2*Rf^2))
                    FilterMatrix(i, j, k) = exp(-0.5 * ((1. + SmoothBuffer - i)*(1. + SmoothBuffer - i)&
                    &+(1. + SmoothBuffer - j)*(1. + SmoothBuffer - j)&
                    &+(1. + SmoothBuffer - k)*(1. + SmoothBuffer - k))/(filterScale * filterScale))
                    seuil = seuil + FilterMatrix(i, j, k)
                enddo
            enddo
        enddo
        !normalisation
        FilterMatrix = FilterMatrix/seuil
        
!        do i = 0, 12
!            write (*,*) i, SmoothBuffer(0, 0, i)
!        enddo
!            
!        call exit(0)

        !saving the filtering scale
        open(unit = 10, file = trim('data/info.txt'), status = 'unknown')
        write(10, *) '#gaussian filtering scale in Ramses unit'
        write(10, *) filterScale/real(nx)
        close(10)
        !---------------------------------------------

        !Read one time each file and communicates so that every processor gets particles with zmin zmax

        call ramses_lecture(xmin, xmax, ymin, ymax, zmincube, zmaxcube, zminwithdzcoarse, zmaxwithdzcoarse, nx, ny, nz, local_nz)
        allocate(x(mynpart, 3))
        do ipart = 1, mynpart !SWITCH ORDER BY HAND NOT OPTIMAL FOR MEMORY
            x(ipart, 1) = xx(1, ipart)
            x(ipart, 2) = xx(2, ipart)
            x(ipart, 3) = xx(3, ipart)
        end do
        deallocate(xx)

        !Allocate cube to compute density
        allocate(cube(0:nx, 0:ny, 0:local_nz))
        cube = 0.

        !-----------------------------------------------
        ! Compute projected mass using CIC smoothing
        !----------------------------------------------
        
        npart2 = mynpart

        !Compute density
        mtot = 0.0d0
        npartin = 0
        
        l=0
        
        nbA = 0
        nbB = 0
        
        do i = 1, npart2
            ok_part = (x(i, 1) >= xmin.and.x(i, 1) <= xmax.and. &
            & x(i, 2) >= ymin.and.x(i, 2) <= ymax.and. &
            & ((x(i, 3) >= zmin .and.x(i, 3) <= zmax).or.&
            & (x(i, 3) >= zmin + 1.and.x(i, 3) <= 1.).or.&
            & (x(i, 3) >= 0. .and.x(i, 3) <= zmax - 1.)))

            ok_part2 = (x(i, 1) >= xxmin.and.x(i, 1) < xxmax.and. &
            & x(i, 2) >= yymin.and.x(i, 2) < yymax.and. &
            & x(i, 3) >= zzmin.and.x(i, 3) < zzmax)
            if (ok_part2)then
                npartin = npartin + 1
                if (x(i, 1) == 1.)write(*, *) 'WARNING x=1'
                if (x(i, 2) == 1.)write(*, *) 'WARNING y=1'
                if (x(i, 3) == 1.)write(*, *) 'WARNING z=1'
            endif

            if (ok_part)then

                ddx = (x(i, idim) - xxmin)/dx
                ddy = (x(i, jdim) - yymin)/dy
                ddz = (x(i, kdim) - zzmin)/dz
                ix = ddx
                if (ddx < 0) ix = ix - 1 !because conversion to integer is different when going to negative numbers. We want a floor.
                iy = ddy
                if (ddy < 0) iy = iy - 1
                iz = ddz
                if (ddz < 0) iz = iz - 1
                ddx = ddx - ix
                ddy = ddy - iy
                ddz = ddz - iz
                dex = 1.0 - ddx
                dey = 1.0 - ddy
                dez = 1.0 - ddz

                ixp1 = ix + 1
                iyp1 = iy + 1
                izp1 = iz + 1
                if (ix == nx) ix = 0
                if (iy == ny) iy = 0
                if (ixp1 == nx) ixp1 = 0
                if (iyp1 == ny) iyp1 = 0

                if (iz + myid * usable_local_nz == 0.and.myid == nproc - 1) then
                    iz = local_nz - totalBuffer -1
                    nbA = nbA+1
                endif
                
                if (izp1 == (nz + totalBuffer).and.myid == 0) then
                    izp1 = totalBuffer
                    nbB = nbB+1
                endif

                if (iz >= 0.and.iz < local_nz) then 
                    cube(ix, iy, iz) = cube(ix, iy, iz) + mp * dex * dey * dez
                    cube(ix, iyp1, iz) = cube(ix, iyp1, iz) + mp * dex * ddy * dez
                    cube(ixp1, iy, iz) = cube(ixp1, iy, iz) + mp * ddx * dey * dez
                    cube(ixp1, iyp1, iz) = cube(ixp1, iyp1, iz) + mp * ddx * ddy * dez
                endif
                if (izp1 >= 0.and.izp1 < local_nz) then
                    cube(ix, iy, izp1) = cube(ix, iy, izp1) + mp * dex * dey * ddz
                    cube(ix, iyp1, izp1) = cube(ix, iyp1, izp1) + mp * dex * ddy * ddz
                    cube(ixp1, iy, izp1) = cube(ixp1, iy, izp1) + mp * ddx * dey * ddz
                    cube(ixp1, iyp1, izp1) = cube(ixp1, iyp1, izp1) + mp * ddx * ddy * ddz
                endif
                mtot = mtot + mp
            end if
        end do
                
        !modif CIC^-1 !moved from wrong place RY 13/12/11
        !cube(nx, 0:ny - 1,:) = cube(0, 0:ny - 1,:)
        !cube(:, ny,:) = cube(:, 0,:)

        deallocate(x)

        ! exchange boundaries from extreme procs
        Call Mpi_Barrier(MPI_comm_world,mpierr)
        if(nproc == 1) then            
            allocate(tempBufferCube(0:nx * ny * totalBuffer))
            do ix = 0, nx - 1
                do iy = 0, ny - 1
                    do iz = 0, totalBuffer - 1
                        tempBufferCube(ix * totalBuffer * ny + iy * totalBuffer + iz) = cube(ix, iy, totalBuffer + iz)
                    enddo
                enddo
            enddo
            do ix = 0, nx - 1
                do iy = 0, ny - 1
                    do iz = 0, totalBuffer - 1
                        cube(ix, iy, local_nz - totalBuffer + iz) = tempBufferCube(ix * totalBuffer * ny + iy * totalBuffer + iz)
                    enddo
                enddo
            enddo
            do ix = 0, nx - 1
                do iy = 0, ny - 1
                    do iz = 0, totalBuffer - 1
                        tempBufferCube(ix * totalBuffer * ny + iy * totalBuffer + iz) = cube(ix, iy, local_nz - totalBuffer * 2 + iz)
                    enddo
                enddo
            enddo
            do ix = 0, nx - 1
                do iy = 0, ny - 1
                    do iz = 0, totalBuffer - 1
                        cube(ix, iy, iz) = tempBufferCube(ix * totalBuffer * ny + iy * totalBuffer + iz)
                    enddo
                enddo
            enddo
            deallocate(tempBufferCube)
            
        else

            if (myid == 0) then                
                allocate(tempBufferCube(0:nx*ny*totalBuffer))
                do ix = 0, nx - 1
                    do iy = 0, ny - 1
                        do iz = 0, totalBuffer - 1
                            tempBufferCube(ix*totalBuffer*ny+iy*totalBuffer+iz) = cube(ix, iy, totalBuffer+iz)
                        enddo
                    enddo
                enddo
                Call Mpi_Send(tempBufferCube, nx*ny*totalBuffer, MPI_Type, nproc -1, 0, MPI_comm_world,mpierr)            
                Call Mpi_Recv(tempBufferCube, nx*ny*totalBuffer, MPI_Type, nproc -1, 0, MPI_comm_world,status,mpierr)
                do ix = 0, nx - 1
                    do iy = 0, ny - 1
                        do iz = 0, totalBuffer - 1
                            cube(ix, iy, iz) = tempBufferCube(ix*totalBuffer*ny+iy*totalBuffer+iz)
                        enddo
                    enddo
                enddo            
                deallocate(tempBufferCube)            
            endif
            if(myid == nproc -1) then                
                allocate(tempBufferCube(0:nx*ny*totalBuffer))
                Call Mpi_Recv(tempBufferCube, nx*ny*totalBuffer, MPI_Type, 0, 0, MPI_comm_world,status,mpierr)
                do ix = 0, nx - 1
                    do iy = 0, ny - 1
                        do iz = 0, totalBuffer - 1
                            cube(ix,iy,local_nz-totalBuffer+iz) = tempBufferCube(ix*totalBuffer*ny+iy*totalBuffer+iz)
                        enddo
                    enddo
                enddo            
                do ix = 0, nx - 1
                    do iy = 0, ny - 1
                        do iz = 0, totalBuffer - 1
                            tempBufferCube(ix*totalBuffer*ny+iy*totalBuffer+iz) = cube(ix,iy,local_nz-totalBuffer*2+iz)
                        enddo
                    enddo
                enddo
                Call Mpi_Send(tempBufferCube, nx*ny*totalBuffer, MPI_Type, 0, 0, MPI_comm_world,mpierr)            
                deallocate(tempBufferCube)            
            endif
        endif
      
        !-----------------------------------------------------
        ! Compute Minimums in the density field and save them
        !-----------------------------------------------------

        !smoothing the spectrum with a gaussian window
        if (myid == 0)write(*, *) 'Smoothing the field with a gaussian window function'
        
        do ix = 0, nx - 1
            do i = 1, 1 + 2 * SmoothBuffer
                shiftx(i) = ix + i - 1 - SmoothBuffer
                if (shiftx(i) < 0) shiftx(i) = shiftx(i) + nx
                if (shiftx(i) >= nx) shiftx(i) = shiftx(i) - nx
            enddo
            do iy = 0, ny - 1
                do i = 1, 1 + 2 * SmoothBuffer
                    shifty(i) = iy + i - 1 - SmoothBuffer
                    if (shifty(i) < 0) shifty(i) = shifty(i) + ny
                    if (shifty(i) >= ny) shifty(i) = shifty(i) - ny
                enddo
                do iz = 0, local_nz - 1 - 2 * (SmoothBuffer)
                    do i = 1, 1 + 2 * SmoothBuffer
                        shiftz(i) = iz + i - 1
                    enddo
                    !Computing the smoothed field
                    if (SmoothBuffer == 0) then
                        Smoothedcube(ix, iy, iz) = cube(ix, iy, iz + SmoothBuffer)
                    else
                        Smoothedcube(ix, iy, iz) = 0.0
                        do i = 1, 1 + 2 * SmoothBuffer
                            do j = 1, 1 + 2 * SmoothBuffer
                                do k = 1, 1 + 2 * SmoothBuffer
                                    Smoothedcube(ix, iy, iz) = Smoothedcube(ix, iy, iz) + cube(shiftx(i), shifty(j), shiftz(k)) * FilterMatrix(i, j, k)
                                enddo
                            enddo
                        enddo
                    endif
                enddo
            enddo
        enddo
        
        !Deallocate
        deallocate(cube)
        
        ! ----------------------
        ! GENERATE DENSITY HISTO
        ! ----------------------
        nb_histo = 1000
        
        allocate(density_histo(0:nb_histo-1))
        
        density_histo=0
        
        do iz = 1, usable_local_nz            
            do ix = 0, nx - 1                
                do iy = 0, ny - 1
                    rank_histo = floor(nb_histo*min(Smoothedcube(ix, iy, iz),3.0)/3.0)
                    density_histo(rank_histo) = density_histo(rank_histo) + 1
                enddo
            enddo
        enddo
        
        density_histo = density_histo * nb_histo / sum(density_histo(0:nb_histo-1)) / nproc
        
        allocate(all_density_histo(0:nb_histo-1, 0:nproc-1))
        call MPI_Gather(density_histo, nb_histo, MPI_Type, all_density_histo, nb_histo, MPI_Type, 0, MPI_comm_world, mpierr)
        
        if(myid == 0) then
            density_histo = sum(all_density_histo,2)
            nomfich = 'data/' // trim(outputfile) // '/' // trim(outputfile) // '_all.deus_histo.txt'        
            open (unit = 3, file = nomfich)
            do iz = 0, nb_histo-1
                write(3,*) density_histo(iz)
            enddo
            close(unit=3)
        endif
        
        
        ! ----------------------
        ! COMPUTE MINIMA
        ! ----------------------
                       
        if (myid == 0)write(*, *) 'Computing minimums in the smoothed density field ...'
                
        !k = number of miniums
        !i = number of miniums such as seuil = 0
        nb_minima = 0
        nb_seuil0 = 0
        
        allocate(list)
        list%next => list        
        curr => list
        
        do iz = 1, usable_local_nz
            izm1 = iz - 1
            izp1 = iz + 1
            
            do ix = 0, nx - 1
                !gestion des CL
                ixm1 = ix - 1
                ixp1 = ix + 1
                if (ix == 0) ixm1 = nx - 1
                if (ix == nx - 1) ixp1 = 0
                
                do iy = 0, ny - 1
                    iym1 = iy - 1
                    iyp1 = iy + 1
                    if (iy == 0) iym1 = ny - 1
                    if (iy == ny - 1) iyp1 = 0
                    !using the +/- 1 buffer region for the local extremum computation
                    
                    if (Smoothedcube(ix, iy, iz) > 0.0 .and. Smoothedcube(ix, iy, iz) < 1.0) then

                        !computing the minimum of the neighbour treshold
                        seuil = 100.0 - Smoothedcube(ix, iy, iz)

                        seuil_moyen = 0.0

                        neighbor_density(1) = Smoothedcube(ixm1, iym1, izm1)
                        neighbor_density(2) = Smoothedcube(ixm1, iym1, iz)
                        neighbor_density(3) = Smoothedcube(ixm1, iym1, izp1)
                        neighbor_density(4) = Smoothedcube(ixm1, iy, izm1)
                        neighbor_density(5) = Smoothedcube(ixm1, iy, iz)
                        neighbor_density(6) = Smoothedcube(ixm1, iy, izp1)
                        neighbor_density(7) = Smoothedcube(ixm1, iyp1, izm1)
                        neighbor_density(8) = Smoothedcube(ixm1, iyp1, iz)
                        neighbor_density(9) = Smoothedcube(ixm1, iyp1, izp1)

                        neighbor_density(10) = Smoothedcube(ix, iym1, izm1)
                        neighbor_density(11) = Smoothedcube(ix, iym1, iz)
                        neighbor_density(12) = Smoothedcube(ix, iym1, izp1)
                        neighbor_density(13) = Smoothedcube(ix, iy, izm1)
                        !neighbor_density(1) = Smoothedcube(ix, iy, iz)
                        neighbor_density(14) = Smoothedcube(ix, iy, izp1)
                        neighbor_density(15) = Smoothedcube(ix, iyp1, izm1)
                        neighbor_density(16) = Smoothedcube(ix, iyp1, iz)
                        neighbor_density(17) = Smoothedcube(ix, iyp1, izp1)

                        neighbor_density(18) = Smoothedcube(ixp1, iym1, izm1)
                        neighbor_density(19) = Smoothedcube(ixp1, iym1, iz)
                        neighbor_density(20) = Smoothedcube(ixp1, iym1, izp1)
                        neighbor_density(21) = Smoothedcube(ixp1, iy, izm1)
                        neighbor_density(22) = Smoothedcube(ixp1, iy, iz)
                        neighbor_density(23) = Smoothedcube(ixp1, iy, izp1)
                        neighbor_density(24) = Smoothedcube(ixp1, iyp1, izm1)
                        neighbor_density(25) = Smoothedcube(ixp1, iyp1, iz)
                        neighbor_density(26) = Smoothedcube(ixp1, iyp1, izp1)

                        do i = 1, 26
                            local_seuil = neighbor_density(i) - Smoothedcube(ix, iy, iz)
                            if (local_seuil < seuil) seuil = local_seuil
                            seuil_moyen = seuil_moyen + local_seuil
                        enddo

                        if (seuil > 0.)then
                            allocate(curr%next)
                            curr => curr%next
                            curr%next => null()
                            curr%extremum(1) = (real(ix) + 0.5)/real(nx)
                            curr%extremum(2) = (real(iy) + 0.5)/real(ny)
                            curr%extremum(3) = (real(iz - 1 + usable_local_nz * myid) + 0.5)/real(nz)
                            curr%extremum(4) = Smoothedcube(ix, iy, iz)
                            curr%extremum(5) = seuil
                            curr%extremum(6) = seuil_moyen

                            nb_minima = nb_minima + 1
                                                        
                        elseif (seuil == 0.) then
                            nb_seuil0 = nb_seuil0 + 1                                                                
                        endif
                        
                    endif
                enddo
            enddo
        enddo
                
        allocate(curr%next)
        curr => curr%next
        curr%next => null()
        
        if (myid == 0)write(*, *) 'Saving files...'        
        nomfich = 'data/' // trim(outputfile) // '/' // trim(outputfile) // '_' // ncharcpu // '.deus_extrema'        
        open(unit = 2, file = trim(nomfich), status = 'unknown', form = 'unformatted')
        
        curr = list        
        allocate(extremumToSave(nb_minima,6))
        do ix = 1, nb_minima            
            curr => curr%next            
            extremumToSave(ix,:) = curr%extremum(:)            
        enddo
        
        write(2) nb_minima
        write(2) extremumToSave(:,:)
        
        close(2)
        
        
        !write(*, *) 'for proc', myid, 'total minimums', nb_minima, 'with seuil=0', nb_seuil0        
        !write(*, *) 'Real mass proc', myid, sum(Smoothedcube(0:nx - 1, 0:ny - 1, 1:usable_local_nz - 1))
        
        allocate(all_seuil0(0:nproc-1))
        
        call MPI_Gather(nb_minima, 1, MPI_INTEGER, all_minima, 1, MPI_INTEGER, 0, MPI_comm_world, mpierr)
        call MPI_Gather(nb_seuil0, 1, MPI_INTEGER, all_seuil0, 1, MPI_INTEGER, 0, MPI_comm_world, mpierr)
        
        total_minima = sum(all_minima)
                
        if (myid == 0) then            
            write(*, *) 'TOTAL MINIMA', total_minima, 'WITH SEUIL=0', sum(all_seuil0)
            
            nomfich = 'data/TOTAL_MIN_COUNT.txt'
            open (unit = 3, file = nomfich, position = 'APPEND')
            write(3,*) trim(outputfile), total_minima, sum(all_seuil0)            
            close(unit=3)
            
        endif
        
        Call Mpi_Barrier(MPI_comm_world,mpierr)
        
        ! ----------------------
        ! GENERATE MINIMA DENSITY HISTO
        ! ----------------------
                
        density_histo=0
        all_density_histo=0
        
        do ix = 1, nb_minima
            rank_histo = floor(nb_histo*extremumToSave(ix,4))            
            density_histo(rank_histo) = density_histo(rank_histo) + 1            
        enddo
        
        !density_histo = density_histo / total_minima
        
        call MPI_Gather(density_histo, nb_histo, MPI_Type, all_density_histo, nb_histo, MPI_Type, 0, MPI_comm_world, mpierr)
        
        if(myid == 0) then
            density_histo = sum(all_density_histo,2)
            nomfich = 'data/' // trim(outputfile) // '/' // trim(outputfile) // '_min.deus_histo.txt'        
            open (unit = 3, file = nomfich)
            do iz = 0, nb_histo -1
                write(3,*) nb_histo * density_histo(iz) / total_minima
            enddo
            close(unit=3)
        endif
                                
        !deallocate
        deallocate(all_minima)
        deallocate(all_seuil0)
        deallocate(Smoothedcube)
        deallocate(FilterMatrix)
        deallocate(shiftx)
        deallocate(shifty)
        deallocate(shiftz)

        !CUSTOM to manage Smooth + minimum search
        local_nz = local_nz - (SmoothBuffer + 1) * 2
        zminwithdzcoarse = zminwithdzcoarseold
        zmaxwithdzcoarse = zmaxwithdzcoarseold

    end subroutine part2extrema

end module modpart2extrema
