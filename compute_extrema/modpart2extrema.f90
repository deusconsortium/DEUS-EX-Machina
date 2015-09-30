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
        
        type(extremaList), pointer :: list, curr
        real(sp), dimension(:,:), allocatable :: extremumToSave
        
#ifdef DOUB    
      MPI_Type = MPI_Double
#else
      MPI_Type = MPI_Float
#endif
        
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
        mp = 1.d0/2.d0**(3. * info % levelmin)

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

        !computing the minimums ...
        if (myid == 0)write(*, *) 'Computing minimums in the smoothed density field ...'
        call title(myid, ncharcpu)
        nomfich = 'data/extrema_' // ncharcpu // '.dat'        
        open(unit = 2, file = trim(nomfich), status = 'unknown', form = 'unformatted')

        !k = number of miniums
        !i = number of miniums such as seuil = 0
        k = 0
        i = 0
        
        allocate(list)
        list%next => list        
        curr => list
        
        do iz = 1, usable_local_nz - 2
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

                    !computing the minimum of the neighbour treshold
                    seuil = 100.0/(real(nx)*real(ny)*real(nz)) - Smoothedcube(ix, iy, iz)
                    seuil_moyen = 0.0
                    if (Smoothedcube(ixp1, iy, iz) - Smoothedcube(ix, iy, iz) < seuil) seuil = Smoothedcube(ixp1, iy, iz) - Smoothedcube(ix, iy, iz)
                    seuil_moyen = seuil_moyen + Smoothedcube(ixp1, iy, iz) - Smoothedcube(ix, iy, iz)
                    if (Smoothedcube(ixm1, iy, iz) - Smoothedcube(ix, iy, iz) < seuil) seuil = Smoothedcube(ixm1, iy, iz) - Smoothedcube(ix, iy, iz)
                    seuil_moyen = seuil_moyen + Smoothedcube(ixm1, iy, iz) - Smoothedcube(ix, iy, iz)
                    if (Smoothedcube(ix, iyp1, iz) - Smoothedcube(ix, iy, iz) < seuil) seuil = Smoothedcube(ix, iyp1, iz) - Smoothedcube(ix, iy, iz)
                    seuil_moyen = seuil_moyen + Smoothedcube(ix, iyp1, iz) - Smoothedcube(ix, iy, iz)
                    if (Smoothedcube(ix, iym1, iz) - Smoothedcube(ix, iy, iz) < seuil) seuil = Smoothedcube(ix, iym1, iz) - Smoothedcube(ix, iy, iz)
                    seuil_moyen = seuil_moyen + Smoothedcube(ix, iym1, iz) - Smoothedcube(ix, iy, iz)
                    if (Smoothedcube(ix, iy, izp1) - Smoothedcube(ix, iy, iz) < seuil) seuil = Smoothedcube(ix, iy, izp1) - Smoothedcube(ix, iy, iz)
                    seuil_moyen = seuil_moyen + Smoothedcube(ix, iy, izp1) - Smoothedcube(ix, iy, iz)
                    if (Smoothedcube(ix, iy, izm1) - Smoothedcube(ix, iy, iz) < seuil) seuil = Smoothedcube(ix, iy, izm1) - Smoothedcube(ix, iy, iz)
                    seuil_moyen = seuil_moyen + Smoothedcube(ix, iy, izm1) - Smoothedcube(ix, iy, iz)

                    seuil_moyen = seuil_moyen/6.0
                    !we need the treshold to be positive

                    !				write(2) (real(ix)+0.5)/real(nx),(real(iy)+0.5)/real(ny),(real(iz + local_nz*myid) + 0.5)/real(nz)&
                    !					        &,Smoothedcube(ix,iy,iz)*real(nx)*real(ny)*real(nz),seuil*real(nx)*real(ny)*real(nz)&
                    !					        &,seuil_moyen*real(nx)*real(ny)*real(nz)
                    
                    
                    if (seuil >= 0.)then
                        
                        allocate(curr%next)                        
                        curr => curr%next
                        curr%next => null()
                        curr%extremum(1) = (real(ix) + 0.5)/real(nx)
                        curr%extremum(2) = (real(iy) + 0.5)/real(ny)
                        curr%extremum(3) = (real(iz - 1 + usable_local_nz * myid) + 0.5)/real(nz)
                        curr%extremum(4) = Smoothedcube(ix, iy, iz)*real(nx)*real(ny)*real(nz)
                        curr%extremum(5) = seuil*real(nx)*real(ny)*real(nz)
                        curr%extremum(6) = seuil_moyen*real(nx)*real(ny)*real(nz)
                        
!                        write(2) &
!                            &(real(ix) + 0.5)/real(nx), &
!                            &(real(iy) + 0.5)/real(ny), &
!                            &(real(iz - totalBuffer + usable_local_nz * myid) + 0.5)/real(nz), &
!                            &Smoothedcube(ix, iy, iz)*real(nx)*real(ny)*real(nz), &
!                            &seuil*real(nx)*real(ny)*real(nz), &
!                            &seuil_moyen*real(nx)*real(ny)*real(nz)
                        
                        k = k + 1
                        
                        if (k==53) then
                            write (*,*) curr%extremum
                        endif

                        if (seuil == 0.) i = i + 1
                    endif
                enddo
            enddo
        enddo
        
        
        if (myid == 0)write(*, *) 'Saving files...'        
        curr = list        
        allocate(extremumToSave(k,6))
        do ix = 1, k
            curr => curr%next
            extremumToSave(ix,:) = curr%extremum(:)
            
            if (ix==53) then
                write (*,*) extremumToSave(ix,:)
            endif
        enddo
        
        write(2) k
        write(2) extremumToSave(:,:)
        
        close(2)
        
        write(*, *) 'for proc', myid, 'total minimums', k, 'with seuil=0', i        
        write(*, *) 'Real mass proc', myid, sum(Smoothedcube(0:nx - 1, 0:ny - 1, 1:usable_local_nz - 1))
       
        !deallocate
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