Module modpart2extrema

    Use modconst
    use common_var_ramses
    use modsmooth_field

    type extremaList
        type(extremaList), pointer :: next
#ifdef DOUB          
        real(dp), dimension(6) :: extremum
#else
        real(sp), dimension(6) :: extremum
#endif      
    end type extremaList

contains

    function getcube(i, j, k, cube, zp, zm, local_nz) result(density)
        integer, intent(in) :: i, j, k, local_nz
#ifdef DOUB    
        real(dp), dimension(:,:,:), allocatable :: cube
        real(dp), dimension(:,:), allocatable :: zp, zm
        real(dp) :: density
#else   
        real(sp), dimension(:,:,:), allocatable :: cube
        real(sp), dimension(:,:), allocatable :: zp, zm
        real(sp) :: density
#endif
        !write(*, *)"getcube with ",i,j,k

        if (k < 0) then
            density = zm(i, j)
        elseif (k >= local_nz) then
            density = zp(i, j)
        else
            density = cube(i, j, k)
        endif

        !write(*,*)'ok'

    end function getcube

    subroutine part2extrema(myid, xmin, xmax, ymin, ymax, zmincube, zmaxcube, zminwithdzcoarse, zmaxwithdzcoarse, nx, ny, nz, local_nz, nproc, filterScale, plan, backplan, total_local_size, local_z_start)
        !  use common_var_ramses,only:cpu_list,cpu_read,ncpu_read,nchar,rep,cubeout
        use modvariable
        use modio, only:ramses_lecture
        implicit none
        integer(i8b) :: plan, backplan
        integer(i4b) :: total_local_size, local_z_start
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
        integer :: imin, imax, jmin, jmax, kmin, kmax, lmin
        real(8) :: xxmin, xxmax, yymin, yymax, zzmin, zzmax, dx, dy, dz, deltax
#ifdef DOUB    
        real(dp), dimension(:,:,:), allocatable :: cube
        real(dp), dimension(:,:), allocatable :: zplus, zmoins
        real(dp), dimension(26) :: neighbor_density
        real(dp) local_seuil
        real(dp), dimension(:,:), allocatable :: extremumToSave
        ! Density Histogram variables
        real(dp), dimension(:), allocatable :: density_histo, all_mass
        real(dp), dimension(:,:), allocatable :: all_density_histo
#else
        real(sp), dimension(:,:,:), allocatable :: cube
        real(sp), dimension(:,:), allocatable :: zplus, zmoins
        real(sp), dimension(26) :: neighbor_density
        real(sp) local_seuil
        real(sp), dimension(:,:), allocatable :: extremumToSave
        ! Density Histogram variables
        real(sp), dimension(:), allocatable :: density_histo, all_mass
        real(sp), dimension(:,:), allocatable :: all_density_histo
#endif
        
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
        integer :: local_nz
        integer :: npartin
        real(8), dimension(:), allocatable :: rho
        integer(8), dimension(:), allocatable :: idp, idptmp
        !seuil pour le calcul du minimum
        real :: filterScale
        real :: xminproc, xmaxproc, yminproc, ymaxproc, zminproc, zmaxproc, seuil, seuil_moyen
        integer :: procid
        real(KIND = 4), dimension(:,:), allocatable :: tmp
        integer :: ierr
        integer :: ipart
        integer :: status(MPI_STATUS_SIZE)
        integer :: MPI_Type

        integer :: nb_histo, rank_histo        
        type(extremaList), pointer :: list, curr        

#ifdef DOUB
        MPI_Type = MPI_Double        
#else
        MPI_Type = MPI_Float
#endif

        ! to gather results in multiprocs
        allocate(all_minima(0:nproc - 1))
        allocate(all_seuil0(0:nproc - 1))
        allocate(all_mass(0:nproc - 1))
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
        mp = 1.0 !1.d0/2.d0**(3. * info % levelmin)

        !The region to compute the cube is different from the region to select particule!!!!
        zmin = zminwithdzcoarse
        zmax = zmaxwithdzcoarse

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
        allocate(cube(0:nx - 1, 0:ny - 1, 0:local_nz - 1))
        
        allocate(zplus(0:nx - 1, 0:ny - 1))
        allocate(zmoins(0:nx - 1, 0:ny - 1))

        cube = 0.

        !-----------------------------------------------
        ! Compute projected mass using CIC smoothing
        !----------------------------------------------

        npart2 = mynpart

        !Compute density
        mtot = 0.0d0
        npartin = 0

        l = 0

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

                if (iz + myid * local_nz == 0.and.myid == nproc - 1) then
                    iz = local_nz
                    nbA = nbA + 1
                endif

                if (izp1 == nz.and.myid == 0) then
                    izp1 = 0
                    nbB = nbB + 1
                endif                
                if (iz == nz.and.myid == 0) then
                    iz = 0                    
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

        call MPI_Gather(sum(cube(0:nx - 1, 0:ny - 1, 0:local_nz - 1)), 1, MPI_Type, all_mass, 1, MPI_Type, 0, MPI_comm_world, mpierr)                
        
        if(myid == 0) write(*, *) 'Real mass proc CIC', myid, sum(all_mass)
        
        !-----------------------------------------------------
        ! Compute Minimums in the density field and save them
        !-----------------------------------------------------

        !smoothing the spectrum with a gaussian window
        if (myid == 0)write(*, *) 'Smoothing the field with a gaussian window function'
        call smooth_field(plan, backplan,total_local_size, filterScale,cube,nx,ny,nz,local_nz,local_z_start)

        ! ----------------------
        ! GENERATE DENSITY HISTO
        ! ----------------------
        nb_histo = 1000

        allocate(density_histo(0:nb_histo))

        density_histo = 0

        k = 0
        do iz = 0, local_nz - 1
            do ix = 0, nx - 1
                do iy = 0, ny - 1
                    rank_histo = floor(nb_histo * min(cube(ix, iy, iz), 3.0)/3.0)
                    if (rank_histo < 0) then
                        write(*, *) "error : rank_histo = ", rank_histo
                    else
                        density_histo(rank_histo) = density_histo(rank_histo) + 1
                    endif
                enddo
            enddo
        enddo

        density_histo = density_histo * nb_histo / sum(density_histo(0:nb_histo)) / nproc

        allocate(all_density_histo(0:nb_histo, 0:nproc - 1))
        all_density_histo = 0
        call MPI_Gather(density_histo, nb_histo + 1, MPI_Type, all_density_histo, nb_histo + 1, MPI_Type, 0, MPI_comm_world, mpierr)

        if (myid == 0) then
            density_histo = sum(all_density_histo, 2)
            nomfich = 'data/' // trim(outputfile) // '/' // trim(outputfile) // '_all_after.deus_histo.txt'
            open (unit = 3, file = nomfich)
            do iz = 0, nb_histo
                write(3, *) density_histo(iz)
            enddo
            close(unit = 3)
        endif


        ! ----------------------
        ! EXCHANGING Zp/m BUFFERS
        ! ----------------------

        if (nproc == 1)then
            zmoins = cube(:,:, local_nz - 1)
            zplus = cube(:,:, 0)
        else
            write(*, *) 'pas content, pas content, pas content. Pas content'
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
        list % next => list
        curr => list

        do iz = 0, local_nz - 1
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

                    if (cube(ix, iy, iz) > 0.0 .and. cube(ix, iy, iz) < 1.0) then

                        !computing the minimum of the neighbour treshold
                        seuil = 100.0 - cube(ix, iy, iz)

                        !write(*, *)ix,iy,iz

                        seuil_moyen = 0.0

                        neighbor_density(1) = getcube(ixm1, iym1, izm1, cube, zplus, zmoins, local_nz)
                        neighbor_density(2) = getcube(ixm1, iym1, iz, cube, zplus, zmoins, local_nz)
                        neighbor_density(3) = getcube(ixm1, iym1, izp1, cube, zplus, zmoins, local_nz)
                        neighbor_density(4) = getcube(ixm1, iy, izm1, cube, zplus, zmoins, local_nz)
                        neighbor_density(5) = getcube(ixm1, iy, iz, cube, zplus, zmoins, local_nz)
                        neighbor_density(6) = getcube(ixm1, iy, izp1, cube, zplus, zmoins, local_nz)
                        neighbor_density(7) = getcube(ixm1, iyp1, izm1, cube, zplus, zmoins, local_nz)
                        neighbor_density(8) = getcube(ixm1, iyp1, iz, cube, zplus, zmoins, local_nz)
                        neighbor_density(9) = getcube(ixm1, iyp1, izp1, cube, zplus, zmoins, local_nz)

                        neighbor_density(10) = getcube(ix, iym1, izm1, cube, zplus, zmoins, local_nz)
                        neighbor_density(11) = getcube(ix, iym1, iz, cube, zplus, zmoins, local_nz)
                        neighbor_density(12) = getcube(ix, iym1, izp1, cube, zplus, zmoins, local_nz)
                        neighbor_density(13) = getcube(ix, iy, izm1, cube, zplus, zmoins, local_nz)
                        neighbor_density(14) = getcube(ix, iy, izp1, cube, zplus, zmoins, local_nz)
                        neighbor_density(15) = getcube(ix, iyp1, izm1, cube, zplus, zmoins, local_nz)
                        neighbor_density(16) = getcube(ix, iyp1, iz, cube, zplus, zmoins, local_nz)
                        neighbor_density(17) = getcube(ix, iyp1, izp1, cube, zplus, zmoins, local_nz)

                        neighbor_density(18) = getcube(ixp1, iym1, izm1, cube, zplus, zmoins, local_nz)
                        neighbor_density(19) = getcube(ixp1, iym1, iz, cube, zplus, zmoins, local_nz)
                        neighbor_density(20) = getcube(ixp1, iym1, izp1, cube, zplus, zmoins, local_nz)
                        neighbor_density(21) = getcube(ixp1, iy, izm1, cube, zplus, zmoins, local_nz)
                        neighbor_density(22) = getcube(ixp1, iy, iz, cube, zplus, zmoins, local_nz)
                        neighbor_density(23) = getcube(ixp1, iy, izp1, cube, zplus, zmoins, local_nz)
                        neighbor_density(24) = getcube(ixp1, iyp1, izm1, cube, zplus, zmoins, local_nz)
                        neighbor_density(25) = getcube(ixp1, iyp1, iz, cube, zplus, zmoins, local_nz)
                        neighbor_density(26) = getcube(ixp1, iyp1, izp1, cube, zplus, zmoins, local_nz)

                        !write(*, *) "done"
                        do i = 1, 26
                            local_seuil = neighbor_density(i) - cube(ix, iy, iz)
                            if (local_seuil < seuil) seuil = local_seuil
                            seuil_moyen = seuil_moyen + local_seuil
                        enddo

                        if (seuil > 0.)then
                            allocate(curr % next)
                            curr => curr % next
                            curr % next => null()
                            curr % extremum(1) = (real(ix) + 0.5)/real(nx)
                            curr % extremum(2) = (real(iy) + 0.5)/real(ny)
                            curr % extremum(3) = (real(iz - 1 + local_nz * myid) + 0.5)/real(nz)
                            curr % extremum(4) = cube(ix, iy, iz)
                            curr % extremum(5) = seuil
                            curr % extremum(6) = seuil_moyen

                            nb_minima = nb_minima + 1

                        elseif (seuil == 0.) then
                            nb_seuil0 = nb_seuil0 + 1
                        endif

                    endif
                enddo
            enddo
        enddo

        allocate(curr % next)
        curr => curr % next
        curr % next => null()

        if (myid == 0)write(*, *) 'Saving files...'
        nomfich = 'data/' // trim(outputfile) // '/' // trim(outputfile) // '_' // ncharcpu // '.deus_extrema'
        open(unit = 2, file = trim(nomfich), status = 'unknown', form = 'unformatted')

        curr = list
        allocate(extremumToSave(nb_minima, 6))
        do ix = 1, nb_minima
            curr => curr % next
            extremumToSave(ix,:) = curr % extremum(:)
        enddo

        write(2) nb_minima
        write(2) extremumToSave(:,:)

        close(2)

        write(*, *) 'for proc', myid, 'total minimums', nb_minima, 'with seuil=0', nb_seuil0        
        write(*, *) 'Real mass proc', myid, sum(cube(0:nx - 1, 0:ny - 1, 0:local_nz - 1))/(nx * ny * local_nz)


        call MPI_Gather(nb_minima, 1, MPI_INTEGER, all_minima, 1, MPI_INTEGER, 0, MPI_comm_world, mpierr)
        call MPI_Gather(nb_seuil0, 1, MPI_INTEGER, all_seuil0, 1, MPI_INTEGER, 0, MPI_comm_world, mpierr)

        total_minima = sum(all_minima)

        if (myid == 0) then
            write(*, *) 'TOTAL MINIMA', total_minima, 'WITH SEUIL=0', sum(all_seuil0)

            nomfich = 'data/TOTAL_MIN_COUNT.txt'
            open (unit = 3, file = nomfich, position = 'APPEND')
            write(3, *) trim(outputfile), total_minima, sum(all_seuil0)
            close(unit = 3)

        endif

        Call Mpi_Barrier(MPI_comm_world, mpierr)

        ! ----------------------
        ! GENERATE MINIMA DENSITY HISTO
        ! ----------------------

        density_histo = 0
        all_density_histo = 0

        do ix = 1, nb_minima
            rank_histo = floor(nb_histo * extremumToSave(ix, 4))
            density_histo(rank_histo) = density_histo(rank_histo) + 1
        enddo

        !density_histo = density_histo / total_minima        
        call MPI_Gather(density_histo, nb_histo, MPI_Type, all_density_histo, nb_histo, MPI_Type, 0, MPI_comm_world, mpierr)

        if (myid == 0) then
            density_histo = sum(all_density_histo, 2)
            nomfich = 'data/' // trim(outputfile) // '/' // trim(outputfile) // '_min_after.deus_histo.txt'
            open (unit = 3, file = nomfich)
            do iz = 0, nb_histo
                write(3, *) nb_histo * density_histo(iz) / total_minima
            enddo
            close(unit = 3)
        endif

        !deallocate
        deallocate(all_minima)
        deallocate(all_seuil0)
        deallocate(cube)

    end subroutine part2extrema

end module modpart2extrema
