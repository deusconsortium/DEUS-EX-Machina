module modsmooth_field

    use grafic_types
    use transform
    use grafic_io
    implicit none

contains

    subroutine smooth_field(plan, backplan, total_local_size, filterScale, input, nx, ny, nz, local_nz&
        &, local_z_start)

        ! Arguments
#ifdef DOUB
        real(dp), dimension(:,:,:), allocatable :: input
        real(dp), dimension(local_nz * ny * 2 * (nx/2 + 1)) :: buffer
#else
        real(sp), dimension(:,:,:), allocatable :: input
        real(sp), dimension(local_nz * ny * 2 * (nx/2 + 1)) :: buffer
#endif
        
        integer(i4b), intent(in) :: local_nz, local_z_start, nx, ny, nz, total_local_size
        integer(i8b), intent(in) :: plan, backplan
        real, intent(in) :: filterScale
        ! Local stuff

        integer(i8b) :: index
        real(dp) :: fact, Rs
        integer(i4b) :: k1, k2, k3, i, j, k, n2x
        Integer(kind = 4) :: myid, ierr
        real(dp) :: dx, dk1, dk2, dk3, d3k
        integer :: n12, n22, n32, n23, n2p1
        real(dp) :: akmax, ak1, ak2, ak3, ak23, ak33, akk, ak
        
        call mpi_comm_rank(mpi_comm_world, myid, ierr)
                
        n2x = 2 * (nx/2 + 1)
        do k = 1, local_nz
            do j = 1, ny
                do i = 1, nx
                    index = ((k - 1) * ny + j - 1) * n2x + i                    
                    !buffer(index) = real(input(i, j, k), kind = 4)
                    buffer(index) = input(i-1,j-1,k-1)                    
                enddo
            enddo
        end do
        
        fact = (1.0_dp * nx) * ny * nz
        buffer = buffer/fact

        !Synchronize after reading
        
        call mpi_barrier(mpi_comm_world, ierr)

        
        call fft_mpi(plan, buffer, total_local_size)
        
        n12 = nx/2
        n22 = ny/2
        n32 = nz/2
        n23 = ny * nz
        n2p1 = 2 * (n12 + 1)

        dx = 1.d0/real(nx, kind = 8)
        dk1 = 2.0 * PI/(nx * dx)
        dk2 = 2.0 * PI/(ny * dx)
        dk3 = 2.0 * PI/(nz * dx)

        d3k = dk1 * dk2 * dk3
        akmax = 2.0 * PI/dx

        !guess ...
        Rs = filterScale/nx
        !~ 
        do k3 = 1, local_nz
            ak3 = (k3 + local_z_start - 1) * dk3
            if (k3 + local_z_start > n32) ak3 = ak3 - akmax
            ak33 = ak3 * ak3
            do k2 = 1, ny
                ak2 = (k2 - 1) * dk2
                if (k2 > n22) ak2 = ak2 - akmax
                ak23 = ak2 * ak2 + ak33
                do k1 = 1, n12 ! Complex field, treat Nyquist (n12+1) separately
                    ak1 = (k1 - 1) * dk1
                    akk = ak1 * ak1 + ak23
                    ak = sqrt(akk)

                    index = int((k3 - 1) * ny + k2 - 1, kind = i8b) * n2p1 + 2 * k1 - 1

                    !smoothing here
                    if (akk > 0.0) then
                        buffer(index) = buffer(index) * exp(-0.5 * (akk * Rs * Rs))
                        buffer(index + 1) = buffer(index + 1) * exp(-0.5 * (akk * Rs * Rs))
                    endif
                    !~              buffer(index) = buffer(index)*1.
                    !~              buffer(index + 1) = buffer(index + 1)*1.

                enddo
                ak1 = 0.5 * akmax
                akk = ak1 * ak1 + ak23
                ak = sqrt(akk)

                index = int((k3 - 1) * ny + k2 - 1, kind = i8b) * n2p1 + nx + 1

                !~           !and here
                if (akk > 0.0) then
                    buffer(index) = buffer(index) * exp(-0.5 * (akk * Rs * Rs))
                    buffer(index + 1) = buffer(index + 1) * exp(-0.5 * (akk * Rs * filterScale))
                endif
                !~           buffer(index) = buffer(index)*1.
                !~           buffer(index + 1) = buffer(index + 1)*1.

            enddo
        enddo

        !turning into real space

        call mpi_barrier(mpi_comm_world, ierr)
        call fft_mpi(backplan, buffer, total_local_size)
        call mpi_barrier(mpi_comm_world, ierr)

        !transforming buffer back to input tab
        n2x = 2 * (nx/2 + 1)
        do k = 1, local_nz
            do j = 1, ny
                do i = 1, nx
                    index = ((k - 1) * ny + j - 1) * n2x + i
                    input(i-1, j-1, k-1) = buffer(index)
                enddo
            enddo
        end do

    end subroutine smooth_field
end module modsmooth_field
