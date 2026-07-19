! wave_mpi.f90
!
! solves the one-dimensional linear advection equation
!     dphi/dt = -c0*dphi/dx
! on a periodic domain using mpi domain decomposition, a fourth-order
! central spatial derivative, and third-order runge-kutta time integration.
!
! this mpi baseline reports timing and numerical error against the
! analytical sine-wave solution. rank 0 gathers and writes the comparison
! result after the parallel solve completes.
!
! output: wave_solution_mpi_fortran.h5
! author: dekeract01, 2023
! university of southampton

program wave_equation_1d_mpi
    use mpi
    use hdf5
    implicit none
    
    ! Simulation parameters
    real(8), parameter :: c0 = 0.5d0          ! Wave speed
    real(8), parameter :: dt = 0.0000001d0    ! Time step
    integer, parameter :: niter = 1000        ! Number of iterations
    integer, parameter :: nx_global = 2000000 ! Total grid points
    real(8), parameter :: dx = 1.0d0/nx_global
    real(8), parameter :: pi = 4.0d0*atan(1.0d0)
    
    ! MPI variables
    integer :: ierr, rank, nprocs
    integer :: nx_local, i_start, i_end
    integer :: left_rank, right_rank
    integer :: status(MPI_STATUS_SIZE)
    
    ! Arrays
    real(8), allocatable :: x0(:), x_global(:), phi(:), phi_global(:)
    real(8), allocatable :: phi1(:), phi2(:), k1(:), k2(:), k3(:), dphi_dx(:)
    real(8) :: phi_left, phi_right  ! Ghost cells for boundary exchange
    
    ! Variables
    integer :: i, n, nx_per_proc, remainder
    real(8) :: t_final, phi_min, phi_max, max_error, l2_error
    real(8) :: start_time, end_time, total_time
    real(8), allocatable :: phi_exact(:), error(:)
    logical :: nan_detected
    
    ! Initialize MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    
    ! Domain decomposition
    nx_per_proc = nx_global / nprocs
    remainder = mod(nx_global, nprocs)
    
    ! Distribute points (give remainder to last process)
    if (rank < nprocs - 1) then
        nx_local = nx_per_proc
        i_start = rank * nx_per_proc + 1
    else
        nx_local = nx_per_proc + remainder
        i_start = rank * nx_per_proc + 1
    end if
    i_end = i_start + nx_local - 1
    
    ! Set neighbor ranks for periodic boundaries
    left_rank = rank - 1
    right_rank = rank + 1
    if (rank == 0) left_rank = nprocs - 1
    if (rank == nprocs - 1) right_rank = 0
    
    ! Allocate local arrays
    allocate(x0(nx_local))
    allocate(phi(nx_local))
    allocate(phi1(nx_local))
    allocate(phi2(nx_local))
    allocate(k1(nx_local))
    allocate(k2(nx_local))
    allocate(k3(nx_local))
    allocate(dphi_dx(nx_local))
    
    ! Print header (only rank 0)
    if (rank == 0) then
        print *, '1D Wave Equation Simulation (Fortran 90 + MPI)'
        print *, repeat('=', 50)
        print '(A,I0)', 'Total grid points: ', nx_global
        print '(A,I0)', 'Number of MPI processes: ', nprocs
        print '(A,I0)', 'Points per process: ', nx_per_proc
        print '(A,F10.6)', 'dx: ', dx
        print '(A,F10.6)', 'dt: ', dt
        print '(A,F10.6)', 'Wave speed c0: ', c0
        print '(A,F10.8)', 'CFL number: ', c0*dt/dx
        print '(A,I0)', 'Number of iterations: ', niter
        print '(A,F10.8)', 'Total time: ', niter*dt
        print *, repeat('=', 50)
    end if
    
    ! Grid setup (local portion)
    do i = 1, nx_local
        x0(i) = (i_start + i - 2) * dx
    end do
    
    ! Initial condition: phi = sin(2*pi*x)
    do i = 1, nx_local
        phi(i) = sin(2.0d0 * pi * x0(i))
    end do
    
    ! Synchronize before starting
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    
    ! Time stepping loop
    if (rank == 0) then
        print *, ''
        print *, 'Starting simulation...'
    end if
    
    start_time = MPI_Wtime()
    
    do n = 1, niter
        ! Exchange boundary data with neighbors
        call exchange_boundaries(phi, phi_left, phi_right, nx_local, &
                                 left_rank, right_rank, rank, nprocs)
        
        ! RK3 time step
        call rk3_step_mpi(phi, phi1, phi2, k1, k2, k3, dphi_dx, &
                         phi_left, phi_right, c0, dx, dt, nx_local, &
                         left_rank, right_rank, rank, nprocs)
        
        ! NaN check and progress report (every 100 steps)
        if (mod(n, 100) == 0) then
            nan_detected = .false.
            do i = 1, nx_local
                if (phi(i) /= phi(i)) then
                    nan_detected = .true.
                    exit
                end if
            end do
            
            if (nan_detected) then
                print *, 'NaN detected at iteration ', n, ' on rank ', rank
                call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            
            if (rank == 0) then
                end_time = MPI_Wtime()
                phi_min = minval(phi)
                phi_max = maxval(phi)
                print '(A,I4,A,I4,A,F6.2,A,F9.6,A,F9.6,A)', &
                    'Iter ', n, '/', niter, ' | Time: ', end_time - start_time, &
                    's | phi: [', phi_min, ', ', phi_max, ']'
            end if
        end if
    end do
    
    end_time = MPI_Wtime()
    total_time = end_time - start_time
    
    ! Gather results to rank 0
    allocate(phi_global(nx_global))
    call gather_solution(phi, phi_global, nx_local, nx_global, rank, nprocs)
    
    if (rank == 0) then
        print *, repeat('=', 50)
        print *, 'Simulation complete!'
        print '(A,F10.8,A)', 'Total time: ', total_time, ' seconds'
        print '(A,F10.8,A)', 'Time per iteration: ', total_time/niter*1000.0d0, ' ms'
        print '(A,F10.1)', 'Iterations per second: ', niter/total_time
        print *, repeat('=', 50)
        
        ! Analytical solution
        allocate(x_global(nx_global))
        allocate(phi_exact(nx_global))
        allocate(error(nx_global))
        
        t_final = niter * dt
        do i = 1, nx_global
            x_global(i) = (i-1)*dx
            phi_exact(i) = sin(2.0d0 * pi * (x_global(i) - c0*t_final))
            error(i) = abs(phi_global(i) - phi_exact(i))
        end do
        
        max_error = maxval(error)
        l2_error = sqrt(sum(error**2) / nx_global)
        
        print *, ''
        print *, 'Numerical vs Analytical:'
        print '(A,ES12.6)', 'Max error: ', max_error
        print '(A,ES12.6)', 'L2 error: ', l2_error
        
        call write_hdf5_solution(x_global, phi_global, phi_exact, error)
        print *, ''
        print *, 'results written to wave_solution_mpi_fortran.h5'

        deallocate(x_global, phi_exact, error)
    end if
    
    ! Cleanup
    deallocate(x0, phi, phi1, phi2, k1, k2, k3, dphi_dx, phi_global)
    
    call MPI_Finalize(ierr)
    
contains

    subroutine write_hdf5_solution(x, phi_numerical, phi_exact, error)
        implicit none
        real(8), intent(in) :: x(:), phi_numerical(:), phi_exact(:), error(:)
        integer(hid_t) :: file_id
        integer :: hdf5_error

        call h5open_f(hdf5_error)
        call h5fcreate_f('wave_solution_mpi_fortran.h5', H5F_ACC_TRUNC_F, &
                         file_id, hdf5_error)
        call check_hdf5_error(hdf5_error, 'could not create output file')

        call write_dataset(file_id, 'x', x)
        call write_dataset(file_id, 'phi_numerical', phi_numerical)
        call write_dataset(file_id, 'phi_exact', phi_exact)
        call write_dataset(file_id, 'error', error)

        call h5fclose_f(file_id, hdf5_error)
        call check_hdf5_error(hdf5_error, 'could not close output file')
        call h5close_f(hdf5_error)
        call check_hdf5_error(hdf5_error, 'could not close hdf5 library')
    end subroutine write_hdf5_solution

    subroutine write_dataset(file_id, name, values)
        implicit none
        integer(hid_t), intent(in) :: file_id
        character(len=*), intent(in) :: name
        real(8), intent(in) :: values(:)
        integer(hid_t) :: dataspace_id, dataset_id
        integer(hsize_t) :: dimensions(1)
        integer :: hdf5_error

        dimensions(1) = size(values, kind=hsize_t)
        call h5screate_simple_f(1, dimensions, dataspace_id, hdf5_error)
        call check_hdf5_error(hdf5_error, 'could not create dataspace')

        call h5dcreate_f(file_id, name, H5T_IEEE_F64LE, dataspace_id, &
                         dataset_id, hdf5_error)
        call check_hdf5_error(hdf5_error, 'could not create dataset')

        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, values, dimensions, hdf5_error)
        call check_hdf5_error(hdf5_error, 'could not write dataset')
        call h5dclose_f(dataset_id, hdf5_error)
        call check_hdf5_error(hdf5_error, 'could not close dataset')
        call h5sclose_f(dataspace_id, hdf5_error)
        call check_hdf5_error(hdf5_error, 'could not close dataspace')
    end subroutine write_dataset

    subroutine check_hdf5_error(hdf5_error, message)
        implicit none
        integer, intent(in) :: hdf5_error
        character(len=*), intent(in) :: message

        if (hdf5_error < 0) then
            print *, 'hdf5 error: ', message
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
    end subroutine check_hdf5_error

    ! Exchange boundary data with neighbors (periodic)
    subroutine exchange_boundaries(phi, phi_left, phi_right, nx_local, &
                                   left_rank, right_rank, rank, nprocs)
        implicit none
        integer, intent(in) :: nx_local, left_rank, right_rank, rank, nprocs
        real(8), intent(in) :: phi(nx_local)
        real(8), intent(out) :: phi_left, phi_right
        integer :: ierr, status(MPI_STATUS_SIZE)
        
        ! Send right boundary to right neighbor, receive from left neighbor
        call MPI_Sendrecv(phi(nx_local), 1, MPI_DOUBLE_PRECISION, right_rank, 0, &
                         phi_left, 1, MPI_DOUBLE_PRECISION, left_rank, 0, &
                         MPI_COMM_WORLD, status, ierr)
        
        ! Send left boundary to left neighbor, receive from right neighbor
        call MPI_Sendrecv(phi(1), 1, MPI_DOUBLE_PRECISION, left_rank, 1, &
                         phi_right, 1, MPI_DOUBLE_PRECISION, right_rank, 1, &
                         MPI_COMM_WORLD, status, ierr)
    end subroutine exchange_boundaries
    
    ! 4th order central difference (with ghost cells)
    subroutine spatial_derivative_4th_mpi(df, f, phi_left, phi_right, dx, n)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: f(n), dx, phi_left, phi_right
        real(8), intent(out) :: df(n)
        integer :: i
        real(8) :: inv_12dx
        real(8) :: fm2, fm1, fp1, fp2
        
        inv_12dx = 1.0d0 / (12.0d0 * dx)
        
        do i = 1, n
            ! Handle boundary points with ghost cells
            if (i == 1) then
                fm2 = phi_left
                fm1 = phi_left
                if (n == 1) then
                    fp1 = phi_right
                    fp2 = phi_right
                else if (i+1 <= n) then
                    fp1 = f(i+1)
                    if (i+2 <= n) then
                        fp2 = f(i+2)
                    else
                        fp2 = phi_right
                    end if
                end if
            else if (i == 2) then
                fm2 = phi_left
                fm1 = f(i-1)
                fp1 = f(i+1)
                if (i+2 <= n) then
                    fp2 = f(i+2)
                else
                    fp2 = phi_right
                end if
            else if (i == n-1) then
                fm2 = f(i-2)
                fm1 = f(i-1)
                fp1 = f(i+1)
                fp2 = phi_right
            else if (i == n) then
                fm2 = f(i-2)
                fm1 = f(i-1)
                fp1 = phi_right
                fp2 = phi_right
            else
                fm2 = f(i-2)
                fm1 = f(i-1)
                fp1 = f(i+1)
                fp2 = f(i+2)
            end if
            
            df(i) = (-fp2 + 8.0d0*fp1 - 8.0d0*fm1 + fm2) * inv_12dx
        end do
    end subroutine spatial_derivative_4th_mpi
    
    ! RHS with MPI boundaries
    subroutine rhs_mpi(dphi, phi, dphi_dx, phi_left, phi_right, c0, dx, n)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: phi(n), c0, dx, phi_left, phi_right
        real(8), intent(out) :: dphi(n), dphi_dx(n)
        integer :: i
        
        call spatial_derivative_4th_mpi(dphi_dx, phi, phi_left, phi_right, dx, n)
        
        do i = 1, n
            dphi(i) = -c0 * dphi_dx(i)
        end do
    end subroutine rhs_mpi
    
    ! RK3 time step with MPI
    subroutine rk3_step_mpi(phi, phi1, phi2, k1, k2, k3, dphi_dx, &
                           phi_left, phi_right, c0, dx, dt, n, &
                           left_rank, right_rank, rank, nprocs)
        implicit none
        integer, intent(in) :: n, left_rank, right_rank, rank, nprocs
        real(8), intent(in) :: c0, dx, dt
        real(8), intent(inout) :: phi(n)
        real(8), intent(out) :: phi1(n), phi2(n), k1(n), k2(n), k3(n), dphi_dx(n)
        real(8), intent(inout) :: phi_left, phi_right
        integer :: i
        real(8) :: phi1_left, phi1_right, phi2_left, phi2_right
        
        ! Stage 1
        call rhs_mpi(k1, phi, dphi_dx, phi_left, phi_right, c0, dx, n)
        do i = 1, n
            phi1(i) = phi(i) + dt * k1(i)
        end do
        call exchange_boundaries(phi1, phi1_left, phi1_right, n, left_rank, right_rank, rank, nprocs)
        
        ! Stage 2
        call rhs_mpi(k2, phi1, dphi_dx, phi1_left, phi1_right, c0, dx, n)
        do i = 1, n
            phi2(i) = phi(i) + dt * (0.25d0*k1(i) + 0.25d0*k2(i))
        end do
        call exchange_boundaries(phi2, phi2_left, phi2_right, n, left_rank, right_rank, rank, nprocs)
        
        ! Stage 3
        call rhs_mpi(k3, phi2, dphi_dx, phi2_left, phi2_right, c0, dx, n)
        do i = 1, n
            phi(i) = phi(i) + dt * (k1(i)/6.0d0 + k2(i)/6.0d0 + 2.0d0*k3(i)/3.0d0)
        end do
    end subroutine rk3_step_mpi
    
    ! Gather solution from all processes to rank 0
    subroutine gather_solution(phi_local, phi_global, nx_local, nx_global, rank, nprocs)
        implicit none
        integer, intent(in) :: nx_local, nx_global, rank, nprocs
        real(8), intent(in) :: phi_local(nx_local)
        real(8), intent(out) :: phi_global(nx_global)
        integer :: ierr
        integer, allocatable :: recvcounts(:), displs(:)
        integer :: i, nx_per_proc, remainder
        
        allocate(recvcounts(nprocs))
        allocate(displs(nprocs))
        
        nx_per_proc = nx_global / nprocs
        remainder = mod(nx_global, nprocs)
        
        do i = 1, nprocs
            if (i < nprocs) then
                recvcounts(i) = nx_per_proc
            else
                recvcounts(i) = nx_per_proc + remainder
            end if
            if (i == 1) then
                displs(i) = 0
            else
                displs(i) = displs(i-1) + recvcounts(i-1)
            end if
        end do
        
        call MPI_Gatherv(phi_local, nx_local, MPI_DOUBLE_PRECISION, &
                        phi_global, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                        0, MPI_COMM_WORLD, ierr)
        
        deallocate(recvcounts, displs)
    end subroutine gather_solution

end program wave_equation_1d_mpi