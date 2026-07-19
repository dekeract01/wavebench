! wave.f90
!
! solves the one-dimensional linear advection equation
!     dphi/dt = -c0*dphi/dx
! on a periodic domain using a fourth-order central spatial derivative
! and third-order runge-kutta time integration.
!
! this serial baseline reports timing and numerical error against the
! analytical sine-wave solution. it is intended as a reference for
! comparing parallel implementations and optimisation strategies.
!
! output: wave_solution_fortran.bin
! author: dekeract01, 2022
! university of southampton

program wave_equation_1d
    implicit none
    
    ! simulation parameters
    real(8), parameter :: c0 = 0.5d0          ! wave speed
    real(8), parameter :: dt = 0.0000001d0    ! time step
    integer, parameter :: niter = 1000        ! number of iterations
    integer, parameter :: nx = 2000000        ! grid points
    real(8), parameter :: dx = 1.0d0/nx       ! grid spacing
    real(8), parameter :: pi = 4.0d0*atan(1.0d0)
    
    ! arrays
    real(8), dimension(nx) :: x0, phi, phi_exact, error
    real(8), dimension(nx) :: phi1, phi2, k1, k2, k3, dphi_dx
    
    ! variables
    integer :: i, n
    real(8) :: t_final, phi_min, phi_max, max_error, l2_error
    real(8) :: start_time, end_time, total_time
    logical :: nan_detected
    
    ! print header
    print *, '1d wave equation simulation (fortran)'
    print *, repeat('=', 50)
    print '(A,I0)', 'grid points: ', nx
    print '(A,F10.6)', 'dx: ', dx
    print '(A,F10.6)', 'dt: ', dt
    print '(A,F10.6)', 'wave speed c0: ', c0
    print '(A,F10.3)', 'cfl number: ', c0*dt/dx
    print '(A,I0)', 'number of iterations: ', niter
    print '(A,F12.8)', 'total time: ', niter*dt
    print *, repeat('=', 50)
    
    ! grid setup (periodic domain from 0 to 1)
    do i = 1, nx
        x0(i) = (i-1)*dx
    end do
    
    ! initial condition: phi = sin(2*pi*x)
    do i = 1, nx
        phi(i) = sin(2.0d0*pi*x0(i))
    end do
    
    ! time stepping loop
    print *, ''
    print *, 'starting simulation...'
    call cpu_time(start_time)
    
    do n = 1, niter
        call rk3_step(phi, phi1, phi2, k1, k2, k3, dphi_dx, c0, dx, dt, nx)
        
        ! nan check and progress report
        if (mod(n, 100) == 0) then
            nan_detected = .false.
            do i = 1, nx
                if (phi(i) /= phi(i)) then  ! nan check
                    nan_detected = .true.
                    exit
                end if
            end do
            
            if (nan_detected) then
                print *, 'nan detected at iteration ', n, '!'
                stop
            end if
            
            call cpu_time(end_time)
            phi_min = minval(phi)
            phi_max = maxval(phi)
            print '(A,I4,A,I4,A,F6.2,A,F9.6,A,F9.6,A)', &
                'iter ', n, '/', niter, ' | time: ', end_time-start_time, &
                's | phi: [', phi_min, ', ', phi_max, ']'
        end if
    end do
    
    call cpu_time(end_time)
    total_time = end_time - start_time
    
    print *, repeat('-', 50)
    print *, 'simulation complete!'
    print '(A,F12.8,A)', 'total time: ', total_time, ' seconds'
    print '(A,F12.8,A)', 'time per iteration: ', total_time/niter*1000.0d0, ' ms'
    print '(A,F10.1)', 'iterations per second: ', niter/total_time
    print *, repeat('-', 50)
    
    ! analytical solution: phi = sin(2*pi*(x-c*t))
    t_final = niter * dt
    do i = 1, nx
        phi_exact(i) = sin(2.0d0*pi*(x0(i)-c0*t_final))
        error(i) = abs(phi(i) - phi_exact(i))
    end do
    
    ! calculate errors
    max_error = maxval(error)
    l2_error = sqrt(sum(error**2) / nx)
    
    print *, ''
    print *, 'numerical vs analytical:'
    print '(A,ES12.6)', 'max error: ', max_error
    print '(A,ES12.6)', 'l2 error: ', l2_error
    
    ! write native real(8) values as x, phi_numerical, phi_exact, error.
    open(unit=10, file='wave_solution_fortran.bin', status='replace', &
         access='stream', form='unformatted', action='write')
    do i = 1, nx
        write(10) x0(i), phi(i), phi_exact(i), error(i)
    end do
    close(10)
    print *, ''
    print *, 'results written to wave_solution_fortran.bin'
    
contains

    ! fourth-order central difference for spatial derivative (periodic)
    subroutine spatial_derivative_fourth(df, f, dx, n)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: f(n), dx
        real(8), intent(out) :: df(n)
        integer :: i
        real(8) :: inv_12dx

        inv_12dx = 1.0d0 / (12.0d0 * dx)

        ! interior: clean contiguous stencil, vectorizes
        do i = 3, n-2
            df(i) = (-f(i+2)+8.0d0*f(i+1)-8.0d0*f(i-1)+f(i-2))*inv_12dx
        end do

        ! periodic wrap: the four edge points
        df(1)   = (-f(3)  +8.0d0*f(2)  -8.0d0*f(n)  +f(n-1))*inv_12dx
        df(2)   = (-f(4)  +8.0d0*f(3)  -8.0d0*f(1)  +f(n)  )*inv_12dx
        df(n-1) = (-f(1)  +8.0d0*f(n)  -8.0d0*f(n-2)+f(n-3))*inv_12dx
        df(n)   = (-f(2)  +8.0d0*f(1)  -8.0d0*f(n-1)+f(n-2))*inv_12dx
    end subroutine
    
    ! rhs function: dphi/dt = -c*dphi/dx
    subroutine rhs(dphi, phi, dphi_dx, c0, dx, n)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: phi(n), c0, dx
        real(8), intent(out) :: dphi(n), dphi_dx(n)
        integer :: i
        
        call spatial_derivative_fourth(dphi_dx, phi, dx, n)
        
        do i = 1, n
            dphi(i) = -c0*dphi_dx(i)
        end do
    end subroutine rhs
    
    ! runge-kutta 3 time integration
    subroutine rk3_step(phi, phi1, phi2, k1, k2, k3, dphi_dx, c0, dx, dt, n)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: c0, dx, dt
        real(8), intent(inout) :: phi(n)
        real(8), intent(out) :: phi1(n), phi2(n), k1(n), k2(n), k3(n), dphi_dx(n)
        integer :: i
        
        ! stage 1
        call rhs(k1, phi, dphi_dx, c0, dx, n)
        do i = 1, n
            phi1(i) = phi(i) + dt*k1(i)
        end do
        
        ! stage 2
        call rhs(k2, phi1, dphi_dx, c0, dx, n)
        do i = 1, n
            phi2(i) = phi(i) + dt*(0.25d0*k1(i) + 0.25d0*k2(i))
        end do
        
        ! stage 3
        call rhs(k3, phi2, dphi_dx, c0, dx, n)
        do i = 1, n
            phi(i) = phi(i) + dt*(k1(i)/6.0d0 + k2(i)/6.0d0 + 2.0d0*k3(i)/3.0d0)
        end do
    end subroutine rk3_step

end program wave_equation_1d