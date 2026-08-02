! wave_stdpar.f90
!
! solves the one-dimensional linear advection equation
!     dphi/dt = -c0*dphi/dx
! on a periodic domain using a fourth-order central spatial derivative
! and third-order runge-kutta time integration.
!
! this standard-parallelism variant reports timing and numerical error
! against the analytical sine-wave solution. the do concurrent loops
! compile to serial, multicore, or GPU from one source, as a reference
! for comparing parallel backends and optimisation strategies.
!
! output: wave_fortran_stdparr.h5
! author: dekeract01, 2026
! university of southampton

program wave_equation_1d_stdparr
    implicit none

    ! assign double precision numbers where apppropriate
    real(8), parameter :: c0 = 0.5d0
    real(8), parameter :: dt = 0.0000001d0
    integer, parameter :: niter = 1000
    integer, parameter :: nx = 2000000
    real(8), parameter :: dx = 1.0d0/nx
    real(8), parameter :: pi = 4.0d0*atan(1.0d0)

    real(8), allocatable :: x0(:), phi(:), phi_exact(:), error(:)   ! allocatable -> managed mem on gpu
    real(8), allocatable :: phi1(:), phi2(:), k1(:), k2(:), k3(:), dphi_dx(:)
 
    integer :: i, n
    real(8) :: t_final, phi_min, phi_max, max_error, l2_error
    integer(8) :: t_start, t_end, clock_rate           ! wall-clock timing
    real(8) :: total_time, elapsed

    allocate(x0(nx), phi(nx), phi_exact(nx), error(nx))
    allocate(phi1(nx), phi2(nx), k1(nx), k2(nx), k3(nx), dphi_dx(nx))

    print *, '1d wave equation simulation (fortran, stdpar / do concurrent)'
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
    ! this is where it differs, the concurrent part. 
    do concurrent (i = 1:nx)
        x0(i) = (i-1)*dx
    end do

    ! initial condition: phi = sin(2*pi*x)
    do concurrent (i = 1:nx)
        phi(i) = sin(2.0d0*pi*x0(i))
    end do

    ! time stepping loop
    print *, ''
    print *, 'starting simulation...'
    call system_clock(t_start, clock_rate)

    do n = 1, niter
        call rk3_step(phi, phi1, phi2, k1, k2, k3, dphi_dx, c0, dx, dt, nx)

        if (mod(n, 100) == 0) then
            if (any(phi /= phi)) then          ! NaN check
                print *, 'nan detected at iteration ', n, '!'
                stop
            end if
            call system_clock(t_end)
            elapsed = real(t_end - t_start, 8)/real(clock_rate, 8)
            phi_min = minval(phi)
            phi_max = maxval(phi)
            print '(A,I4,A,I4,A,F6.2,A,F9.6,A,F9.6,A)', &
                'iter ', n, '/', niter, ' | time: ', elapsed, &
                's | phi: [', phi_min, ', ', phi_max, ']'
        end if
 
    end do

    call system_clock(t_end)
    total_time = real(t_end - t_start, 8)/real(clock_rate, 8)

     print *, repeat('-', 50)
    print *, 'simulation complete!'
    print '(A,F12.8,A)', 'total time: ', total_time, ' seconds'
    print '(A,F12.8,A)', 'time per iteration: ', total_time/niter*1000.0d0, ' ms'
    print '(A,F10.1)', 'iterations per second: ', niter/total_time
    print *, repeat('-', 50)

    ! analytical solution + error
    t_final = niter*dt
    do concurrent (i = 1:nx)
        phi_exact(i) = sin(2.0d0*pi*(x0(i)-c0*t_final))
        error(i) = abs(phi(i) - phi_exact(i))
    end do

    max_error = maxval(error)
    l2_error = sqrt(sum(error**2)/nx)

    print *, ''
    print *, 'numerical vs analytical:'
    print '(A,ES12.6)', 'max error: ', max_error
    print '(A,ES12.6)', 'l2 error: ', l2_error

    open(unit=10, file='wave_solution_stdpar.bin', status='replace', &
         access='stream', form='unformatted', action='write')
    do i = 1, nx
        write(10) x0(i), phi(i), phi_exact(i), error(i)
    end do
    close(10)
    print *, ''
    print *, 'results written to wave_solution_stdpar.bin'

    deallocate(x0, phi, phi_exact, error)
    deallocate(phi1, phi2, k1, k2, k3, dphi_dx)

contains

    ! 4th-order central diff, periodic via modular indexing (one on-device kernel).
    subroutine spatial_derivative_fourth(df, f, dx, n)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: f(n), dx
        real(8), intent(out) :: df(n)
        integer :: i
        real(8) :: inv_12dx

        inv_12dx = 1.0d0/(12.0d0*dx)

        do concurrent (i = 1:n)
            df(i) = ( -f(modulo(i+1,n)+1) + 8.0d0*f(modulo(i,n)+1) &
                      -8.0d0*f(modulo(i-2,n)+1) + f(modulo(i-3,n)+1) )*inv_12dx
        end do
    end subroutine spatial_derivative_fourth

    ! rhs: dphi/dt = -c*dphi/dx
    subroutine rhs(dphi, phi, dphi_dx, c0, dx, n)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: phi(n), c0, dx
        real(8), intent(out) :: dphi(n), dphi_dx(n)
        integer :: i

        call spatial_derivative_fourth(dphi_dx, phi, dx, n)

        do concurrent (i = 1:n)
            dphi(i) = -c0*dphi_dx(i)
        end do
    end subroutine rhs

    ! rk3 step changes with concurrent loop
    subroutine rk3_step(phi, phi1, phi2, k1, k2, k3, dphi_dx, c0, dx, dt, n)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: c0, dx, dt
        real(8), intent(inout) :: phi(n)
        real(8), intent(out) :: phi1(n), phi2(n), k1(n), k2(n), k3(n), dphi_dx(n)
        integer :: i

        call rhs(k1, phi, dphi_dx, c0, dx, n)
        do concurrent (i = 1:n)
            phi1(i) = phi(i) + dt*k1(i)
        end do

        call rhs(k2, phi1, dphi_dx, c0, dx, n)
        do concurrent (i = 1:n)
            phi2(i) = phi(i) + dt*(0.25d0*k1(i) + 0.25d0*k2(i))
        end do

        call rhs(k3, phi2, dphi_dx, c0, dx, n)
        do concurrent (i = 1:n)
            phi(i) = phi(i) + dt*(k1(i)/6.0d0 + k2(i)/6.0d0 + 2.0d0*k3(i)/3.0d0)
        end do
    end subroutine rk3_step

end program wave_equation_1d_stdparr