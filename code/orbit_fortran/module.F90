        module system_env
            #ifdef _OPENACC
            use openacc
            use cudafor
            #endif
        end module system_env

        module precision
            use system_env
            integer, parameter :: doubleprec=selected_real_kind(12),&
                                  singleprec=selected_real_kind(6),defaultprec=kind(0.0)
            #ifdef DOUBLE_PRECISION
                integer, parameter :: wp=doubleprec
            #else
                integer, parameter :: wp=singleprec
            #endif
        end module precision


        module global_parameters
            use precision

! control parameters
            integer, parameter ::  c = 299792458, np = 10
            real(wp), parameter ::  e = 1.602176565e-19, m_pr = 1.672621777e-27, m_el = 9.10938291e-31

            end module global_parameters

        module equilibrium
            real sibry, simag, m, q, dt, b0
            real, dimension(:), allocatable :: r, z
            real, dimension(:) :: vpa(np), vper(np), rn1(np),rn2(np), rn3(np), rn4(np), &
                                  x(np), y(np), vx(np), vy(np), vz(np)
            real, dimension(:,:), allocatable :: Brxy, Btxy, Bzxy, Rxy, Zxy

            call random_number(r1)
            call random_number(r2)
            call random_number(r3)
            call random_number(r4)
       end module equilibrium

