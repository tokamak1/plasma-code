!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!                           general orbit                                    !
!                          Version 1, 2019                                   !
!                             Chen Zhao                                      !
!                                                                            !
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        program orbit
            use global_parameters
            use equilibrium
            implicit none

            integer i, fid
            integer, parameter :: mstap = 1000


!load parameters from gfile and setup the particle parameters 
            CALL SETUP
! main time loop
            fid = 444
            open(fid, file = 'orbitm.dat', status = 'replace')
            do istep=1,mstep
                call rk4ob
                write(fid, *)x, y, z
            enddo
            close(fid)
        end program orbit 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine rk4ob
            use global_parameters
            use equilibrium
            implicit none
            real(wp) br, bt, bz, bx, by, ax0, ay0, az0, ax1, ax2, ax3,&
                     ay1, ay2, ay3, az1, az2, az3, dvx, dvy, dvz, dx, dy, dz 
            real(wp), dimension(:) :: x0(np), y0(np), z0(np), r0(np), phi0(np), vpa(np), vper(np), &
                                      vx0(np), vy0(np), vz0(np), vx1(np), vy1(np), r1(np), r2(np), &
                                      vz1(np), x1(np), y1(np), z1(np), x2(np), y2(np), r3(np), &
                                      z2(np), vx2(np), vy2(np), vz2(np), x3(np), &
                                      y3(np), z3(np), vx3(np), vy3(np), vz3(np)
            integer ip

            x0 = x
            y0 = y
            z0 = z
            r0 = r
            phi0 = phi
            vx0 = vx
            vy0 = vy
            vz0 = vz

            do ip = 1, np

!       br = interp2(r,z,fr,r0,z0)
!       bt = interp2(r,z,ft,r0,z0)
!       bz = interp2(r,z,fz,r0,z0)

                b = sqrt(br^2 + bt^2 + bz^2)
                bx = br * x0(ip) / r0(ip) - bt * y0(ip) / r0(ip)
                by = br * y0(ip) / r0(ip) + bt * x0(ip) / r0(ip)
                ax0 = q * (vy0(ip) * bz - vz0(ip) * by) / m
                ay0 = q * (vz0(ip) * bx - vx0(ip) * bz) / m;
                az0 = q * (vx0(ip) * by - vy0(ip) * bx) / m;
                vx1(ip) = vx0(ip) + ax0 * dt / 2
                vy1(ip) = vy0(ip) + ay0 * dt / 2
                vz1(ip) = vz0(ip) + az0 * dt / 2
                x1(ip) = x0(ip) + vx1(ip) * dt / 2
                y1(ip) = y0(ip) + vy1(ip) * dt / 2
                z1(ip) = z0(ip) + vz1(ip) * dt / 2
                r1(ip) = sqrt(x1(ip)^2 + y1(ip)^2)
!br = interp2(r,z,fr,r1,z1)
!bt = interp2(r,z,ft,r1,z1);
!bz = interp2(r,z,fz,r1,z1);
                bx = br * x1(ip) / r1(ip) - bt * y1(ip) / r1(ip)
                by = br * y1(ip) / r1(ip) + bt * x1(ip) / r1(ip)
                ax1 = q * (vy1(ip) * bz - vz1(ip) * by) / m
                ay1 = q * (vz1(ip) * bx - vx1(ip) * bz) / m
                a1z = q * (vx(ip) * by - vy(ip) * bx) / m
                vx2(ip) = vx0(ip) + ax1 * dt / 2
                vy2(ip) = vy0(ip) + ay1 * dt / 2
                vz2(ip) = vz0(ip) + az1 * dt / 2
                x2(ip) = x0(ip) + vx2(ip) * dt / 2
                y2(ip) = y0(ip) + vy2(ip) * dt / 2
                z2(ip) = z0(ip) + vz2(ip) * dt / 2
                r2(ip) = sqrt(x2(ip)^2 + y2(ip)^2)
!br = interp2(r,z,fr,r2,z2)
!bt = interp2(r,z,ft,r2,z2)
!bz = interp2(r,z,fz,r2,z2)
                bx = br * x2(ip) / r2(ip) - bt * y2(ip) / r2(ip)
                by = br * y2(ip) / r2(ip) + bt * x2(ip) / r2(ip)
                ax2 = q * (vy2(ip) * bz - vz2(ip) * by) / m
                ay2 = q * (vz2(ip) * bx - vx2(ip) * bz) / m
                az2 = q * (vx2(ip) * by - vy2(ip) * bx) / m
                vx3(ip) = vx0(ip) + ax2 * dt
                vy3(ip) = vy0(ip) + ay2 * dt
                vz3(ip) = vz0(ip) + az2 * dt
                x3(ip) = x0(ip) + vx3 * dt
                y3(ip) = y0(ip) + vy3 * dt
                z3(ip) = z0(ip) + vz3 * dt
                r3(ip) = sqrt(x3(ip)^2 + y3(ip)^2)
!                br = interp2(r,z,fr,r3,z3);
!bt = interp2(r,z,ft,r3,z3);
!bz = interp2(r,z,fz,r3,z3);
                bx = br * x3(ip) / r3(ip) - bt * y3(ip) / r3(ip)
                by = br * y3(ip) / r3(ip) + bt * x3(ip) / r3(ip)
                ax3 = q * (vy3(ip) * bz - vz3(ip) * by) / m
                ay3 = q * (vz3(ip) * bx - vx3(ip) * bz) / m
                az3 = q * (vx3(ip) * by - vy(ip) * bx) / m
                dvx = (ax0 + 2 * ax1 + 2 * ax2 + ax3) * dt / 6
                dvy = (ay0 + 2 * ay1 + 2 * ay2 + ay3) * dt / 6
                dvz = (az0 + 2 * az1 + 2 * az2 + az3) * dt / 6
                dx = (vx0(ip) + 2 * vx1(ip) + 2 * vx2(ip) + vx3(ip)) * dt / 6
                dy = (vy0(ip) + 2 * vy1(ip) + 2 * vy2(ip) + vy3(ip)) * dt / 6
                dz = (vz0(ip) + 2 * vz1(ip) + 2 * vz2(ip) + vz3(ip)) * dt / 6
                vx(ip) = vx(ip) + dvx
                vy(ip) = vy(ip) + dvy
                vz(ip) = vz(ip) + dvz
                x(ip) = x(ip) + dx
                y(ip) = y(ip) + dy
                z(ip) = z(ip) + dz
            enddo         
