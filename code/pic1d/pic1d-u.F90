!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!                                pic1d                                       !
!                          Version 1, 2019                                   !
!                             Chen Zhao                                      !
!                                                                            !
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program pic1d
  implicit none

  integer i, it, nt, np, ng, gf, gb, flag
  real q, m, l, vt, T, e, mp, me, k, dt, rho0, vmax, dx, Ep, Eg0, &
       ap0, ap1, ap2, ap3, vp1, vp2, vp3, dxp, dvp, omega_p, wb, wf
  real, dimension(:), allocatable :: xp, xp0, xp1, xp2, xp3, vp, vp0, &
                                     Eg, rhog, time, r1, r2
  double precision, parameter :: pi = 3.141592653589793239
                                     

! initialize paricles

  np = 1e6
  ng = 64 
  k = 0.5
  l = 2 * pi / k
  mp = 1836 
  me = 1
  m = me
  e = -1
  T = 1
  vt = sqrt(T / m)
  vmax = 20 * vt
  omega_p = 1
  q = omega_p ** 2 / (e * np / l)
  rho0 = - q * np / l
  dx = l / ng
  allocate(xp(np), xp1(np), xp2(np), xp3(np), vp(np), r1(np), r2(np))
  allocate(Eg(ng), rhog(ng), xp0(np), vp0(np))
  allocate(time(nt))
  do i = 1, np
     xp(i) = l / (np - 1) * (i - 1)
     xp(i) = xp(i) + 0.01 * cos(k * xp(i))
  enddo
  call random_number(r1)
  call random_number(r2)
  flag = 0
  do i = 1, np
     if (flag == 0) then
        vp(i) = vt * sqrt(-2 * log(r1(i))) * cos(2.0 * pi * r2(i))
        flag = 1
     endif
     if (flag == 1) then
        vp(i) = vt * sqrt(-2 * log(r1(i))) * sin(2.0 * pi * r2(i))
        flag = 0
     endif
     if (vp(i) > vmax) vp(i) = vmax
     if (vp(i) < -vmax) vp(i) = -vmax
  enddo
 
  nt = 1000
  dt = 0.05
  open(unit=233, file='pic1d.dat',form='binary')
  write(233)nt, np, ng
  write(233)q, dx, dt
! particle push(RK4) 
  do it = 1, nt + 1
     
     rhog = 0.0
     Eg = 0.0
     Eg0 = 0.0
!     write(*,*)Eg,rhog
     do i = 1, np
        gb = floor(xp(i) / dx - 0.5) + 1
        gf = gb + 1
        wb = 1 - abs(xp(i) / dx + 0.5 - gb)
        wf = 1 - wb
        if (gb < 1) gb = gb + ng
        if (gf < 1) gf = gf + ng
        if (gb > ng) gb = gb - ng
        if (gf > ng) gf = gf - ng
        rhog(gb) = rhog(gb) + q * wb / dx
        rhog(gf) = rhog(gf) + q * wf / dx
     enddo
     do i = 1, ng - 1
        Eg(i+1) = Eg(i) + (rhog(i) + rhog(i + 1) + 2 * rho0) / 2 * dx
     enddo
     Eg(1) = Eg(ng) + (rhog(ng) + rho0) * dx
     Eg0 = sum(Eg) / ng 
     do i = 1, ng
        Eg(i) = Eg(i) - Eg0
     enddo
     
!     write(*,*)Eg,rhog
     write(233)xp, vp, Eg
      
     do i = 1, np
        
        vp0(i) = vp(i)
        xp0(i) = xp(i)
     
        gb = floor(xp0(i) / dx - 0.5) + 1
        gf = gb + 1
        wb = 1 - abs(xp0(i) / dx + 0.5 - gb)
        wf = 1 - wb
        if (gb < 1) gb = gb + ng
        if (gf < 1) gf = gf + ng
        if (gb > ng) gb = gb - ng
        if (gf > ng) gf = gf - ng
        Ep = Eg(gb) * wb + Eg(gf) * wf
        ap0 = e * Ep / m
     
        vp1 = vp0(i) + ap0 * dt / 2
        xp1(i) = xp0(i) + vp1 * dt / 2
        vp(i) = vp(i) + ap0 * dt / 6
        xp(i) = xp(i) + (vp0(i) + 2 * vp1) * dt /6 
     enddo
     Eg = 0.0
     rhog = 0.0
     do i = 1, np
        gb = floor(xp1(i) / dx - 0.5) + 1
        gf = gb + 1
        wb = 1 - abs(xp1(i) / dx + 0.5 - gb)
        wf = 1 - wb
        if (gb < 1) gb = gb + ng
        if (gf < 1) gf = gf + ng
        if (gb > ng) gb = gb - ng
        if (gf > ng) gf = gf - ng
        rhog(gb) = rhog(gb) + q * wb / dx
        rhog(gf) = rhog(gf) + q * wf / dx
     enddo
     do i = 1, ng - 1
        Eg(i+1) = Eg(i) + (rhog(i) + rhog(i + 1) + 2 * rho0) / 2 * dx
     enddo
     Eg(1) = Eg(ng) + (rhog(ng) + rho0) * dx
     Eg0 = sum(Eg) / ng
     do i = 1, ng
        Eg(i) = Eg(i) - Eg0
     enddo
     do i = 1, np
        gb = floor(xp1(i) / dx - 0.5) + 1
        gf = gb + 1
        wb = 1 - abs(xp1(i) / dx + 0.5 - gb)
        wf = 1 - wb
        if (gb < 1) gb = gb + ng
        if (gf < 1) gf = gf + ng
        if (gb > ng) gb = gb - ng
        if (gf > ng) gf = gf - ng
        Ep = Eg(gb) * wb + Eg(gf) * wf
        ap1 = e * Ep / m
     
        vp2 = vp0(i) + ap1 * dt / 2
        xp2(i) = xp0(i) + vp2 * dt / 2
        vp(i) = vp(i) + ap1 * dt / 3
        xp(i) = xp(i) + vp2 * dt / 3
     enddo
     Eg = 0.0
     rhog = 0.0
     do i = 1, np
        gb = floor(xp2(i) / dx - 0.5) + 1
        gf = gb + 1
        wb = 1 - abs(xp2(i) / dx + 0.5 - gb)
        wf = 1 - wb
        if (gb < 1) gb = gb + ng
        if (gf < 1) gf = gf + ng
        if (gb > ng) gb = gb - ng
        if (gf > ng) gf = gf - ng
        rhog(gb) = rhog(gb) + q * wb / dx
        rhog(gf) = rhog(gf) + q * wf / dx
     enddo
     do i = 1, ng - 1
        Eg(i+1) = Eg(i) + (rhog(i) + rhog(i + 1) + 2 * rho0) / 2 * dx
     enddo
     Eg(1) = Eg(ng) + (rhog(ng) + rho0) * dx
     Eg0 = sum(Eg) / ng
     do i = 1, ng
        Eg(i) = Eg(i) - Eg0
     enddo
     do i = 1, np
        gb = floor(xp2(i) / dx - 0.5) + 1
        gf = gb + 1
        wb = 1 - abs(xp2(i) / dx + 0.5 - gb)
        wf = 1 - wb
        if (gb < 1) gb = gb + ng
        if (gf < 1) gf = gf + ng
        if (gb > ng) gb = gb - ng
        if (gf > ng) gf = gf - ng
        Ep = Eg(gb) * wb + Eg(gf) * wf
        ap2 = e * Ep / m
     
        vp3 = vp0(i) + ap2 *dt
        xp3(i) = xp0(i) + vp3 * dt
        vp(i) = vp(i) + ap2 * dt / 3
        xp(i) = xp(i) + vp3 * dt / 6
        
    enddo
    Eg = 0.0
    rhog = 0.0
    do i = 1, np
        gb = floor(xp3(i) / dx - 0.5) + 1
        gf = gb + 1
        wb = 1 - abs(xp3(i) / dx + 0.5 - gb)
        wf = 1 - wb
        if (gb < 1) gb = gb + ng
        if (gf < 1) gf = gf + ng
        if (gb > ng) gb = gb - ng
        if (gf > ng) gf = gf - ng
        rhog(gb) = rhog(gb) + q * wb / dx
        rhog(gf) = rhog(gf) + q * wf / dx
     enddo
     do i = 1, ng - 1
        Eg(i+1) = Eg(i) + (rhog(i) + rhog(i + 1) + 2 * rho0) / 2 * dx
     enddo
     Eg(1) = Eg(ng) + (rhog(ng) + rho0) * dx
     Eg0 = sum(Eg) / ng
     do i = 1, ng
        Eg(i) = Eg(i) - Eg0
     enddo
     do i = 1, np
        gb = floor(xp3(i) / dx - 0.5) + 1
        gf = gb + 1
        wb = 1 - abs(xp3(i) / dx + 0.5 - gb)
        wf = 1 - wb
        if (gb < 1) gb = gb + ng
        if (gf < 1) gf = gf + ng
        if (gb > ng) gb = gb - ng
        if (gf > ng) gf = gf - ng
        Ep = Eg(gb) * wb + Eg(gf) * wf
        ap3 = e * Ep / m
     
        vp(i) = vp(i) + ap3 * dt / 6
        xp(i) = xp(i) / l + 20.0
        xp(i) = l * (xp(i) - floor(xp(i)))
    enddo

    write(*,*)it
  enddo

  close(233)

end program pic1d 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
