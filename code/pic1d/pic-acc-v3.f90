!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           pic1d code with openAcc                    !                   
!          by Chen Zhao^1 & Mauricio Ferrato^2         !
!                   1. PPPL                            !
!                   2. University of Delaware          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program pic1d
  implicit none

  integer i, it, nt, np, ng, gf, gb, flag
  real q, m, l, vt, T, e, mp, me, k, dt, rho0, vmax, dx, Ep, Eg0, &
       ap0, ap1, ap2, ap3, vp1, vp2, vp3, dxp, dvp, &
       omega_p, wb, wf, Ef0, Ek0, efg, ekg, beginning, rate
  real, dimension(:), allocatable :: xp, xp1, xp2, xp3, vp, Eg, rhog, &
                                     time, r1, r2, xp0, vp0, Ek, Ef
  double precision, parameter :: pi = 3.141592653589793239

! initialize paricles
  call cpu_time(beginning)
  np = 3e7
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
  nt = 1000
  dt = 0.005
  allocate(xp(np), xp1(np), xp2(np), xp3(np), vp(np), r1(np), r2(np))
  allocate(Eg(ng), rhog(ng), xp0(np), vp0(np), Ef(nt), Ek(nt))
  !allocate(time(nt))
  
  !$acc enter  data create(xp(np), xp1(np), xp2(np), xp3(np), vp(np), Eg(ng), rhog(ng))
  !$acc enter data create(xp0(np), vp0(np))
  !$acc parallel loop present(xp)
  do i = 1, np
     xp(i) = l / (np - 1) * (i - 1)
     xp(i) = xp(i) + 0.03 * cos(k * xp(i))
  enddo
  call random_number(r1)
  call random_number(r2)
  !$acc enter data copyin(r1, r2)
  !$acc parallel loop present(vp, r2, r1)
  do i = 1, np
     if (mod(np, 2) == 0) then
        vp(i) = vt * sqrt(-2 * log(r1(i))) * cos(2.0 * pi * r2(i))
     endif
     if (mod(np, 2) == 1) then
        vp(i) = vt * sqrt(-2 * log(r1(i))) * sin(2.0 * pi * r2(i))
     endif
     if (vp(i) > vmax) vp(i) = vmax
     if (vp(i) < -vmax) vp(i) = -vmax
  enddo
  Ef = 0.0
  Ek = 0.0
  open(unit=233, file='/gpfs/wolf/gen127/scratch/zcshr/pic1d-gpu.dat',form = 'binary')
  write(233)nt, np, ng
  write(233)q, dx, dt
! particle push(RK4) 
  do it = 1, nt + 1
     
     !$acc parallel loop present(rhog, Eg)
     do i = 1, ng
        rhog(i) = 0.0
        Eg(i) = 0.0
     enddo

     Eg0 = 0.0
     Ef0 = 0.0
     Ek0 = 0.0
    
    !$acc parallel loop private(gb, gf, wb, wf) present(xp, rhog)
     do i = 1, np
        gb = floor(xp(i) / dx - 0.5) + 1
        gf = gb + 1
        wb = 1 - abs(xp(i) / dx + 0.5 - gb)
        wf = 1 - wb
        if (gb < 1) gb = gb + ng
        if (gf < 1) gf = gf + ng
        if (gb > ng) gb = gb - ng
        if (gf > ng) gf = gf - ng
        !$acc atomic update
        rhog(gb) = rhog(gb) + q * wb / dx
        !$acc atomic update
        rhog(gf) = rhog(gf) + q * wf / dx
     enddo
     !$acc update self(rhog)
     !$acc update self(Eg)
     do i = 1, ng - 1
        Eg(i+1) = Eg(i) + (rhog(i) + rhog(i + 1) + 2 * rho0) / 2 * dx
     enddo
     Eg(1) = Eg(ng) + (rhog(ng) + rho0) * dx
     Eg0 = sum(Eg) / ng
     !$acc update device(Eg)
     !$acc parallel loop independent reduction(+:Ef0) present(Eg) 
     do i = 1, ng
        Eg(i) = Eg(i) - Eg0
        Ef0 = Ef0 + Eg(i) * Eg(i) * dx / 2
     enddo
     !$acc parallel loop independent private(gb, gf, wb, wf, ap0, Ep, vp1) reduction(+:Ek0) present(vp, xp1, vp0, xp0, xp, Eg)
     do i = 1, np
        
        vp0(i) = vp(i)
        xp0(i) = xp(i)
        Ek0 = Ek0 + abs(q) * vp(i) * vp(i) / 2
     
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
     !$acc parallel loop independent private(gb, gf, wb, wf, vp2, ap1, Ep) &
     !$acc present(vp, xp1, xp0, xp2, vp0, xp, Eg)
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
     !$acc parallel loop independent private(gb, gf, wb, wf, vp3, ap2, Ep) &
     !$acc present(vp, xp, xp3, Eg, vp0, xp0, xp2)
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
     !$acc parallel loop independent private(gb, gf, wb, wf, ap3, Ep) &
     !$acc present(Eg, vp, xp, xp3)
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
    Ek(it) = Ek0
    Ef(it) = Ef0
    write(*,*)it
  enddo

  write(233)Ef, Ek
  close(233)
  !$acc exit data delete(xp, xp1, xp2, xp3, vp, r1, r2)
  !$acc exit data delete(Eg, rhog, xp0, vp0)
  
  deallocate(xp, xp1, xp2, xp3, vp, r1, r2)
  deallocate(Eg, rhog, xp0, vp0, Ef, Ek)
  !deallocate(time)
  call cpu_time(rate)
  print *, "elapsed time: ", rate - beginning

end program pic1d
