subroutine setup
    use global_parameters
    use equilibrium
    implicit none
    integer Bid, RZid, nh, nw
    integer ip, nt
    real T, vmax, vt

    CHARACTER(*), PARAMETER :: bc = 'natural'
    
    nw = 511 
    nh = 511
    
    allocate(Btxy(nw, nh), Brxy(nw, nh), Bzxy(nw, nh))
    allocate(Rxy(nw, nh), Zxy(nw, nh))
    allocate(rsp(nw), zsp(nh))
    allocate(Btxy(16, nw, nh), Brsp(16, nw, nh), Bzsp(16, nw, nh))

    Bid = 233
    open(Bid, file = 'Brz.dat')
    read(Bid, *)Brxy
    read(Bid, *)Btxy
    read(Bid, *)Bzxy
    close(Bid)
    RZid = 433
    open(RZid , file = 'RZ.dat')
    read(RZid, *)Rxy
    read(RZid, *)Zxy
    close(RZid)

    ! construct spline for magnetic fields
    rsp = Rxy(:,1)
    zsp = Zxy(1,:)
    Btsp(0, :, :) = Btxy
    Brsp(0, :, :) = Brxy
    Bzsp(0, :, :) = Bzxy
    call construct_spline2d(nw, nh, rsp, zsp, Btsp, bc, bc)
    call construct_spline2d(nw, nh, rsp, zsp, Brsp, bc, bc)
    call construct_spline2d(nw, nh, rsp, zsp, Bzsp, bc, bc)

    m = m_pr
    q = e
    T = 1e5 * e
    vt = sqrt(2 * T / m)
    vmax = 20 * vt
    nt = 50
    b0 = 1.7
    dt = abs(2 * pi * m / q / b0) / nt 
    allocate(r(np), z(np))
    call random_number(rnd1)
    call random_number(rnd2)
    call random_number(rnd3)
    call random_number(rnd4)
    do ip = 1, np
        if (mod(np, 2) == 0) then
            vpa(ip) = vt * sqrt(-2 * log(rnd1(i))) * cos(2.0 * pi * rnd2(i))
            vper(ip) = abs(vt * sqrt(-2 * log(rnd3(i))) * &
                        cos(2.0 * pi * rnd4(i)))
        endif
        if (mod(np, 2) == 1) then
            vpa(ip) = vt * sqrt(-2 * log(rnd1(i))) * sin(2.0 * pi * rnd2(i))
            vper(ip) = abs(vt * sqrt(-2 * log(rnd3(i))) * &
                        cos(2.0 * pi * rnd4(i)))
        endif
        if (vpa(ip) > vmax) vpa(ip) = vmax
        if (vpa(ip) < -vmax) vpa(ip) = -vmax
        if (vper(ip) > vmax) vper(ip) = vmax
        r(ip) = Rxy(233,233)
        z(ip) = Zxy(233,233)
        phi(ip) = 0.0
        x(ip) = r(ip) * cos(phi(ip))
        y(ip) = r(ip) * sin(phi(ip))
        br = spline2d(r, z, nw, nh, r1, z1, Brsp)
        bt = spline2d(r, z, nw, nh, r1, z1, Btsp)
        bz = spline2d(r, z, nw, nh, r1, z1, Bzsp)
        b = sqrt(br^2 + bt^2 + bz^2)
        bx = br * x0(ip) / r0(ip) - bt * y0(ip) / r0(ip)
        by = br * y0(ip) / r0(ip) + bt * x0(ip) / r0(ip)
        vx(ip) = (vpa(ip) * bx + vper(ip) * bx * bz / sqrt(bx^2 + by^2)) / b
        vy(ip) = (vpa(ip) * by + vper(ip) * by * bz / sqrt(bx^2 + by^2)) / b
        vz(ip) = (vpa(ip) * bz - vper(ip) * sqrt(bx^2 + by^2)) / b
    enddo
end subroutine setup
