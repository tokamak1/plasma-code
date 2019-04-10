subroutine construct_spline(nsp, xsp, ysp, bc)
    ! bc is a string setting boundary conditions,
    !   it can be 'natural' (2nd derivative = 0),
    !   or 'periodic' (match 1st & 2nd derivatives of two ends)
    use precision
    implicit none

    integer, intent(in) :: nsp
    real(wp), intent(in) :: xsp(0:nsp-1)
    real(wp), intent(inout) :: ysp(0:3, 0:nsp-1)
    character(*), intent(in) :: bc

    integer i, info
    integer, allocatable :: pivot(:)
    real(wp) mu, h(nsp-1), m(0:nsp-1)
    real(wp), allocatable :: coef_mat(:, :), g(:)

    if(bc == 'natural')then
        allocate(coef_mat(nsp-2, nsp-2), g(nsp-2), pivot(nsp-2))

        h(1) = xsp(1) - xsp(0)
        do i = 1, nsp-2
            h(i+1) = xsp(i+1) - xsp(i)
            mu = h(i)/(h(i)+h(i+1))

            if(i > 1)coef_mat(i, i-1) = mu
            coef_mat(i, i) = 2.0
            if(i < nsp-1)coef_mat(i, i+1) = 1 - mu

            g(i) = 6.0/(h(i)+h(i+1))*((ysp(0, i+1)-ysp(0, i))/h(i+1)&
                -(ysp(0, i)-ysp(0, i-1))/h(i))
        enddo
        call dgesv(nsp-2, 1, coef_mat, nsp-2, pivot, g, nsp-2, info)

        m(0) = 0.0
        m(nsp-1) = 0.0
        m(1:nsp-2) = g
    elseif(bc == 'periodic')then
        allocate(coef_mat(nsp-1, nsp-1), g(nsp-1), pivot(nsp-1))

        h = xsp(1:nsp-1) - xsp(0:nsp-2)
        do i = 1, nsp-1
            mu = h(i)/(h(i)+h(mod(i, nsp-1)+1))

            coef_mat(i, modulo(i-2, nsp-1)+1) = mu
            coef_mat(i, i) = 2.0
            coef_mat(i, mod(i, nsp-1)+1) = 1 - mu

            g(i) = 6.0/(h(i)+h(mod(i, nsp-1)+1)) *&
                ((ysp(0, mod(i, nsp-1)+1)-ysp(0, mod(i, nsp-1))/h(mod(i, nsp-1)+1))&
                -(ysp(0, i)-ysp(0, i-1))/h(i))
        enddo

        call dgesv(nsp-1, 1, coef_mat, nsp-1, pivot, g, nsp-1, info)

        m(0) = g(nsp-1)
        m(1:nsp-1) = g
    endif

    do i = 1, nsp-1
        ysp(1, i-1) = (ysp(0, i)-ysp(0, i-1))/h(i) &
            - h(i)*(2.0*m(i-1)+m(i))/6.0
        ysp(2, i-1) = m(i-1)/2.0
        ysp(3, i-1) = (m(i)-m(i-1))/(6.0*h(i))
    enddo

end subroutine

subroutine cs_sg(nsp, dx, ysp, bc)
    ! construct spline on structrued 1d grid
    use precision
    implicit none

    integer, intent(in) :: nsp
    real(wp), intent(in) :: dx
    real(wp), intent(inout) :: ysp(0:3, 0:nsp-1)
    character(*), intent(in) :: bc

    integer i
    real(wp) xsp(0:nsp-1)

    xsp = (/(real(i, wp)*dx, i = 0, nsp-1)/)

    call construct_spline(nsp, xsp, ysp, bc)

end subroutine cs_sg

subroutine construct_spline2d(nx, ny, xsp, ysp, zsp, xbc, ybc)
    ! construct spline on 2d grid
    use precision
    implicit none

    integer, intent(in)  :: nx, ny
    real(wp), intent(in) :: xsp(0:nx-1), ysp(0:ny-1)
    real(wp), intent(inout) :: zsp(0:15, 0:nx-1, 0:ny-1)
    character(*), intent(in):: xbc, ybc

    integer i, j
    real(wp) tmp_x(0:3, 0:nx-1), tmp_y(0:3, 0:ny-1)
    !TODO: check it
    do i = 0, ny-1
        tmp_x(0, :) = zsp(0, :, i)
        call construct_spline(nx, xsp, tmp_x, xbc)
        zsp(1, :, i) = tmp_x(1, :)
        zsp(2, :, i) = tmp_x(2, :)
        zsp(3, :, i) = tmp_x(3, :)
    enddo

    do i = 0, nx-1
        do j = 0, 3
            tmp_y(0, :) = zsp(j, i, :)
            call construct_spline(ny, ysp, tmp_y, ybc)
            zsp(j+4, i, :) = tmp_y(1, :)
            zsp(j+8, i, :) = tmp_y(2, :)
            zsp(j+12, i, :) = tmp_y(3, :)
        enddo
    enddo

end subroutine construct_spline2d

function spline(x0, nsp, xsp, ysp) result(y0)
    use precision
    implicit none

    integer i, nsp
    real(wp) x0, y0, xsp(0:nsp-1), ysp(0:3, 0:nsp-1), dx

    do i = 0, nsp-2
        if(x0 >= xsp(i) .and. x0 < xsp(i+1))then
            dx = x0 - xsp(i)
            y0 = ysp(0, i) + ysp(1, i)*dx + ysp(2, i)*dx*dx +&
                ysp(3, i)*dx*dx*dx
            exit
        endif
    enddo

end function spline

function dspline(x0, nsp, xsp, ysp) result(y0)
    use precision
    implicit none

    integer i, nsp
    real(wp) x0, y0, xsp(0:nsp-1), ysp(0:3, 0:nsp-1), dx

    do i = 0, nsp-2
        if(x0 >= xsp(i) .and. x0 < xsp(i+1))then
            dx = x0 - xsp(i)
            y0 = ysp(1, i) + ysp(2, i)*2.0*dx + ysp(3, i)*3.0*dx*dx
            exit
        endif
    enddo

end function dspline

function spline2d(x0, y0, nx, ny, xsp, ysp, zsp) result(z0)
    use precision
    implicit none

    integer i, j, k, nx, ny
    real(wp) x0, y0, z0, dx, dy, xsp(0:nx-1), ysp(0:ny-1),&
        zsp(0:15, 0:nx-1, 0:ny-1), ds(0:15)
    !TODO: check it
    do i = 0, nx-2
        if(x0 >= xsp(i) .and. x0 < xsp(i+1))dx = x0 - xsp(i)
    enddo
    do j = 0, ny-2
        if(y0 >= ysp(i) .and. y0 < ysp(i+1))dy = y0 - ysp(i)
    enddo

    ds = (/(dx**(mod(k, 4))*dy**(k/4), k=0,15)/)
    z0 = dot_product(zsp(:, i, j), ds)

end function spline2d

function dsplin2d(xy, x0, y0, nx, ny, xsp, ysp, zsp) result(z0)
    ! xy variable indicates derivative w.r.t. x(xy==1) or y(xy==2)
    use precision
    implicit none

    integer i, j, k, nx, ny, xy
    real(wp) x0, y0, z0, dx, dy, xsp(0:nx-1), ysp(0:ny-1),&
        zsp(0:15, 0:nx-1, 0:ny-1), ds(0:15)
    do i = 0, nx-2
        if(x0 >= xsp(i) .and. x0 < xsp(i+1))dx = x0 - xsp(i)
    enddo
    do j = 0, ny-2
        if(y0 >= ysp(i) .and. y0 < ysp(i+1))dy = y0 - ysp(i)
    enddo

    if(xy == 1)then
        ds(0:3) = (/0.0, 1.0, 2.0*dx, 3.0*dx*dx/)
        ds(4:7) = ds(0:3)*dy
        ds(8:11) = ds(4:7)*dy
        ds(12:15) = ds(8:11)*dy
    elseif(xy == 2)then
        ds(0:3) = 0.0
        ds(4:7) = (/(dx**k, k = 0,3)/)
        ds(8:11) = ds(4:7)*2.0*dy
        ds(12:15) = ds(8:11)*1.5*dy
    else
        STOP 'dspline2d parameter error'
    endif

    z0 = dot_product(zsp(:, i, j), ds)
    !TODO: check it
    
end function