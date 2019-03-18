subroutine setup

    use global_parameters
    use equilibrium
    implicit none

    integer nw, nh, edge_r, edge_z, fid, i, j, iiii, limitr, maxis_r, &
            maxis_z, nbbbs, xdum, zmid, nw1, nh1, rl

    real(wp) current, rcentr, bcentr, rdim, rleft, rmaxis, zmaxis, zdim, minfpol
    character(len = 40) string
             
  
    real(wp), dimension(:), allocatable :: norm_flux, norm_flux0, flux_profile, &
                                         norm_flux1, rbbbs, rbbbs_d, temp0, &
                                         zbbbs, zbbbs_d, rlim, zlim, temp1, &
                                         psirz1, range_efit, temp, x0, y0

    real(wp), dimension(:,:), allocatable :: flux_efit11, psirz, Rxy, Zxy, & 
                                           fpol, fpol_tmp, fpol_rz

    real(wp), dimension(:,:,:), allocatable :: flux_efit0 
    real(wp), external :: spline, spline2d

    character(len = 7), parameter :: bc = 'natural'

    allocate(range_efit(5))
 ! load gfile
    fid = 233
    allocate(temp0(23))
    open(fid, file = "g096000.00100")
    read(fid,*)temp0
    nw = temp0(2)
    nh = temp0(3)
    range_efit = temp0(4 : 8)
    rdim = temp0(4)
    zdim = temp0(5)
    rcentr = temp0(6)
    rleft = temp0(7)
    zmid = temp0(8)
    rmaxis = temp0(9)
    zmaxis = temp0(10)
    simag = temp0(11)
    sibry = temp0(12)
    bcentr = temp0(13)
    current = temp0(14)
    simag = temp0(15)
    xdum = temp0(16)
    rmaxis = temp0(17)
    xdum = temp0(18)
    zmaxis = temp0(19)
    xdum = temp0(20)
    sibry = temp0(21)
    xdum = temp0(22)
    xdum = temp0(23)
    allocate(temp(4 * nw + nw * nh + nw + 2))
    read(fid, *)temp
    allocate(psirz1(nh * nw), fpol_tmp(0:3, nw))
    fpol_tmp(0, :) = temp(1 : nw);
    psirz1 = temp(1 + 4 * nw : 4 * nw + nw * nh)
    nbbbs = temp(1 + 4 * nw + nw * nh + nw)
    limitr = temp(2 + 4 * nw + nw * nh + nw)
    allocate(temp1(2 * (nbbbs + limitr)))
    read(fid, *)temp1
    allocate(rbbbs(nbbbs), zbbbs(nbbbs), rlim(limitr), zlim(limitr))
    rbbbs = temp1(1 : nbbbs)
    zbbbs = temp1(1 + nbbbs : 2 * nbbbs)
    rlim = temp1(1 + 2 * nbbbs : 2 * nbbbs + limitr)
    zlim = temp1(1 + 2 * nbbbs + limitr : 2 * nbbbs + 2 * limitr) 
    close(fid)
    allocate(flux_efit0(0:15, nh, nw))
    do i = 1, nw
       flux_efit0(0, :, i) = psirz1((i - 1) * nh + 1 : i * nh)
    enddo
    allocate(psirz(nw,nh))
    psirz = flux_efit0(0, :, :)

    nw1 = 511
    nh1 = 511

! interpolation
    allocate(x0(nw), y0(nh), x(nw1), y(nh1), fpol(0:3, nw1))
   !x0 = range_efit(4) : range_efit(1) / (nw-1) : range_efit(4) + &
   !    range_efit(1)
   !y0 = -range_efit(2) / 2 + range_efit(5):range_efit(2)/(nh-1) : &
   !    range_efit(2)/2+range_efit(5)
    do i = 1, nw
       x0(i) = range_efit(4) + (i - 1) * range_efit(1) / (nw - 1)
    enddo
    do i = 1, nh 
       y0(i) = -range_efit(2) / 2 + range_efit(5) + (i - 1) * range_efit(2) / (nh - 1)
    enddo
   call construct_spline2d(nh, nw, y0, x0, flux_efit0, bc, bc)

   !x = range_efit(4) : range_efit(1) / (nw-1) : range_efit(4) + &
   !    range_efit(1)
   !y = -range_efit(2) / 2 + range_efit(5):range_efit(2)/(nh-1) : &
   !    range_efit(2)/2+range_efit(5)
    do i = 1, nw1
       x(i) = range_efit(4) + (i - 1) * range_efit(1) / (nw1 - 1)
    enddo
    do i = 1, nh1 
       y(i) = -range_efit(2) / 2 + range_efit(5) + (i - 1) * range_efit(2) / (nh1 - 1)
    enddo
    allocate(Rxy(nw1, nh1), Zxy(nw1, nh1), flux_efit(nw1,nh1))
    do i = 1, nw1
        do j = 1, nh1
 !           call construt_spline2d(nh, nw, y0, x0, flux_efit0, bc, bc)
            flux_efit(i, j) = spline2d(y(i), x(j), nh1, nw1, y0, x0, flux_efit0)
            Rxy(i, j) = x(i)
            Zxy(i, j) = y(j)
        enddo
    enddo
!   flux_efit = interp(flux_efit0,x0,y0,x,y)
    allocate(rbbbs_d(nbbbs),zbbbs_d(nbbbs))
    do i = 1, nbbbs
       rbbbs_d(i) = anint((rbbbs(i) - range_efit(4)) * (nw1-1) / range_efit(1)) + 1
       zbbbs_d(i) = anint((zbbbs(i) + (range_efit(2) - range_efit(5)) / 2 ) * &
                   (nh1-1)/range_efit(2)) + 1
    enddo

    maxis_r = anint((rmaxis - range_efit(4)) * (nh1 - 1) / range_efit(1)) + 1
    maxis_z = anint((zmaxis + (range_efit(2) - range_efit(5)) / 2) * & 
                   (nh1-1) / range_efit(2)) + 1

    do iiii = anint(real(size(zbbbs) / 2)) - 10, size(zbbbs)
       if ((zbbbs(iiii) - zmaxis) * (zmaxis - zbbbs(iiii - 1)) >= 0) then
          exit
       endif
    enddo
    iiii=iiii-1;
    edge_r = anint((rbbbs(iiii) - range_efit(4)) * (nh-1) / range_efit(1)) + 1
    edge_z = anint((zbbbs(iiii) + (range_efit(2) - range_efit(5)) / 2) * &
            (nh-1) / range_efit(2)) + 1

    psirz = transpose(flux_efit)

    allocate(fpol(0:3, maxis_r:edge_r))
    call construct_spline(nw-1, x0, fpol_tmp, bc)
    minfpol = fpol(0,maxis_r)
    do i = maxis_r, edge_r
!        call construt_spline(nw-1, x0, fpol_tmp, bc)
        fpol(0,i) = spline(x(i), x0, nw, fpol_tmp)
        if (i > maxis_r .and. fpol(0,i) < fpol(0,i - 1)) &
           minfpol = fpol(0,i)
    enddo
   !fpol_tmp=interp(fpol,x0,x)
    
    rl = edge_r - maxis_r
    allocate(flux_profile(rl + 1))
    flux_profile = psirz(maxis_r : edge_r, maxis_z)
   call construct_spline(rl, flux_profile, fpol, bc)
   !p = polyfit(flux_profile,fpol_tmp,4)
    allocate(Btxy(nw1, nh1), Brxy(nw1, nh1), Bzxy(nw1, nh1), fpol_rz(nw1,nh1))
    do i = 1, nw1
        do j = 1, nh1
!            call construct_spline(rl, flux_profile, fpol, bc)
            fpol_rz(i, j) = max(minfpol / 100, &
                            spline(psirz(i, j), rl, flux_profile, fpol))
            Btxy(i, j) = fpol_rz(i, j) / Rxy(i, j)
            !fpol_rz(i,j)=polyval(p,psirz(i,j));
        enddo
    enddo

    !do i = 1, nw
    !   do j = 1, nh
    !      if (fpol_rz(i,j) < min(fpol) / 100) then
    !          fpol_rz(i,j)=min(fpol)/100
    !      endif
    !   enddo
    !enddo



! Brxy = diff(psirz,y)
! Bzxy = diff(psirz,x)

    end subroutine setup










