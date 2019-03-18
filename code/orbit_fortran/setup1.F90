subroutine setup

  use global_parameters
  use equilibrium
  implicit none

  integer nw, nh, fid, edge_r, edge_z, fid, i, iiii, limitr, maxis_r, &
          maxis_z, nbbbs, xdum, zmid

  real(wp) current, rcentr, rdim, rleft, rmaxis, sibry, simag, zmaxis, &
           zdim
  
  real(wp), dimension(:), allocatable :: fpol, norm_flux, norm_flux0, &
                                         norm_flux1, rbbbs, rbbbs_d, rlim, &
                                         x, y, zbbbs, zbbbs_d, rlim, zlim, &
                                         psirz1, range_efit, temp, x0, y0

  real(wp), dimension(:,:), allocatable :: flux_efit, flux_efit0, &
                                           flux_efit11, psirz, Rxy, Zxy, &

 ! load gfile
   fid = 233
   open(fid, file = "g096000.00100", status = "replace")
   read(fid,*)temp
   nw = temp(42)
   nh = temp(43)
   range_efit = temp(44 : 48)
   rdim = temp(44)
   zdim = temp(45)
   rcentr = temp(46)
   rleft = temp(47)
   zmid = temp(48)
   rmaxis = temp(49)
   zmaxis = temp(50)
   simag = temp(51)
   sibry = temp(52)
   bcentr = temp(53)
   current = temp(54)
   simag = temp(55)
   xdum = temp(56)
   rmaxis = temp(57)
   xdum = temp(58)
   zmaxis = temp(59)
   xdum = temp(60)
   sibry = temp(61)
   xdum = temp(62)
   xdum = temp(63)
   fpol = temp(64 : 64 + nw - 1);
   psirz1 = temp(64 + 4 * nw : 63 + 4 * nw + nw * nh)
   nbbbs = temp(64 + 4 * nw + nw * nh + nw)
   limitr = temp(65 + 4 * nw + nw * nh + nw)
   rbbbs = temp(66 + 4 * nw + nw * nh + nw : & 
           65 + 4 * nw + nw * nh + nw + nbbbs)
   zbbbs = temp(66 + 4 * nw + nw * nh + nw + nbbbs : & 
           65 + 4 * nw + nw * nh + nw + 2 * nbbbs)
   rlim = temp(66 + 4 * nw + nw * nh + nw + 2 * nbbbs: & 
           65 + 4 * nw + nw * nh + nw + 2 * nbbbs + limitr)
   zlim = temp(66 + 4 * nw + nw * nh + nw + 2 * nbbbs + limitr: & 
           65 + 4 * nw + nw * nh + nw + 2 * nbbbs + 2 * limitr) 
   do i = 1, nw
      flux_efit0(:,i) = psirz1((i - 1) * nh + 1 : i * nh)
   enddo
   psirz = flux_efit0

! interpolation

   x0 = range_efit(4) : range_efit(1) / (nw-1) : range_efit(4) + &
       range_efit(1)
   y0 = -range_efit(2) / 2 + range_efit(5):range_efit(2)/(nh-1) : &
       range_efit(2)/2+range_efit(5)
  
   nw = 511
   nh = 511

   x = range_efit(4) : range_efit(1) / (nw-1) : range_efit(4) + &
       range_efit(1)
   y = -range_efit(2) / 2 + range_efit(5):range_efit(2)/(nh-1) : &
       range_efit(2)/2+range_efit(5)

!   flux_efit = interp(flux_efit0,x0,y0,x,y)


   rbbbs_d = round((rbbbs - range_efit(4)) * (nw-1) / range_efit(1)) + 1
   zbbbs_d = round((zbbbs + (range_efit(2) - range_efit(5)) / 2 ) * &
             (nh-1)/range_efit(2)) + 1

   maxis_r = round((rmaxis - range_efit(4)) * (nh - 1) / range_efit(1)) + 1
   maxis_z = round((zmaxis + (range_efit(2) - range_efit(5)) / 2) * $
             (nh-1) / range_efit(2)) + 1

   do iiii = round(length(zbbbs) / 2) - 10, length(zbbbs)
      if ((zbbbs(iiii) - zmaxis) * (zmaxis - zbbbs(iiii - 1)) >= 0) then
         exit
      endif
   enddo
   iiii=iiii-1;
   edg_r = round((rbbbs(iiii) - range_efit(4)) * (nh-1) / range_efit(1)) + 1
   edg_z = round((zbbbs(iiii) + (range_efit(2) - range_efit(5)) / 2) * &
           (nh-1) / range_efit(2)) + 1

   psirz = transpose(flux_efit)
   !fpol_tmp=interp(fpol,x0,x)
    
   flux_profile1 = transpose(psirz(maxis_z, maxis_r : edg_r))
   !flux_profile = interp(flux_profile1,x0,x)
   !p = polyfit(flux_profile,fpol_tmp,4)

    do i = 1, nw
        do j = 1, nh
            !fpol_rz(i,j)=polyval(p,psirz(i,j));
        enddo
    enddo

    do i = 1, nw
       do j = 1, nh
          if (fpol_rz(i,j) < min(fpol) / 100) then
              fpol_rz(i,j)=min(fpol)/100
          endif
       enddo
    enddo

    Btxy=fpol_rz./Rxy;


! Brxy = diff(psirz,y)
! Bzxy = diff(psirz,x)












