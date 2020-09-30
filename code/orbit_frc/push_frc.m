

for im = 1 : nw
    x0 = xp(im + (i-1)*nw);
    y0 = yp(im + (i-1)*nw);
    z0 = zp(im + (i-1)*nw);
    v0x = vxp(im + (i-1)*nw);
    v0y = vyp(im + (i-1)*nw);
    v0z = vzp(im + (i-1)*nw);
    
    RK4_frc;
    
    xp(im+1 + (i-1)*nw) = x0 + dx;
    yp(im+1 + (i-1)*nw) = y0 + dy;
    zp(im+1 + (i-1)*nw) = z0 + dz;
    vxp(im+1 + (i-1)*nw) = v0x + dvx;
    vyp(im+1 + (i-1)*nw) = v0y + dvy;
    vzp(im+1 + (i-1)*nw) = v0z + dvz;
end
    
    