

for m = 1 : nw
    x0 = xp(m + (i-1)*nw);
    y0 = yp(m + (i-1)*nw);
    z0 = zp(m + (i-1)*nw);
    v0x = vxp(m + (i-1)*nw);
    v0y = vyp(m + (i-1)*nw);
    v0z = vzp(m + (i-1)*nw);
    
    RK4_frc;
    
    xp(m+1 + (i-1)*nw) = x0 + dx;
    yp(m+1 + (i-1)*nw) = y0 + dy;
    zp(m+1 + (i-1)*nw) = z0 + dz;
    vxp(m+1 + (i-1)*nw) = v0x + dvx;
    vyp(m+1 + (i-1)*nw) = v0y + dvy;
    vzp(m+1 + (i-1)*nw) = v0z + dvz;
end
    
    