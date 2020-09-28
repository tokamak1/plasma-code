bx = 0; by = sqrt(x0^2+y0^2+z0^2) - 0.2; bz = 0;
h = dtp;
for k = 1 : 8
    [bx0, by0, bz0] = Bcoil(k, x0, y0, z0, I0, omega, t);
    bx = bx + bx0;
    by = by + by0;
    bz = bz + bz0;
end
a0x = q * (v0y * bz - v0z * by) / m;
a0y = q * (v0z * bx - v0x * bz) / m;
a0z = q * (v0x * by - v0y * bx) / m;
v1x = v0x + a0x * h / 2;
v1y = v0y + a0y * h / 2;
v1z = v0z + a0z * h / 2;
xp1 = x0 + v1x * h / 2;
yp1 = y0 + v1y * h / 2;
zp1 = z0 + v1z * h / 2;
bx = 0; by = 0; bz = 0;
for k = 1 : 8
    [bx0, by0, bz0] = Bcoil(k, xp1, yp1, zp1, I0, omega, t);
    bx = bx + bx0;
    by = by + by0;
    bz = bz + bz0;
end
a1x = q * (v1y * bz - v1z * by) / m;
a1y = q * (v1z * bx - v1x * bz) / m;
a1z = q * (v1x * by - v1y * bx) / m;
v2x = v0x + a1x * h / 2;
v2y = v0y + a1y * h / 2;
v2z = v0z + a1z * h / 2;
xp2 = x0 + v2x * h / 2;
yp2 = y0 + v2y * h / 2;
zp2 = z0 + v2z * h / 2;
bx = 0; by = 0; bz = 0;
for k = 1 : 8
    [bx0, by0, bz0] = Bcoil(k, xp2, yp2, zp2, I0, omega, t);
    bx = bx + bx0;
    by = by + by0;
    bz = bz + bz0;
end
a2x = q * (v2y * bz - v2z * by) / m;
a2y = q * (v2z * bx - v2x * bz) / m;
a2z = q * (v2x * by - v2y * bx) / m;
v3x = v0x + a2x * h;
v3y = v0y + a2y * h;
v3z = v0z + a2z * h;
xp3 = x0 + v3x * h;
yp3 = y0 + v3y * h;
zp3 = z0 + v3z * h;
bx = 0; by = 0; bz = 0;
for k = 1 : 8
    [bx0, by0, bz0] = Bcoil(k, xp3, yp3, zp3, I0, omega, t);
    bx = bx + bx0;
    by = by + by0;
    bz = bz + bz0;
end
a3x = q * (v3y * bz - v3z * by) / m;
a3y = q * (v3z * bx - v3x * bz) / m;
a3z = q * (v3x * by - v3y * bx) / m;
dvx = (a0x + 2 * a1x + 2 * a2x + a3x) * h / 6;
dvy = (a0y + 2 * a1y + 2 * a2y + a3y) * h / 6;
dvz = (a0z + 2 * a1z + 2 * a2z + a3z) * h / 6;
dx = (v0x + 2 * v1x + 2 * v2x + v3x) * h / 6;
dy = (v0y + 2 * v1y + 2 * v2y + v3y) * h / 6;
dz = (v0z + 2 * v1z + 2 * v2z + v3z) * h / 6 ;
