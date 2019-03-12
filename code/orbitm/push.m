mu = []; E = [];
E(1) = K;
mu(1) = m * vper^2 / 2 / b0;
b = B(Brxy,Btxy,Bzxy,fr,fz,Rxy,Zxy);
br = interp2(Rxy,Zxy,Brxy,fr,fz);
bt = interp2(Rxy,Zxy,Btxy,fr,fz);
bz = interp2(Rxy,Zxy,Bzxy,fr,fz);
bx = br * fx / fr - bt * fy / fr;
by = br * fy / fr + bt * fx / fr;
vx = (vpara * bx + vper * bx * bz / sqrt(bx^2 + by^2)) / b;
vy = (vpara * by + vper * by * bz / sqrt(bx^2 + by^2)) / b;
vz = (vpara * bz - vper * sqrt(bx^2 + by^2)) / b;

Xf = [];
Zf = [];
Yf = [];
Vx = [];
Vy = [];
Vz = [];

for i = 1 : nt
    Xf(i) = fx;
    Yf(i) = fy;
    Zf(i) = fz;
    Vx(i) = vx;
    Vy(i) = vy;
    Vz(i) = vz;
    
    [dvx,dvy,dvz,dx,dy,dz] = RK4ob(fx,fy,fz,vx,vy,vz,Brxy,Btxy,Bzxy,Rxy,Zxy,dt,q,m);
    
    vx = vx + dvx;
    vy = vy + dvy;
    vz = vz + dvz;

    fx = fx + dx;
    fy = fy + dy;
    fz = fz + dz;

    b = B(Brxy,Btxy,Bzxy,fx,fz,Rxy,Zxy);
    mu(i + 1) = m * (vx^2 + vz^2)/ 2 / b;
    E(i + 1) = m * (vx^2 + vz^2 + vy^2);
end
figure;
contour(xx,yy,flux_efit',100);
hold on;shading interp;axis equal;
Rf = sqrt(Xf .* Xf + Yf .* Yf);
plot(Rf,Zf,'r','LineWidth',2);hold on;
plot(Rf(1),Zf(1),'ob');
if ilim
plot(rlim,zlimter,'k','linewidth',2);
end
saveas(gcf,'./output/orbit_fig/2D_orbit.jpg');

figure;
nfig = 1;
CS = RSmin ./ RS(nsa(nf):nsa(nf + 1) - 1,1 : end / 1.5);
h = surf(XS(nsa(nf):nsa(nf + 1) - 1,1 : end / 1.5),YS(nsa(nf):nsa(nf + 1) - 1,1 : end / 1.5),...
ZS(nsa(nf):nsa(nf + 1) - 1,1 : end / 1.5),- CS); 
alpha (h,0.8);shading interp;hold on;
if ilim
surf(xlimi,ylimi,zlimi,- RSmin ./ rlimi, 'LineStyle', ':','facealpha',0);hold on;
end
plot3(Xf(1 : nfig : end),Yf(1 : nfig : end),Zf(1 : nfig : end),'-b','LineWidth',2);
hold on;axis equal;axis off;
plot3(Xf(1),Yf(1),Zf(1),'or','linewidth',2);hold on;axis equal;
colormap(jet);

ORBITxyz = [Xf;Yf;Zf];
saveas(gcf,'./output/orbit_fig/3D_orbit.jpg');
save('./output/ORBITxyz.dat','ORBITxyz','-ascii','-double','-tabs');
