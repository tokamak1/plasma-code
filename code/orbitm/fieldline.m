Xl = [];
Yl = [];
Zl = [];

for i = 1 : nl
    Xl(i) = lx;
    Yl(i) = ly;
    Zl(i) = lz;

    [dx,dy,dz] = RK4(lx,ly,lz,Brxy,Btxy,Bzxy,Rxy,Zxy,dh);

    lx = lx + dx;
    ly = ly + dy;
    lz = lz + dz;

end
FLxyz = [Xl;Yl;Zl];

figure;
contour(xx,yy,flux_efit',100);hold on;shading interp;
hold on;axis equal;axis off;
Rl = sqrt(Xl .* Xl + Yl .* Yl);
plot(Rl,Zl,'r','LineWidth',2);hold on;
if ilim
plot(rlim,zlimter,'k','linewidth',2);
end
hold on;axis equal;axis off;
saveas(gcf,'./output/surface_fig/2D_fieldline.jpg');
nfig = 1;

figure;
CS = 0 * RS(nsa(nf) : nsa(nf + 1) - 1,1 : end / 1.5) + psis_n(nf);
h = surf(XS(nsa(nf) : nsa(nf + 1) - 1,1 : end / 1.5), YS(nsa(nf) : nsa(nf + 1) - 1,1 : end / 1.5),...
ZS(nsa(nf) : nsa(nf + 1) - 1,1 : end / 1.5),- CS); 
alpha (h,0.5);shading interp;hold on;axis equal;
colormap(spring);
plot3(Xl,Yl,Zl,'LineWidth',2);hold on;
if ilim
surf(xlimi,ylimi,zlimi,- RSmin ./ rlimi, 'LineStyle', ':','facealpha',0);
end
saveas(gcf,'./output/surface_fig/3D_fieldline.jpg');
save('./output/FLxyz.dat','FLxyz','-ascii','-double','-tabs');
