
 RSmin = min(min(RS));
 contour(xx,yy,flux_efit',100);hold on;
 if ilim
 plot(rlim,zlimter,'r','linewidth',2);
 end
 
 hold on;shading interp;axis equal;
 saveas(gcf,'./output/surface_fig/2D_flux.jpg');
 figure;
for i = 1 : npsi
     CS = RSmin ./ RS(nsa(i) : nsa(i + 1) - 1, 1 : end / 1.5);
     h = surf(XS(nsa(i) : nsa(i + 1) - 1,1 : end / 1.5),YS(nsa(i) : nsa(i + 1) - 1, 1 : end / 1.5),...
         ZS(nsa(i) : nsa(i + 1) - 1, 1 : end / 1.5),- CS); 
     alpha (h,0.3 + 0.05 * i);
     shading interp;hold on;axis equal;
     colormap(jet);
end
if ilim
hl = surf(xlimi,ylimi,zlimi,- RSmin ./ rlimi, 'LineStyle', ':');
alpha (hl,0);
end
saveas(gcf,'./output/surface_fig/3D_surface.jpg'); 
