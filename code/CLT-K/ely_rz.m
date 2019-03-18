figure;
ng = 128;
jp = 1;
dr = 3 / (ng - 1);
dz = 3 / (ng - 1);
r0 = -1.5 : dr : 1.5;
z0 = 2.5 : dr : 5.5;
[r, z] = meshgrid(r0, z0);
for it = 2 : 50 : nstop / ndiag + 1
%     figure; 
    ely = zeros(ng);
    n0 = ones(ng);
%     set(gca,'nextplot','replacechildren');
    set(gcf,'DefaultAxesFontSize',15);
    set(gcf,'Position',get(0,'ScreenSize'));
    for i = 1 : 1 : length(tag)
        eval(['R = pdata',num2str(tag(i)),'(1, :);']);
        eval(['Z = pdata',num2str(tag(i)),'(3, :);']);

        if it > length(R)
            continue;
        else
            gbr = floor(abs(R(it - 1)) / dr) + 1; 
            gfr = gbr + 1;
            wbr = 1-abs(abs(R(it - 1)) / dr - gbr);
            wfr = 1 - wbr;
            gbz = floor(Z(it - 1) / dz) + 1; 
            gfz = gbl + 1;
            wbz = 1-abs(Z(it - 1) / dz - gbz);
            wfz = 1 - wbz;
            rr = R(it) - R(it-1);
            rz = Z(it) - Z(it - 1);
            r0 = sqrt(rr .* rr + rz .* rz);
            ely0 = max(r0);
            ely(gbr,gbz) = ely(gbr,gbz) + ely0 * wbr * wbz;
            ely(gbr,gfz) = ely(gbr,gfz) + ely0 * wbr * wfz;
            ely(gfp,gbl) = ely(gfl,gbp) + ely0 * wfr * wbz;
            ely(gfp,gfl) = ely(gfp,gfl) + ely0 * wfr * wfz;
            n0(gbr,gbz) = n0(gbr,gbz) + 1;
            n0(gbr,gfz) = n0(gbr,gfz) + 1;
            n0(gfr,gbz) = n0(gbr,gbz) + 1;
            n0(gfr,gfz) = n0(gbr,gbz) + 1;
        end
        pcolor(r,z,ely ./ n0);shading interp;
        xlim([0, 2 * pi]);ylim([0, 15]);
        xlabel('R','fontsize',18);
        ylabel('Z','fontsize',18);
        title(['Lyapunov exponent t = ',num2str(timep(it - 1))]);
%         caxis([-2,2]);
        colorbar();
%         drawnow;
    end
    Frz(jp) = getframe(gcf);
    jp = jp + 1;
    

end