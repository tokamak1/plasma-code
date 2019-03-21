figure;
ng = 128;
jp = 1;
dp = 15 / (ng - 1);
dl = 2.2 * pi / (ng - 1);
P0 = 0 : dp : 15;
L0 = 0 : dl : 2.2 * pi;
[P, L] = meshgrid(P0, L0);
for it = 2 : 50 : nstop / ndiag + 1
%     figure; 
    ely = zeros(ng);
    n0 = ones(ng);
%     set(gca,'nextplot','replacechildren');
    set(gcf,'DefaultAxesFontSize',15);
    set(gcf,'Position',get(0,'ScreenSize'));
    for i = 1 : 1 : length(tag)
        eval(['Pphi = pdata',num2str(tag(i)),'(5, :);']);
        eval(['Phi = pdata',num2str(tag(i)),'(2, :);']);
        eval(['vpara = pdata',num2str(tag(i)),'(4, :);']);

        if it > length(Phi)
            continue;
        else
            gbp = floor(abs(Pphi(it - 1)) / dp) + 1; 
            gfp = gbp + 1;
            wbp = 1-abs(abs(Pphi(it - 1)) / dp - gbp);
            wfp = 1 - wbp;
            gbl = floor(Phi(it - 1) / dl) + 1; 
            gfl = gbl + 1;
            wbl = 1-abs(Phi(it - 1) / dl - gbl);
            wfl = 1 - wbl;
            rp = Pphi(it) - Pphi(it-1);
            rl = Phi(it) - Phi(it - 1) - vpara(it) * timep(it) + vpara(it - 1) * timep(it - 1);
            rl = 0;
            r0 = sqrt(rp .* rp + rl .* rl);
            ely0 = max(r0);
            ely(gbp,gbl) = ely(gbp,gbl) + ely0 * wbp * wbl;
            ely(gbp,gfl) = ely(gbp,gfl) + ely0 * wbp * wfl;
            ely(gfp,gbl) = ely(gfp,gbl) + ely0 * wfp * wbl;
            ely(gfp,gfl) = ely(gfp,gfl) + ely0 * wfp * wfl;
            n0(gbp,gbl) = n0(gbp,gbl) + 1;
            n0(gbp,gfl) = n0(gbp,gfl) + 1;
            n0(gfp,gbl) = n0(gbp,gbl) + 1;
            n0(gfp,gfl) = n0(gbp,gbl) + 1;
        end
        pcolor(L,P,ely ./ n0);shading interp;
        xlim([0, 2 * pi]);ylim([0, 15]);
        xlabel('\Phi|','fontsize',18);
        ylabel('P_\phi','fontsize',18);
        title(['Lyapunov exponent t = ',num2str(timep(it - 1))]);
%         caxis([-2,2]);
        colorbar();
%         drawnow;
    end
    F3(jp) = getframe(gcf);
    jp = jp + 1;
    

end