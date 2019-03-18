figure;
jp = 1;
for it = 1 : 100 : nstop / ndiag
%     figure; 
    set(gca,'nextplot','replacechildren');
    set(gcf,'DefaultAxesFontSize',15);
    set(gcf,'Position',get(0,'ScreenSize'));
    for i = 1 : length(tag)
        if mod(i,50) == 0
        eval(['Pphi = pdata',num2str(tag(i)),'(5, :);']);
        eval(['Phi = pdata',num2str(tag(i)),'(2, :);']);

            if it >= length(Pphi)
                plot(Phi(1 : 10 : end),abs(Pphi(1 : 10 : end)),'o','linewidth',2);hold on;
                plot(Phi(end),abs(Pphi(end)),'*k','linewidth',2);hold on;
%                 pause(0.1);
            else
                plot(Phi(1 : 10 : it),abs(Pphi(1 : 10 : it)),'o','linewidth',2);hold on;
                plot(Phi(it),abs(Pphi(it)),'*k','linewidth',2);hold on;
            end
            xlim([0, 2 * pi]);ylim([0, 15]);
            ylabel('|P_\phi|','fontsize',18);
            xlabel('\Phi','fontsize',18);

        end
            title(['particle phase motion t = ',num2str(timep(it))]);
    end
    F1(jp) = getframe(gcf);
    jp = jp + 1;

end