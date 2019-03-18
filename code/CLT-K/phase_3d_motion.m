figure;
jp = 1;
for it = 1 : 100 : nstop / ndiag
%     figure; 
    set(gca,'nextplot','replacechildren');
    set(gcf,'DefaultAxesFontSize',15);
    set(gcf,'Position',get(0,'ScreenSize'));
    for i = 1 : length(tag)
        if mod(i,15) == 0
        eval(['Pphi = pdata',num2str(tag(i)),'(5, :);']);
        eval(['E = pdata',num2str(tag(i)),'(6, :);']);
        eval(['lamda = pdata',num2str(tag(i)),'(7, :);']);

            if it >= length(E)
                plot3(Pphi,E,lamda,'-o','linewidth',2);hold on;
%                 pause(0.1);
            else
                plot3(Pphi(1:it),E(1:it),lamda(1:it),'-o','linewidth',2);hold on;
            end
            xlim([-15, 15]);ylim([0, 3.5]);zlim([0, 1.4]);
            xlabel('P_\phi','fontsize',18);
            ylabel('E','fontsize',18);
            zlabel('\Lambda','fontsize',18);

        end
            title(['particle phase motion t = ',num2str(timep(it))]);
    end
    Fp(jp) = getframe(gcf);
    jp = jp + 1;

end  