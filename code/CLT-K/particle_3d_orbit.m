figure;
set(gca,'nextplot','replacechildren');
set(gcf,'DefaultAxesFontSize',15);
set(gcf,'Position',get(0,'ScreenSize'));
for i = 1 : length(tag)
    if mod(i,150) == 0
    eval(['R = pdata',num2str(tag(i)),'(1,:);']);
    eval(['Z = pdata',num2str(tag(i)),'(3,:);']);
    eval(['Phi = pdata',num2str(tag(i)),'(2,:);']);
    X = R .* cos(Phi);
    Y = R .* sin(Phi);
    plot3(X,Y,Z,'linewidth',2);
    eval(['time = pdata',num2str(tag(i)),'(8,:);']);
    if length(R) == nstop / ndiag + 1
        timep = time;
    end
    xlabel('X','fontsize',18);ylabel('Y','fontsize',18);zlabel('Z','fontsize',18);
    title(['particle ',num2str(tag(i)),' t = 0 ~ ',num2str(time(end))]);
    xlim([-6, 6]);ylim([-6, 6]);zlim([-4, 4]);
    drawnow;
    pause(0.1);
    F=getframe(gcf);
    saveas(gcf, ['p',num2str(tag(i))], 'jpg');
    end
end