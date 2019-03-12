clear;
close all;
clc;
partdat = './particle_orbit.dat';
fid = fopen(partdat,'rb');
nstop = 400000;
ndiag = 40;
ntrack = 4000;
itag = 1;
for i = 1 : (nstop / ndiag + 1) * ntrack
    idata1 = fread(fid,8,'double');
    idata2 = fread(fid,2,'int');    
    if sum(idata2) == 0
        break;
    end
    if itag
        inp = idata2(2);
        tag(inp) = idata2(1);
        eval(['pdata',num2str(idata2(1)),'= idata1;']);
        if inp == 1 && i == 1
           tag0 = idata2(1);
           pdata0 = idata1;
        end
            
        if inp == 1 && i > 1
            itag = 0;
            tag(1) = tag0;
            eval(['pdata',num2str(idata2(1)),'= [pdata0 , idata1];']);
        end
    else
        eval(['pdata',num2str(idata2(1)),'= [pdata',num2str(idata2(1)),...
             ' , idata1];']);
    end      
end
figure;
set(gca,'nextplot','replacechildren');
set(gcf,'DefaultAxesFontSize',15);
set(gcf,'Position',get(0,'ScreenSize'));
jp = 1;
for i = 1 : length(tag)
    if mod(i,15) == 0
%     if i == 300
    eval(['R = pdata',num2str(tag(i)),'(1,:);']);
    eval(['Z = pdata',num2str(tag(i)),'(3,:);']);
    eval(['Phi = pdata',num2str(tag(i)),'(2,:);']);
    X = R .* cos(Phi);
    Y = R .* sin(Phi);
    plot3(X,Y,Z,'linewidth',2);
    xlabel('X','fontsize',18);ylabel('Y','fontsize',18);zlabel('Z','fontsize',18);
    title(['particle ',num2str(tag(i)),' t = 0 ~ ',num2str((length(R) - 1) * ndiag)]);
    xlim([-6, 6]);ylim([-6, 6]);zlim([-4, 4]);
    drawnow;
    pause(0.1);
    F=getframe(gcf);
    saveas(gcf, ['p',num2str(tag(i))], 'jpg');
    end
end
figure;
for it = 1 : 100 : nstop / ndiag
%     figure; 
    set(gca,'nextplot','replacechildren');
    set(gcf,'DefaultAxesFontSize',15);
    set(gcf,'Position',get(0,'ScreenSize'));
    jp = 1;
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
            title(['phase space motion t = ' num2str((it - 1) * ndiag)]);
        end
    end
    F(jp) = getframe(gcf);
    jp = jp + 1;

end  
    writegif(['phase_movie_','t = 0 ~ ',num2str(nstop),'.gif'],F,0.1);
        