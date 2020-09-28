plot3(x1,y1,z1,'-r',x2,y2,z2,'-r',x3,y3,z3,'-r',x4,y4,z4,'-r',...
    coil1x,coil1y,coil1z,'--k',coil2x,coil2y,coil2z,'--k',...
    coil3x,coil3y,coil3z,'--k',coil4x,coil4y,coil4z,'--k',...
    coil5x,coil5y,coil5z,'--k',coil6x,coil6y,coil6z,'--k',...
    coil7x,coil7y,coil7z,'--k',coil8x,coil8y,coil8z,'--k');
xlabel('X','fontsize',18);
ylabel('Y','fontsize',18);
zlabel('Z','fontsize',18);
axis equal;
%     xlim([-1.5 1.5]);
%     ylim([-1.5 1.5]);
%     zlim([-1.5 1.5]);
drawnow;
% pause(0.5);
% frame = getframe(fig);
% im{i} = frame2im(frame);
% [A,map] = rgb2ind(im{i},256);
% if i==1
%   imwrite(A,map,[name0, '.gif'],'gif','LoopCount',Inf,'DelayTime',0.1);
% else
%   imwrite(A,map,[name0, '.gif'],'gif','WriteMode','append','DelayTime',0.1);      
% end