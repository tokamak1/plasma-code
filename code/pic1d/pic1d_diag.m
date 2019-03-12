clear;
close all;
clc;
data = [0.1   1.0152   -4.75613E-15
0.2   1.06398   -5.51074E-05
0.3   1.15985   -0.0126204
0.4   1.28506   -0.066128
0.5   1.41566   -0.153359
0.6   1.54571   -0.26411
0.7   1.67387   -0.392401
0.8   1.7999   -0.534552
0.9   1.92387   -0.688109
1.0   2.0459   -0.85133
1.5   2.63233   -1.77571
2   3.18914   -2.8272];
id=2;
k=data(id,1);
wr=data(id,2); wi=data(id,3);
picdat = './pic1d.dat';
fid = fopen(picdat,'rb');
nt = fread(fid,1,'int');
np = fread(fid,1,'int');
ng = fread(fid,1,'int');
q = fread(fid,1,'float');
dx = fread(fid,1,'float');
dt = fread(fid,1,'float');
for i = 1 : nt + 1
    xp = fread(fid,np,'float');
    vp = fread(fid,np,'float');
    Eg = fread(fid,ng,'float');
    Ek(i) = 0.5 * abs(q) * sum(vp.^2); % kinetic energy  
    Ef(i) = 0.5 * sum(Eg .* Eg) * dx; % potential energy
    t(i) = (i - 1) * dt;
    if mod((i - 1), 50) == 0
        plot(xp,vp,'b.','Markersize',2);
        xlim([0,ng * dx]);ylim([-6, 6]);
        title(['Phase space plotting, vp-x, t=',num2str(t(i))]);
        xlabel('xp');ylabel('vp');
        drawnow;
    end
end
%%
h = figure('Unit','Normalized','position',...
    [0.02 0.4 0.6 0.3],'DefaultAxesFontSize',15);
subplot(121);plot(t,Ek,t,Ef,t,Ek+Ef,'r:','LineWidth',2);
title(['(a) k=',num2str(k),', \omega_{theory}=',...
    num2str(wr+1i*wi)],'fontsize',15);
xlabel('t'); ylabel('Energy');legend('E_k','E_e','E_{tot}');
legend('boxoff');

subplot(122);

% Find the corresponding indexes of the extreme max values
lndE=log(sqrt(real((Ef)))); % EEf=E^2=[exp(gam*t)]^2=exp(2*gam*t)
it0=floor(nt*1/20); it1=floor(nt*17/20);
yy=lndE(it0:it1);
extrMaxIndex = find(diff(sign(diff(yy)))==-2)+1;
t1=t(it0+extrMaxIndex(1));t2=t(it0+extrMaxIndex(end));
y1=yy(extrMaxIndex(1));y2=yy(extrMaxIndex(end));
plot(t,lndE,[t1,t2],[y1,y2],'r*--','LineWidth',2);
omega=pi/((t2-t1)/(length(extrMaxIndex)-1));
gammas=(real(y2)-real(y1))/(t2-t1);
plot(t,lndE);
title(['(b) \omega^S=',num2str(omega),', \gamma^S=',num2str(gammas)]);
axis tight;