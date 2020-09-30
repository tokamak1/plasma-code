close all;
clc;
clear;

t = 0;
nt = 10000;
nw = 10;
ntp = nw * nt;
dl = 2e-3;
nl = 1;
omega = 5e5; % frequency of rotating field

I0 = 1e6;
B0 = 1;
init_fl;
init_p; % initial particle parameter 
omega_c = q*B0/m; % cyclotron frequency
dt = 2*pi/omega/50;
dtp = 2*pi/omega_c/50;

Coil;
name = 'FRC_B';
name0 = [name, '_t=', num2str(nt*dt)]; 
fig = figure;

for i = 1 : nt
    
%     rotate_B;
    
%     show_fig;
    
    push_frc;


    t = t + dt;    
end
rotate_B;
show_fig;
hold on;
plot3(xp,yp,zp);



