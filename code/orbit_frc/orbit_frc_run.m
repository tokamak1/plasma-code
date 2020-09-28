close all;
clc;
clear;

t = 0;
nt = 50000;
nw = 10;
ntp = nw * nt;
dl = 2e-3;
nl = 1;
dt = 5e-3;
dtp = dt/nw;
omega = 2*pi;
I0 = 1;
x1 = zeros(1,nl+1);
y1 = zeros(1,nl+1);
z1 = zeros(1,nl+1);
y1(1) = -0.95;
x1(1) = 0.05;
z1(1) = 0.05;
x2 = zeros(1,nl+1);
y2 = zeros(1,nl+1);
z2 = zeros(1,nl+1);
y2(1) = -0.95;
x2(1) = -0.05;
z2(1) = -0.05;
x3 = zeros(1,nl+1);
y3 = zeros(1,nl+1);
z3 = zeros(1,nl+1);
y3(1) = -0.95;
x3(1) = 0.05;
z3(1) = -0.05;
x4 = zeros(1,nl+1);
y4 = zeros(1,nl+1);
z4 = zeros(1,nl+1);
y4(1) = -0.95;
x4(1) = -0.05;
z4(1) = 0.05;
q = 20; m = 1;
xp = zeros(1, ntp+1);
yp = zeros(1, ntp+1);
zp = zeros(1, ntp+1);
yp(1) = 0.0;
xp(1) = 0.2;
zp(1) = 0.2;
vxp = zeros(1, ntp+1);
vyp = zeros(1, ntp+1);
vzp = zeros(1, ntp+1);
vyp(1) = 0;
vxp(1) = 0.0;
vzp(1) = 0.03;

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



