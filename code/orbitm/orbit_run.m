% zhao chen 2019.2.3

clc;
clear;
close all;

load_para;

nw = 511;
nh = 511;
psis_n = [0.2,0.5,0.8,1];  % initialize psi of flux surface
npsi = 4;
cal_flux;

ilim = 1;
if ilim
limiter;
end

magsurface;

cal_Bfield;

nf = 3; % choose flux surface;

% start point of fieldline
lr = RS(nsa(nf) + 1000,1);
lz = ZS(nsa(nf) + 1000,1);
phi0 = phi(1);
lx = lr * cos(phi0);
ly = lr * sin(phi0);
% RK4 parameters : dl & number of steps
nl = 1000;
dh = 0.2;

fieldline;

e = 1.602176565e-19; % Elementary charge (Coulomb)
m_pr = 1.672621777e-27; % Proton mass (kg)
m_el = 9.10938291e-31; % Electron mass (kg)
c = 299792458; % speed of light (m/s)

% particle paramaters
m = 4 * m_pr;
q = 2 * e;
K = 3.5e6 * e;
v = sqrt(2 * K / m);
sign = 1;
pitch = 0.2; % particle pitch angle

vper = sqrt(pitch * 2 * K / m); % perpendicular velocity
vpara = sign * sqrt(v^2 - vper^2); % parallel velocity

% particle initial position
phi0 = phi(1);
fr = RS(nsa(nf) + 1000,1);
fz = ZS(nsa(nf) + 1000,1);
fx = fr * cos(phi0);
fy = fr * sin(phi0);

% time parameters
b0 = B(Brxy,Btxy,Bzxy,rmaxis,zmaxis,Rxy,Zxy);
T = abs(2 * pi * m / q / b0);
dt = T / 50 ;
nt = 2000 * T / dt;

push;
