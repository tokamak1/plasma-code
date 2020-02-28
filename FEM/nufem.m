close all;
clc;
clear;

n1 = 50; % # of 1st kind of grid
n2 = 40; % # of 2nd kind of grid
w1 = 0.5; % weight of 1st kind of grid
w2 = 1 - w1; % weight of 2nd kind of grid

l = 2*pi;

% non-uniform grid
h0 = [0,ones(1,n1)*l/n1*w1,ones(1,n2)*l/n2*w2];

% total # of grid
nh = n1+n2;

% build cubic spline basis using C2 contium and boundary condition
nubspline;

plot_basis;

% build eigenfunction matrix
build_mat;
% I = Iv2\Ivp2; % u'=0 bc
I = Iv2(2:nh+2,2:nh+2)\Ivp2(2:nh+2,2:nh+2); % u=0 bc

% solve eigenvalue
[v,D] = eig(I);

N = 100;
n = 76;

% interpolation
x = 0 : 2*pi/N : 2*pi;

% V = v(:,n)'; % u'=0 bc
V = [0,v(:,n)',0]; % u=0 bc
y = u(x,V,hs,C);

figure;
plot(x,y);
axis tight;
    
    
    
    
    
    


