close all;
clc;
clear;

n1 = 50; % # of 1st kind of grid
n2 = 50; % # of 2nd kind of grid
w1 = 0.5; % weight of 1st kind of grid
w2 = 1 - w1; % weight of 2nd kind of grid

l = 2*pi;
h1 = w1*l/n1;
h2 = w2*l/n2;

% non-uniform grid
h0 = [ones(1,n1)*h1,ones(1,n2)*h2];

% total # of grid
nh = n1+n2;

% x coordinate of grid
hs = zeros(1,nh+1);
for i = 2 : nh+1
    hs(i) = sum(h0(1:i-1));
end
hs(end) = l;

% build cubic spline basis using C2 contium and boundary condition
% nubspline;

local_bspline;

plot_basis;

% build eigenfunction matrix
build_mat;
I = Iv2\Ivp2; % u'=0 bc
% I = Iv2(2:nh+2,2:nh+2)\Ivp2(2:nh+2,2:nh+2); % u=0 bc

% solve eigenvalue
[V,D] = eig(I);

N = 1000;
n = 90;

% interpolation
x0 = 0 : 2*pi/N : 2*pi;
y = zeros(1,N+1);
y(1) = u(0,V(1:4,n)',c_l1);
% V = [0,v(:,n)',0]; % u=0 bc

for i = 2 : N+1
    
    ix = sum(x0(i)>hs);
    v = V(ix:ix+3, n)';
    x = x0(i) - hs(ix);
    
    if ix == 1
        
        y(i) = u(x,v,c_l1);
        
    elseif ix == 2
        
        y(i) = u(x,v,c_l2);
        
    elseif ix == 3

        y(i) = u(x,v,c_l3);
        
    elseif ix < n1-2
        
        y(i) = u(x,v,flipud(c_center1));
        
    elseif ix == n1-2
        
        y(i) = u(x,v,c_x11);
        
    elseif ix == n1-1
        
        y(i) = u(x,v,c_x12);
        
    elseif ix == n1
        
        y(i) = u(x,v,c_x13);
        
    elseif ix == n1+1
        
        y(i) = u(x,v,c_x21);
        
    elseif ix == n1+2
        
        y(i) = u(x,v,c_x22);
        
    elseif ix == n1+3
        
        y(i) = u(x,v,c_x23);
        
    elseif ix < nh-2
        
        y(i) = u(x,v,flipud(c_center2));
        
    elseif ix == nh-2
        
        y(i) = u(x,v,c_r1);
        
    elseif ix == nh-1
        
        y(i) = u(x,v,c_r2);
        
    elseif ix == nh
        
        y(i) = u(x,v,c_r3);
        
    end
    
end


% figure;
plot(x0,y);
axis tight;

figure;
plot(sqrt(abs(real(diag(D)))),imag(diag(D)),'o');
    
    
    
    
    
    


