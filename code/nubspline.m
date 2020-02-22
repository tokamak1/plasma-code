close all;
clc;
% clear;
b = zeros(10,1);
b(end) = 1;
n1 = 10; 
n2 = 5;
w1 = 0.5;
w2 = 1 - w1;
l = 2*pi;
h0 = [0,ones(1,n1)*l/n1*w1,ones(1,n2)*l/n2*w2];  % non-uniform grid
nh = n1+n2;
l0 = 0;
for n = 3 : nh-1
    l0 = l0 + h0(n-2);
    h = h0(n-1:n+2);
    hs = zeros(1,4);
    for i = 1 : 4
        hs(i) = sum(h(1:i));
    end
    A = [hs(1)^3,-1,-hs(1),-hs(1)^2,-hs(1)^3,0,0,0,0,0;... 
        % 0th order continum at 1st grid point
        3*hs(1)^2,0,-1,-2*hs(1),-3*hs(1)^2,0,0,0,0,0;...   
        % 1st order continum at 1st grid poin
        6*hs(1),0,0,-2,-6*hs(1),0,0,0,0,0;...  
        % 2st order continum at 1st grid poin
        0,1,hs(2),hs(2)^2,hs(2)^3,-1,-hs(2),-hs(2)^2,-hs(2)^3,0;... 
        %... at 2nd grid point
        0,0,1,2*hs(2),3*hs(2)^2,0,-1,-2*hs(2),-3*hs(2)^2,0;... 
        %... at 2nd grid poin
        0,0,0,2,6*hs(2),0,0,-2,-6*hs(2),0;... 
        %... at 2nd grid poin
        0,0,0,0,0,1,hs(3),hs(3)^2,hs(3)^3,h(4)^3;... 
        %... at 3rd grid poin
        0,0,0,0,0,0,1,2*hs(3),3*hs(3)^2,-3*h(4)^2;... 
        %... at 3rd grid poin
        0,0,0,0,0,0,0,2,6*hs(3),6*h(4)]; 
        %... at 3rd grid poin
        
    if n1 > n2 
        % h1 < h2
        A = [A;0,1,hs(4)/4,(hs(4)/4)^2,(hs(4)/4)^3,...
            1,hs(4)/2,(hs(4)/2)^2,(hs(4)/2)^3,-(hs(4)/4)^3]; % sum = 1
    else
        % h1 > h2
        A = [A;(hs(4)/4)^3,1,hs(4)/2,(hs(4)/2)^2,(hs(4)/2)^3,...
            1,3*hs(4)/4,(3*hs(4)/4)^2,(3*hs(4)/4)^3,0]; % sum = 1
    end
    
    c = inv(A)*b;
    x1 = 0:0.01*h(1):hs(1);
    phi1 = c(1)*x1.^3;  
    x2 = hs(1):0.01*h(2):hs(2);
    phi2= c(2) + c(3)*x2 + c(4)*x2.^2 + c(5)*x2.^3;
    x3 = hs(2):0.01*h(3):hs(3);
    phi3= c(6) + c(7)*x3 + c(8)*x3.^2 + c(9)*x3.^3;
    x4 = hs(3):0.01*h(4):hs(4);
    phi4= c(10)*(x4-hs(4)).^3;

    plot(x1+l0,phi1,x2+l0,phi2,x3+l0,phi3,x4+l0,phi4);
    hold on;
end
h = h0(2:4);
A = [h(1)^2,h(1)^3,-1,-h(1),-h(1)^2,-h(1)^3,0;... 
    % 0th order continum at 1st grid point
    2*h(1),3*h(1)^2,0,-1,-2*h(1),-3*h(1)^2,0;...   
    % 1st order continum at 1st grid poin
    2,6*h(1),0,0,-2,-6*h(1),0;...  
    % 2st order continum at 1st grid poin
    0,0,1,h(1)+h(2),(h(1)+h(2))^2,(h(1)+h(2))^3,h(3)^3;... 
    %... at 3rd grid poin
    0,0,0,1,2*(h(1)+h(2)),3*(h(1)+h(2))^2,-3*(h(3))^2;... 
    %... at 3rd grid poin
    0,0,0,0,2,6*(h(1)+h(2)),6*h(3)]; 
    %... at 3rd grid poin
    
% A = [A;(sum(h)/4)^2,(sum(h)/4)^3,1,sum(h)/2,(sum(h)/2)^2,(sum(h)/2)^3,...
%     -(sum(h)/4)^3]; % sum = 1
A = [A;0,0,0,0,0,0,-h(3)^3]; 


b = [0,0,0,0,0,0,1/6]';

c = inv(A)*b;

x1 = 0:0.01*h(1):h(1);
phi1 = c(1)*x1.^2 + c(2)*x1.^3;  
x2 = h(1):0.01*h(2):h(1)+h(2);
phi2= c(3) + c(4)*x2 + c(5)*x2.^2 + c(6)*x2.^3;
x3 = h(1)+h(2):0.01*h(3):sum(h);
phi3= c(7)*(x3-sum(h)).^3;

plot(x1,phi1,x2,phi2,x3,phi3);
hold on;

h = h0(2:3);
A = [h(1),h(1)^2,h(1)^3,h(2)^3;... 
    % 0th order continum at 1st grid point
    1,2*h(1),3*h(1)^2,-3*h(2)^2;...   
    % 1st order continum at 1st grid poin
    0,2,6*h(1),6*h(2)];  
    % 2st order continum at 1st grid poin
    
% A = [A;h(1)/2+h(1),(h(1)/2)^2+h(1)^2,(h(1)/2)^3+h(1)^3,-(h(2)/2)^3]; % sum = 1
A = [A;0,0,0,-h(2)^3]; 

b = [0,0,0,1/4]';

c = inv(A)*b;

phi1 = c(1)*x1 + c(2)*x1.^2 + c(3)*x1.^3;  
phi2= c(4)*(x2-sum(h)).^3;

plot(x1,phi1,x2,phi2);
hold on;

plot(x1,-((x1-h(1))/h(1)).^3);

h = h0(nh-1:nh+1);
l0 = sum(h0(1:nh-2));
A = [h(1)^3,-1,-h(1),-h(1)^2,-h(1)^3,0,0;... 
    % 0th order continum at 1st grid point
    3*h(1)^2,0,-1,-2*h(1),-3*h(1)^2,0,0;...   
    % 1st order continum at 1st grid poin
    6*h(1),0,0,-2,-6*h(1),0,0;...  
    % 2st order continum at 1st grid poin
    0,1,h(1)+h(2),(h(1)+h(2))^2,(h(1)+h(2))^3,-h(3)^2,h(3)^3;... 
    %... at 3rd grid poin
    0,0,1,2*(h(1)+h(2)),3*(h(1)+h(2))^2,2*h(3),-3*h(3)^2;... 
    %... at 3rd grid poin
    0,0,0,2,6*(h(1)+h(2)),-2,6*h(3)]; 
    %... at 3rd grid poin
    
% A = [A;(sum(h)/4)^2,(sum(h)/4)^3,1,sum(h)/2,(sum(h)/2)^2,(sum(h)/2)^3,...
%     -(sum(h)/4)^3]; % sum = 1
A = [A;h(1)^3,0,0,0,0,0,0]; 


b = [0,0,0,0,0,0,1/6]';

c = inv(A)*b;

x1 = 0:0.01*h(1):h(1);
phi1 = c(1)*x1.^3;  
x2 = h(1):0.01*h(2):h(1)+h(2);
phi2= c(2) + c(3)*x2 + c(4)*x2.^2 + c(5)*x2.^3;
x3 = h(1)+h(2):0.01*h(3):sum(h);
phi3= c(6)*(x3-sum(h)).^2 + c(7)*(x3-sum(h)).^3;

plot(x1+l0,phi1,x2+l0,phi2,x3+l0,phi3);
hold on;

h = h0(nh:nh+1);
l0 = sum(h0(1:nh-1));
A = [h(1)^3,h(2),-h(2)^2,h(2)^3;... 
    % 0th order continum at 1st grid point
    3*h(1)^2,-1,2*h(2),-3*h(2)^2;...   
    % 1st order continum at 1st grid poin
    6*h(1),0,-2,6*h(2)];  
    % 2st order continum at 1st grid poin
    
% A = [A;h(1)/2+h(1),(h(1)/2)^2+h(1)^2,(h(1)/2)^3+h(1)^3,-(h(2)/2)^3]; % sum = 1
A = [A;h(1)^3,0,0,0]; 

b = [0,0,0,1/4]';

c = inv(A)*b;

phi1 =c(1)*x1.^3;  
phi2= c(2)*(x2-sum(h)) + c(3)*(x2-sum(h)).^2 + c(4)*(x2-sum(h)).^3;

plot(x1+l0,phi1,x2+l0,phi2);
hold on;

l0 = sum(h0(1:nh));
plot(x1+l0,(x1/h(1)).^3);
