for i = 1 : nh+1
    hs(i) = sum(h0(1:i));
end

C = [];
b = [zeros(15,1); 1];
for n = 3 : nh-1
    

    A = [phi(0,hs(n-2)), phi(-1,0), phi(-1,0), phi(-1,0);...
        phi(1,hs(n-2)), phi(-1,0), phi(-1,0), phi(-1,0);...
        phi(2,hs(n-2)), phi(-1,0), phi(-1,0), phi(-1,0);
        phi(0,hs(n-1)), -phi(0,hs(n-1)), phi(-1,0), phi(-1,0);... 
        phi(1,hs(n-1)), -phi(1,hs(n-1)), phi(-1,0), phi(-1,0);...   
        phi(2,hs(n-1)), -phi(2,hs(n-1)), phi(-1,0), phi(-1,0);...  
        phi(-1,0), phi(0,hs(n)), -phi(0,hs(n)), phi(-1,0);... 
        phi(-1,0), phi(1,hs(n)), -phi(1,hs(n)), phi(-1,0);... 
        phi(-1,0), phi(2,hs(n)), -phi(2,hs(n)), phi(-1,0);... 
        phi(-1,0), phi(-1,0), phi(0,hs(n+1)), -phi(0,hs(n+1));... 
        phi(-1,0), phi(-1,0), phi(1,hs(n+1)), -phi(1,hs(n+1));... 
        phi(-1,0), phi(-1,0), phi(2,hs(n+1)), -phi(2,hs(n+1));...
        phi(-1,0), phi(-1,0), phi(-1,0), phi(0,hs(n+2));...
        phi(-1,0), phi(-1,0), phi(-1,0), phi(1,hs(n+2));...
        phi(-1,0), phi(-1,0), phi(-1,0), phi(2,hs(n+2));...
        zeros(1,16)]; 
    

    xp1 = hs(n-2) + sum(h0(n-1:n+2))/4;
    xp2 = hs(n-2) + sum(h0(n-1:n+2))/2;
    xp3 = hs(n-2) + sum(h0(n-1:n+2))*3/4;
    p1 = sum(xp1>hs(n-1:n+2));
    p2 = sum(xp2>hs(n-1:n+2));
    p3 = sum(xp3>hs(n-1:n+2));
    
    A(16, p1*4+1) = A(end, p1*4+1) + 1;
    A(16, p1*4+2) = A(end, p1*4+2) + xp1;
    A(16, p1*4+3) = A(end, p1*4+3) + xp1^2;
    A(16, p1*4+4) = A(end, p1*4+4) + xp1^3;
    A(16, p2*4+1) = A(end, p2*4+1) + 1;
    A(16, p2*4+2) = A(end, p2*4+2) + xp2;
    A(16, p2*4+3) = A(end, p2*4+3) + xp2^2;
    A(16, p2*4+4) = A(end, p2*4+4) + xp2^3;
    A(16, p3*4+1) = A(end, p3*4+1) + 1;
    A(16, p3*4+2) = A(end, p3*4+2) + xp3;
    A(16, p3*4+3) = A(end, p3*4+3) + xp3^2;
    A(16, p3*4+4) = A(end, p3*4+4) + xp3^3;
    
    c = inv(A)*b;
    C = [C, c];
    
end

A = [phi(0,0), phi(-1,0), phi(-1,0);... 
    phi(1,0), phi(-1,0), phi(-1,0);...   
    phi(0,hs(2)), -phi(0,hs(2)), phi(-1,0);...  
    phi(1,hs(2)), -phi(1,hs(2)), phi(-1,0);... 
    phi(2,hs(2)), -phi(2,hs(2)), phi(-1,0);...
    phi(-1,0), phi(0,hs(3)), -phi(0,hs(3));...  
    phi(-1,0), phi(1,hs(3)), -phi(1,hs(3));... 
    phi(-1,0), phi(2,hs(3)), -phi(2,hs(3));...
    phi(-1,0), phi(-1,0), phi(0,hs(4));...  
    phi(-1,0), phi(-1,0), phi(1,hs(4));... 
    phi(-1,0), phi(-1,0), phi(2,hs(4));...
    phi(-1,0), phi(-1,0), phi(0,hs(3))]; 
     


b = [zeros(1,11), 1/6]';

c = inv(A)*b;
C = [[c;phi(-1,0)'], C];

A = [phi(0,0), phi(-1,0);...   
    phi(0,hs(2)), -phi(0,hs(2));...  
    phi(1,hs(2)), -phi(1,hs(2));... 
    phi(2,hs(2)), -phi(2,hs(2));...
    phi(-1,0), phi(0,hs(3));...  
    phi(-1,0), phi(1,hs(3));... 
    phi(-1,0), phi(2,hs(3));...
    phi(-1,0), phi(0,hs(2))];

b = [zeros(1,7), 1/4]';

c = inv(A)*b;
C = [[c;phi(-1,0)';phi(-1,0)'], C];

% C = [[1;-3/h0(2);3/h0(2)^2;-1/h0(2)^3;...
%     phi(-1,0)';phi(-1,0)';phi(-1,0)'], C];

C = [[phi(-1,0)';phi(-1,0)';phi(-1,0)';phi(-1,0)'], C];

A = [phi(0,hs(nh-2)), phi(-1,0), phi(-1,0);... 
    phi(1,hs(nh-2)), phi(-1,0), phi(-1,0);... 
    phi(2,hs(nh-2)), phi(-1,0), phi(-1,0);...   
    phi(0,hs(nh-1)), -phi(0,hs(nh-1)), phi(-1,0);...  
    phi(1,hs(nh-1)), -phi(1,hs(nh-1)), phi(-1,0);... 
    phi(2,hs(nh-1)), -phi(2,hs(nh-1)), phi(-1,0);...
    phi(-1,0), phi(0,hs(nh)), -phi(0,hs(nh));...  
    phi(-1,0), phi(1,hs(nh)), -phi(1,hs(nh));... 
    phi(-1,0), phi(2,hs(nh)), -phi(2,hs(nh));...
    phi(-1,0), phi(-1,0), phi(0,hs(nh+1));...  
    phi(-1,0), phi(-1,0), phi(1,hs(nh+1));...
    phi(-1,0), phi(0,hs(nh-1)), phi(-1,0)]; 
     


b = [zeros(1,11),1/6]';

c = inv(A)*b;
C = [C,[c;phi(-1,0)']];

A = [phi(0,hs(nh-1)), phi(-1,0);...  
    phi(1,hs(nh-1)), phi(-1,0);... 
    phi(2,hs(nh-1)), phi(-1,0);...   
    phi(0,hs(nh)), -phi(0,hs(nh));...  
    phi(1,hs(nh)), -phi(1,hs(nh));... 
    phi(2,hs(nh)), -phi(2,hs(nh));...
    phi(-1,0), phi(0,hs(nh+1));...
    phi(-1,0), phi(0,hs(nh))];

b = [zeros(1,7),1/4]';

c = inv(A)*b;
C = [C,[c;phi(-1,0)';phi(-1,0)']];

% C = [C,[...
%     -hs(nh)^3/h0(nh+1)^3;3*hs(nh)^2/h0(nh+1)^3;-3*hs(nh)/h0(nh+1)^3;1/h0(nh+1)^3;...
%     phi(-1,0)';phi(-1,0)';phi(-1,0)']];


C = [C,[phi(-1,0)';phi(-1,0)';phi(-1,0)';phi(-1,0)']];

    
