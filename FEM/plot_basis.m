x1 = 0:0.01*h1:h1;
x2 = 0:0.01*h2:h2;

x_h1 = phi(0,x1);
x_h2 = phi(0,x2);

y_center1 = c_center1 * x_h1;
y_center2 = c_center2 * x_h2;

for i = 1 : n1-3
    plot(x1+hs(i),y_center1(1,:),x1+hs(i+1),y_center1(2,:),...
        x1+hs(i+2),y_center1(3,:),x1+hs(i+3),y_center1(4,:)); 
    hold on;
end

for i = n1+1 : nh-3
    plot(x2+hs(i),y_center2(1,:),x2+hs(i+1),y_center2(2,:),...
        x2+hs(i+2),y_center1(3,:),x2+hs(i+3),y_center2(4,:));
end
y_bl1 = c_bl1 * x_h1;

plot(x1,y_bl1(1,:));

y_bl2 = c_bl2 * x_h1;

plot(x1,y_bl2(1,:),x1+h1,y_bl2(2,:));

y_bl3 = c_bl3 * x_h1;

plot(x1,y_bl3(1,:),x1+h1,y_bl3(2,:),x1+2*h1,y_bl3(3,:));

y_br1 = c_br1 * x_h2;

plot(x2+hs(nh-2),y_br1(1,:), ... 
    x2+hs(nh-1),y_br1(2,:),x2+hs(nh),y_br1(3,:));

y_br2 = c_br2 * x_h2;

plot(x2+hs(nh-1),y_br2(1,:),x2+hs(nh),y_br2(2,:));

y_br3 = c_br3 * x_h2;

plot(x2+hs(nh),y_br3(1,:));

y_x11 = c_x1(1,:) * phi(0,x1);
y_x12 = c_x1(2,:) * phi(0,x1);
y_x13 = c_x1(3,:) * phi(0,x1);
y_x14 = c_x1(4,:) * phi(0,x2);

plot(x1+hs(n1-2),y_x11,x1+hs(n1-1),y_x12,x1+hs(n1),y_x13,x2+hs(n1+1),y_x14);

y_x21 = c_x2(1,:) * phi(0,x1);
y_x22 = c_x2(2,:) * phi(0,x1);
y_x23 = c_x2(3,:) * phi(0,x2);
y_x24 = c_x2(4,:) * phi(0,x2);

plot(x1+hs(n1-1),y_x21,x1+hs(n1),y_x22,x2+hs(n1+1),y_x23,x2+hs(n1+2),y_x24);

y_x31 = c_x3(1,:) * phi(0,x1);
y_x32 = c_x3(2,:) * phi(0,x2);
y_x33 = c_x3(3,:) * phi(0,x2);
y_x34 = c_x3(4,:) * phi(0,x2);

plot(x1+hs(n1),y_x31,x2+hs(n1+1),y_x32,...
    x2+hs(n1+2),y_x33,x2+hs(n1+3),y_x34);
axis tight;





                
    