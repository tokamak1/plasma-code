Iv2 = zeros(nh+3);
Ivp2 = zeros(nh+3);

c_l1 = [c_bl1;c_bl2(1,:);c_bl3(1,:);c_center1(1,:)];
mat_v2_1 = c_l1 * y2(h1) * c_l1';
mat_vp2_1 = c_l1 * yp2(h1) * c_l1';

Iv2(1:4,1:4) = mat_v2_1;
Ivp2(1:4,1:4) = mat_vp2_1;

c_l2 = [c_bl2(2,:);c_bl3(2,:);c_center1(2,:);c_center1(1,:)];
mat_v2_2 = c_l2 * y2(h1) * c_l2';
mat_vp2_2 = c_l2 * yp2(h1) * c_l2';

Iv2(2:5,2:5) = Iv2(2:5,2:5) + mat_v2_2;
Ivp2(2:5,2:5) = Ivp2(2:5,2:5) + mat_vp2_2;

c_l3 = [c_bl3(3,:);c_center1(3,:);c_center1(2,:);c_center1(1,:)];
mat_v2_3 = c_l3 * y2(h1) * c_l3';
mat_vp2_3 = c_l3 * yp2(h1) * c_l3';

Iv2(3:6,3:6) = Iv2(3:6,3:6) + mat_v2_3;
Ivp2(3:6,3:6) = Ivp2(3:6,3:6) + mat_vp2_3;

mat_v2_center1 = flipud(c_center1) * y2(h1) * flipud(c_center1)';
mat_vp2_center1 = flipud(c_center1) * yp2(h1) * flipud(c_center1)';

for i = 4 : n1-3
    
    Iv2(i:i+3,i:i+3) = Iv2(i:i+3,i:i+3) + mat_v2_center1;
    Ivp2(i:i+3,i:i+3) = Ivp2(i:i+3,i:i+3) + mat_vp2_center1;
    
end

c_x11 = [c_center1(4,:);c_center1(3,:);c_center1(2,:);c_x1(1,:)];
mat_v2_x11 = c_x11 * y2(h1) * c_x11';
mat_vp2_x11 = c_x11 * yp2(h1) * c_x11';

Iv2(n1-2:n1+1,n1-2:n1+1) = Iv2(n1-2:n1+1,n1-2:n1+1) + mat_v2_x11;
Ivp2(n1-2:n1+1,n1-2:n1+1) = Ivp2(n1-2:n1+1,n1-2:n1+1) + mat_vp2_x11;

c_x12 = [c_center1(4,:);c_center1(3,:);c_x1(2,:);c_x2(1,:)];
mat_v2_x12 = c_x12 * y2(h1) * c_x12';
mat_vp2_x12 = c_x12 * yp2(h1) * c_x12';

Iv2(n1-1:n1+2,n1-1:n1+2) = Iv2(n1-1:n1+2,n1-1:n1+2) + mat_v2_x12;
Ivp2(n1-1:n1+2,n1-1:n1+2) = Ivp2(n1-1:n1+2,n1-1:n1+2) + mat_vp2_x12;

c_x13 = [c_center1(4,:);c_x1(3,:);c_x2(2,:);c_x3(1,:)];
mat_v2_x13 = c_x13 * y2(h1) * c_x13';
mat_vp2_x13 = c_x13 * yp2(h1) * c_x13';

Iv2(n1:n1+3,n1:n1+3) = Iv2(n1:n1+3,n1:n1+3) + mat_v2_x13;
Ivp2(n1:n1+3,n1:n1+3) = Ivp2(n1:n1+3,n1:n1+3) + mat_vp2_x13;

c_x21 = [c_x1(4,:);c_x2(3,:);c_x3(2,:);c_center2(1,:)];
mat_v2_x21 = c_x21 * y2(h2) * c_x21';
mat_vp2_x21 = c_x21 * yp2(h2) * c_x21';

Iv2(n1+1:n1+4,n1+1:n1+4) = Iv2(n1+1:n1+4,n1+1:n1+4) + mat_v2_x21;
Ivp2(n1+1:n1+4,n1+1:n1+4) = Ivp2(n1+1:n1+4,n1+1:n1+4) + mat_vp2_x21;

c_x22 = [c_x2(4,:);c_x3(3,:);c_center2(2,:);c_center2(1,:)];
mat_v2_x22 = c_x22 * y2(h2) * c_x22';
mat_vp2_x22 = c_x22 * yp2(h2) * c_x22';

Iv2(n1+2:n1+5,n1+2:n1+5) = Iv2(n1+2:n1+5,n1+2:n1+5) + mat_v2_x22;
Ivp2(n1+2:n1+5,n1+2:n1+5) = Ivp2(n1+2:n1+5,n1+2:n1+5) + mat_vp2_x22;

c_x23 = [c_x3(4,:);c_center2(3,:);c_center2(2,:);c_center2(1,:)];
mat_v2_x23 = c_x23 * y2(h2) * c_x23';
mat_vp2_x23 = c_x23 * yp2(h2) * c_x23';

Iv2(n1+3:n1+6,n1+3:n1+6) = Iv2(n1+3:n1+6,n1+3:n1+6) + mat_v2_x23;
Ivp2(n1+3:n1+6,n1+3:n1+6) = Ivp2(n1+3:n1+6,n1+3:n1+6) + mat_vp2_x23;


mat_v2_center2 = flipud(c_center2) * y2(h2) * flipud(c_center2)';
mat_vp2_center2 = flipud(c_center2) * yp2(h2) * flipud(c_center2)';

for i = n1+4 : nh-3
    
    Iv2(i:i+3,i:i+3) = Iv2(i:i+3,i:i+3) + mat_v2_center2;
    Ivp2(i:i+3,i:i+3) = Ivp2(i:i+3,i:i+3) + mat_vp2_center2;
    
end

c_r1 = [c_center2(4,:);c_center2(3,:);c_center2(2,:);c_br1(1,:)];
mat_v2_e3 = c_r1 * y2(h2) * c_r1';
mat_vp2_e3 = c_r1 * yp2(h2) * c_r1';

Iv2(nh-2:nh+1,nh-2:nh+1) = Iv2(nh-2:nh+1,nh-2:nh+1) + mat_v2_e3;
Ivp2(nh-2:nh+1,nh-2:nh+1) = Ivp2(nh-2:nh+1,nh-2:nh+1) + mat_vp2_e3;

c_r2 = [c_center2(4,:);c_center2(3,:);c_br1(2,:);c_br2(1,:)];
mat_v2_e2 = c_r2 * y2(h2) * c_r2';
mat_vp2_e2 = c_r2 * yp2(h2) * c_r2';

Iv2(nh-1:nh+2,nh-1:nh+2) = Iv2(nh-1:nh+2,nh-1:nh+2) + mat_v2_e2;
Ivp2(nh-1:nh+2,nh-1:nh+2) = Ivp2(nh-1:nh+2,nh-1:nh+2) + mat_vp2_e2;

c_r3 = [c_center2(4,:);c_br1(3,:);c_br2(2,:);c_br3];
mat_v2_e1 = c_r3 * y2(h2) * c_r3';
mat_vp2_e1 = c_r3 * yp2(h2) * c_r3';

Iv2(nh:nh+3,nh:nh+3) = Iv2(nh:nh+3,nh:nh+3) + mat_v2_e1;
Ivp2(nh:nh+3,nh:nh+3) = Ivp2(nh:nh+3,nh:nh+3) + mat_vp2_e1;

