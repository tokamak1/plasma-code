n = 100;

for i = 1 : nh + 3
    if i < 4
        bx1 = 0 : h0(2)/n : hs(2);
        bx2 = hs(2) : h0(3)/n : hs(3);
        bx3 = hs(3) : h0(4)/n : hs(4);
        by1 = C(1,i) + C(2,i)*bx1 + C(3,i)*bx1.^2 + C(4,i)*bx1.^3;
        by2 = C(5,i) + C(6,i)*bx2 + C(7,i)*bx2.^2 + C(8,i)*bx2.^3;
        by3 = C(9,i) + C(10,i)*bx3 + C(11,i)*bx3.^2 + C(12,i)*bx3.^3;
        plot(bx1,by1, bx2,by2, bx3,by3); hold on;
    elseif i < nh+1
        bx1 = hs(i-3) : h0(i-2)/n : hs(i-2);
        bx2 = hs(i-2) : h0(i-1)/n : hs(i-1);
        bx3 = hs(i-1) : h0(i)/n : hs(i);
        bx4 = hs(i) : h0(i+1)/n : hs(i+1);
        by1 = C(1,i) + C(2,i)*bx1 + C(3,i)*bx1.^2 + C(4,i)*bx1.^3;
        by2 = C(5,i) + C(6,i)*bx2 + C(7,i)*bx2.^2 + C(8,i)*bx2.^3;
        by3 = C(9,i) + C(10,i)*bx3 + C(11,i)*bx3.^2 + C(12,i)*bx3.^3;
        by4 = C(13,i) + C(14,i)*bx4 + C(15,i)*bx4.^2 + C(16,i)*bx4.^3;
        plot(bx1,by1, bx2,by2, bx3,by3, bx4,by4); hold on;
    else
        bx = hs(nh-2) : h0(nh-1)/n : hs(nh-1);
        bx1 = bx + (i-nh-1)*h0(nh-1);
        bx2 = bx + (i-nh)*h0(nh-1);
        bx3 = bx + (i-nh+1)*h0(nh-1);        
        by1 = C(1,i) + C(2,i)*bx1 + C(3,i)*bx1.^2 + C(4,i)*bx1.^3;
        by2 = C(5,i) + C(6,i)*bx2 + C(7,i)*bx2.^2 + C(8,i)*bx2.^3;
        by3 = C(9,i) + C(10,i)*bx3 + C(11,i)*bx3.^2 + C(12,i)*bx3.^3;
        plot(bx1,by1, bx2,by2, bx3,by3, bx4,by4); hold on;
    end
end

axis tight;
                
    