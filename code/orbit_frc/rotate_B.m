

Bz = zeros(1,nl);
Bx = zeros(1,nl);
By =0.0*ones(1,nl);



for j = 1 : nl
    Bx1 = Bx(j);
    By1 = By(j);
    Bz1 = Bz(j);
    Bx2 = Bx(j);
    By2 = By(j);
    Bz2 = Bz(j);
    Bx3 = Bx(j);
    By3 = By(j);
    Bz3 = Bz(j);
    Bx4 = Bx(j);
    By4 = By(j);
    Bz4 = Bz(j);
    for k = 1 : 8
        [bx1, by1, bz1] = Bcoil(k, x1(j), y1(j), z1(j), I0, omega, t);
        Bx1 = Bx1 + bx1;
        By1 = By1 + by1;
        Bz1 = Bz1 + bz1;
        [bx2, by2, bz2] = Bcoil(k, x2(j), y2(j), z2(j), I0, omega, t);
        Bx2 = Bx2 + bx2;
        By2 = By2 + by2;
        Bz2 = Bz2 + bz2;
        [bx3, by3, bz3] = Bcoil(k, x3(j), y3(j), z3(j), I0, omega, t);
        Bx3 = Bx3 + bx3;
        By3 = By3 + by3;
        Bz3 = Bz3 + bz3;
        [bx4, by4, bz4] = Bcoil(k, x4(j), y4(j), z4(j), I0, omega, t);
        Bx4 = Bx4 + bx4;
        By4 = By4 + by4;
        Bz4 = Bz4 + bz4;
    end

    B1 = sqrt(Bx1^2 + By1^2 + Bz1^2);
    x1(j+1) = x1(j) + Bx1/B1*dl;
    y1(j+1) = y1(j) + By1/B1*dl;
    z1(j+1) = z1(j) + Bz1/B1*dl;
    B2 = sqrt(Bx2^2 + By2^2 + Bz2^2);
    x2(j+1) = x2(j) + Bx2/B2*dl;
    y2(j+1) = y2(j) + By2/B2*dl;
    z2(j+1) = z2(j) + Bz2/B2*dl;
    B3 = sqrt(Bx3^2 + By3^2 + Bz3^2);
    x3(j+1) = x3(j) + Bx3/B3*dl;
    y3(j+1) = y3(j) + By3/B3*dl;
    z3(j+1) = z3(j) + Bz3/B3*dl;
    B4 = sqrt(Bx4^2 + By4^2 + Bz4^2);
    x4(j+1) = x4(j) + Bx4/B4*dl;
    y4(j+1) = y4(j) + By4/B4*dl;
    z4(j+1) = z4(j) + Bz4/B4*dl;
end   




    
    
    
    
    
    
    
