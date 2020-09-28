ng = 10;

x11 = 0.55; y11 = -1.05; z11 = 0.5;
x12 = x11; y12 = y11; z12 = -0.5;
x13 = x11; y13 = -0.05; z13 = z12;
x14 = x11; y14 = y13; z14 = z11;
dz11 = (z12 - z11)/ng;
dy12 = (y13 - y12)/ng;
dz13 = (z14 - z13)/ng;
dy14 = (y11 - y14)/ng;

coil1x = x11*ones(1,4*(ng+1));
coil1y = [y11*ones(1,ng+1), (y12:dy12:y13), y13*ones(1,ng+1), (y14:dy14:y11)];
coil1z = [(z11:dz11:z12), z12*ones(1,ng+1), (z13:dz13:z14), z11*ones(1,ng+1)];


x21 = 0.55; y21 = 0.05; z21 = 0.5;
x22 = x11; y22 = y21; z22 = -0.5;
x23 = x11; y23 = 1.05; z23 = z22;
x24 = x11; y24 = y23; z24 = z21;
dz21 = (z22 - z21)/ng;
dy22 = (y23 - y22)/ng;
dz23 = (z24 - z23)/ng;
dy24 = (y21 - y24)/ng;

coil2x = x21*ones(1,4*(ng+1));
coil2y = [y21*ones(1,ng+1), (y22:dy22:y23), y23*ones(1,ng+1), (y24:dy24:y21)];
coil2z = [(z21:dz21:z22), z22*ones(1,ng+1), (z23:dz23:z24), z21*ones(1,ng+1)];


x31 = -0.55; y31 = 0.05; z31 = 0.5;
x32 = x31; y32 = y31; z32 = -0.5;
x33 = x31; y33 = 1.05; z33 = z32;
x34 = x31; y34 = y33; z34 = z31;
dz31 = (z32 - z31)/ng;
dy32 = (y33 - y32)/ng;
dz33 = (z34 - z33)/ng;
dy34 = (y31 - y34)/ng;

coil3x = x31*ones(1,4*(ng+1));
coil3y = [y31*ones(1,ng+1), (y32:dy32:y33), y33*ones(1,ng+1), (y34:dy34:y31)];
coil3z = [(z31:dz31:z32), z32*ones(1,ng+1), (z33:dz33:z34), z31*ones(1,ng+1)];


x41 = -0.55; y41 = -1.05; z41 = 0.5;
x42 = x41; y42 = y41; z42 = -0.5;
x43 = x41; y43 = -0.05; z43 = z42;
x44 = x41; y44 = y43; z44 = z41;
dz41 = (z42 - z41)/ng;
dy42 = (y43 - y42)/ng;
dz43 = (z44 - z43)/ng;
dy44 = (y41 - y44)/ng;

coil4x = x41*ones(1,4*(ng+1));
coil4y = [y41*ones(1,ng+1), (y42:dy42:y43), y43*ones(1,ng+1), (y44:dy44:y41)];
coil4z = [(z41:dz41:z42), z42*ones(1,ng+1), (z43:dz43:z44), z41*ones(1,ng+1)];


x51 = 0.5; y51 = -1.05; z51 = 0.55;
x52 = x51; y52 = -0.05; z52 = z51;
x53 = -0.5; y53 = y52; z53 = z51;
x54 = x53; y54 = y51; z54 = z51;
dy51 = (y52 - y51)/ng;
dx52 = (x53 - x52)/ng;
dy53 = (y54 - y53)/ng;
dx54 = (x51 - x54)/ng;

coil5z = z51*ones(1,4*(ng+1));
coil5x = [x51*ones(1,ng+1), (x52:dx52:x53), x53*ones(1,ng+1), (x54:dx54:x51)];
coil5y = [(y51:dy51:y52), y52*ones(1,ng+1), (y53:dy53:y54), y51*ones(1,ng+1)];


x61 = 0.5; y61 = 0.05; z61 = 0.55;
x62 = x61; y62 = 1.05; z62 = z61;
x63 = -0.5; y63 = y62; z63 = z61;
x64 = x63; y64 = y61; z64 = z61;
dy61 = (y62 - y61)/ng;
dx62 = (x63 - x62)/ng;
dy63 = (y64 - y63)/ng;
dx64 = (x61 - x64)/ng;

coil6z = z61*ones(1,4*(ng+1));
coil6x = [x61*ones(1,ng+1), (x62:dx62:x63), x63*ones(1,ng+1), (x64:dx64:x61)];
coil6y = [(y61:dy61:y62), y62*ones(1,ng+1), (y63:dy63:y64), y61*ones(1,ng+1)];



x71 = 0.5; y71 = 0.05; z71 = -0.55;
x72 = x71; y72 = 1.05; z72 = z71;
x73 = -0.5; y73 = y72; z73 = z71;
x74 = x73; y74 = y71; z74 = z71;
dy71 = (y72 - y71)/ng;
dx72 = (x73 - x72)/ng;
dy73 = (y74 - y73)/ng;
dx74 = (x71 - x74)/ng;

coil7z = z71*ones(1,4*(ng+1));
coil7x = [x71*ones(1,ng+1), (x72:dx72:x73), x73*ones(1,ng+1), (x74:dx74:x71)];
coil7y = [(y71:dy71:y72), y72*ones(1,ng+1), (y73:dy73:y74), y71*ones(1,ng+1)];


x81 = 0.5; y81 = -1.05; z81 = -0.55;
x82 = x81; y82 = -0.05; z82 = z81;
x83 = -0.5; y83 = y82; z83 = z81;
x84 = x83; y84 = y81; z84 = z81;
dy81 = (y82 - y81)/ng;
dx82 = (x83 - x82)/ng;
dy83 = (y84 - y83)/ng;
dx84 = (x81 - x84)/ng;

coil8z = z81*ones(1,4*(ng+1));
coil8x = [x81*ones(1,ng+1), (x82:dx82:x83), x83*ones(1,ng+1), (x84:dx84:x81)];
coil8y = [(y81:dy81:y82), y82*ones(1,ng+1), (y83:dy83:y84), y81*ones(1,ng+1)];


