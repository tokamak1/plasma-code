coil = 0;
if coil
position_coil = [
    0.963	1.908	0.144	0.215	0 4 6 24 738
    2.1065	1.335	0.144	0.215	0 4 6 24 738
    2.1065	0.445	0.144	0.215	0 4 6 24 738
    0.963	-1.908	0.144	0.215	0 4 6 24 738
    2.1065	-1.335	0.144	0.215	0 4 6 24 738
    2.1065	-0.445	0.144	0.215	0 4 6 24 738
    ];
end

rlim(length(rlim)+1)=rlim(1);
zlimter(length(zlimter)+1)=zlimter(1);
RZ_limiter = [rlim;zlimter];
save('./output/RZ_limiter.dat','RZ_limiter','-ascii','-double','-tabs');