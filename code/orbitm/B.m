function y = B(bx,by,bz,x0,z0,R,Z)

yx = interp2(R,Z,bx,x0,z0);
yy = interp2(R,Z,by,x0,z0);
yz = interp2(R,Z,bz,x0,z0);

y = sqrt(yx^2 + yy^2 + yz^2);
 
end

