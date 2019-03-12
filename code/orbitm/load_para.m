pass_gfile='./gfile/gfile-cfetr5.7-baseline-129';
fid=fopen(pass_gfile);
temp=fread(fid,40);
temp=fscanf(fid,'%g',3);
nw=temp(2);
nh=temp(3);

temp=fscanf(fid,'%g',4*5+nw*4+nw*nh+nw+2);
range_efit=temp(1:5);
rdim=temp(1);
zdim=temp(2);
rcentr=temp(3);
rleft=temp(4);
zmid=temp(5);
rmaxis=temp(6);
zmaxis=temp(7);
simag=temp(8);
sibry=temp(9);
bcentr=temp(10);
current=temp(11);
simag=temp(12);
xdum=temp(13);
rmaxis=temp(14);
xdum=temp(15);
zmaxis=temp(16);
xdum=temp(17);
sibry=temp(18);
xdum=temp(19);
xdum=temp(20);
fpol=temp(21:21+nw-1);
psirz1=temp(21+4*nw:20+4*nw+nw*nh);
for i=1:nw
    flux_efit0(:,i)=psirz1((i-1)*nh+1:i*nh);
end
psirz=flux_efit0;

nbbbs=temp(20+4*nw+nw*nh+nw+1);
limitr=temp(20+4*nw+nw*nh+nw+2);

temp=fscanf(fid,'%g %g',[2,nbbbs]);
rbbbs=temp(1,:);
zbbbs=temp(2,:);
temp=fscanf(fid,'%g %g',[2,limitr]);
rlim=temp(1,:);
zlimter=temp(2,:);


