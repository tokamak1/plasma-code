Rxy=xx;
Zxy=yy;
  
fpol_tmp=imresize(fpol,[nw 1],'bicubic');
    
flux_profile1 = psirz(maxis_z,maxis_r:edg_r);
flux_profile=imresize(flux_profile1,[1 nw],'bicubic');
p = polyfit(flux_profile',fpol_tmp,4);

for i=1:nw
    for j=1:nh
        fpol_rz(i,j)=polyval(p,psirz(i,j));
    end
end

for i=1:nw
    for j=1:nh
        if fpol_rz(i,j)<min(fpol)/100
            fpol_rz(i,j)=min(fpol)/100;
        end
    end
end

Btxy=fpol_rz./Rxy;

dr0=range_efit(1)/(nw-1);
dz0=range_efit(2)/(nh-1);

for i=1:nw
    for j=1:nh
        
         if j==1
           dpsir(i,j)=psirz(i,j+1)-psirz(i,j);
           Bzxy(i,j)=1/Rxy(i,j)*dpsir(i,j)/dr0;
         elseif j==nh
           dpsir(i,j)=psirz(i,j)-psirz(i,j-1);
           Bzxy(i,j)=1/Rxy(i,j)*dpsir(i,j)/dr0;
         else
           dpsir(i,j)=psirz(i,j+1)-psirz(i,j-1);
           Bzxy(i,j)=1/Rxy(i,j)*dpsir(i,j)/2/dr0;
         end

         if i==1
            dpsiz(i,j)=psirz(i+1,j)-psirz(i,j);
            Brxy(i,j)=-1/Rxy(i,j)*dpsiz(i,j)/dz0;
         elseif i==nw
            dpsiz(i,j)=psirz(i,j)-psirz(i-1,j);
            Brxy(i,j)=-1/Rxy(i,j)*dpsiz(i,j)/dz0;
         else
            dpsiz(i,j)=psirz(i+1,j)-psirz(i-1,j);
            Brxy(i,j)=-1/Rxy(i,j)*dpsiz(i,j)/2/dz0;
         end
     
    end
end

Brz = [Brxy;Btxy;Bzxy];
save('./output/Brz.dat','Brz','-ascii','-double','-tabs');