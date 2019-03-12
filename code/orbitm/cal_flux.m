x=range_efit(4):range_efit(1)/(nw-1):range_efit(4)+range_efit(1);
y=-range_efit(2)/2+range_efit(5):range_efit(2)/(nh-1):range_efit(2)/2+range_efit(5);
flux_efit=imresize(flux_efit0,[nw nh],'bicubic');


rbbbs_d=round((rbbbs-range_efit(4))*(nw-1)/range_efit(1))+1;
zbbbs_d=round((zbbbs+(range_efit(2)-range_efit(5))/2)*(nh-1)/range_efit(2))+1;

maxis_r=round((rmaxis-range_efit(4))*(nh-1)/range_efit(1))+1;
maxis_z=round((zmaxis+(range_efit(2)-range_efit(5))/2)*(nh-1)/range_efit(2))+1;

for iiii=round(length(zbbbs)/2)-10:length(zbbbs)
    if (zbbbs(iiii)-zmaxis)*(zmaxis-zbbbs(iiii-1))>=0
        break
    end
end
iiii=iiii-1;
edg_r=round((rbbbs(iiii)-range_efit(4))*(nh-1)/range_efit(1))+1;
edg_z=round((zbbbs(iiii)+(range_efit(2)-range_efit(5))/2)*(nh-1)/range_efit(2))+1;

psirz = flux_efit';



[xx,yy]=meshgrid(x,y);
XS = []; YS = []; ZS = []; RS = [];
nst = 0; 
nsa = []; nsa(1) = 1;
ns = 0;

 for k = 1 : npsi
     psis = simag + psis_n(k) * abs(sibry - simag);
     [cc,hh] = contour(xx,yy,flux_efit',[psis,psis]);
     close;
     n0 = length(cc(1,:));
     n1 = cc(2,1);
     sr1 = cc(1,2 : 1 + n1);
     sz1 = cc(2,2 : 1 + n1);
     if abs(max(sz1) + min(sz1)) < 0.0001
         ns = n1;
         sr = sr1;
         sz = sz1;
     end
     if abs(sum(sz1)/length(sz1)) < 1
         ns = n1;
         sr = sr1;
         sz = sz1;
     end
     if n1 < n0 - 1
         n2 = cc(2,2 + n1);
         sr2 = cc(1,3 + n1 : 2 + n1 + n2);
         sz2 = cc(2,3 + n1 : 2 + n1 + n2);
         if abs(max(sz2) + min(sz2)) < 0.0001
             ns = n2;
             sr = sr2;
             sz = sz2;
         end
         if abs(sum(sz2)/length(sz2)) < 1
             ns = n2;
             sr = sr2;
             sz = sz2;
         end
          if n1 + n2 < n0 - 2
             n3 = cc(2,3 + n1 + n2);
             sr3 = cc(1,4 + n1 + n2 : 3 + n1 + n2 + n3);
             sz3 = cc(2,4 + n1 + n2 : 3 + n1 + n2 + n3);
             if abs(max(sz3) + min(sz3)) < 0.0001
                 ns = n3;
                 sr = sr3;
                 sz = sz3;
             end
             if ns == 0
                 ns = n3;
                 sr = sr3;
                 sz = sz3;
             end
              if abs(sum(sz3)/length(sz3)) < 1
                 ns = n3;
                 sr = sr3;
                 sz = sz3;
              end
          end
     end

    
    nt = 50;
    
    phi = (- pi / 3 : 2 * pi / (nt - 1) : 5 * pi / 3);
    for i = 1 : ns
        for j = 1 : nt
            XS(i + nst,j) = sr(i) * cos(phi(j));
            ZS(i + nst,j) = sz(i);
            YS(i + nst,j) = sr(i) * sin(phi(j));
            RS(i + nst,j) = sr(i);
        end
    end
    nst = nst + ns;
    nsa(k + 1) = nst + 1;
 end
 
 nlimi = length(rlim);
 
 xlimi = []; ylimi = []; zlimi = []; rlimi = [];
 for i = 1 : nlimi
    for j = 1 : nt 
        xlimi(i, j) = rlim(i) * cos(phi(j));
        zlimi(i, j) = zlimter(i);
        ylimi(i, j) = rlim(i) * sin(phi(j));
        rlimi(i, j) = rlim(i);
    end
end

RZ = [xx;yy];
PSIrz = flux_efit';
save('./output/RZ.dat','RZ','-ascii','-double','-tabs');
save('./output/PSIrz.dat','PSIrz','-ascii','-double','-tabs');