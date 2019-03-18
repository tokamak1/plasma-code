clear;
close all;
clc;
partdat = './particle_orbit.dat';
fid = fopen(partdat,'rb');
nstop = 400000;
ndiag = 40;
ntrack = 40000;
itag = 1;
for i = 1 : (nstop / ndiag + 1) * ntrack
    idata1 = fread(fid,8,'double');
    idata2 = fread(fid,2,'int');    
    if sum(idata2) == 0
        break;
    end
    if itag
        inp = idata2(2);
        tag(inp) = idata2(1);
        eval(['pdata',num2str(idata2(1)),'= idata1;']);
        if inp == 1 && i == 1
           tag0 = idata2(1);
           pdata0 = idata1;
        end
            
        if inp == 1 && i > 1
            itag = 0;
            tag(1) = tag0;
            eval(['pdata',num2str(idata2(1)),'= [pdata0 , idata1];']);
        end
    else
        eval(['pdata',num2str(idata2(1)),'= [pdata',num2str(idata2(1)),...
             ' , idata1];']);
    end      
end
save;