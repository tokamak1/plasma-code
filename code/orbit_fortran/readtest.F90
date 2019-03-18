        program readfile
            implicit none
            integer fid
            real, dimension(:, :) :: Br(511, 511), Bt(511, 511), &
                                     Bz(511,511)
         
            fid = 233
            open(fid, file = 'Brz.dat')
            read(fid, *)Br
            read(fid, *)Bt
            read(fid, *)Bz
            close(fid)
            open(fid + 1, file = 'B.dat', status = 'replace')
            write(fid + 1, *)Br
            write(fid + 1, *)Bt
            write(fid + 1, *)Bz
            close(fid + 1)
        end program readfile    
