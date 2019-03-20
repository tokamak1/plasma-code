!************************************************************************************************
!                                                                                               *
!                                diagnosis #0 particle modified by C. Zhao                      *
!                                                                                               *
!************************************************************************************************
    subroutine diagnosis_orbit_zhao
        use var
        use distribution 
        implicit none
        include 'mpif.h'
        
        integer :: STATUS(MPI_STATUS_SIZE), itrack, pfid, nl
        integer, parameter :: ntrack = 40000, ndiag = 40 
        integer, dimension(ntrack) :: trackid
        integer :: np0, np
        real*8,dimension(8) :: data_diag
!        real*8,allocatable :: data_diag_total(:,:)
        character(len = 64) pf, pf0
!        logical :: flag_diag
!        logical,allocatable :: flag_diag_total(:)
!        allocate(data_diag_total(8, ntrack))
!        allocate(flag_diag_total(nsize))
        if (nstep < 2) then
            do i = 1, ntrack
               trackid(i) = nsize * N_initia / ntrack * i - 10
            enddo
            if (nrank == 0 .and. nstep == 0) then
                write(*, *)233333
                open(433, file = './pdata/tag.dat',status = 'replace',form = 'binary')
                write(*, *)233333
                write(433)trackid
                write(*, *)233333
                write(*, *)size(trackid)
                close(433)
                write(*, *)233333
            endif
        endif

        if(nstep < 2) then
           do p=1,N
                marker(p)%TRACK=.false.
           enddo
           do itrack = 1, ntrack
                do p=1,N
                        if(marker(p)%id == trackid(itrack)) then
                              marker(p)%TRACK=.true.
                        endif
                enddo
           enddo
        endif



        if (mod(nstep, 40) == 0) then
           nl = floor(1.0 * nstep / ndiag) + 1 
           np0 = 0
!           do itrack = 1, ntrack
!                flag_diag = .false.
                do p=1,N
!                        if(marker(p)%id == trackid(itrack)) then
                        if(marker(p)%TRACK) then

                        i = floor((marker(p)%X(1) - myxmin)/dxx) + 3
                        j = floor((marker(p)%X(2) - myymin)/dyy) + 3
                        k = floor((marker(p)%X(3) - myzmin)/dzz) + 3
             
                        R   = xx(i)
                        phi = yy(j)
                        Z   = zz(k)
             
                        dR2 = xx(i+1)**2 - xx(i)**2
             
                        weight_line(1,1) = (marker(p)%X(1)**2 - R**2)/dR2
                        weight_line(2,1) = (marker(p)%X(2) - phi)/dyy
                        weight_line(3,1) = (marker(p)%X(3) - Z)/dzz
             
                        weight_line(1,2) = 1.0 - weight_line(1,1)
                        weight_line(2,2) = 1.0 - weight_line(2,1)
                        weight_line(3,2) = 1.0 - weight_line(3,1)
             
                        do ii=1,2
                            do jj=1,2
                                do kk=1,2
                                    weight_cubic(ii,jj,kk) = weight_line(1,3-ii)*weight_line(2,3-jj)*weight_line(3,3-kk)
                                end do
                            end do
                        end do

                        weight_line(1,1) = (marker(p)%X(1) - R)/dxx
                        weight_line(2,1) = (marker(p)%X(3) - Z)/dzz
                     
                        weight_line(1,2) = 1.0 - weight_line(1,1)
                        weight_line(2,2) = 1.0 - weight_line(2,1)
                        
             
                        do ii=1,2
                            do jj=1,2
                                weight_square(ii,jj) = weight_line(1,3-ii)*weight_line(2,3-jj)
                            end do
                        end do
 
                        B   = 0
                        do ii=1,2
                            do jj=1,2
                                do kk=1,2!note that B in CLT code is Bx(R,Z,phi)
                                    B = B + x(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
                                end do
                            end do
                        end do   
                
                        abs_B = sqrt(dot_product(B,B))
                
                        psi_par = 0
                        do ii=1,2
                            do kk=1,2
                                psi_par = psi_par  + psi(i-1+ii,k-1+kk)*weight_square(ii,kk)
                            end do
                        end do  
                        marker(p)%v_perp = sqrt(2*marker(p)%mu*abs_B/m)
                        marker(p)%P_phi  = m*marker(p)%v_para*marker(p)%X(1)*(B(2)/abs_B) - Zh*psi_par
                        E      = 0.5*m*marker(p)%v_para**2 + marker(p)%mu*abs_B
                        Lambda = marker(p)%mu*B0/E
                        data_diag(1) = marker(p)%X(1)
                        data_diag(2) = marker(p)%X(2)
                        data_diag(3) = marker(p)%X(3)
                        data_diag(4) = marker(P)%v_para
                        data_diag(5) = marker(P)%P_phi
                        data_diag(6) = E
                        data_diag(7) = Lambda
                        data_diag(8) = time
!                        data_diag(9, itrack, nl) = trackid(itrack)
!                        flag_diag    = .true.
!                        if (flag_diag) then
                        np0 = np0 + 1
!                        pfid = trackid(itrack)
                        pfid = marker(p)%id
                        write(pf0,'(i7.7,".dat")')pfid
                        pf = './pdata/p'//trim(pf0)
                        if (nl > 1) then
                            open(unit = pfid,file = pf,status = 'old', &
                                 position='Append',form = 'binary')
                        else
                            open(unit = pfid,file = pf,status = 'replace',form = 'binary')
                        endif
                        write(pfid)data_diag
                        close(pfid)
!                        endif
!                        if (nrank == 0) then
!                            write(*,*)data_diag(:,itrack,nl),itrack,trackid(itrack),nl
!                        endif
                      end if
                   end do

!                   call MPI_REDUCE(data_diag, 9 * ntrack * (nstop / ndiag + 1), MPI_DOUBLE_PRECISION, data_diag_total, 9 * ntrack * (nstop / ndiag + 1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)
!                   call MPI_GATHER(flag_diag, 1, MPI_LOGICAL, flag_diag_total, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, IERROR)

!                   if(nrank == 0 .and. nstep == nstop) then
!                     open(unit=32,file='particle_orbit.dat',form='binary')
!                     do i=1,nsize
!                        if(flag_diag_total(i)) then
!                           write(32)data_diag_total
!                           write(*, *)data_diag_total(:,433,:)
!                        endif
!                     end do            
!                   end if
!                enddo
                call MPI_REDUCE(np0, np, 1 ,MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, IERROR)

                if(nrank.eq.0) then
                    open(unit=32,file='Pnumber.dat')
                    write(32,*)'np=',np
                    close(32)
                endif
!                if (nrank == 0) then
!                    if (nl > 1) then 
!                        do i=1,ntrack
!                            pfid = trackid(i)
!                            write(pf0,'(i7.7,".dat")')pfid
!                            pf = './pdata/p'//trim(pf0)
!                            open(unit = pfid,file = pf,status = 'old', position='Append',form = 'binary')
!                            write(pfid)data_diag_total(:,i)
!                            close(pfid)
!                         enddo
!                     else
!                        do i=1,ntrack
!                            pfid = trackid(i)
!                            write(pf0,'(i7.7,".dat")')pfid
!                            pf = './pdata/p'//trim(pf0)
!                            open(unit = pfid,file = pf,status = 'replace',form = 'binary')
!                            write(pfid)data_diag_total(:,i)
!                            close(pfid)
!                         enddo
!                     endif
                  
!                end if
            endif

    end subroutine diagnosis_orbit_zhao
