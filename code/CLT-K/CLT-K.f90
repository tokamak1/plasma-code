!  CLTK.f90 
!
!  FUNCTIONS:
!  CLTK - Entry point of console application.
!
    
    module distribution
        implicit none
        type particle
            real*8  :: X(3),v_para,v_perp,mu
            real*8  :: f_over_g,g,w
            real*8  :: P_phi
            integer :: id
            logical :: OUT_OF_BOUNDARY
            logical :: TRACK
        end type particle
        integer, parameter :: N_initia = 120000, N_buf_min = 500
        integer :: N, deltaN, N_total
        type(particle), allocatable :: marker(:), marker_initia(:)
        real*8 :: c_f, c_f_center
        real*8 :: beta_f
        real*8 :: mybeta_f
        real*8 :: beta_f_center
        real*8 :: mybeta_f_center
        real*8,parameter :: beta=0.008, beta_center=0.060
        integer :: p
        real*8  :: psmin_clt_k, psmax_clt_k
    end module distribution
    
    module var
        use var_CLT
        implicit none
        
        real*8, parameter :: alpha = 20!omega_c/omega_A        
        real*8, parameter :: Zh=1*alpha,m=1 !coefficient alpha is absorb into charge number Zh!!!!!
        real*8, parameter :: m_div_Zh = m/Zh, Zh_div_m = Zh/m !specific charge :  e/m
        real*8, parameter :: v_A = 1.0, v_0=1.7*v_A, v_c = 0.58*v_0, c_1 = 0.37, Delta_Lambda = 0.10, Lambda_0 = 0.10, B0 = 1.0, a = 1.0, PI = 3.1415926
        real*8, parameter :: v_para_max = 1.5*v_0, v_para_min = -1.5*v_0, v_perp_max= 1.5*v_0, v_perp_min= 0.0*v_0, v_max = 1.5*v_0
        real*8, parameter :: Delta_v = 0.20*v_A
!parameters for the distribution in real space        
        real*8, parameter :: psi_n = 0.25, Delta_n = 0.15, L_n = 0.05
        real*8 :: bracket_psi_norm
        
        real*8 ::  v_parallel,w,delta_f
        integer :: RK4
        real*8 :: Xold(3),Vold,Wold
        real*8 :: Xa(3),Va,Wa
        real*8 :: coeff,timestep,told
        real*8 :: B(3),B_star(3),B_star_parallel,abs_B,b_unit(3),curl_b(3),grad_B(3),delta_E(3),grad_RB_phi(3)
        real*8 :: B_eq(3),B_star_eq(3),B_star_parallel_eq,abs_B_eq,b_unit_eq(3),curl_b_eq(3),grad_B_eq(3),psi_par
!polarization drift related variables
        real*8 :: E_star(3), E_star_eq(3), v_E(3), curl_v_E(3), grad_v_E_2(3), pv_Ept(3), pb_unitpt(3), v_banos, v_banos_eq
        real*8, dimension(3) :: v_d, curl_M
        real*8 :: R2,R,phi,Z
        real*8 :: RHS_X(3),RHS_V,RHS_W,RHS_X_eq(3),RHS_V_eq
        real*8 :: mu,E,v,Lambda,P_phi
        real*8 :: f0,dP_phidt,dEdt,dLambdadt,df0dP_phi,df0dE,df0dLambda,df0dt
        real*8 :: grad_P_phi(3),dP_phidv_parallel,dXdt1(3),dv_paralleldt1,pBpt(3),pabs_Bpt
        real*8 :: bracket_psi
        integer :: sgn_v_parallel

        real*8           :: Delta_psi, R0
        integer :: i, j, k, ii, jj, kk, s
        real*8 :: dR2, weight_line(3,2), weight_square(2,2), weight_cubic(2,2,2)
        real*8, dimension(mx,mz,my)   :: B_grid
        real*8, dimension(mx,mz,my,3) :: grad_B_grid, curl_b_grid, b_unit_grid, pb_unitpt_grid
        real*8, dimension(mx,mz)      :: B_grid_eq
        real*8, dimension(mx,mz,3)    :: grad_B_grid_eq, curl_b_grid_eq, b_unit_grid_eq
        real*8, dimension(mx,mz,my)   :: B_grid_eq_3D 
        real*8, dimension(mx,mz,my,3) :: b_unit_grid_eq_3D, curl_b_grid_eq_3D 
        real*8  :: myxmin, myxmax, myymin, myymax, myzmin, myzmax
!MPI variables
        integer :: IERROR
    contains
        function cross_product(x,y)
            implicit none
            real*8, dimension(3), intent(in) :: x,y
            real*8, dimension(3) :: cross_product
            cross_product(1) = x(2)*y(3) - x(3)*y(2)
            cross_product(2) = x(3)*y(1) - x(1)*y(3)
            cross_product(3) = x(1)*y(2) - x(2)*y(1)
            return          
        end function cross_product
        
        function sech(x)
            implicit none
            real*8, intent(in) :: x
            real*8 :: sech
            sech = 1.0/cosh(x)
            return          
        end function sech
    end module var

    module diganostic
!diganostic variables for f vs. P_phi and E 2d-plot
        real*8, dimension(5) :: P_phi_max, P_phi_min, E_max, E_min
        real*8, dimension(5) :: dP_phi, dE
        integer, parameter :: GRID_P_PHI = 100, GRID_E = 100        
        integer, parameter :: NSTEP_AVG = 20 ! time average  for deltaf and f
        integer, parameter :: NSTEP_INT = 1000 ! time interval for deltaf and f
        integer, parameter :: NSTEP_START = 50*NSTEP_INT ! time to start recording
        real*8, dimension(GRID_P_PHI,GRID_E)   ::   f_vs_P_phi_and_E,   deltaf_vs_P_phi_and_E,   num_2d
        real*8, dimension(GRID_P_PHI,GRID_E,5) :: myf_vs_P_phi_and_E, mydeltaf_vs_P_phi_and_E, mynum_2d
        real*8, dimension(GRID_P_PHI,GRID_E)   :: df0dP_phi_vs_P_phi_and_E, df0dE_vs_P_phi_and_E
        real*8, dimension(GRID_P_PHI,GRID_E,5) :: mydf0dP_phi_vs_P_phi_and_E, mydf0dE_vs_P_phi_and_E
        integer :: num_plot, mynum_plot(5)
        integer :: nstep_local
    end module diganostic

    !num2str

    module strings

        ! GLOBAL FUNCTIONS
        public :: num2str

        ! Everything else is private
        private

        interface num2str
            module procedure num2str_int
            module procedure num2str_real
            module procedure num2str_double_precision
            module procedure num2str_double_precision_array
        end interface 

    contains
        function num2str_int(number)
            implicit none
            integer,intent(in) :: number
            character(len=10)   :: num2str_int
            character(len=10)   :: tmp

            write(tmp,'(I10)')number
            num2str_int = adjustl(tmp)
        end function

        function num2str_real(number)
            implicit none
            real,intent(in)    :: number
            character(len=10)   :: num2str_real
            character(len=10)   :: tmp

            write(tmp,'(F10.2)')number
            num2str_real = adjustl(tmp)
        end function

        function num2str_double_precision(number)
            implicit none
            double precision,intent(in)    :: number
            character(len=10)   :: num2str_double_precision
            character(len=10)   :: tmp

            if((number>1e-20 .and. number<1e-1) .or. number>1e+3) then
                write(tmp,'(E10.2)')number
            else
                write(tmp,'(F10.2)')number
            endif
            num2str_double_precision = adjustl(tmp)
        end function
        
        function num2str_double_precision_array(number,dim)
            implicit none
            integer,intent(in)   :: dim
            double precision,intent(in)    :: number(dim)
            character(len=(6+1)*dim)       :: num2str_double_precision_array
            character(len=(6+1)*dim)       :: tmp

            write(tmp,'(<dim>(F6.2))')number
            num2str_double_precision_array = adjustl(tmp)
        end function
    
    end module

!diganostic variables for Calculating Orbit Frequency
    module COF_var
        implicit none
        real*8  :: dZ_of_step1
        logical :: Z_reverse_sign ! Z form "-" to "+" or form "+" to "-"
        real*8  :: R_initia, phi_initia, Z_initia, E_initia, P_phi_initia, mu_initia
        real*8, parameter :: eps_COF = 1e-2
        real*8  :: omega_phi, omega_theta
        logical :: FINISH_COF
        integer, parameter :: NUM_E=100, NUM_PPHI=100
        integer :: TIMESTEP_COF
    end module COF_var
    
    
!****************************************************************************
!
!  PROGRAM: CLTK
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    subroutine CLTK
        use var
        use distribution
        use diganostic
        implicit none
        include 'mpif.h'

        ! Variables

        ! Body of CLTK
        if(nrank == 0) print *, 'Hello World'

        call calc_grad_B        
        call diagnosis_orbit_zhao
        call diagnosis_field

        if(nrank == 0) print *, '111111111111111'
        call calc_polarization_drift_var
        call particle_orbit_polarization
!        call particle_orbit
        call calc_current
        !diagnostic for distrubution
        if(FEP .and. nstep == 0) then
            do s=1,5
                call free_energy_vs_P_phi_and_E(Lambda_DEP(s), dLambda_DEP, s)
            end do
        end if

        if(DEP .or. DEB) then
            if(nstep >= NSTEP_START) then
                if(mod(nstep,NSTEP_INT) == 0) nstep_local = 0
                if(nstep_local<=NSTEP_AVG) then
                    do s=1,5
                        call diagnostics_output_f_vs_P_phi_and_E(Lambda_DEP(s), dLambda_DEP, s)
                    end do
                    nstep_local = nstep_local + 1
                end if
            end if
        end if


    end subroutine CLTK
 
    
!************************************************************************************************
!                                                                                               *
!                       push particle using 4th Runge-Kutta method                              * 
!                                                                                               *
!************************************************************************************************     
    subroutine particle_orbit
        use distribution
        use var
        implicit none
        include 'mpif.h'

        
        do p=1,N
            do RK4=1,4
                select case(RK4)
                case(1)
                    Xold     = marker(p)%X
                    Vold     = marker(p)%v_para
                    Wold     = marker(p)%w
                    Xa       = marker(p)%X
                    Va       = marker(p)%v_para
                    Wa       = marker(p)%w
                    told     = time
                    timestep = 0.5*dt
                    coeff    = 1.0/6.0
                case(2)
                    timestep = 0.5*dt
                    coeff    = 1.0/3.0
                case(3)
                    timestep = dt
                    coeff    = 1.0/3.0
                case(4)
                    timestep = dt
                    coeff    = 1.0/6.0
                end select
                
                if(GYRO) then
                    call gyro_average_xyz(N_GY)
                else
                    
                    i = floor((marker(p)%X(1) - myxmin)/dxx) + 3
                    j = floor((marker(p)%X(2) - myymin)/dyy) + 3
                    k = floor((marker(p)%X(3) - myzmin)/dzz) + 3
                
                    if(i>=mx .or. i<1 .or. j>=my .or. j<1 .or. k>=mz .or. k<1) then
                        write(*,*)"zzzzzzzzzzzzzzzzzaaaaaaaaaaaaaa : ",i,j,k,marker(p)%X,marker(p)%id,nrank,time
                        write(*,*)"nnnnnnnnnnnnnnnnnbbbbbbbbbbbbbb : ",i,j,k,marker_initia(p)%X,marker_initia(p)%v_para,marker_initia(p)%v_perp,marker_initia(p)%id,nrank,time
                    end if
                
             
                    R   = xx(i)
                    phi = yy(j)
                    Z   = zz(k)
             
                    dR2 = xx(i+1)**2 - xx(i)**2
             
                    weight_line(1,1) = (marker(p)%X(1)**2 - R**2)/dR2
                    weight_line(2,1) = (marker(p)%X(2) - phi)/dyy
                    weight_line(3,1) = (marker(p)%X(3) - Z)/dzz
                
                    if(psi(i,k)<psmax .and. (psi(i+1,k)>psmax .or. psi(i,k+1)>psmax .or. psi(i+1,k+1)>psmax)) then
                        weight_line(1,1) = 0.0
                        weight_line(3,1) = 0.0
                    endif
             
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
                
                    if(psi(i,k)<psmax .and. (psi(i+1,k)>psmax .or. psi(i,k+1)>psmax .or. psi(i+1,k+1)>psmax)) then
                        weight_line(1,1) = 0.0
                        weight_line(2,1) = 0.0
                    endif
             
                    weight_line(1,2) = 1.0 - weight_line(1,1)
                    weight_line(2,2) = 1.0 - weight_line(2,1)

             
                    do ii=1,2
                        do jj=1,2
                            weight_square(ii,jj) = weight_line(1,3-ii)*weight_line(2,3-jj)
                        end do
                    end do                
                
                    B = 0
                    grad_B  = 0
                    curl_b  = 0
                    delta_E = 0
                    grad_RB_phi = 0
                    pBpt    = 0

                    do ii=1,2
                        do jj=1,2
                            do kk=1,2!note that B in CLT code is Bx(R,Z,phi)
                                B       = B + x(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
                                grad_B  = grad_B  + grad_B_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                                curl_b  = curl_b  + curl_b_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                                delta_E = delta_E + Ef(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                                pBpt    = pBpt    + xdif(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
                                grad_RB_phi(1) = grad_RB_phi(1) + (x(i-1+ii,k-1+kk,j-1+jj,7)+xx(i-1+ii)*xr(i-1+ii,k-1+kk,j-1+jj,7))*weight_cubic(ii,jj,kk)
                                grad_RB_phi(2) = grad_RB_phi(2) + (                                     xy(i-1+ii,k-1+kk,j-1+jj,7))*weight_cubic(ii,jj,kk)
                                grad_RB_phi(3) = grad_RB_phi(3) + (                          xx(i-1+ii)*xz(i-1+ii,k-1+kk,j-1+jj,7))*weight_cubic(ii,jj,kk)
                            end do
                        end do
                    end do
                    abs_B    = sqrt(dot_product(B,B))
                    b_unit   = B/abs_B
                    pabs_Bpt = dot_product(b_unit,pBpt)
                
                
                
                    B_eq = 0
                    b_unit_eq = 0
                    grad_B_eq = 0
                    curl_b_eq = 0
                    psi_par   = 0
                    do ii=1,2
                        do kk=1,2
                            B_eq      = B_eq + xint(i-1+ii,k-1+kk,6:8)*weight_square(ii,kk)
                            b_unit_eq = b_unit_eq + b_unit_grid_eq(i-1+ii,k-1+kk,:)*weight_square(ii,kk)
                            grad_B_eq = grad_B_eq + grad_B_grid_eq(i-1+ii,k-1+kk,:)*weight_square(ii,kk)
                            curl_b_eq = curl_b_eq + curl_b_grid_eq(i-1+ii,k-1+kk,:)*weight_square(ii,kk)
                            psi_par   = psi_par  + psi(i-1+ii,k-1+kk)*weight_square(ii,kk)
                        end do
                    end do    
                    
                end if
                

                R          = marker(p)%X(1)
                v_parallel = marker(p)%v_para
                mu         = marker(p)%mu
                E          = 0.5*m*v_parallel**2 + mu*abs_B
                v          = sqrt(2*E/m)
                P_phi      = m*v_parallel*R*(B(2)/abs_B) - Zh*psi_par
                Lambda     = mu*B0/E
                
            
                B_star          = B + m_div_Zh*v_parallel*curl_b
                B_star_parallel = dot_product(B_star,b_unit)
        
                B_star_eq          = B_eq + m_div_Zh*v_parallel*curl_b_eq
                B_star_parallel_eq = dot_product(B_star_eq,b_unit_eq)
 
                
                RHS_X    =  1.0/B_star_parallel*(v_parallel*B_star + cross_product(b_unit,mu/Zh*grad_B-delta_E))
                RHS_V    =  1.0/m/B_star_parallel*dot_product(B_star,Zh*delta_E-mu*grad_B)
                RHS_X_eq =  1.0/B_star_parallel_eq*(v_parallel*B_star_eq + cross_product(b_unit_eq,mu/Zh*grad_B_eq))
                RHS_V_eq =  1.0/m/B_star_parallel_eq*dot_product(B_star_eq,-mu*grad_B_eq)
    !delta-f part
!                dXdt1          = 1.0/B_star_parallel*(v_parallel*(B-B_eq) - cross_product(b_unit_eq,delta_E))
!                dv_paralleldt1 = 1.0/m/B_star_parallel*(dot_product(B_star,Zh*delta_E) - dot_product(mu*(B-B_eq),grad_B))
                dXdt1          = RHS_X - RHS_X_eq 
                dv_paralleldt1 = RHS_V - RHS_V_eq
                
                
                v_d = 1.0/(Zh*B_star_parallel)*(m*v_parallel**2*curl_b + mu*cross_product(b_unit,grad_B))
        
    
                grad_P_phi(1) =   m*v_parallel/abs_B*grad_RB_phi(1) - m*v_parallel*R*B(2)/abs_B**2*grad_B(1) + Zh*B(3)*R
                grad_P_phi(2) =   m*v_parallel/abs_B*grad_RB_phi(2) - m*v_parallel*R*B(2)/abs_B**2*grad_B(2)
                grad_P_phi(3) =   m*v_parallel/abs_B*grad_RB_phi(3) - m*v_parallel*R*B(2)/abs_B**2*grad_B(3) - Zh*B(1)*R
                dP_phidv_parallel = m*R*B(2)/abs_B
        
        
                dP_phidt = dot_product(dXdt1,grad_P_phi) + dv_paralleldt1*dP_phidv_parallel
                dEdt     = Zh*dot_product(v_d + v_parallel*B/B_star_parallel,delta_E) + mu*pabs_Bpt
                
                if(FGP) then
                    dLambdadt = -Lambda**2*dEdt
                end if                
        
                if(v_parallel>0) then
                    sgn_v_parallel = +1
                else
                    sgn_v_parallel = -1
                end if
      
                if(ORBIT_AVERAGE_METHOD == 1) then
                    if((1-Lambda) > 0) then
                        bracket_psi = - P_phi/Zh + m_div_Zh*sgn_v_parallel*v*R0*sqrt(1-Lambda)
                    else
                        bracket_psi = - P_phi/Zh
                    end if
                else if(ORBIT_AVERAGE_METHOD == 2) then
                    bracket_psi = - P_phi/Zh
                end if
                
                    
    
        !Finite width of Gaussian distribution of Pitch angle  : f(Lambda) ~ exp[(Lambda-Lambda_0)^2/Delta_Lambda^2]
    
    
                if(DFR_type == 1) then
                    if(FGP) then
                        f0 = c_f*1.0/(v**3+v_c**3)*(1+erf((v_0-v)/Delta_v))*exp(-(bracket_psi/(c_1*Delta_psi)))*exp(-((Lambda-Lambda_0)/Delta_Lambda)**2)
                    else
                        f0 = c_f*1.0/(v**3+v_c**3)*(1+erf((v_0-v)/Delta_v))*exp(-(bracket_psi/(c_1*Delta_psi)))
                    end if
                else if(DFR_type == 2) then
                    bracket_psi_norm = (bracket_psi - psmin)/(psmax-psmin)
                    f0 = c_f*1.0/(v**3+v_c**3)*(1+erf((v_0-v)/Delta_v))*exp(-Delta_n/L_n*tanh((bracket_psi_norm - psi_n)/Delta_n))
                end if

                if(DFR_type == 1) then
                    df0dP_phi =  1.0/(Zh*c_1*Delta_psi)*f0
                else if(DFR_type == 2) then
                    df0dP_phi =  1.0/(Zh*L_n*Delta_psi)*sech((bracket_psi_norm - psi_n)/Delta_n)**2*f0
                end if
        
                if(ORBIT_AVERAGE_METHOD == 1) then
                    if((1-Lambda) > 0) then
                        if(DFR_type == 1) then
                            df0dE = -3/m*v/(v**3+ v_c**3)*f0 - exp(-((v_0-v)/Delta_v)**2)/(Delta_v*sqrt(PI*E*m/2)*(1+erf((v_0-v)/Delta_v)))*f0 - sgn_v_parallel*R0/(Zh*c_1*Delta_psi*sqrt(v*v-2*mu*B0/m))*f0
                        else if(DFR_type == 2) then
                            df0dE = -3/m*v/(v**3+ v_c**3)*f0 - exp(-((v_0-v)/Delta_v)**2)/(Delta_v*sqrt(PI*E*m/2)*(1+erf((v_0-v)/Delta_v)))*f0 - sgn_v_parallel*R0/(Zh*c_1*Delta_psi*sqrt(v*v-2*mu*B0/m))*sech((bracket_psi_norm - psi_n)/Delta_n)**2*f0
                        end if
                    else
                        df0dE = -3/m*v/(v**3+ v_c**3)*f0 - exp(-((v_0-v)/Delta_v)**2)/(Delta_v*sqrt(PI*E*m/2)*(1+erf((v_0-v)/Delta_v)))*f0
                    end if
                else if(ORBIT_AVERAGE_METHOD == 2) then
                    df0dE = -3/m*v/(v**3+ v_c**3)*f0 - exp(-((v_0-v)/Delta_v)**2)/(Delta_v*sqrt(PI*E*m/2)*(1+erf((v_0-v)/Delta_v)))*f0
                end if
                
                    
    
                if(FGP) then
                    df0dLambda =   -2*(Lambda-Lambda_0)/Delta_Lambda**2*f0
                end if
    
                if(FGP) then
                    df0dt = dP_phidt*df0dP_phi + dEdt*df0dE + dLambdadt*df0dLambda
                else
                    df0dt = dP_phidt*df0dP_phi + dEdt*df0dE
                end if
                
                
                
        
                RHS_W    = -1.0/marker(p)%g*df0dt
                
                if(abs(RHS_W) > 1E5) then
!                    write(*,*)"ahahahahaha",dP_phidt,df0dP_phi,dEdt,df0dE,dLambdadt,df0dLambda,bracket_psi,f0,v,c_f*1.0/(v**3+v_c*v_c*v_c),exp(-(bracket_psi/(c_1*Delta_psi))),Delta_psi,bracket_psi/(c_1*Delta_psi)
                    write(*,*)"ahahahahaha",f0,v,exp(-(bracket_psi/(c_1*Delta_psi))),bracket_psi,P_phi/Zh, m_div_Zh*sgn_v_parallel*v*R0*sqrt(1-Lambda),m*v_parallel*R*(B(2)/abs_B),Zh*psi_par
                    write(*,*)"gagagagga",i,j,k,marker(p)%X(1),marker(p)%X(2),marker(p)%X(3)

                    exit 
                endif
                
!                RHS_W    = m*v_parallel*RHS_V + mu*dot_product(RHS_X,grad_B)
        
        
                RHS_X(2) = RHS_X(2)/R
    
                marker(p)%X          = Xold + timestep*RHS_X
                marker(p)%v_para     = Vold + timestep*RHS_V
                marker(p)%w          = Wold + timestep*RHS_W
                Xa    = Xa + coeff*dt*RHS_X
                Va    = Va + coeff*dt*RHS_V
                Wa    = Wa + coeff*dt*RHS_W
                time  = told + timestep
                
                if(psi_par >= psmax_clt_k) then
                    marker(p)%OUT_OF_BOUNDARY = .true.
                    exit
                end if
                
            end do
            
            
            if(marker(p)%OUT_OF_BOUNDARY) then
                marker(p)%X          = Xold
                marker(p)%v_para     = Vold
                marker(p)%w          = Wold
                marker(p)%X(3)       =  - marker(p)%X(3)
            else
                marker(p)%X          = Xa
                marker(p)%v_para     = Va
                marker(p)%w          = Wa
            end if

!boundary condition
!		    do while(marker(p)%X(2) >= 2*PI)
!                marker(p)%X(2) = marker(p)%X(2) - 2*PI
!            end do
!            do while(marker(p)%X(2) < 0)
!                marker(p)%X(2) = marker(p)%X(2) + 2*PI
!            end do
            
            time = told
        end do
        
!        do p=1,N
!            if(marker(p)%id >N_initia*64) write(*,*)"333333hehehehehehehehe : ",marker(p)%id,time
!        end do
        
        if(RECYCLE_METHOD == 1) call recycle_particle
        if(RECYCLE_METHOD == 2) call remove_particle
!        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!        if(nrank == 0) write(*,*)"222222222222222222"
!        do p=1,N
!            if(marker(p)%id >N_initia*64) write(*,*)"hehehehehehehehehehehe : ",marker(p)%id,time
!        end do
        
        call update
!        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!        if(nrank == 0) write(*,*)"333333333333333333"
!        do p=1,N
!            if(marker(p)%id >N_initia*64) write(*,*)"121212121heheheheheheh : ",marker(p)%id,time,nrank
!        end do
        
        !output number of particle vs time
        call MPI_REDUCE(N,N_total,1,MPI_INTEGER,MPI_SUM, 0, MPI_COMM_WORLD, IERROR)
        if(nrank == 0) write(*,*)"number of particle : ",N_total," at time = ",time
        if(N > N_initia+deltaN) write(*,*)"more than maxima of N in per process!!!!!   N: ",N," N_initia+deltaN :", N_initia+deltaN," nrank :",nrank
        

                    
    end subroutine particle_orbit
 
!************************************************************************************************
!                                                                                               *
!                       recycle particles between 2 processes (from -Z(Z) to Z(-Z))             *
!                                                                                               *
!************************************************************************************************      
    subroutine recycle_particle_back
        use var
        use distribution
        implicit none
        include 'mpif.h'
        integer :: N_in,N_out
        integer :: dest,source, dest_nrkz
        integer, dimension(deltaN)    :: index_out
        type(particle), dimension(deltaN) :: buffer_out,buffer_in
        
        integer :: STATUS(MPI_STATUS_SIZE)
        
        N_out = 0
        N_in  = 0
        do p=1,N
            if(marker(p)%OUT_OF_BOUNDARY) then
                N_out = N_out  + 1
                marker(p)%OUT_OF_BOUNDARY = .false.
                buffer_out(N_out) = marker(p)
                marker(p)%id = -1
                index_out(N_out) = p
            end if
        end do
        
        if(N_out > deltaN) write(*,*)"aahahahahyuahauaiaao: ",N_out, N_in, deltaN
        

        !this is on vaild for even value of nrkz
        dest_nrkz = (nprz - 1) - nrkz(nrank)
        dest = nrky(nrank)*nprx*nprz + dest_nrkz*nprx + nrkx(nrank)
        
        source = dest
        call MPI_SENDRECV(N_out, 1, MPI_INTEGER, dest, 0, N_in, 1, MPI_INTEGER, source, 0, MPI_COMM_WORLD, STATUS, IERROR)

        call MPI_SENDRECV(buffer_out(1:N_out), sizeof(buffer_out(1))*N_out, MPI_BYTE, dest, 0, buffer_in(1:N_in), sizeof(buffer_in(1))*N_in, MPI_BYTE, source, 0, MPI_COMM_WORLD, STATUS, IERROR)

        if(N_in /= 0 .or. N_out /=0) then
        
            if(N_in >= N_out) then
                do i=1,N_out
                    marker(index_out(i)) = buffer_in(i)
                end do
                do i=1,N_in-N_out
                    marker(i+N) = buffer_in(i+N_out)
                end do
            else
                do i=1,N_in
                    marker(index_out(i)) = buffer_in(i)
                end do
                p = N
                k = 0!hole
                i = 0!number of moved marker
                j = 0
                do while(i+k < N_out - N_in)
                    j = j + 1
                    if(j > N) then
                        write(*,*)"cuocuocuocuo!!!!",nrank
                        exit
                    end if
                    
                    do while(marker(p)%id == -1 .and. p>0)
                        p = p - 1
                        k = k + 1
                    end do
                    if(k == N_out - N_in .or. i+k == N_out - N_in) then
                        exit
                    else
                        i = i + 1
                        marker(index_out(i+N_in)) = marker(p)
                        marker(p)%id = -1
                        if(p < index_out(i+N_in)) write(*,*)"yyyyyyyyyyyccccccccccccccc",p,nrank,N,N_in,N_out,index_out(i+N_in),index_out(1+N_in:i+N_out)
                        if(p < 1) write(*,*)"bbbbbbbbbbbbb",nrank,p,N_out,N_in,N,time,index_out
                        p = p -1
                    end if            
                end do
            end if
           
	        N = N + (N_in - N_out)
        end if
        
        do p=1,N
            if(marker(p)%id <0) write(*,*)"lalalalalalalalala: ",marker(p)%id,time,N,p,nrank
        end do
        
        if(N > N_initia+deltaN) write(*,*)"!!!!!N > N_initia+deltaN in recycle_particle()!!!! nrank : ",nrank
    
    end subroutine recycle_particle_back
    
    
!************************************************************************************************
!                                                                                               *
!                                  remove particles at boundary                                 *
!                                                                                               *
!************************************************************************************************  
    subroutine remove_particle_back
        use var
        use distribution
        implicit none
        include 'mpif.h'
        integer :: N_loss
        integer, dimension(deltaN)    :: index_loss
        
        
        N_loss = 0
        do p=1,N
            if(marker(p)%OUT_OF_BOUNDARY) then
                N_loss = N_loss  + 1
                marker(p)%id = -1
                index_loss(N_loss) = p
            end if
        end do
        

 
        if(N_loss /= 0) then
        
            p = N
            k = 0!hole
            i = 0!number of moved marker
            j = 0
            do while(i+k < N_loss)
                j = j + 1
                if(j > N) then
                    write(*,*)"cuocuocuocuo!!!!",nrank
                    exit
                end if
                    
                do while(marker(p)%id == -1)
                    p = p - 1
                    k = k + 1
                    if(p <= 0) exit
                end do
                if(k == N_loss .or. i+k == N_loss) then
                    exit
                else
                    i = i + 1
                    marker(index_loss(i)) = marker(p)
                    marker(p)%id = -1
                    p = p -1
                end if            
            end do
           
	        N = N - N_loss
        end if
    
    end subroutine remove_particle_back
    
!************************************************************************************************
!                                                                                               *
!                       transfer particles between adjacent 26 processes                        *
!                                                                                               *
!************************************************************************************************  
    subroutine update_back()
    use var
    use distribution
    implicit none
    include 'mpif.h'
    
    integer, dimension(3)         :: OUT_OF_PROCESS 
    integer, dimension(27)        :: N_in,N_out
    integer, dimension(deltaN)    :: index_out 
    integer, dimension(27)        :: dest, source
    integer :: dest_nrkx, dest_nrky, dest_nrkz
    integer :: index
    integer :: N_in_total, N_out_total
    type(particle), dimension(27,deltaN) :: buffer_out, buffer_in
    type(particle), dimension(27*deltaN) :: buffer_in_total
    
    integer :: STATUS(MPI_STATUS_SIZE)
    
    N_out = 0
    N_in  = 0
    i = 1
    
    do p=1,N
        if(marker(p)%X(1)< myxmin) then
            OUT_OF_PROCESS(1) = -1
        else if(marker(p)%X(1) > myxmax) then
            OUT_OF_PROCESS(1) = +1
        else
            OUT_OF_PROCESS(1) = 0
        end if  
        
        if(marker(p)%X(2)< myymin) then
            OUT_OF_PROCESS(2) = -1
        else if(marker(p)%X(2) > myymax) then
            OUT_OF_PROCESS(2) = +1
        else
            OUT_OF_PROCESS(2) = 0
        end if
        
        if(marker(p)%X(3)< myzmin) then
            OUT_OF_PROCESS(3) = -1
        else if(marker(p)%X(3) > myzmax) then
            OUT_OF_PROCESS(3) = +1
        else
            OUT_OF_PROCESS(3) = 0
        end if    
        
!        write(*,*)"hahahahah OUT_OF_PROCESS :",OUT_OF_PROCESS

        !boundary condition
        if(nrky(nrank) == npry-1 .and. OUT_OF_PROCESS(2) == +1) then!if particle is located at phi = PI, OUT_OF_PROCESS(2) = +1 should be changed to -(npry-1)
            marker(p)%X(2) = marker(p)%X(2) - 2*PI
        else if (nrky(nrank) == 0 .and. OUT_OF_PROCESS(2) == -1) then!if particle is located at phi = 0, OUT_OF_PROCESS(2) = -1 should be changed  to +(npry-1)
            marker(p)%X(2) = marker(p)%X(2) + 2*PI
        end if
        
        
        !N_out(14) = 0
        if(.not.(OUT_OF_PROCESS(1) == 0 .and. OUT_OF_PROCESS(2) == 0 .and. OUT_OF_PROCESS(3) == 0)) then
            index = (OUT_OF_PROCESS(1)+1)*9 + (OUT_OF_PROCESS(2)+1)*3 + (OUT_OF_PROCESS(3)+1) + 1
            N_out(index) = N_out(index) + 1
            index_out(i) = p
            buffer_out(index,N_out(index)) = marker(p)
			marker(p)%id = -1
            i = i + 1
!            if(nrank == 47) write(*,*)N_out
        end if
        
        
    end do
    
    do i=-1,1
        do j=-1,1
            do k=-1,1
                OUT_OF_PROCESS(1) = i
                OUT_OF_PROCESS(2) = j
                OUT_OF_PROCESS(3) = k
                index     = (OUT_OF_PROCESS(1)+1)*9 + (OUT_OF_PROCESS(2)+1)*3 + (OUT_OF_PROCESS(3)+1) + 1
                dest_nrkx = nrkx(nrank) + OUT_OF_PROCESS(1)
                dest_nrky = nrky(nrank) + OUT_OF_PROCESS(2)
                dest_nrkz = nrkz(nrank) + OUT_OF_PROCESS(3)
                
                !boundary condition
                if(dest_nrky == npry) dest_nrky = 0
                if(dest_nrky == -1)   dest_nrky = npry-1
                
                if(dest_nrkx > nprx - 1 .or. dest_nrkx < 0 .or. dest_nrkz > nprz - 1 .or. dest_nrkz < 0) then
                    dest(index) = -1
                else
                    dest(index) = dest_nrky*nprx*nprz + dest_nrkz*nprx + dest_nrkx
                end if
            end do
        end do
    end do
    
    source = dest  

    do i=1,27
        if(i /= 14 .and. dest(i) /= -1) then
            call MPI_SENDRECV(N_out(i), 1, MPI_INTEGER, dest(i), 0, N_in(i), 1, MPI_INTEGER, source(i), 0, MPI_COMM_WORLD, STATUS, IERROR)
         end if
    end do

    do i=1,27
        if(i /= 14 .and. dest(i) /= -1) call MPI_SENDRECV(buffer_out(i,1:N_out(i)), sizeof(buffer_out(i,1))*N_out(i), MPI_BYTE, dest(i), 0, buffer_in(i,1:N_in(i)), sizeof(buffer_in(i,1))*N_in(i), MPI_BYTE, source(i), 0, MPI_COMM_WORLD, STATUS, IERROR)
    end do

    N_in_total  = sum(N_in)
    N_out_total = sum(N_out)
    
!    do i=1,27
!        do p=1,N_out(i)
!            if(buffer_out(i,p)%id >N_initia*64) write(*,*)"wwuwuwuwuwuw : ",buffer_out(i,p)%id,buffer_out(i,p)%X,time,nrank
!        end do
        
!        do p=1,N_in(i)
!            if(buffer_in(i,p)%id >N_initia*64) write(*,*)"lalalalalwuwuwuw : ",buffer_in(i,p)%id,buffer_in(i,p)%X,time,nrank
!        end do
!    end do
    
    i = 1
    do j=1,27
        do k=1,N_in(j)
            buffer_in_total(i) = buffer_in(j,k) 
            i = i + 1
        end do
    end do
    
!    do i=1,N_in_total
!        if(buffer_in_total(i)%id >N_initia*64) write(*,*)"henghenghengwuwuwuw : ",buffer_in_total(i)%id,buffer_in_total(i)%X,time,nrank
!    end do
    
  
    if(N_in_total >= N_out_total) then
        do i=1,N_out_total
            marker(index_out(i)) = buffer_in_total(i)
        end do
        do i=1,N_in_total-N_out_total
            marker(i+N) = buffer_in_total(i+N_out_total)
        end do
    else
        do i=1,N_in_total
            marker(index_out(i)) = buffer_in_total(i)
        end do
        p = N
        k = 0!hole
        i = 0!number of moved marker
        do while(i+k < N_out_total - N_in_total)
            if(p > 0) then
                do while(marker(p)%id == -1)
                    p = p - 1
                    k = k + 1
                    if(p <= 0) exit
                end do
            end if
            if(k == N_out_total - N_in_total .or. i+k == N_out_total - N_in_total) then
                exit
            else
                i = i + 1
                marker(index_out(i+N_in_total)) = marker(p)
                marker(p)%id = -1
                if(p < index_out(i+N_in_total)) write(*,*)"yyyyyyyyyyyaaaaaaaaaaaaaaaaa",p,nrank,N,N_in_total,N_out_total,index_out(i+N_in_total),index_out(1+N_in_total:i+N_out_total)
                if(p < 1) write(*,*)"jiumingjiiiiiiiiiiiiiiiii",nrank,p,N_out_total,N_in_total,N,time,index_out
                p = p -1
            end if            
        end do
    end if
    
        
	N = N + (N_in_total - N_out_total)
    
    if(N > N_initia+deltaN) write(*,*)"!!!!!N > N_initia+deltaN in update()!!!! nrank : ",nrank
    
    
    end subroutine update_back

!************************************************************************************************
!                                                                                               *
!                                   loading particles                                           *
!               uniform loading in (R^2, phi, Z, v_para, v_perp^2) 5d phase space               *
!                                                                                               *
!************************************************************************************************
    subroutine particle_load_back
        use distribution
        use var
        implicit none
        include 'mpif.h'

        real*8 :: rn(3)!random number
        real*8, dimension(N_initia) :: f,g,f_over_g
        real*8 :: Volume, myVolume, V_cell(2)
        real*8, dimension(mx,mz,my) :: n_R_Z
        integer :: N_OUT_OF_BOUNDARY
        integer, dimension(N) :: index_hole


!initialize 
        mybeta_f = 0
        mybeta_f_center = 0
        
        
!  particle's position loading
!-------------------------------------------------------------------------------|
!                                                                               |
!                   Random loading in real*8 space (R^2,phi,Z)                  |  
!                                                                               |
!-------------------------------------------------------------------------------|   
        call init_random_seed(nrank)
        N_OUT_OF_BOUNDARY = 0
        do p=1,N
            call RANDOM_NUMBER(rn)
            R2  = myxmin**2 + (myxmax**2 - myxmin**2)*rn(1)
            phi = myymin + (myymax - myymin)*rn(2)
            Z   = myzmin + (myzmax - myzmin)*rn(3)
            marker(p)%X(1) = sqrt(R2)
            marker(p)%X(2) = phi
            marker(p)%X(3) = Z
                
            i = floor((marker(p)%X(1) - myxmin)/dxx) + 3
            j = floor((marker(p)%X(2) - myymin)/dyy) + 3
            k = floor((marker(p)%X(3) - myzmin)/dzz) + 3
             
            R   = xx(i)
            phi = yy(j)
            Z   = zz(k)
             
            weight_line(1,1) = (marker(p)%X(1) - R)/dxx
            weight_line(2,1) = (marker(p)%X(3) - Z)/dzz
             
            weight_line(1,2) = 1.0 - weight_line(1,1)
            weight_line(2,2) = 1.0 - weight_line(2,1)
                
             
            do ii=1,2
                do jj=1,2
                    weight_square(ii,jj) = weight_line(1,3-ii)*weight_line(2,3-jj)
                end do
            end do
                 
            psi_par = 0
            do ii=1,2
                do kk=1,2
                    psi_par = psi_par  + psi(i-1+ii,k-1+kk)*weight_square(ii,kk)
                end do
            end do 
                
            if( psi_par < psmax_clt_k ) then
                marker(p)%OUT_OF_BOUNDARY = .false.
                marker(p)%id   = nrank*N + p
            else
                marker(p)%OUT_OF_BOUNDARY = .false.
                marker(p)%id   = -1
                N_OUT_OF_BOUNDARY = N_OUT_OF_BOUNDARY + 1
                index_hole(N_OUT_OF_BOUNDARY) = p
            end if                
        end do
       ! sorting particles
        if(N_OUT_OF_BOUNDARY /= 0) then
            p = N
            k = 0!hole
            i = 0!number of moved marker
            do while(i+k < N_OUT_OF_BOUNDARY)
                do while(marker(p)%id == -1)
                    p = p - 1
                    k = k + 1
                    if(p < 1) exit
                end do
                if(k == N_OUT_OF_BOUNDARY .or. i+k == N_OUT_OF_BOUNDARY) then
                    exit
                else
                    i = i + 1
                    marker(index_hole(i)) = marker(p)
                    marker(p)%id = -1
                    if(p < index_hole(i)) write(*,*)"yyyyyyyyyyyccccccccccccccc"
                    if(p < 1) write(*,*)"bbbbbbbbbbbbb"
                    p = p -1
                end if            
            end do
        end if
        

        !re-define value of N
        N = N - N_OUT_OF_BOUNDARY
        
        
        
!-------------------------------------------------------------------------------|
!                                                                               |
!   Random loading particle in velocity space (v_parallel,v_perpendicular^2)    |
!                                                                               |
!-------------------------------------------------------------------------------|    
        do p=1,N
            marker(p)%v_para = v_max
            marker(p)%v_perp = v_max
            do while(sqrt(marker(p)%v_para**2 + marker(p)%v_perp**2) > v_max)
                call RANDOM_NUMBER(rn)
                marker(p)%v_para = v_para_min + (v_para_max  - v_para_min)*rn(1)
                !for mu = 0 case!!!
                if((.not. IIV) .and. (.not. FGP)) then
                    marker(p)%v_perp = 0.0
                else
                    marker(p)%v_perp = sqrt(v_perp_min**2 + (v_perp_max**2-v_perp_min**2)*rn(2))
                end if
            end do
        end do
 
        
!initial value of weight w is zero
        do p=1,N
            marker(p)%w  = 0
        end do
        
                    
!-------------------------------------------------------------------------------|
!                                                                               |
!                       Compute Distrubution Function                           |
!                                                                               |
!-------------------------------------------------------------------------------|
!        if(nrank == 47) N = N_initia
        do p=1,N
!            if(nrank == 47) then
!                marker(p)%X(1) = myxmin + (myxmax - myxmin)*0.9
!                marker(p)%X(2) = myymin + (myymax - myymin)*0.5
!                marker(p)%X(3) = myzmin + (myzmax - myzmin)*0.5
!                marker(p)%v_para = 0.3*v_max
!                marker(p)%v_perp = 0.8*v_max
!                marker(p)%X(1) = 3.63358396135151
!                marker(p)%X(2) = 1.05892179798067
!                marker(p)%X(3) = 0.454638965415636
!                marker(p)%v_para = 0.580036127656559
!                marker(p)%v_perp = 0.751244187653075
!                marker(p)%X(1) = 3.93358396135151
!                marker(p)%X(2) = 1.05892179798067
!                marker(p)%X(3) = 0.254638965415636
!                marker(p)%v_para = -0.380036127656559
!                marker(p)%v_perp = 0.851244187653075
!                marker(p)%id   = nrank*N + p
!                write(*,*)"marker(p)%X : ",marker(p)%X,"marker(p)%v_para",marker(p)%v_para,"marker(p)%v_perp",marker(p)%v_perp," nrank :",nrank
!            end if
            
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
             
            B = 0
            b_unit  = 0
            curl_b  = 0
            do ii=1,2
                do jj=1,2
                    do kk=1,2!note that B in CLT code is Bx(R,Z,phi)
                        B      = B + x(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
                        b_unit = b_unit + b_unit_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                        curl_b = curl_b + curl_b_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                    end do
                end do
            end do
            abs_B = sqrt(dot_product(B,B))
            
            psi_par = 0
            do ii=1,2
                do kk=1,2
                    psi_par       = psi_par  + psi(i-1+ii,k-1+kk)*weight_square(ii,kk)
                end do
            end do
            
             
            marker(p)%mu    = 0.5*m*marker(p)%v_perp**2/abs_B
            marker(p)%P_phi = m*marker(p)%v_para*marker(p)%X(1)*(B(2)/abs_B) - Zh*psi_par
             
            E      = 0.5*m*marker(p)%v_para**2 + marker(p)%mu*abs_B
            v      = sqrt(2*E/m)
            mu     = marker(p)%mu
            P_phi  = marker(p)%P_phi
            Lambda = marker(p)%mu*B0/E
            v_parallel = marker(p)%v_para
             
            if(v_parallel>0) then
                sgn_v_parallel = +1
            else
                sgn_v_parallel = -1
            end if
      
            if(ORBIT_AVERAGE_METHOD == 1) then
                if((1-Lambda)>0) then
                    bracket_psi = - P_phi/Zh + m_div_Zh*sgn_v_parallel*v*R0*sqrt(1-Lambda)
                else
                    bracket_psi = - P_phi/Zh
                end if
            else if(ORBIT_AVERAGE_METHOD == 2) then
                bracket_psi = - P_phi/Zh
            end if
            
                

            
!		uniform g in guiding center coordinates is : B/B^*_||
            B_star          = B + m_div_Zh*v_parallel*curl_b
            B_star_parallel = dot_product(B_star,b_unit)
            g(p) = abs_B/B_star_parallel
            
!-------------------------------------------------------------------------------|
!                                                                               |
!       f(P_phi,E,Lambda) ~   1/(v^3+v_c^3)                                     |
!                            *exp(-(bracket_psi/(c_1*Delta psi_par)))           |
!                            *exp(-((Lambda-Lambda_0)/Delta_Lambda)**2)         |
!                                                                               |
!-------------------------------------------------------------------------------|
            if(DFR_type == 1) then
                if(FGP) then
                    f(p) = 1.0/(v**3+v_c**3)*(1+erf((v_0-v)/Delta_v))*exp(-(bracket_psi/(c_1*Delta_psi)))*exp(-((Lambda-Lambda_0)/Delta_Lambda)**2)
                else
                    f(p) = 1.0/(v**3+v_c**3)*(1+erf((v_0-v)/Delta_v))*exp(-(bracket_psi/(c_1*Delta_psi)))
                end if
            else if(DFR_type == 2) then
                bracket_psi_norm = (bracket_psi - psmin)/(psmax-psmin)
                f(p) = 1.0/(v**3+v_c**3)*(1+erf((v_0-v)/Delta_v))*exp(-Delta_n/L_n*tanh((bracket_psi_norm - psi_n)/Delta_n))
            end if
                

!-------------------------------------------------------------------------------|
!                                                                               |
!                calculate volume average fast ion beta in code                 |
!                                                                               |
!-------------------------------------------------------------------------------|

            mybeta_f = mybeta_f + (m*marker(p)%v_para**2 + marker(p)%mu*abs_B)/abs_B/abs_B*f(p)/g(p)
            
!-------------------------------------------------------------------------------|
!                                                                               |
!                     calculate ceneter fast ion beta in code                   |
!                                                                               |
!-------------------------------------------------------------------------------|

            if(psi_par < (psmin + 0.1*(psmax - psmin))) then
                mybeta_f_center = mybeta_f_center + (m*marker(p)%v_para**2 + marker(p)%mu*abs_B)/abs_B/abs_B*f(p)/g(p)
            end if
                 
        end do
        
        call MPI_REDUCE(mybeta_f, beta_f, 1, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, IERROR)
        call MPI_BCAST(beta_f, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)
        
        call MPI_REDUCE(mybeta_f_center, beta_f_center, 1, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, IERROR)
        call MPI_BCAST(beta_f_center, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)
        
        beta_f_center = beta_f_center/(PI*(0.1*a)**2*(2*PI*R0))
        !normalize by center beta
        c_f_center = beta_center/beta_f_center

        
        myVolume = 0.0
        do i=3,mx-2
            do j = 3,my-2
                do k = 3,mz-2
                    if(psi(i,k) <= psmax_clt_k) then
                        myVolume = myVolume + dyy*dzz*dR2/2
                    end if
                end do
            end do
        end do
        
        call MPI_REDUCE(myVolume,Volume,1,MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, IERROR)
        call MPI_BCAST(Volume, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)        
        
        beta_f = beta_f/Volume

!normalize by volume averaged beta
        c_f = beta/beta_f
        
        f = c_f*f
        do p=1,N
            marker(p)%f_over_g = f(p)/g(p)
            marker(p)%g        = g(p)
            marker_initia(p)   = marker(p)
        end do
        

        
        write(*,*)"Volume =",Volume,"V_0 = ",(PI*(0.1*a)**2*(2*PI*R0))," c_f =",c_f," c_f_center =",c_f_center," mybeta_f = ",mybeta_f," N = ",N,"myid = ",nrank
        
        
        call MPI_REDUCE(N,N_total,1,MPI_INTEGER,MPI_SUM, 0, MPI_COMM_WORLD, IERROR)
        if(nrank == 0) write(*,*)"initial number of particle : ",N_total
        
        
!------------------------------------------------------------------------------|
!                                                                              |
!           output initial particle distrbution                                |
!                                                                              |
!------------------------------------------------------------------------------|
        !n_R_Z = 0.0
        !open(unit=31,file='g_R_Z.dat')
        !do p=1,N
        !    i = floor((marker(p)%X(1) - myxmin)/dxx) + 3
        !    j = floor((marker(p)%X(2) - myymin)/dyy) + 3
        !    k = floor((marker(p)%X(3) - myzmin)/dzz) + 3
        !     
        !    R   = xx(i)
        !    phi = yy(j)
        !    Z   = zz(k)
        !     
        !    dR2 = xx(i+1)**2 - xx(i)**2
        !     
        !    weight_line(1,1) = (marker(p)%X(1)**2 - R**2)/dR2
        !    weight_line(2,1) = (marker(p)%X(2) - phi)/dyy
        !    weight_line(3,1) = (marker(p)%X(3) - Z)/dzz
        !     
        !    weight_line(1,2) = 1.0 - weight_line(1,1)
        !    weight_line(2,2) = 1.0 - weight_line(2,1)
        !    weight_line(3,2) = 1.0 - weight_line(3,1)
        !     
        !    do ii=1,2
        !        do jj=1,2
        !            do kk=1,2
        !                weight_cubic(ii,jj,kk) = weight_line(1,3-ii)*weight_line(2,3-jj)*weight_line(3,3-kk)
        !            end do
        !        end do
        !    end do
        !    
        !    V_cell = dyy*dzz*dR2/2
        !    do ii=1,2
        !        do jj=1,2
        !            do kk=1,2!note that field in CLT code is Bx(R,Z,phi)
        !                if(j == my .and. jj == 2) then
        !                    n_R_Z(i-1+ii,k-1+kk,j-1+jj) = n_R_Z(i-1+ii,k-1+kk,1) + 1.0/V_cell*weight_cubic(ii,jj,kk)
        !                else
        !                    n_R_Z(i-1+ii,k-1+kk,j-1+jj) = n_R_Z(i-1+ii,k-1+kk,j-1+jj) + 1.0/V_cell*weight_cubic(ii,jj,kk)
        !                end if
        !            end do
        !        end do
        !    end do
        !end do        
        !do i=1,mx
        !    do k=1,mz
        !        write(31,"(1x,e12.5)")n_R_Z(i,k,1)
        !    end do
        !end do
        !
         
        
    end subroutine particle_load_back
 
    
!************************************************************************************************
!                                                                                               *
!                       calculate current of enegertic particles : J_h                          *
!                                                                                               *
!************************************************************************************************      
    subroutine calc_current
        use distribution
        use var
        implicit none
        include 'mpif.h'
        real*8, dimension(2) :: V_cell
        real*8               :: d1f2, d1fc
        real*8 f_jm1,f_jm2,f_j,f_jp1,f_jp2,x_jm1,x_j,x_jp1,coeff_a,coeff_b,coeff_c,coeff_d
        integer ix_first, ix_last,iz_first, iz_last,iy_first,iy_last
        integer jx,jy,jz
        real*8 boundary_damping

!  define statement functions
!  d1f2= d f / dx  with second-order accuracy central difference
        
        d1f2(f_jm1,f_j,f_jp1,x_jm1,x_j,x_jp1)= &
        ((x_jm1-x_j)/(x_jp1-x_j)*(f_jp1-f_j) &
        -(x_jp1-x_j)/(x_jm1-x_j)*(f_jm1-f_j))/(x_jm1-x_jp1) 
        
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(f_jm2,f_jm1,f_j,f_jp1,f_jp2,coeff_a,coeff_b,coeff_c,coeff_d)= &
       coeff_a*(f_jp1-f_j)+coeff_b*(f_j-f_jm1)+coeff_c*(f_jp2-f_j)+coeff_d*(f_j-f_jm2)
        
        
        ix_first=1
        ix_last=mx
        iz_first=1
        iz_last=mz
        iy_first=1
        iy_last=my
        
        p_h_perp = 0
        p_h_para = 0
        nv_para  = 0
        n_h      = 0
        
        do p=1,N
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
            
            if(psi(i,k)<psmax .and. (psi(i+1,k)>psmax .or. psi(i,k+1)>psmax .or. psi(i+1,k+1)>psmax)) then
                weight_line(1,1) = 0.0
                weight_line(3,1) = 0.0
            endif
             
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
            
            B = 0
            do ii=1,2
                do jj=1,2
                    do kk=1,2!note that B in CLT code is Bx(R,Z,phi)
                        B      = B + x(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
                    end do
                end do
            end do
            abs_B = sqrt(dot_product(B,B))
                
            R          = marker(p)%X(1)
            v_parallel = marker(p)%v_para
            mu         = marker(p)%mu
                
            delta_f = marker(p)%w*marker(p)%g
            
            dR2 = ((xx(i)+dxx+xx(i))/2)**2 - ((xx(i)+xx(i)-dxx)/2)**2
            V_cell(1) = dyy*dzz*dR2/2
            dR2 = ((xx(i+1)+dxx+xx(i+1))/2)**2 - ((xx(i+1)+xx(i+1)-dxx)/2)**2
            V_cell(2) = dyy*dzz*dR2/2
            
            do ii=1,2
                do jj=1,2
                    do kk=1,2!note that field in CLT code is Bx(R,Z,phi)
                       if(psi(i-1+ii,k-1+kk) <= (psmin_clt_k + 1.0*(psmax_clt_k - psmin_clt_k))) then
                            p_h_perp(i-1+ii,k-1+kk,j-1+jj) = p_h_perp(i-1+ii,k-1+kk,j-1+jj) + mu*abs_B*delta_f/V_cell(ii)*weight_cubic(ii,jj,kk)
                            p_h_para(i-1+ii,k-1+kk,j-1+jj) = p_h_para(i-1+ii,k-1+kk,j-1+jj) + m*v_parallel**2*delta_f/V_cell(ii)*weight_cubic(ii,jj,kk)
                            nv_para(i-1+ii,k-1+kk,j-1+jj)  = nv_para(i-1+ii,k-1+kk,j-1+jj)  + v_parallel*delta_f/V_cell(ii)*weight_cubic(ii,jj,kk)
                    if(VEB) n_h(i-1+ii,k-1+kk,j-1+jj)      = n_h(i-1+ii,k-1+kk,j-1+jj)      + delta_f/V_cell(ii)*weight_cubic(ii,jj,kk)
                       endif
                    end do
                end do
            end do
        end do
        
!update p_h_perp,p_h_para,nv_para,n_h on boundary of MPI sub-domian
        call update_scalar_sub_domain_bounnary(p_h_perp)
        call update_scalar_sub_domain_bounnary(p_h_para)
        call update_scalar_sub_domain_bounnary(nv_para)
        if(VEB) call update_scalar_sub_domain_bounnary(n_h)
        

!first smooth!!!!!
!        do s=1,10
!            call valbm_atlastgrid_v1(p_h_perp,1,0)
!            call smooth26(p_h_perp, coeff_smooth, 1)
!        end do
            
!possion smooth!!!!!
        call possion_solver_3D(p_h_para,J_h_old(:,:,:,1))
        call possion_solver_3D(p_h_perp,J_h_old(:,:,:,2))
        call possion_solver_3D(nv_para,J_h_old(:,:,:,3))
        if(VEB) call possion_solver_3D(n_h,n_h_old)
        
        J_h_old(:,:,:,1) = p_h_para
        J_h_old(:,:,:,2) = p_h_perp
        J_h_old(:,:,:,3) = nv_para
        if(VEB) n_h_old = n_h
        
        
        if(ADI) then
            !substract adiabatic term : xi dot grad f0!!!!
            !calc plasma displayment 
            xi = xi + x(:,:,:,3:5)*dt
        
            do i=3,mx-2
                do j=3,mz-2
                    do k=3,my-2
                        p_h_para(i,j,k) = p_h_para(i,j,k) + dot_product(xi(i,j,k,:),grad_p_h_para0(i,j,:))
                        p_h_perp(i,j,k) = p_h_perp(i,j,k) + dot_product(xi(i,j,k,:),grad_p_h_perp0(i,j,:))
                        nv_para(i,j,k)  = nv_para(i,j,k)  + dot_product(xi(i,j,k,:),grad_nv_para0(i,j,:))
                if(VEB) n_h(i,j,k)      = n_h(i,j,k)      + dot_product(xi(i,j,k,:),grad_n_h0(i,j,:))
                    end do
                end do
            end do
        end if
        
!mpi_transfersm to perpare grad_p_h_perp       
        call mpi_transfersm(p_h_perp,1)
        call valbm_atlastgrid_v1(p_h_perp,1,0)
        
        do jy=iy_first,iy_last
            do jz=iz_first+2,iz_last-2
                do jx=ix_first+2,ix_last-2
                    if(FOUR_TH) then
                        grad_p_h_perp(jx,jz,jy,1) = d1fc(p_h_perp(jx-2,jz,jy),p_h_perp(jx-1,jz,jy),p_h_perp(jx,jz,jy),p_h_perp(jx+1,jz,jy),p_h_perp(jx+2,jz,jy),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
                        grad_p_h_perp(jx,jz,jy,3) = d1fc(p_h_perp(jx,jz-2,jy),p_h_perp(jx,jz-1,jy),p_h_perp(jx,jz,jy),p_h_perp(jx,jz+1,jy),p_h_perp(jx,jz+2,jy),az1(jz),bz1(jz),cz1(jz),dz1(jz))
                    else
                        grad_p_h_perp(jx,jz,jy,1) = d1f2(p_h_perp(jx-1,jz,jy),p_h_perp(jx,jz,jy),p_h_perp(jx+1,jz,jy),xx(jx-1),xx(jx),xx(jx+1))
                        grad_p_h_perp(jx,jz,jy,3) = d1f2(p_h_perp(jx,jz-1,jy),p_h_perp(jx,jz,jy),p_h_perp(jx,jz+1,jy),zz(jz-1),zz(jz),zz(jz+1))
                    end if
                end do
            end do
        end do
        
        do jz=iz_first,iz_last
            do jx=ix_first,ix_last
                do jy=iy_first+2,iy_last-2
                    if(FOUR_TH) then
                        grad_p_h_perp(jx,jz,jy,2) = d1fc(p_h_perp(jx,jz,jy-2),p_h_perp(jx,jz,jy-1),p_h_perp(jx,jz,jy),p_h_perp(jx,jz,jy+1),p_h_perp(jx,jz,jy+2),ay1(jy),by1(jy),cy1(jy),dy1(jy))/xx(jx)
                    else
                        grad_p_h_perp(jx,jz,jy,2) = d1f2(p_h_perp(jx,jz,jy-1),p_h_perp(jx,jz,jy),p_h_perp(jx,jz,jy+1),yy(jy-1),yy(jy),yy(jy+1))/xx(jx)
                    end if
                end do
            end do
        end do
        
        
        
!guiding centercurrent + magnetization current
        J_h(:,:,:,1)  = Zh*nv_para*b_unit_grid(:,:,:,1) + (p_h_para-p_h_perp)*curl_b_grid(:,:,:,1)/B_grid + (b_unit_grid(:,:,:,2)*grad_p_h_perp(:,:,:,3) - b_unit_grid(:,:,:,3)*grad_p_h_perp(:,:,:,2))/B_grid
        J_h(:,:,:,2)  = Zh*nv_para*b_unit_grid(:,:,:,2) + (p_h_para-p_h_perp)*curl_b_grid(:,:,:,2)/B_grid + (b_unit_grid(:,:,:,3)*grad_p_h_perp(:,:,:,1) - b_unit_grid(:,:,:,1)*grad_p_h_perp(:,:,:,3))/B_grid 
        J_h(:,:,:,3)  = Zh*nv_para*b_unit_grid(:,:,:,3) + (p_h_para-p_h_perp)*curl_b_grid(:,:,:,3)/B_grid + (b_unit_grid(:,:,:,1)*grad_p_h_perp(:,:,:,2) - b_unit_grid(:,:,:,2)*grad_p_h_perp(:,:,:,1))/B_grid
        
        if(Varying_b_unit) then
            J_h(:,:,:,1) = J_h(:,:,:,1) + Zh*nv_para0_3D*(b_unit_grid(:,:,:,1)-b_unit_grid_eq_3D(:,:,:,1)) + &
            ( B_grid_eq_3D*(p_h_para0_3D-p_h_perp0_3D)*(curl_b_grid(:,:,:,1)-curl_b_grid_eq_3D(:,:,:,1)) - &
              (B_grid-B_grid_eq_3D)*(p_h_para0_3D-p_h_perp0_3D)*curl_b_grid_eq_3D(:,:,:,1) )/( B_grid*B_grid_eq_3D ) + &
            ( B_grid_eq_3D* ( (b_unit_grid(:,:,:,2)-b_unit_grid_eq_3D(:,:,:,2))*grad_p_h_perp0_3D(:,:,:,3) - &
                              (b_unit_grid(:,:,:,3)-b_unit_grid_eq_3D(:,:,:,3))*grad_p_h_perp0_3D(:,:,:,2) ) - &
              (B_grid-B_grid_eq_3D)*( b_unit_grid_eq_3D(:,:,:,2)*grad_p_h_perp0_3D(:,:,:,3) - &
                                      b_unit_grid_eq_3D(:,:,:,3)*grad_p_h_perp0_3D(:,:,:,2) ) )/( B_grid*B_grid_eq_3D )

            J_h(:,:,:,2) = J_h(:,:,:,2) + Zh*nv_para0_3D*(b_unit_grid(:,:,:,2)-b_unit_grid_eq_3D(:,:,:,2)) + &
            ( B_grid_eq_3D*(p_h_para0_3D-p_h_perp0_3D)*(curl_b_grid(:,:,:,2)-curl_b_grid_eq_3D(:,:,:,2)) - &
              (B_grid-B_grid_eq_3D)*(p_h_para0_3D-p_h_perp0_3D)*curl_b_grid_eq_3D(:,:,:,2) )/( B_grid*B_grid_eq_3D ) + &
            ( B_grid_eq_3D* ( (b_unit_grid(:,:,:,3)-b_unit_grid_eq_3D(:,:,:,3))*grad_p_h_perp0_3D(:,:,:,1) - &
                              (b_unit_grid(:,:,:,1)-b_unit_grid_eq_3D(:,:,:,1))*grad_p_h_perp0_3D(:,:,:,3) ) - &
              (B_grid-B_grid_eq_3D)*( b_unit_grid_eq_3D(:,:,:,3)*grad_p_h_perp0_3D(:,:,:,1) - &
                                      b_unit_grid_eq_3D(:,:,:,1)*grad_p_h_perp0_3D(:,:,:,3) ) )/( B_grid*B_grid_eq_3D )

            J_h(:,:,:,3) = J_h(:,:,:,3) + Zh*nv_para0_3D*(b_unit_grid(:,:,:,3)-b_unit_grid_eq_3D(:,:,:,3)) + &
            ( B_grid_eq_3D*(p_h_para0_3D-p_h_perp0_3D)*(curl_b_grid(:,:,:,3)-curl_b_grid_eq_3D(:,:,:,3)) - &
              (B_grid-B_grid_eq_3D)*(p_h_para0_3D-p_h_perp0_3D)*curl_b_grid_eq_3D(:,:,:,3) )/( B_grid*B_grid_eq_3D ) + &
            ( B_grid_eq_3D* ( (b_unit_grid(:,:,:,1)-b_unit_grid_eq_3D(:,:,:,1))*grad_p_h_perp0_3D(:,:,:,2) - &
                              (b_unit_grid(:,:,:,2)-b_unit_grid_eq_3D(:,:,:,2))*grad_p_h_perp0_3D(:,:,:,1) ) - &
              (B_grid-B_grid_eq_3D)*( b_unit_grid_eq_3D(:,:,:,1)*grad_p_h_perp0_3D(:,:,:,2) - &
                                      b_unit_grid_eq_3D(:,:,:,2)*grad_p_h_perp0_3D(:,:,:,1) ) )/( B_grid*B_grid_eq_3D )
        endif


        if(VEB) then
        if(.not.ADD_NH0) then
            J_h(:,:,:,1) = J_h(:,:,:,1) + m*(nv_para*curl_v_E_grid(:,:,:,1) + 0.5*n_h*(b_unit_grid(:,:,:,2)*grad_v_E_2_grid(:,:,:,3) - b_unit_grid(:,:,:,3)*grad_v_E_2_grid(:,:,:,2)) + nv_para*(b_unit_grid(:,:,:,2)*pb_unitpt_grid(:,:,:,3) - b_unit_grid(:,:,:,3)*pb_unitpt_grid(:,:,:,2)))/B_grid
            J_h(:,:,:,2) = J_h(:,:,:,2) + m*(nv_para*curl_v_E_grid(:,:,:,2) + 0.5*n_h*(b_unit_grid(:,:,:,3)*grad_v_E_2_grid(:,:,:,1) - b_unit_grid(:,:,:,1)*grad_v_E_2_grid(:,:,:,3)) + nv_para*(b_unit_grid(:,:,:,3)*pb_unitpt_grid(:,:,:,1) - b_unit_grid(:,:,:,1)*pb_unitpt_grid(:,:,:,3)))/B_grid
            J_h(:,:,:,3) = J_h(:,:,:,3) + m*(nv_para*curl_v_E_grid(:,:,:,3) + 0.5*n_h*(b_unit_grid(:,:,:,1)*grad_v_E_2_grid(:,:,:,2) - b_unit_grid(:,:,:,2)*grad_v_E_2_grid(:,:,:,1)) + nv_para*(b_unit_grid(:,:,:,1)*pb_unitpt_grid(:,:,:,2) - b_unit_grid(:,:,:,2)*pb_unitpt_grid(:,:,:,1)))/B_grid
 !polarization current
        else
            J_h(:,:,:,1) = J_h(:,:,:,1) + m*((nv_para+nv_para0_3D)*curl_v_E_grid(:,:,:,1) + 0.5*(n_h+n_h0_3D)*(b_unit_grid(:,:,:,2)*grad_v_E_2_grid(:,:,:,3) - b_unit_grid(:,:,:,3)*grad_v_E_2_grid(:,:,:,2)) + (nv_para+nv_para0_3D)*(b_unit_grid(:,:,:,2)*pb_unitpt_grid(:,:,:,3) - b_unit_grid(:,:,:,3)*pb_unitpt_grid(:,:,:,2)))/B_grid
            J_h(:,:,:,2) = J_h(:,:,:,2) + m*((nv_para+nv_para0_3D)*curl_v_E_grid(:,:,:,2) + 0.5*(n_h+n_h0_3D)*(b_unit_grid(:,:,:,3)*grad_v_E_2_grid(:,:,:,1) - b_unit_grid(:,:,:,1)*grad_v_E_2_grid(:,:,:,3)) + (nv_para+nv_para0_3D)*(b_unit_grid(:,:,:,3)*pb_unitpt_grid(:,:,:,1) - b_unit_grid(:,:,:,1)*pb_unitpt_grid(:,:,:,3)))/B_grid
            J_h(:,:,:,3) = J_h(:,:,:,3) + m*((nv_para+nv_para0_3D)*curl_v_E_grid(:,:,:,3) + 0.5*(n_h+n_h0_3D)*(b_unit_grid(:,:,:,1)*grad_v_E_2_grid(:,:,:,2) - b_unit_grid(:,:,:,2)*grad_v_E_2_grid(:,:,:,1)) + (nv_para+nv_para0_3D)*(b_unit_grid(:,:,:,1)*pb_unitpt_grid(:,:,:,2) - b_unit_grid(:,:,:,2)*pb_unitpt_grid(:,:,:,1)))/B_grid
        endif
            if(POL) then
            if(.not.ADD_NH0) then
                J_h(:,:,:,1) = J_h(:,:,:,1) + m*n_h*(b_unit_grid(:,:,:,2)*pv_Ept_grid(:,:,:,3) - b_unit_grid(:,:,:,3)*pv_Ept_grid(:,:,:,2))/B_grid
                J_h(:,:,:,2) = J_h(:,:,:,2) + m*n_h*(b_unit_grid(:,:,:,3)*pv_Ept_grid(:,:,:,1) - b_unit_grid(:,:,:,1)*pv_Ept_grid(:,:,:,3))/B_grid
                J_h(:,:,:,3) = J_h(:,:,:,3) + m*n_h*(b_unit_grid(:,:,:,1)*pv_Ept_grid(:,:,:,2) - b_unit_grid(:,:,:,2)*pv_Ept_grid(:,:,:,1))/B_grid
            else
                J_h(:,:,:,1) = J_h(:,:,:,1) + m*(n_h+n_h0_3D)*(b_unit_grid(:,:,:,2)*pv_Ept_grid(:,:,:,3) - b_unit_grid(:,:,:,3)*pv_Ept_grid(:,:,:,2))/B_grid
                J_h(:,:,:,2) = J_h(:,:,:,2) + m*(n_h+n_h0_3D)*(b_unit_grid(:,:,:,3)*pv_Ept_grid(:,:,:,1) - b_unit_grid(:,:,:,1)*pv_Ept_grid(:,:,:,3))/B_grid
                J_h(:,:,:,3) = J_h(:,:,:,3) + m*(n_h+n_h0_3D)*(b_unit_grid(:,:,:,1)*pv_Ept_grid(:,:,:,2) - b_unit_grid(:,:,:,2)*pv_Ept_grid(:,:,:,1))/B_grid
            end if
            end if
        end if

!for diagnostic
        M_h(:,:,:,1) = p_h_perp
        M_h(:,:,:,2) = p_h_para
        M_h(:,:,:,3) = nv_para
        J_M = grad_p_h_perp   

!3 times smooth!
!        do s=1,3
!            call bndry3_ex(J_h,0)
!            call smooth26(J_h, coeff_smooth, 3)
!        end do
        
    end subroutine calc_current

    subroutine calc_grad_B_initia
        use var
        implicit none
        include 'mpif.h'

        do i=1,mx
            do j=1,mz
                do k=1,my!note that field in CLT code is Bx(R,Z,phi)
                    B_grid(i,j,k)    = sqrt(dot_product(x(i,j,k,6:8),x(i,j,k,6:8)))
                end do
            end do
        end do 
        grad_B_grid(:,:,:,1) = xr(:,:,:,6)*x(:,:,:,6)/B_grid + xr(:,:,:,7)*x(:,:,:,7)/B_grid + xr(:,:,:,8)*x(:,:,:,8)/B_grid
        grad_B_grid(:,:,:,2) = xy(:,:,:,6)*x(:,:,:,6)/B_grid + xy(:,:,:,7)*x(:,:,:,7)/B_grid + xy(:,:,:,8)*x(:,:,:,8)/B_grid
        grad_B_grid(:,:,:,3) = xz(:,:,:,6)*x(:,:,:,6)/B_grid + xz(:,:,:,7)*x(:,:,:,7)/B_grid + xz(:,:,:,8)*x(:,:,:,8)/B_grid        
        !dBdphi/R
        do i=1,mx
            grad_B_grid(i,:,:,2) = grad_B_grid(i,:,:,2)/xx(i)
        end do
        
                
        b_unit_grid(:,:,:,1) = x(:,:,:,6)/B_grid
        b_unit_grid(:,:,:,2) = x(:,:,:,7)/B_grid
        b_unit_grid(:,:,:,3) = x(:,:,:,8)/B_grid
        
        do i=1,mx
            do j=1,mz
                do k=1,my!note that field in CLT code is Bx(R,Z,phi)
                    curl_b_grid(i,j,k,:)    = (cur(i,j,k,:)+cint(i,j,:))/B_grid(i,j,k) - cross_product(grad_B_grid(i,j,k,:),x(i,j,k,6:8))/(B_grid(i,j,k)**2)
                end do
            end do
        end do     

        
        do i=1,mx
            do j=1,mz
                B_grid_eq(i,j) = sqrt(dot_product(xint(i,j,6:8),xint(i,j,6:8)))
            end do
        end do 
        b_unit_grid_eq(:,:,1) = xint(:,:,6)/B_grid_eq
        b_unit_grid_eq(:,:,2) = xint(:,:,7)/B_grid_eq
        b_unit_grid_eq(:,:,3) = xint(:,:,8)/B_grid_eq
        
        grad_B_grid_eq(:,:,1) = xint_dx(:,:,6)*xint(:,:,6)/B_grid_eq + xint_dx(:,:,7)*xint(:,:,7)/B_grid_eq + xint_dx(:,:,8)*xint(:,:,8)/B_grid_eq
        grad_B_grid_eq(:,:,2) = 0
        grad_B_grid_eq(:,:,3) = xint_dz(:,:,6)*xint(:,:,6)/B_grid_eq + xint_dz(:,:,7)*xint(:,:,7)/B_grid_eq + xint_dz(:,:,8)*xint(:,:,8)/B_grid_eq
        
        do i=1,mx
            do j=1,mz
                curl_b_grid_eq(i,j,:)   = cint(i,j,:)/B_grid_eq(i,j) - cross_product(grad_B_grid_eq(i,j,:),xint(i,j,6:8))/(B_grid_eq(i,j)**2)
            end do
        end do      

        do j=1,my
                b_unit_grid_eq_3D(:,:,j,:) = b_unit_grid_eq
                curl_b_grid_eq_3D(:,:,j,:) = curl_b_grid_eq
                B_grid_eq_3D(:,:,j)        = B_grid_eq
        enddo
        
    end subroutine calc_grad_B_initia
    
    
    subroutine calc_grad_B
        use var
        implicit none
        include 'mpif.h'

        do i=1,mx
            do j=1,mz
                do k=1,my!note that field in CLT code is Bx(R,Z,phi)
                    B_grid(i,j,k)    = sqrt(dot_product(x(i,j,k,6:8),x(i,j,k,6:8)))
                end do
            end do
        end do 
        grad_B_grid(:,:,:,1) = xr(:,:,:,6)*x(:,:,:,6)/B_grid + xr(:,:,:,7)*x(:,:,:,7)/B_grid + xr(:,:,:,8)*x(:,:,:,8)/B_grid
        grad_B_grid(:,:,:,2) = xy(:,:,:,6)*x(:,:,:,6)/B_grid + xy(:,:,:,7)*x(:,:,:,7)/B_grid + xy(:,:,:,8)*x(:,:,:,8)/B_grid
        grad_B_grid(:,:,:,3) = xz(:,:,:,6)*x(:,:,:,6)/B_grid + xz(:,:,:,7)*x(:,:,:,7)/B_grid + xz(:,:,:,8)*x(:,:,:,8)/B_grid        
        !dBdphi/R
        do i=1,mx
            grad_B_grid(i,:,:,2) = grad_B_grid(i,:,:,2)/xx(i)
        end do
        
                
        b_unit_grid(:,:,:,1) = x(:,:,:,6)/B_grid
        b_unit_grid(:,:,:,2) = x(:,:,:,7)/B_grid
        b_unit_grid(:,:,:,3) = x(:,:,:,8)/B_grid
        
        do i=1,mx
            do j=1,mz
                do k=1,my!note that field in CLT code is Bx(R,Z,phi)
                    curl_b_grid(i,j,k,:)    = (cur(i,j,k,:)+cint(i,j,:))/B_grid(i,j,k) - cross_product(grad_B_grid(i,j,k,:),x(i,j,k,6:8))/(B_grid(i,j,k)**2)
                end do
            end do
        end do
        
!        pb_unitpt = (pBpt - pabs_Bpt*b_unit)/abs_B
        do i=1,mx
            do j=1,mz
                do k=1,my!note that field in CLT code is Bx(R,Z,phi)
                    pb_unitpt_grid(i,j,k,1:3) = (xdif(i,j,k,6:8) - dot_product(b_unit_grid(i,j,k,1:3),xdif(i,j,k,6:8))*b_unit_grid(i,j,k,1:3))/B_grid(i,j,k)
                end do
            end do
        end do 
        
    end subroutine calc_grad_B
  
    
!************************************************************************************************
!                                                                                               *
!                                diagnosis #0 particle                                          *
!                                                                                               *
!************************************************************************************************
    subroutine diagnosis_orbit
        use var
        use distribution 
        implicit none
        include 'mpif.h'
        
        integer :: STATUS(MPI_STATUS_SIZE)
        real*8,dimension(8) :: data_diag
        real*8,allocatable :: data_diag_total(:,:)
        logical :: flag_diag
        logical,allocatable :: flag_diag_total(:)
        allocate(data_diag_total(8,nsize))
        allocate(flag_diag_total(nsize))
        flag_diag = .false.
        do p=1,N
            if(marker(p)%id == 47*N_initia+1) then
!!update P_phi and v_perp
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
                            B      = B + x(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
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
                flag_diag    = .true.
            end if
        end do
        

        call MPI_GATHER(data_diag, 8, MPI_DOUBLE_PRECISION, data_diag_total(1:8,:), 8, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)
        call MPI_GATHER(flag_diag, 1, MPI_LOGICAL, flag_diag_total, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, IERROR)
        if(nrank == 0) then
            open(unit=32,file='particle_orbit.dat')
            do i=1,nsize
                if(flag_diag_total(i)) write(32,"(8(1x,e12.5),1x,i0)")data_diag_total(:,i),i-1
            end do            
        end if

    end subroutine diagnosis_orbit
    
    subroutine initia_CLTK
        use var
        use distribution
        implicit none
        include 'mpif.h'
        
!        N = 0
!        if(nrank == 47)    N = N_initia
        N = N_initia
        deltaN = floor(0.5*N) + N_buf_min
        allocate(marker(N_initia+deltaN))
        allocate(marker_initia(N_initia+deltaN))
        psmin_clt_k = psmin
        psmax_clt_k = psmin + 0.9*(psmax - psmin)
!        Delta_psi = psmax_clt_k - psmin_clt_k
        Delta_psi = psmax - psmin  !Delta_psi should be same with MHD part!!!
        R0        = xzero
        myxmin    = xx(3) 
        if(nrkx(nrank) == nprx-1) then
            myxmax    = xx(mx-2)
        else
            myxmax    = xx(mx-1)
        end if
        myymin    = yy(3) 
        myymax    = yy(my-1)
        myzmin    = zz(3) 
        if(nrkz(nrank) == nprz-1) then
            myzmax    = zz(mz-2)
        else
            myzmax    = zz(mz-1)
        end if
        J_h = 0
        n_h = 0
        write(*,*)"R0 :",R0,"myzmin",myzmin,"myzmax",myzmax,"myxmin",myxmin,"myxmax",myxmax," nrank:",nrank
    
    end subroutine initia_CLTK
    
    
    subroutine initia_parameter_output
        use DECLARE 
        use var
        use distribution
        use diganostic
        use strings
        implicit none
    
        if(nrank == 0) then
            open(unit=31,file='parameter_out',RECL=256)
            write(31,*)"number of particle               : ",num2str(N_total)
            write(31,*)"volume average fast ion beta     : ",num2str(beta)           
            write(31,*)"coefficient of smooth            : ",num2str(coeff_smooth)
            write(31,*)"viscosity                        : ",num2str(fmu0)
            write(31,*)"diffusion                        : ",num2str(pmu0)
            write(31,*)"resistivity                      : ",num2str(eta0)
            if(RECYCLE_METHOD == 1) then
            write(31,*)"treat bounary particles:         : recycle"
            else if(RECYCLE_METHOD == 2) then
            write(31,*)"treat bounary particles:         : remove" 
            end if
            if(ORBIT_AVERAGE_METHOD == 1) then
            write(31,*)"Orbit averaging method:          : <psi> = -P_phi/e + m/e*sgn(v_||)v*R0*sqrt(1-mu*B0/E) for passing and -P_phi/e for trapped particles"
            else if(ORBIT_AVERAGE_METHOD == 2) then
            write(31,*)"Orbit averaging method:          : <psi> = -P_phi/e for both passing and -P_phi/e for trapped particles"
            end if
            if(VEB) then
            write(31,*)"Lagrangian including E*B flow    : A_star = A + m/e*v_para*b + m/e*v_E"
            else
            write(31,*)"Lagrangian without E*B flow      : A_star = A + m/e*v_para*b"
            end if
            if(POL) then
            write(31,*)"Including polarization drift     : pv_Ept"
            else
            write(31,*)"Including polarization drift     : no"
            end if
            if(BAS) then 
            write(31,*)"Including banos drift            : mu/2/e*b*curl_b"
            else
            write(31,*)"Including banos drift            : no"
            end if
            if(FLT) then
            write(31,*)"single mode                      : n=",num2str(n_filt)
            if(NMC) then
            write(31,*)"Nonlinear MHD coupling           : Yes"
            else
            write(31,*)"Nonlinear MHD coupling           : No"
            end if
            else
            write(31,*)"multiple mode                    : no filter mode"
            end if
            write(31,*)"coefficient of possion smooth    : ",num2str(c_possion)
            if(DFR_type == 1) then
            write(31,*)"distribution in real space       : exp(-<psi>/L_n)"
            else if(DFR_type == 2) then
            write(31,*)"distribution in real space       : exp(-Delta_n/L_n*tanh((<psi>-psi_n)/Delta_n)) with psi_n = ",num2str(psi_n),", L_n = ",num2str(L_n),", Delta_n = ",num2str(Delta_n)
            end if
            if(IIV) then
            write(31,*)"distrubution in velocity         : isotropy"
            else if(FGP) then
            write(31,*)"distrubution in velocity         : Finite width of Gaussian distrbution of Pitch angle:",num2str(Lambda_0)
            else
            write(31,*)"distrubution in velocity         : anisotropy, mu=0"
            end if
            if(LINEAR_MHD) then
            write(31,*)"MHD reponse                      : linear"
            else
            write(31,*)"MHD reponse                      : nonlinear"
            end if
            if(EQ_WITH_EP) then
            write(31,*)"Equilibrium                      : with fast ions"
            else
            write(31,*)"Equilibrium                      : without fast ions"
            end if
            if(FOUR_TH) then
            write(31,*)"Accuracy of   finite difference  : 4th"
            else
            write(31,*)"Accuracy of   finite difference  : 2ed"
            end if
            if(FHW) then
            write(31,*)"Toridal mode number limitation   : n<",n_limit
            end if
            if(DEP) then
            write(31,*)"--------------parameter for diagnostics of distribution function--------------"
		    write(31,*)"grid on P_phi                                                   :   ",num2str(GRID_P_PHI)
            write(31,*)"grid on Energy                                                  :   ",num2str(GRID_E)
            write(31,*)"pitch angle variable Lambda0=mu*B0/E for diagnostic             :   ",num2str(Lambda_DEP,5)
            write(31,*)"how many time step to average distribution function data        :   ",num2str(NSTEP_AVG)
            write(31,*)"how many time step to record distribution function data         :   ",num2str(NSTEP_INT)
            write(31,*)"how many time step to  start recording data                     :   ",num2str(NSTEP_START)
            write(31,*)"size of diagnostic windows for fixed Lambda and E               :   ",num2str(dLambda_DEP)
            write(31,*)"output distribution function as function of                     :   ","f(P_phi,E,Lambda_0)"
            end if
            if(DEB) then
            write(31,*)"--------------parameter for diagnostics of distribution function--------------"
		    write(31,*)"grid on <psi>                                                   :   ",num2str(GRID_P_PHI)
            write(31,*)"grid on Energy                                                  :   ",num2str(GRID_E)
            write(31,*)"pitch angle variable Lambda0=mu*B0/E for diagnostic             :   ",num2str(Lambda_DEP,5)
            write(31,*)"how many time step to average distribution function data        :   ",num2str(NSTEP_AVG)
            write(31,*)"how many time step to record distribution function data         :   ",num2str(NSTEP_INT)
            write(31,*)"how many time step to  start recording data                     :   ",num2str(NSTEP_START)
            write(31,*)"size of diagnostic windows for fixed Lambda and E               :   ",num2str(dLambda_DEP)
            if(DEB_type == -1) then
            write(31,*)"output distribution function as function of                     :   ","f(<psi>,E,Lambda_0) with counter-going particles "
            else if(DEB_type == 0) then
            write(31,*)"output distribution function as function of                     :   ","f(<psi>,E,Lambda_0)"
            else if(DEB_type == +1) then
            write(31,*)"output distribution function as function of                     :   ","f(<psi>,E,Lambda_0) with co-going particles "
            end if
            end if
        end if
        
    end subroutine initia_parameter_output
    
    subroutine update_J_h_bounnary
        use var
        use distribution
        implicit none
        include 'mpif.h'
        integer :: up_id, down_id, left_id, right_id, front_id, back_id
        integer :: up_left_id, up_right_id, down_left_id, down_right_id
        integer :: left_front_id, left_back_id, right_front_id, right_back_id
        integer :: up_front_id, up_back_id, down_front_id, down_back_id
        integer :: dest_nrkz_up, dest_nrkz_down, dest_nrkx_left, dest_nrkx_right, dest_nrky_front, dest_nrky_back
        real*8, dimension(mx,my,3) :: J_h_send_down, J_h_send_up,    J_h_recv_down, J_h_recv_up
        real*8, dimension(mz,my,3) :: J_h_send_left, J_h_send_right, J_h_recv_left, J_h_recv_right
        real*8, dimension(mx,mz,3) :: J_h_send_back, J_h_send_front, J_h_recv_back, J_h_recv_front
        real*8, dimension(my,3) :: J_h_send_up_left, J_h_send_up_right, J_h_send_down_left, J_h_send_down_right
        real*8, dimension(my,3) :: J_h_recv_up_left, J_h_recv_up_right, J_h_recv_down_left, J_h_recv_down_right
        real*8, dimension(mz,3) :: J_h_send_left_front, J_h_send_left_back, J_h_send_right_front, J_h_send_right_back
        real*8, dimension(mz,3) :: J_h_recv_left_front, J_h_recv_left_back, J_h_recv_right_front, J_h_recv_right_back
        real*8, dimension(mx,3) :: J_h_send_up_front, J_h_send_up_back, J_h_send_down_front, J_h_send_down_back
        real*8, dimension(mx,3) :: J_h_recv_up_front, J_h_recv_up_back, J_h_recv_down_front, J_h_recv_down_back
        
        integer, dimension(2)       :: dest_rl, dest_ud, dest_fb
        integer, dimension(2,2,2)   :: dest
        real*8,  dimension(2,2,2,3) :: J_h_send_rl_ud_fb, J_h_recv_rl_ud_fb
        
        integer :: STATUS(MPI_STATUS_SIZE)
        
        dest_nrkx_left  = nrkx(nrank) - 1
        dest_nrkx_right = nrkx(nrank) + 1
        dest_nrky_back  = nrky(nrank) - 1
        dest_nrky_front = nrky(nrank) + 1
        dest_nrkz_down  = nrkz(nrank) - 1
        dest_nrkz_up    = nrkz(nrank) + 1
        
        if(dest_nrkx_right > nprx - 1) then 
            right_id = -1
        else
            right_id = nrky(nrank)*nprx*nprz + nrkz(nrank)*nprx + dest_nrkx_right
        endif            
        if(dest_nrkx_left < 0) then
            left_id  = -1
        else
            left_id = nrky(nrank)*nprx*nprz + nrkz(nrank)*nprx + dest_nrkx_left
        endif 
        if(dest_nrkz_up > nprz - 1) then
            up_id    = -1
        else
            up_id = nrky(nrank)*nprx*nprz + dest_nrkz_up*nprx + nrkx(nrank)
        endif 
        if(dest_nrkz_down < 0) then
            down_id  = -1
        else
            down_id = nrky(nrank)*nprx*nprz + dest_nrkz_down*nprx + nrkx(nrank)
        endif 
        
        !boundary condition
        if(dest_nrky_front == npry) dest_nrky_front = 0
        if(dest_nrky_back == -1)    dest_nrky_back = npry-1
        
        front_id = dest_nrky_front*nprx*nprz + nrkz(nrank)*nprx + nrkx(nrank)
        back_id  = dest_nrky_back*nprx*nprz + nrkz(nrank)*nprx + nrkx(nrank)
     
        
        J_h_send_down  = J_h(:,3,:,:)
        J_h_send_up    = J_h(:,mz-1,:,:)
        J_h_send_left  = J_h(3,:,:,:)
        J_h_send_right = J_h(mx-1,:,:,:)
        J_h_send_back  = J_h(:,:,3,:)
        J_h_send_front = J_h(:,:,my-1,:) 
        if(down_id  /= -1) call MPI_SENDRECV(J_h_send_down(1:mx,1:my,1:3),  mx*my*3, MPI_DOUBLE_PRECISION, down_id,  0, J_h_recv_down(1:mx,1:my,1:3),  mx*my*3, MPI_DOUBLE_PRECISION, down_id,  0, MPI_COMM_WORLD, STATUS, IERROR)
        if(up_id    /= -1) call MPI_SENDRECV(J_h_send_up(1:mx,1:my,1:3),    mx*my*3, MPI_DOUBLE_PRECISION, up_id,    0, J_h_recv_up(1:mx,1:my,1:3),    mx*my*3, MPI_DOUBLE_PRECISION, up_id,    0, MPI_COMM_WORLD, STATUS, IERROR)
        if(left_id  /= -1) call MPI_SENDRECV(J_h_send_left(1:mz,1:my,1:3),  my*mz*3, MPI_DOUBLE_PRECISION, left_id,  0, J_h_recv_left(1:mz,1:my,1:3),  my*mz*3, MPI_DOUBLE_PRECISION, left_id,  0, MPI_COMM_WORLD, STATUS, IERROR)
        if(right_id /= -1) call MPI_SENDRECV(J_h_send_right(1:mz,1:my,1:3), my*mz*3, MPI_DOUBLE_PRECISION, right_id, 0, J_h_recv_right(1:mz,1:my,1:3), my*mz*3, MPI_DOUBLE_PRECISION, right_id, 0, MPI_COMM_WORLD, STATUS, IERROR)
        if(front_id /= back_id) then
            if(front_id /= -1) call MPI_SENDRECV(J_h_send_front(1:mx,1:mz,1:3), mx*mz*3, MPI_DOUBLE_PRECISION, front_id, 0, J_h_recv_front(1:mx,1:mz,1:3), mx*mz*3, MPI_DOUBLE_PRECISION, front_id, 0, MPI_COMM_WORLD, STATUS, IERROR)
            if(back_id  /= -1) call MPI_SENDRECV(J_h_send_back(1:mx,1:mz,1:3),  mx*mz*3, MPI_DOUBLE_PRECISION, back_id,  0, J_h_recv_back(1:mx,1:mz,1:3),  mx*mz*3, MPI_DOUBLE_PRECISION, back_id,  0, MPI_COMM_WORLD, STATUS, IERROR)
        else
            J_h_recv_front = J_h_send_back 
            J_h_recv_back  = J_h_send_front
        endif
!        if(nrank == 11 ) write(*,*)"kaiakia",J_h_recv_down(10,3,1),down_id
!        if(nrank == 11 ) write(*,*)"yauayua",J_h_send_down(10,3,1),down_id  


!        if(nrank == 11 ) write(*,*)"wawawaw",J_h(10,3,3,1),down_id       
        
        if(up_id /= -1 .and. left_id /= -1) then 
            up_left_id =  up_id - 1
        else
            up_left_id  = -1
        end if
        
        if(up_id /= -1 .and. right_id /= -1) then 
            up_right_id =  up_id + 1
        else
            up_right_id  = -1
        end if
        
        if(down_id /= -1 .and. left_id /= -1) then 
            down_left_id =  down_id - 1
        else
            down_left_id  = -1
        end if
        
        if(down_id /= -1 .and. right_id /= -1) then 
            down_right_id =  down_id + 1
        else
            down_right_id = -1
        end if
        
        J_h_send_up_left    = J_h(3,mz-1,:,:)
        J_h_send_up_right   = J_h(mx-1,mz-1,:,:)
        J_h_send_down_left  = J_h(3,3,:,:)
        J_h_send_down_right = J_h(mx-1,3,:,:)
 
        if(up_left_id    /= -1) call MPI_SENDRECV(J_h_send_up_left,    my*3, MPI_DOUBLE_PRECISION, up_left_id,    0, J_h_recv_up_left,    my*3, MPI_DOUBLE_PRECISION, up_left_id,    0, MPI_COMM_WORLD, STATUS, IERROR)
        if(down_right_id /= -1) call MPI_SENDRECV(J_h_send_down_right, my*3, MPI_DOUBLE_PRECISION, down_right_id, 0, J_h_recv_down_right, my*3, MPI_DOUBLE_PRECISION, down_right_id, 0, MPI_COMM_WORLD, STATUS, IERROR)
        if(up_right_id   /= -1) call MPI_SENDRECV(J_h_send_up_right,   my*3, MPI_DOUBLE_PRECISION, up_right_id,   0, J_h_recv_up_right,   my*3, MPI_DOUBLE_PRECISION, up_right_id,   0, MPI_COMM_WORLD, STATUS, IERROR)
        if(down_left_id  /= -1) call MPI_SENDRECV(J_h_send_down_left,  my*3, MPI_DOUBLE_PRECISION, down_left_id,  0, J_h_recv_down_left,  my*3, MPI_DOUBLE_PRECISION, down_left_id,  0, MPI_COMM_WORLD, STATUS, IERROR)

        
        
        if(up_id /= -1 .and. front_id /= -1) then 
            up_front_id =  front_id + nprx
        else
            up_front_id  = -1
        end if
        
        if(up_id /= -1 .and. back_id /= -1) then 
            up_back_id =  back_id + nprx
        else
            up_back_id  = -1
        end if
        
        if(down_id /= -1 .and. front_id /= -1) then 
            down_front_id =  front_id - nprx
        else
            down_front_id = -1
        end if
        
        if(down_id /= -1 .and. back_id /= -1) then 
            down_back_id =  back_id - nprx
        else
            down_back_id = -1
        end if
        
        J_h_send_up_front   = J_h(:,mz-1,my-1,:)
        J_h_send_up_back    = J_h(:,mz-1,3,:)
        J_h_send_down_front = J_h(:,3,my-1,:)
        J_h_send_down_back  = J_h(:,3,3,:)
 
        if(up_front_id   /= -1) call MPI_SENDRECV(J_h_send_up_front,   mx*3, MPI_DOUBLE_PRECISION, up_front_id,   0, J_h_recv_up_front,   mx*3, MPI_DOUBLE_PRECISION, up_front_id,   0, MPI_COMM_WORLD, STATUS, IERROR)
        if(up_back_id    /= -1) call MPI_SENDRECV(J_h_send_up_back,    mx*3, MPI_DOUBLE_PRECISION, up_back_id,    0, J_h_recv_up_back,    mx*3, MPI_DOUBLE_PRECISION, up_back_id,    0, MPI_COMM_WORLD, STATUS, IERROR)
        if(down_front_id /= -1) call MPI_SENDRECV(J_h_send_down_front, mx*3, MPI_DOUBLE_PRECISION, down_front_id, 0, J_h_recv_down_front, mx*3, MPI_DOUBLE_PRECISION, down_front_id, 0, MPI_COMM_WORLD, STATUS, IERROR)
        if(down_back_id  /= -1) call MPI_SENDRECV(J_h_send_down_back,  mx*3, MPI_DOUBLE_PRECISION, down_back_id,  0, J_h_recv_down_back,  mx*3, MPI_DOUBLE_PRECISION, down_back_id,  0, MPI_COMM_WORLD, STATUS, IERROR)

        
        if(right_id /= -1 .and. front_id /= -1) then 
            right_front_id =  front_id + 1
        else
            right_front_id = -1
        end if
        
        if(right_id /= -1 .and. back_id /= -1) then 
            right_back_id =  back_id + 1
        else
            right_back_id  = -1
        end if
        
        if(left_id /= -1 .and. front_id /= -1) then 
            left_front_id =  front_id - 1
        else
            left_front_id = -1
        end if
        
        if(left_id /= -1 .and. back_id /= -1) then 
            left_back_id =  back_id - 1
        else
            left_back_id = -1
        end if
        
        J_h_send_right_front = J_h(mx-1,:,my-1,:)
        J_h_send_right_back  = J_h(mx-1,:,3,:)
        J_h_send_left_front  = J_h(3,:,my-1,:)
        J_h_send_left_back   = J_h(3,:,3,:)
 
        if(right_front_id /= -1) call MPI_SENDRECV(J_h_send_right_front, mz*3, MPI_DOUBLE_PRECISION, right_front_id, 0, J_h_recv_right_front, mz*3, MPI_DOUBLE_PRECISION, right_front_id, 0, MPI_COMM_WORLD, STATUS, IERROR)
        if(right_back_id  /= -1) call MPI_SENDRECV(J_h_send_right_back,  mz*3, MPI_DOUBLE_PRECISION, right_back_id,  0, J_h_recv_right_back,  mz*3, MPI_DOUBLE_PRECISION, right_back_id,  0, MPI_COMM_WORLD, STATUS, IERROR)
        if(left_front_id  /= -1) call MPI_SENDRECV(J_h_send_left_front,  mz*3, MPI_DOUBLE_PRECISION, left_front_id,  0, J_h_recv_left_front,  mz*3, MPI_DOUBLE_PRECISION, left_front_id,  0, MPI_COMM_WORLD, STATUS, IERROR)
        if(left_back_id   /= -1) call MPI_SENDRECV(J_h_send_left_back,   mz*3, MPI_DOUBLE_PRECISION, left_back_id,   0, J_h_recv_left_back,   mz*3, MPI_DOUBLE_PRECISION, left_back_id,   0, MPI_COMM_WORLD, STATUS, IERROR)

        
        
        dest_rl(1) = left_id
        dest_rl(2) = right_id
        dest_ud(1) = down_id
        dest_ud(2) = up_id
        dest_fb(1) = back_id
        dest_fb(2) = front_id
        do i=1,2
            do j=1,2
                do k=1,2
                    if(dest_rl(i) /= -1 .and. dest_ud(j) /= -1 .and. dest_fb(k) /= -1) then
                        dest(i,j,k) = dest_fb(k) + (j*2-3)*nprx + (i*2-3)
                        J_h_send_rl_ud_fb(i,j,k,:) = J_h(3+(i-1)*(mx-1-3),3+(j-1)*(mz-1-3),3+(k-1)*(my-1-3),:)
                        call MPI_SENDRECV(J_h_send_rl_ud_fb(i,j,k,:), 3, MPI_DOUBLE_PRECISION, dest(i,j,k), 0, J_h_recv_rl_ud_fb(i,j,k,:), 3, MPI_DOUBLE_PRECISION, dest(i,j,k), 0, MPI_COMM_WORLD, STATUS, IERROR)
                    end if
                end do
            end do
        end do
        
        

        
        
        
        J_h(:,3,:,:)    = J_h(:,3,:,:)    + J_h_recv_down
        J_h(:,mz-1,:,:) = J_h(:,mz-1,:,:) + J_h_recv_up

        J_h(3,:,:,:)    = J_h(3,:,:,:)    + J_h_recv_left
        J_h(mx-1,:,:,:) = J_h(mx-1,:,:,:) + J_h_recv_right
        
        J_h(:,:,3,:)    = J_h(:,:,3,:)    + J_h_recv_back
        J_h(:,:,my-1,:) = J_h(:,:,my-1,:) + J_h_recv_front

        
        J_h(3,mz-1,:,:)    = J_h(3,mz-1,:,:)    + J_h_recv_up_left
        J_h(mx-1,mz-1,:,:) = J_h(mx-1,mz-1,:,:) + J_h_recv_up_right
        J_h(3,3,:,:)       = J_h(3,3,:,:)       + J_h_recv_down_left
        J_h(mx-1,3,:,:)    = J_h(mx-1,3,:,:)    + J_h_recv_down_right
        
        J_h(:,mz-1,my-1,:) = J_h(:,mz-1,my-1,:) + J_h_recv_up_front 
        J_h(:,mz-1,3,:)    = J_h(:,mz-1,3,:)    + J_h_recv_up_back
        J_h(:,3,my-1,:)    = J_h(:,3,my-1,:)    + J_h_recv_down_front 
        J_h(:,3,3,:)       = J_h(:,3,3,:)       + J_h_recv_down_back 
        
        J_h(mx-1,:,my-1,:) = J_h(mx-1,:,my-1,:) + J_h_recv_right_front 
        J_h(mx-1,:,3,:)    = J_h(mx-1,:,3,:)    + J_h_recv_right_back
        J_h(3,:,my-1,:)    = J_h(3,:,my-1,:)    + J_h_recv_left_front 
        J_h(3,:,3,:)       = J_h(3,:,3,:)       + J_h_recv_left_back 
        
        do i=1,2
            do j=1,2
                do k=1,2
                    J_h(3+(i-1)*(mx-1-3),3+(j-1)*(mz-1-3),3+(k-1)*(my-1-3),:)  = J_h(3+(i-1)*(mx-1-3),3+(j-1)*(mz-1-3),3+(k-1)*(my-1-3),:) + J_h_recv_rl_ud_fb(i,j,k,:)
                end do
            end do
        end do
!        J_h(:,3,:,:)    = J_h_send_down  + J_h_recv_down
!        J_h(:,mz-1,:,:) = J_h_send_up    + J_h_recv_up
!        J_h(3,:,:,:)    = J_h_send_left  + J_h_recv_left
!        J_h(mx-1,:,:,:) = J_h_send_right + J_h_recv_right
!        J_h(:,:,3,:)    = J_h_send_back  + J_h_recv_back
!        J_h(:,:,my-1,:) = J_h_send_front + J_h_recv_front

    end subroutine update_J_h_bounnary
    
    
    
    subroutine diagnosis_field
        use var
        use distribution 
        implicit none
        include 'mpif.h'
        
        if(nrank == 21) then
            open(unit=33,file='field_delta_B_R.dat')
            write(33,"(2(1x,e12.5))")x(mx/2,mz/2,1,6)-xint(mx/2,mz/2,6),time
        end if
        
    end subroutine diagnosis_field


    
    subroutine init_random_seed(myid)

        integer :: i, n, clock
        integer :: myid
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)

        allocate(seed(n))

        call system_clock(count=clock)

        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        seed = seed*(myid+1)

        call random_seed(put = seed)

        deallocate(seed)

    end subroutine init_random_seed
    
!************************************************************************************************
!                                                                                               *
!                               7 points smooth method                                          *
!                                                                                               *
!************************************************************************************************      
!coeff_smooth : large coefficient means heavily smooth
    subroutine smooth(data_in, coefficient_smooth)
        use var
        implicit none
        include 'mpif.h'
        
        real*8, dimension(mx,mz,my,3), intent(inout) :: data_in
        real*8, intent(in) :: coefficient_smooth
        
        real*8 :: coeff_left, coeff_right, coeff_back, coeff_front, coeff_down, coeff_up
        real*8 :: coeff_total        
        real*8, dimension(mx,mz,my,3) :: data_smooth     

        call mpi_transfersm(data_in,3)
        
        !points on boundary should not be smoothed!!
        data_smooth = data_in
        
        do i=3,mx-2
            coeff_total = 2.0/dxx + 2.0/(xx(i)*dyy) + 2.0/dzz
            coeff_left  = (1.0/dxx)/coeff_total
            coeff_right = coeff_left
            coeff_back  = (1.0/(xx(i)*dyy))/coeff_total
            coeff_front = coeff_back
            coeff_down  = (1.0/dzz)/coeff_total
            coeff_up    = coeff_down
            do j=3,mz-2
                do k=3,my-2
                    if(psi(i,j)<psmax .and. psi(i+1,j)<psmax .and. psi(i,j+1)<psmax .and. psi(i+1,j+1)<psmax .and. psi(i-1,j)<psmax .and. psi(i,j-1)<psmax .and. psi(i-1,j-1)<psmax) then
                        data_smooth(i,j,k,:) = (1.0 - coefficient_smooth)*data_in(i,j,k,:) &
                                            + coefficient_smooth*(coeff_left*data_in(i-1,j,k,:) + coeff_right*data_in(i+1,j,k,:) &
                                            +               coeff_back*data_in(i,j,k-1,:) + coeff_front*data_in(i,j,k+1,:) &
                                            +               coeff_down*data_in(i,j-1,k,:) +    coeff_up*data_in(i,j+1,k,:))
                    end if
                end do
            end do
        end do
        
        data_in = data_smooth
        
        call mpi_transfersm(data_in,3)
        
    end subroutine smooth
    
!************************************************************************************************
!                                                                                               *
!                               26 points smooth method                                         *
!                                                                                               *
!************************************************************************************************      
!coeff_smooth : large coefficient means heavily smooth
    subroutine smooth26(data_in, coefficient_smooth, dim)
        use var
        implicit none
        include 'mpif.h'
        
        integer, intent(in) :: dim
        real*8, dimension(mx,mz,my,dim), intent(inout) :: data_in
        real*8, intent(in) :: coefficient_smooth
        
        real*8, dimension(3,3,3) :: weight
        real*8 :: distance, weight_total        
        real*8, dimension(mx,mz,my,dim) :: data_smooth

        call mpi_transfersm(data_in,dim)
        
        !points on boundary should not be smoothed!!
        data_smooth = data_in
        
        do i=3,mx-2
            weight_total = 0.0
            do ii=1,3
                do jj=1,3
                    do kk=1,3
                        if(ii /= 2 .or. jj /=2 .or. kk /= 2) then
                            distance = sqrt(abs(ii-2)*dxx**2 + abs(jj-2)*dzz**2 + abs(kk-2)*(xx(i)*dyy)**2)
                            weight(ii,jj,kk) = 1.0/distance
                            weight_total = weight_total + weight(ii,jj,kk)
                        end if
                    end do
                end do
            end do
            
            weight = weight/weight_total*coefficient_smooth
            weight(2,2,2) = (1.0 - coefficient_smooth)
            
            do j=3,mz-2
                do k=3,my-2
                    if(psi(i,j)<psmax .and. psi(i+1,j)<psmax .and. psi(i,j+1)<psmax .and. psi(i+1,j+1)<psmax .and. psi(i-1,j)<psmax .and. psi(i,j-1)<psmax .and. psi(i-1,j-1)<psmax .and. psi(i+1,j-1)<psmax .and. psi(i-1,j+1)<psmax) then
                        data_smooth(i,j,k,1:dim) = 0.0
                        do ii=1,3
                            do jj=1,3
                                do kk=1,3
                                    data_smooth(i,j,k,1:dim) = data_smooth(i,j,k,1:dim) + weight(ii,jj,kk)*data_in(i+ii-2,j+jj-2,k+kk-2,1:dim)
                                end do
                            end do
                        end do
                    end if
                end do
            end do
        end do
        
        data_in = data_smooth
        
        call mpi_transfersm(data_in,dim)
        
    end subroutine smooth26

    
!************************************************************************************************
!                                                                                               *
!                               26 points smooth method for eq field                            *
!                                                                                               *
!************************************************************************************************      
!coeff_smooth : large coefficient means heavily smooth
    subroutine smooth26_eq(data_in, coefficient_smooth)
        use var
        implicit none
        include 'mpif.h'
        
        real*8, dimension(mx,mz,my,3), intent(inout) :: data_in
        real*8, intent(in) :: coefficient_smooth
        
        real*8, dimension(3,3) :: weight
        real*8 :: distance, weight_total        
        real*8, dimension(mx,mz,my,3) :: data_smooth

        call mpi_transfersm(data_in,3)
        
        !points on boundary should not be smoothed!!
        data_smooth = data_in
        
        do i=3,mx-2
            weight_total = 0.0
            do ii=1,3
                do jj=1,3
                    if(ii /= 2 .or. jj /=2) then
                        distance = sqrt(abs(ii-2)*dxx**2 + abs(jj-2)*dzz**2)
                        weight(ii,jj) = 1.0/distance
                        weight_total = weight_total + weight(ii,jj)
                    end if
                end do
            end do
            
            weight = weight/weight_total*coefficient_smooth
            weight(2,2) = (1.0 - coefficient_smooth)
            
!do smooth in R,Z plane!
            
            do j=3,mz-2
                do k=3,my-2
                    if(psi(i,j)<psmax .and. psi(i+1,j)<psmax .and. psi(i,j+1)<psmax .and. psi(i+1,j+1)<psmax .and. psi(i-1,j)<psmax .and. psi(i,j-1)<psmax .and. psi(i-1,j-1)<psmax .and. psi(i+1,j-1)<psmax .and. psi(i-1,j+1)<psmax) then
                        data_smooth(i,j,k,:) = 0.0
                        do ii=1,3
                            do jj=1,3
                               data_smooth(i,j,k,:) = data_smooth(i,j,k,:) + weight(ii,jj)*data_in(i+ii-2,j+jj-2,k,:)
                            end do
                        end do
                    end if
                end do
            end do
        end do

!do smooth in phi direction!!
        do k=3+1,my-2
            data_smooth(:,:,3,:) = data_smooth(:,:,3,:) + data_smooth(:,:,k,:)
        end do
        data_smooth(:,:,3,:) = data_smooth(:,:,3,:)/(my-4)

        do k=3,my-2
            data_in(:,:,k,:) = data_smooth(:,:,3,:)
        end do
        
        call mpi_transfersm(data_in,3)
        
    end subroutine smooth26_eq
    !!!clt-k!!!
!************************************************************************************************
!                                                                                               *
!                         filter wave on toridal direction                                      *
!                                                                                               *
!************************************************************************************************    
    subroutine filter_wave(data_in, mm, mode_filt)
        use DECLARE
        implicit none
        include 'mpif.h'
        include 'fftw_f77.i'
        integer status(mpi_status_size)
        
        integer, intent(in) :: mm, mode_filt
        real*8, dimension(mx,mz,my,mm) :: data_in
        real*8, dimension(myt/2+1) :: coeff_mode
        integer*8 :: plan

        coeff_mode = 0

        coeff_mode(mode_filt+1) = 1 !n=mode_filt
!        coeff_mode(1) = 1 !keep n=0

        
        do jx=ix_first,ix_last
            do jz=iz_first,iz_last
                do m=1,mm
                    call mpi_transfersy1(data_in(jx,jz,:,m),data0)
                    
!                    if(nrank == 0 .and. m ==1) then
!                        do jy=1,myt
!                            data0(jy) = jy
!                            write(*,*)"jiyjjyjuyjuy : ",jy
!                        end do
!                    end if
                    
                    
                    call dfftw_plan_dft_r2c_1d(plan,myt,data0,spec,FFTW_ESTIMATE)
                    call dfftw_execute_dft_r2c(plan,data0,spec)
                    call dfftw_destroy_plan(plan)
                    
                    spec=spec*coeff_mode
                    
                    spec=spec*fac
                    
                    !spec(myt/2+1)=c0
                    
                    call dfftw_plan_dft_c2r_1d(plan,myt,spec,data0,FFTW_ESTIMATE)
                    call dfftw_execute_dft_c2r(plan,spec,data0)
                    call dfftw_destroy_plan(plan)
                    
!                    if(nrank == 0 .and. m ==1) then
!                        do jy=1,myt
!                            write(*,*)"joyojojoyojouoyo : ",data0(jy)
!                        end do
!                    end if
                    
                    do jy=iy_first+2,iy_last-2
                        data_in(jx,jz,jy,m) = data0(nrky(nrank)*mym+jy-2)
                    end do
                    
                end do
            end do
        end do
        
        call mpi_transfersm(data_in,mm)
        
        return
    end subroutine filter_wave
 
    !!!clt-k!!!
!************************************************************************************************
!                                                                                               *
!                         filter high-n wave on toridal direction                               *
!                                                                                               *
!************************************************************************************************  
    subroutine filter_high_mode(data_in, mm, mode_limit)
        use DECLARE
        implicit none
        include 'mpif.h'
        include 'fftw_f77.i'
        integer status(mpi_status_size)
        
        integer, intent(in) :: mm, mode_limit
        real*8, dimension(mx,mz,my,mm) :: data_in
        real*8, dimension(myt/2+1) :: coeff_mode
        integer*8 :: plan
        integer :: i

!       filter high n (n>4)
        coeff_mode = 1
        do i=mode_limit,myt/2
            coeff_mode(i+1) = 0
        end do
            
        
        
        do jx=ix_first,ix_last
            do jz=iz_first,iz_last
                do m=1,mm
                    call mpi_transfersy1(data_in(jx,jz,:,m),data0)
                    
!                    if(nrank == 0 .and. m ==1) then
!                        do jy=1,myt
!                            data0(jy) = jy
!                            write(*,*)"jiyjjyjuyjuy : ",jy
!                        end do
!                    end if
                    
                    
                    call dfftw_plan_dft_r2c_1d(plan,myt,data0,spec,FFTW_ESTIMATE)
                    call dfftw_execute_dft_r2c(plan,data0,spec)
                    call dfftw_destroy_plan(plan)
                    
                    spec=spec*coeff_mode
                    
                    spec=spec*fac
                    
                    !spec(myt/2+1)=c0
                    
                    call dfftw_plan_dft_c2r_1d(plan,myt,spec,data0,FFTW_ESTIMATE)
                    call dfftw_execute_dft_c2r(plan,spec,data0)
                    call dfftw_destroy_plan(plan)
                    
!                    if(nrank == 0 .and. m ==1) then
!                        do jy=1,myt
!                            write(*,*)"joyojojoyojouoyo : ",data0(jy)
!                        end do
!                    end if
                    
                    do jy=iy_first+2,iy_last-2
                        data_in(jx,jz,jy,m) = data0(nrky(nrank)*mym+jy-2)
                    end do
                    
                end do
            end do
        end do
        
        call mpi_transfersm(data_in,mm)
        
        return
    end subroutine filter_high_mode
!************************************************************************************************
!                                                                                               *
!       update boundary of vector sub_domain on 2,3,mx(mz,my)-1,mx(mz,my)-2 grid pionts         *
!             related to 26 sub_domain(8 points + 12 lines + 6 faces)                           *
!                                                                                               *
!************************************************************************************************
    subroutine update_vector_sub_domain_bounnary(data_in)
        use var
        implicit none
        include 'mpif.h'

        real*8, dimension(mx,mz,my,3), intent(inout) :: data_in
        real*8, dimension(mx,mz,my,3) :: data_tmp
        integer, parameter :: num_trans = 1
        real*8, dimension(:,:) ,    allocatable :: data_sent_point, data_recv_point
        real*8, dimension(:,:,:),   allocatable :: data_sent_line, data_recv_line
        real*8, dimension(:,:,:,:), allocatable :: data_sent_face, data_recv_face
        integer :: dest, source , dest_nrkx, dest_nrky, dest_nrkz
        integer :: size_1st, size_2nd, max_length, num_zero
        integer :: s1,s2,s3
        integer, dimension(num_trans) :: start_x, start_y, start_z, end_x, end_y, end_z
        integer, dimension(num_trans) :: iii,jjj,kkk
        
        
        integer :: STATUS(MPI_STATUS_SIZE)
        
        max_length = max(max(mx,my),mz)
        
        allocate(data_sent_point(3,num_trans*num_trans*num_trans))
        allocate(data_recv_point(3,num_trans*num_trans*num_trans))
    
        allocate(data_sent_line(max_length,3,num_trans*num_trans))
        allocate(data_recv_line(max_length,3,num_trans*num_trans))
    
        allocate(data_sent_face(max_length,max_length,3,num_trans))
        allocate(data_recv_face(max_length,max_length,3,num_trans))
        
        
        start_x(1) = 3
        start_y(1) = 3
        start_z(1) = 3
        end_x(1) = mx -1
        end_y(1) = my -1
        end_z(1) = mz -1
!if num_trans = 1
!        start_x(2) = 2
!        start_y(2) = 2
!        start_z(2) = 2
!        end_x(2) = mx -2
!        end_y(2) = my -2
!        end_z(2) = mz -2

    
        data_tmp = 0.0 
        
        do i=-1,1
            do j=-1,1
                do k=-1,1
                    dest_nrkx = nrkx(nrank) + i
                    dest_nrky = nrky(nrank) + j
                    dest_nrkz = nrkz(nrank) + k
                
                    !boundary condition
                    if(dest_nrky == npry) dest_nrky = 0
                    if(dest_nrky == -1)   dest_nrky = npry-1
                
                    if(dest_nrkx > nprx - 1 .or. dest_nrkx < 0 .or. dest_nrkz > nprz - 1 .or. dest_nrkz < 0) then
                        dest = -1
                    else
                        dest = dest_nrky*nprx*nprz + dest_nrkz*nprx + dest_nrkx
                    end if
                
                    if(dest /= -1 .and. dest /= nrank) then
                        source = dest
                        num_zero = 0
                        if(i == 0) num_zero = num_zero + 1
                        if(j == 0) num_zero = num_zero + 1
                        if(k == 0) num_zero = num_zero + 1
        
        
                    !transfer 8 points
                        if(num_zero == 0) then
    
                            do s1 = 1,num_trans
                                if(i == -1) then
                                    iii(s1) = start_x(s1)
                                else
                                    iii(s1) = end_x(s1)
                                end if
                                
                                if(j == -1) then
                                    jjj(s1) = start_y(s1)
                                else
                                    jjj(s1) = end_y(s1)
                                end if
                                
                                if(k == -1) then
                                    kkk(s1) = start_z(s1)
                                else
                                    kkk(s1) = end_z(s1)
                                end if
                            end do                               
                            
                            s = 1
                            do s1 = 1,num_trans
                                do s2 = 1,num_trans
                                    do s3 = 1,num_trans
                                        data_sent_point(:,s) = data_in(iii(s1),kkk(s3),jjj(s2),:)
                                        s = s + 1
                                    end do
                                end do
                            end do
    
                            call MPI_SENDRECV(data_sent_point(:,:), 3*num_trans*num_trans*num_trans, MPI_DOUBLE_PRECISION, dest, 0, data_recv_point(:,:), 3*num_trans*num_trans*num_trans, MPI_DOUBLE_PRECISION, source, 0, MPI_COMM_WORLD, STATUS, IERROR)
    
                            s = 1
                            do s1 = 1,num_trans
                                do s2 = 1,num_trans
                                    do s3 = 1,num_trans
                                        data_tmp(iii(s1),kkk(s3),jjj(s2),:) = data_tmp(iii(s1),kkk(s3),jjj(s2),:) + data_recv_point(:,s)
                                        s = s + 1
                                    end do
                                end do
                            end do
                        !transfer 12 lines
                        else if(num_zero == 1) then    
    
                            do s1 = 1,num_trans
                                if(i == -1) then
                                    iii(s1) = start_x(s1)
                                else
                                    iii(s1) = end_x(s1)
                                end if
                                
                                if(j == -1) then
                                    jjj(s1) = start_y(s1)
                                else
                                    jjj(s1) = end_y(s1)
                                end if
                                
                                if(k == -1) then
                                    kkk(s1) = start_z(s1)
                                else
                                    kkk(s1) = end_z(s1)
                                end if
                            end do
    
                            s = 1
                            do s1 = 1,num_trans  
                                do s2 = 1,num_trans
                                    if(i == 0) then
                                        size_1st = mx
                                        data_sent_line(1:size_1st,:,s) = data_in(1:size_1st,kkk(s1),jjj(s2),:)
                                    end if
                                    if(j == 0) then
                                        size_1st = my
                                        data_sent_line(1:size_1st,:,s) = data_in(iii(s1),kkk(s2),1:size_1st,:)
                                    end if
                                    if(k == 0) then
                                        size_1st= mz 
                                        data_sent_line(1:size_1st,:,s) = data_in(iii(s1),1:size_1st,jjj(s2),:)
                                    end if
                                    s = s + 1
                                end do
                            end do
    
                            call MPI_SENDRECV(data_sent_line(1:size_1st,:,:), 3*size_1st*num_trans*num_trans, MPI_DOUBLE_PRECISION, dest, 0, data_recv_line(1:size_1st,:,:), 3*size_1st*num_trans*num_trans, MPI_DOUBLE_PRECISION, source, 0, MPI_COMM_WORLD, STATUS, IERROR)
    
                            s = 1
                            do s1 = 1,num_trans  
                                do s2 = 1,num_trans
                                    if(i == 0) data_tmp(1:mx,kkk(s1),jjj(s2),:) = data_tmp(1:mx,kkk(s1),jjj(s2),:) + data_recv_line(1:size_1st,:,s)
                                    if(j == 0) data_tmp(iii(s1),kkk(s2),1:my,:) = data_tmp(iii(s1),kkk(s2),1:my,:) + data_recv_line(1:size_1st,:,s)
                                    if(k == 0) data_tmp(iii(s1),1:mz,jjj(s2),:) = data_tmp(iii(s1),1:mz,jjj(s2),:) + data_recv_line(1:size_1st,:,s)
                                    s = s + 1
                                end do
                            end do

                        !transfer 6 faces
                        else if(num_zero == 2) then
                            do s1 = 1,num_trans
                                if(i == -1) then
                                    iii(s1) = start_x(s1)
                                else
                                    iii(s1) = end_x(s1)
                                end if
                                
                                if(j == -1) then
                                    jjj(s1) = start_y(s1)
                                else
                                    jjj(s1) = end_y(s1)
                                end if
                                
                                if(k == -1) then
                                    kkk(s1) = start_z(s1)
                                else
                                    kkk(s1) = end_z(s1)
                                end if
                            end do
    
                            s = 1
                            do s1 = 1,num_trans  
                                if(i /= 0) then
                                    size_1st = mz
                                    size_2nd = my
                                    data_sent_face(1:mz,1:my,:,s) = data_in(iii(s1),1:mz,1:my,:)
                                end if
                                if(j /= 0) then
                                    size_1st = mx
                                    size_2nd = mz
                                    data_sent_face(1:mx,1:mz,:,s) = data_in(1:mx,1:mz,jjj(s1),:)
                                end if
                                if(k /= 0) then
                                    size_1st = mx
                                    size_2nd = my
                                    data_sent_face(1:mx,1:my,:,s) = data_in(1:mx,kkk(s1),1:my,:)
                                end if
                                s = s + 1
                            end do
    
                            call MPI_SENDRECV(data_sent_face(1:size_1st,1:size_2nd,:,:), 3*size_1st*size_2nd*num_trans, MPI_DOUBLE_PRECISION, dest, 0, data_recv_face(1:size_1st,1:size_2nd,:,:), 3*size_1st*size_2nd*num_trans, MPI_DOUBLE_PRECISION, source, 0, MPI_COMM_WORLD, STATUS, IERROR)

                            s = 1
                            do s1 = 1,num_trans  
                                if(i /= 0) data_tmp(iii(s1),1:mz,1:my,:) = data_tmp(iii(s1),1:mz,1:my,:) + data_recv_face(1:size_1st,1:size_2nd,:,s)
                                if(j /= 0) data_tmp(1:mx,1:mz,jjj(s1),:) = data_tmp(1:mx,1:mz,jjj(s1),:) + data_recv_face(1:size_1st,1:size_2nd,:,s)
                                if(k /= 0) data_tmp(1:mx,kkk(s1),1:my,:) = data_tmp(1:mx,kkk(s1),1:my,:) + data_recv_face(1:size_1st,1:size_2nd,:,s)
                                s = s + 1
                            end do
    
                        end if
                    end if
                    
                end do
            end do
        end do
        
        data_in = data_in + data_tmp
    
        deallocate(data_sent_point)
        deallocate(data_recv_point)
        deallocate(data_sent_line)
        deallocate(data_recv_line)
        deallocate(data_sent_face)
        deallocate(data_recv_face)

    end subroutine update_vector_sub_domain_bounnary
    
    
!************************************************************************************************
!                                                                                               *
!       update boundary of scalar sub_domain on 2,3,mx(mz,my)-1,mx(mz,my)-2 grid pionts         *
!               related to 26 sub_domain(8 points + 12 lines + 6 faces)                         *
!                                                                                               *
!************************************************************************************************
    subroutine update_scalar_sub_domain_bounnary(data_in)
        use var
        implicit none
        include 'mpif.h'

        real*8, dimension(mx,mz,my), intent(inout) :: data_in
        real*8, dimension(mx,mz,my) :: data_tmp
        integer, parameter :: num_trans = 1
        real*8, dimension(:) ,    allocatable :: data_sent_point, data_recv_point
        real*8, dimension(:,:),   allocatable :: data_sent_line, data_recv_line
        real*8, dimension(:,:,:), allocatable :: data_sent_face, data_recv_face
        integer :: dest, source , dest_nrkx, dest_nrky, dest_nrkz
        integer :: size_1st, size_2nd, max_length, num_zero
        integer :: s1,s2,s3
        integer, dimension(num_trans) :: start_x, start_y, start_z, end_x, end_y, end_z
        integer, dimension(num_trans) :: iii,jjj,kkk
        
        
        integer :: STATUS(MPI_STATUS_SIZE)
        
        max_length = max(max(mx,my),mz)
        
        allocate(data_sent_point(num_trans*num_trans*num_trans))
        allocate(data_recv_point(num_trans*num_trans*num_trans))
    
        allocate(data_sent_line(max_length,num_trans*num_trans))
        allocate(data_recv_line(max_length,num_trans*num_trans))
    
        allocate(data_sent_face(max_length,max_length,num_trans))
        allocate(data_recv_face(max_length,max_length,num_trans))
        
        
        start_x(1) = 3
        start_y(1) = 3
        start_z(1) = 3
        end_x(1) = mx -1
        end_y(1) = my -1
        end_z(1) = mz -1
!if num_trans = 1
!        start_x(2) = 2
!        start_y(2) = 2
!        start_z(2) = 2
!        end_x(2) = mx -2
!        end_y(2) = my -2
!        end_z(2) = mz -2

    
        data_tmp = 0.0 
        
        do i=-1,1
            do j=-1,1
                do k=-1,1
                    dest_nrkx = nrkx(nrank) + i
                    dest_nrky = nrky(nrank) + j
                    dest_nrkz = nrkz(nrank) + k
                
                    !boundary condition
                    if(dest_nrky == npry) dest_nrky = 0
                    if(dest_nrky == -1)   dest_nrky = npry-1
                
                    if(dest_nrkx > nprx - 1 .or. dest_nrkx < 0 .or. dest_nrkz > nprz - 1 .or. dest_nrkz < 0) then
                        dest = -1
                    else
                        dest = dest_nrky*nprx*nprz + dest_nrkz*nprx + dest_nrkx
                    end if
                
                    if(dest /= -1 .and. dest /= nrank) then
                        source = dest
                        num_zero = 0
                        if(i == 0) num_zero = num_zero + 1
                        if(j == 0) num_zero = num_zero + 1
                        if(k == 0) num_zero = num_zero + 1
        
        
                    !transfer 8 points
                        if(num_zero == 0) then
    
                            do s1 = 1,num_trans
                                if(i == -1) then
                                    iii(s1) = start_x(s1)
                                else
                                    iii(s1) = end_x(s1)
                                end if
                                
                                if(j == -1) then
                                    jjj(s1) = start_y(s1)
                                else
                                    jjj(s1) = end_y(s1)
                                end if
                                
                                if(k == -1) then
                                    kkk(s1) = start_z(s1)
                                else
                                    kkk(s1) = end_z(s1)
                                end if
                            end do                               
                            
                            s = 1
                            do s1 = 1,num_trans
                                do s2 = 1,num_trans
                                    do s3 = 1,num_trans
                                        data_sent_point(s) = data_in(iii(s1),kkk(s3),jjj(s2))
                                        s = s + 1
                                    end do
                                end do
                            end do
    
                            call MPI_SENDRECV(data_sent_point(:), num_trans*num_trans*num_trans, MPI_DOUBLE_PRECISION, dest, 0, data_recv_point(:), num_trans*num_trans*num_trans, MPI_DOUBLE_PRECISION, source, 0, MPI_COMM_WORLD, STATUS, IERROR)
    
                            s = 1
                            do s1 = 1,num_trans
                                do s2 = 1,num_trans
                                    do s3 = 1,num_trans
                                        data_tmp(iii(s1),kkk(s3),jjj(s2)) = data_tmp(iii(s1),kkk(s3),jjj(s2)) + data_recv_point(s)
                                        s = s + 1
                                    end do
                                end do
                            end do
                        !transfer 12 lines
                        else if(num_zero == 1) then    
    
                            do s1 = 1,num_trans
                                if(i == -1) then
                                    iii(s1) = start_x(s1)
                                else
                                    iii(s1) = end_x(s1)
                                end if
                                
                                if(j == -1) then
                                    jjj(s1) = start_y(s1)
                                else
                                    jjj(s1) = end_y(s1)
                                end if
                                
                                if(k == -1) then
                                    kkk(s1) = start_z(s1)
                                else
                                    kkk(s1) = end_z(s1)
                                end if
                            end do
    
                            s = 1
                            do s1 = 1,num_trans  
                                do s2 = 1,num_trans
                                    if(i == 0) then
                                        size_1st = mx
                                        data_sent_line(1:size_1st,s) = data_in(1:size_1st,kkk(s1),jjj(s2))
                                    end if
                                    if(j == 0) then
                                        size_1st = my
                                        data_sent_line(1:size_1st,s) = data_in(iii(s1),kkk(s2),1:size_1st)
                                    end if
                                    if(k == 0) then
                                        size_1st= mz 
                                        data_sent_line(1:size_1st,s) = data_in(iii(s1),1:size_1st,jjj(s2))
                                    end if
                                    s = s + 1
                                end do
                            end do
    
                            call MPI_SENDRECV(data_sent_line(1:size_1st,:), size_1st*num_trans*num_trans, MPI_DOUBLE_PRECISION, dest, 0, data_recv_line(1:size_1st,:), size_1st*num_trans*num_trans, MPI_DOUBLE_PRECISION, source, 0, MPI_COMM_WORLD, STATUS, IERROR)
    
                            s = 1
                            do s1 = 1,num_trans  
                                do s2 = 1,num_trans
                                    if(i == 0) data_tmp(1:mx,kkk(s1),jjj(s2)) = data_tmp(1:mx,kkk(s1),jjj(s2)) + data_recv_line(1:size_1st,s)
                                    if(j == 0) data_tmp(iii(s1),kkk(s2),1:my) = data_tmp(iii(s1),kkk(s2),1:my) + data_recv_line(1:size_1st,s)
                                    if(k == 0) data_tmp(iii(s1),1:mz,jjj(s2)) = data_tmp(iii(s1),1:mz,jjj(s2)) + data_recv_line(1:size_1st,s)
                                    s = s + 1
                                end do
                            end do

                        !transfer 6 faces
                        else if(num_zero == 2) then
                            do s1 = 1,num_trans
                                if(i == -1) then
                                    iii(s1) = start_x(s1)
                                else
                                    iii(s1) = end_x(s1)
                                end if
                                
                                if(j == -1) then
                                    jjj(s1) = start_y(s1)
                                else
                                    jjj(s1) = end_y(s1)
                                end if
                                
                                if(k == -1) then
                                    kkk(s1) = start_z(s1)
                                else
                                    kkk(s1) = end_z(s1)
                                end if
                            end do
    
                            s = 1
                            do s1 = 1,num_trans  
                                if(i /= 0) then
                                    size_1st = mz
                                    size_2nd = my
                                    data_sent_face(1:mz,1:my,s) = data_in(iii(s1),1:mz,1:my)
                                end if
                                if(j /= 0) then
                                    size_1st = mx
                                    size_2nd = mz
                                    data_sent_face(1:mx,1:mz,s) = data_in(1:mx,1:mz,jjj(s1))
                                end if
                                if(k /= 0) then
                                    size_1st = mx
                                    size_2nd = my
                                    data_sent_face(1:mx,1:my,s) = data_in(1:mx,kkk(s1),1:my)
                                end if
                                s = s + 1
                            end do
    
                            call MPI_SENDRECV(data_sent_face(1:size_1st,1:size_2nd,:), size_1st*size_2nd*num_trans, MPI_DOUBLE_PRECISION, dest, 0, data_recv_face(1:size_1st,1:size_2nd,:), size_1st*size_2nd*num_trans, MPI_DOUBLE_PRECISION, source, 0, MPI_COMM_WORLD, STATUS, IERROR)

                            s = 1
                            do s1 = 1,num_trans  
                                if(i /= 0) data_tmp(iii(s1),1:mz,1:my) = data_tmp(iii(s1),1:mz,1:my) + data_recv_face(1:size_1st,1:size_2nd,s)
                                if(j /= 0) data_tmp(1:mx,1:mz,jjj(s1)) = data_tmp(1:mx,1:mz,jjj(s1)) + data_recv_face(1:size_1st,1:size_2nd,s)
                                if(k /= 0) data_tmp(1:mx,kkk(s1),1:my) = data_tmp(1:mx,kkk(s1),1:my) + data_recv_face(1:size_1st,1:size_2nd,s)
                                s = s + 1
                            end do
    
                        end if
                    end if
                    
                end do
            end do
        end do
        
        data_in = data_in + data_tmp
    
        deallocate(data_sent_point)
        deallocate(data_recv_point)
        deallocate(data_sent_line)
        deallocate(data_recv_line)
        deallocate(data_sent_face)
        deallocate(data_recv_face)

    end subroutine update_scalar_sub_domain_bounnary
    
    subroutine calc_analytical_form
        use distribution
        use var
        implicit none
        include 'mpif.h'
        
        integer :: nn
        real*8, dimension(mx,mz)   :: p_h_para0, p_h_perp0, nv_para0
        real*8, dimension(mx,mz)   :: n_h0
        real*8, dimension(mx,mz,3) :: J_h_para0
        real*8, dimension(mx,mz,my):: data_tmp
        real*8 :: mybeta_avg, beta_avg, myVolume, Volume
        real*8 :: v_para, v_perp2, dv_para, dv_perp2
        real*8 :: c_norm
        real*8 :: d1f2, d1fc
        real*8 :: f_jm1,f_jm2,f_j,f_jp1,f_jp2,x_jm1,x_j,x_jp1,coeff_a,coeff_b,coeff_c,coeff_d
        integer ix_first, ix_last,iz_first, iz_last,iy_first,iy_last
        integer jx,jy,jz

!  define statement functions
!  d1f2= d f / dx  with second-order accuracy central difference
        
        d1f2(f_jm1,f_j,f_jp1,x_jm1,x_j,x_jp1)= &
        ((x_jm1-x_j)/(x_jp1-x_j)*(f_jp1-f_j) &
        -(x_jp1-x_j)/(x_jm1-x_j)*(f_jm1-f_j))/(x_jm1-x_jp1)  

!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(f_jm2,f_jm1,f_j,f_jp1,f_jp2,coeff_a,coeff_b,coeff_c,coeff_d)= &
       coeff_a*(f_jp1-f_j)+coeff_b*(f_j-f_jm1)+coeff_c*(f_jp2-f_j)+coeff_d*(f_j-f_jm2)
        
        
        ix_first=1
        ix_last=mx
        iz_first=1
        iz_last=mz
        iy_first=1
        iy_last=my
        
        nn = 1001
        dv_para  = (v_para_max - v_para_min)/(nn-1)
        dv_perp2 = (v_perp_max**2 - v_perp_min**2)/(nn-1)
                        
        mybeta_avg = 0
        myVolume   = 0
        J_h_para0  = 0
        p_h_para0  = 0
        p_h_perp0  = 0
        nv_para0   = 0
        n_h0       = 0
        
        do i=1,mx
            do j=1,mz
                if(psi(i,j) < psmax) then
                    R      = xx(i)
                    B      = xint(i,j,6:8)
                    abs_B  = sqrt(dot_product(B,B))
                    b_unit = B/abs_B
                    !integrate in d^3v
                    do ii=1,nn
                        do jj=1,nn
                            v_para  = v_para_min    + (ii-1)*dv_para
                            v_perp2 = v_perp_min**2 + (jj-1)*dv_perp2
                            
                            !for mu = 0 case!!!
                            if((.not. IIV) .and. (.not. FGP)) v_perp2 = 0.0
                        
                            v      = sqrt(v_para**2 + v_perp2)                      
                            E      = 0.5*m*v**2                        
                            mu     = 0.5*m*v_perp2/abs_B
                            P_phi  = m*v_para*R*B(2)/abs_B - Zh*psi(i,j)
                            if(E > 1e-12) then
                                Lambda = mu*B0/E
                            else
                                Lambda = 0
                            end if

                            if(v_para>0) then
                                sgn_v_parallel = +1
                            else
                                sgn_v_parallel = -1
                            end if
      
                            if(ORBIT_AVERAGE_METHOD == 1) then
                                if((1-Lambda) > 0) then
                                    bracket_psi = - P_phi/Zh + m_div_Zh*sgn_v_parallel*v*R0*sqrt(1-Lambda)
                                else
                                    bracket_psi = - P_phi/Zh
                                end if
                            else if(ORBIT_AVERAGE_METHOD == 2) then
                                bracket_psi = - P_phi/Zh
                            end if
                        
                            
                           if(DFR_type == 1) then
                                if(FGP) then
                                    f0 = 1.0/(v**3+v_c**3)*(1+erf((v_0-v)/Delta_v))*exp(-(bracket_psi/(c_1*Delta_psi)))*exp(-((Lambda-Lambda_0)/Delta_Lambda)**2)
                                else
                                    f0 = 1.0/(v**3+v_c**3)*(1+erf((v_0-v)/Delta_v))*exp(-(bracket_psi/(c_1*Delta_psi)))
                                end if
                            else if(DFR_type == 2) then
                                bracket_psi_norm = (bracket_psi - psmin)/(psmax-psmin)
                                f0 = 1.0/(v**3+v_c**3)*(1+erf((v_0-v)/Delta_v))*exp(-Delta_n/L_n*tanh((bracket_psi_norm - psi_n)/Delta_n))
                            end if
                        
                            J_h_para0(i,j,:) = J_h_para0(i,j,:) + v_para*b_unit*f0*dv_para*dv_perp2
                            p_h_para0(i,j)   = p_h_para0(i,j)   + m*v_para**2*f0*dv_para*dv_perp2
                            p_h_perp0(i,j)   = p_h_perp0(i,j)   + 0.5*m*v_perp2*f0*dv_para*dv_perp2
                            nv_para0(i,j)    = nv_para0(i,j)    + v_para*f0*dv_para*dv_perp2
                            n_h0(i,j)        = n_h0(i,j)        + f0*dv_para*dv_perp2
                        end do
                    end do
                                                
                    mybeta_avg = mybeta_avg + (p_h_para0(i,j) + p_h_perp0(i,j))/abs_B**2*dxx*dzz
                    myVolume   = myVolume + dxx*dzz
                end if
                

            end do
        end do
        
        call MPI_REDUCE(mybeta_avg, beta_avg, 1, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, IERROR)
        call MPI_BCAST(beta_avg, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)
        
        call MPI_REDUCE(myVolume, Volume, 1, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, IERROR)
        call MPI_BCAST(Volume, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)
        
        c_norm = beta/(beta_avg/Volume)
        
        J_h_para0 = c_norm*J_h_para0
        p_h_para0 = c_norm*p_h_para0
        p_h_perp0 = c_norm*p_h_perp0        
        nv_para0  = c_norm*nv_para0
        n_h0      = c_norm*n_h0
        
        do i=1,my
            data_tmp(:,:,i) = p_h_perp0
        end do
        
        call mpi_transfersm(data_tmp,1)
        
        p_h_perp0 = data_tmp(:,:,3)
        
        grad_p_h_para0 = 0
        grad_p_h_perp0 = 0
        grad_nv_para0  = 0
        grad_n_h0      = 0
        do jz=iz_first+2,iz_last-2
            do jx=ix_first+2,ix_last-2
                if(FOUR_TH) then
                    grad_p_h_para0(jx,jz,1) = d1fc(p_h_para0(jx-2,jz),p_h_para0(jx-1,jz),p_h_para0(jx,jz),p_h_para0(jx+1,jz),p_h_para0(jx+2,jz),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
                    grad_p_h_para0(jx,jz,3) = d1fc(p_h_para0(jx,jz-2),p_h_para0(jx,jz-1),p_h_para0(jx,jz),p_h_para0(jx,jz+1),p_h_para0(jx,jz+2),az1(jz),bz1(jz),cz1(jz),dz1(jz))
                    grad_p_h_perp0(jx,jz,1) = d1fc(p_h_perp0(jx-2,jz),p_h_perp0(jx-1,jz),p_h_perp0(jx,jz),p_h_perp0(jx+1,jz),p_h_perp0(jx+2,jz),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
                    grad_p_h_perp0(jx,jz,3) = d1fc(p_h_perp0(jx,jz-2),p_h_perp0(jx,jz-1),p_h_perp0(jx,jz),p_h_perp0(jx,jz+1),p_h_perp0(jx,jz+2),az1(jz),bz1(jz),cz1(jz),dz1(jz))
                    grad_nv_para0(jx,jz,1)  = d1fc(nv_para0(jx-2,jz),nv_para0(jx-1,jz),nv_para0(jx,jz),nv_para0(jx+1,jz),nv_para0(jx+2,jz),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
                    grad_nv_para0(jx,jz,3)  = d1fc(nv_para0(jx,jz-2),nv_para0(jx,jz-1),nv_para0(jx,jz),nv_para0(jx,jz+1),nv_para0(jx,jz+2),az1(jz),bz1(jz),cz1(jz),dz1(jz))
                    grad_n_h0(jx,jz,1)      = d1fc(n_h0(jx-2,jz), n_h0(jx-1,jz), n_h0(jx,jz), n_h0(jx+1,jz), n_h0(jx+2,jz), az1(jx),bz1(jx),cz1(jx),dz1(jx))
                    grad_n_h0(jx,jz,3)      = d1fc(n_h0(jx,jz-2), n_h0(jx,jz-1), n_h0(jx,jz), n_h0(jx,jz+1), n_h0(jx,jz+2), az1(jz),bz1(jz),cz1(jz),dz1(jz))
                else                    
                    grad_p_h_para0(jx,jz,1) = d1f2(p_h_para0(jx-1,jz),p_h_para0(jx,jz),p_h_para0(jx+1,jz),xx(jx-1),xx(jx),xx(jx+1))
                    grad_p_h_para0(jx,jz,3) = d1f2(p_h_para0(jx,jz-1),p_h_para0(jx,jz),p_h_para0(jx,jz+1),zz(jz-1),zz(jz),zz(jz+1))
                    grad_p_h_perp0(jx,jz,1) = d1f2(p_h_perp0(jx-1,jz),p_h_perp0(jx,jz),p_h_perp0(jx+1,jz),xx(jx-1),xx(jx),xx(jx+1))
                    grad_p_h_perp0(jx,jz,3) = d1f2(p_h_perp0(jx,jz-1),p_h_perp0(jx,jz),p_h_perp0(jx,jz+1),zz(jz-1),zz(jz),zz(jz+1))
                    grad_nv_para0(jx,jz,1)  = d1f2(nv_para0(jx-1,jz),nv_para0(jx,jz),nv_para0(jx+1,jz),xx(jx-1),xx(jx),xx(jx+1))
                    grad_nv_para0(jx,jz,3)  = d1f2(nv_para0(jx,jz-1),nv_para0(jx,jz),nv_para0(jx,jz+1),zz(jz-1),zz(jz),zz(jz+1))
                    grad_n_h0(jx,jz,1)      = d1f2(n_h0(jx-1,jz), n_h0(jx,jz), n_h0(jx+1,jz), xx(jx-1),xx(jx),xx(jx+1))
                    grad_n_h0(jx,jz,3)      = d1f2(n_h0(jx,jz-1), n_h0(jx,jz), n_h0(jx,jz+1), zz(jz-1),zz(jz),zz(jz+1))
                end if
            end do
        end do
        grad_p_h_perp0(:,:,2) = 0
        grad_p_h_para0(:,:,2) = 0
        grad_nv_para0(:,:,2)  = 0
        grad_n_h0(:,:,2)      = 0
        
        do i=1,mx
            do j=1,mz
                if(psi(i,j) < psmax) then
                    B      = xint(i,j,6:8)
                    abs_B  = sqrt(dot_product(B,B))
                    b_unit = b_unit_grid_eq(i,j,:)
                    curl_b = curl_b_grid_eq(i,j,:)

                    J_h0_ana(i,j,:) = J_h_para0(i,j,:) + 1.0/abs_B*(p_h_para0(i,j)-p_h_perp0(i,j))*curl_b - 1.0/abs_B*cross_product(grad_p_h_perp0(i,j,:),b_unit)
                end if
            end do
        end do        
        
        do i=1,my
            J_h0(:,:,i,:) = J_h0_ana
            p_h_perp0_3D(:,:,i) = p_h_perp0
            p_h_para0_3D(:,:,i) = p_h_para0
            nv_para0_3D(:,:,i)  = nv_para0
            n_h0_3D(:,:,i) = n_h0 
            grad_p_h_para0_3D(:,:,i,:) = grad_p_h_para0
            grad_p_h_perp0_3D(:,:,i,:) = grad_p_h_perp0
        end do
        call bndry3_ex(J_h0, 0)
        J_h0_ana = J_h0(:,:,3,:)
        
        !define p_h0
        p_h0    = 0.5*(p_h_para0 + p_h_perp0)
        dp_h0dR = 0.5*(grad_p_h_para0(:,:,1) + grad_p_h_perp0(:,:,1))
        dp_h0dZ = 0.5*(grad_p_h_para0(:,:,3) + grad_p_h_perp0(:,:,3))
        
      
        
    end subroutine calc_analytical_form    
    
    subroutine calc_polarization_drift_var
        use distribution
        use var
        implicit none
        include 'mpif.h'    

        real*8, dimension(mx,mz,my,3) :: dv_EdR, dv_Edphi, dv_EdZ
        real*8 :: d1f2, d1fc
        real*8 :: f_jm1,f_jm2,f_j,f_jp1,f_jp2,x_jm1,x_j,x_jp1,coeff_a,coeff_b,coeff_c,coeff_d
        integer ix_first, ix_last,iz_first, iz_last,iy_first,iy_last
        integer jx,jy,jz
        
!  define statement functions
!  d1f2= d f / dx  with second-order accuracy central difference
        
        d1f2(f_jm1,f_j,f_jp1,x_jm1,x_j,x_jp1)= &
        ((x_jm1-x_j)/(x_jp1-x_j)*(f_jp1-f_j) &
        -(x_jp1-x_j)/(x_jm1-x_j)*(f_jm1-f_j))/(x_jm1-x_jp1)
        
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(f_jm2,f_jm1,f_j,f_jp1,f_jp2,coeff_a,coeff_b,coeff_c,coeff_d)= &
       coeff_a*(f_jp1-f_j)+coeff_b*(f_j-f_jm1)+coeff_c*(f_jp2-f_j)+coeff_d*(f_j-f_jm2)
        
        ix_first=1
        ix_last=mx
        iz_first=1
        iz_last=mz
        iy_first=1
        iy_last=my
        
        do i=1,mx
            do j=1,mz
                do k=1,my
                    abs_B = sqrt(dot_product(x(i,j,k,6:8),x(i,j,k,6:8)))
                    v_E_grid(i,j,k,:)  = cross_product(Ef(i,j,k,:), x(i,j,k,6:8))/abs_B
                    v_E2_grid(i,j,k) = dot_product(v_E_grid(i,j,k,:),v_E_grid(i,j,k,:))
                end do
            end do
        end do
        
        call mpi_transfersm(v_E_grid,3)
        call mpi_transfersm(v_E2_grid,1)
        
        
        do jy=iy_first,iy_last
            do jz=iz_first+2,iz_last-2
                do jx=ix_first+2,ix_last-2
                    if(FOUR_TH) then
                        do s=1,3                        
                            dv_EdR(jx,jz,jy,s) = d1fc(v_E_grid(jx-2,jz,jy,s),v_E_grid(jx-1,jz,jy,s),v_E_grid(jx,jz,jy,s),v_E_grid(jx+1,jz,jy,s),v_E_grid(jx+2,jz,jy,s),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
                            dv_EdZ(jx,jz,jy,s) = d1fc(v_E_grid(jx,jz-2,jy,s),v_E_grid(jx,jz-1,jy,s),v_E_grid(jx,jz,jy,s),v_E_grid(jx,jz+1,jy,s),v_E_grid(jx,jz+2,jy,s),az1(jz),bz1(jz),cz1(jz),dz1(jz))
                        end do
                        grad_v_E_2_grid(jx,jz,jy,1) = d1fc(v_E2_grid(jx-2,jz,jy),v_E2_grid(jx-1,jz,jy),v_E2_grid(jx,jz,jy),v_E2_grid(jx+1,jz,jy),v_E2_grid(jx+2,jz,jy),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
                        grad_v_E_2_grid(jx,jz,jy,3) = d1fc(v_E2_grid(jx,jz-2,jy),v_E2_grid(jx,jz-1,jy),v_E2_grid(jx,jz,jy),v_E2_grid(jx,jz+1,jy),v_E2_grid(jx,jz+2,jy),az1(jz),bz1(jz),cz1(jz),dz1(jz))
                    else
                        do s=1,3
                            dv_EdR(jx,jz,jy,s) = d1f2(v_E_grid(jx-1,jz,jy,s),v_E_grid(jx,jz,jy,s),v_E_grid(jx+1,jz,jy,s),xx(jx-1),xx(jx),xx(jx+1))
                            dv_EdZ(jx,jz,jy,s) = d1f2(v_E_grid(jx,jz-1,jy,s),v_E_grid(jx,jz,jy,s),v_E_grid(jx,jz+1,jy,s),zz(jz-1),zz(jz),zz(jz+1))
                        end do
                        grad_v_E_2_grid(jx,jz,jy,1) = d1f2(v_E2_grid(jx-1,jz,jy),v_E2_grid(jx,jz,jy),v_E2_grid(jx+1,jz,jy),xx(jx-1),xx(jx),xx(jx+1))
                        grad_v_E_2_grid(jx,jz,jy,3) = d1f2(v_E2_grid(jx,jz-1,jy),v_E2_grid(jx,jz,jy),v_E2_grid(jx,jz+1,jy),zz(jz-1),zz(jz),zz(jz+1))
                    end if
                end do
            end do
        end do
        
        
        do jz=iz_first,iz_last
            do jx=ix_first,ix_last
                do jy=iy_first+2,iy_last-2
                    if(FOUR_TH) then
                        do s=1,3
                            dv_Edphi(jx,jz,jy,s) = d1fc(v_E_grid(jx,jz,jy-2,s),v_E_grid(jx,jz,jy-1,s),v_E_grid(jx,jz,jy,s),v_E_grid(jx,jz,jy+1,s),v_E_grid(jx,jz,jy+2,s),ay1(jy),by1(jy),cy1(jy),dy1(jy))/xx(jx)
                        end do
                        grad_v_E_2_grid(jx,jz,jy,2) = d1fc(v_E2_grid(jx,jz,jy-2),v_E2_grid(jx,jz,jy-1),v_E2_grid(jx,jz,jy),v_E2_grid(jx,jz,jy+1),v_E2_grid(jx,jz,jy+2),ay1(jy),by1(jy),cy1(jy),dy1(jy))/xx(jx)
                    else
                        do s=1,3
                            dv_Edphi(jx,jz,jy,s) = d1f2(v_E_grid(jx,jz,jy-1,s),v_E_grid(jx,jz,jy,s),v_E_grid(jx,jz,jy+1,s),yy(jy-1),yy(jy),yy(jy+1))/xx(jx)
                        end do
                        grad_v_E_2_grid(jx,jz,jy,2) = d1f2(v_E2_grid(jx,jz,jy-1),v_E2_grid(jx,jz,jy),v_E2_grid(jx,jz,jy+1),yy(jy-1),yy(jy),yy(jy+1))/xx(jx)
                    end if
                end do
            end do
        end do
              
        curl_v_E_grid(:,:,:,1)   = dv_Edphi(:,:,:,3) - dv_EdZ(:,:,:,2)
        curl_v_E_grid(:,:,:,2)   = dv_EdZ(:,:,:,1)   - dv_EdR(:,:,:,3)
        
        do jz=iz_first,iz_last
            do jx=ix_first,ix_last
                do jy=iy_first,iy_last
                    curl_v_E_grid(jx,jz,jy,3)   = v_E_grid(jx,jz,jy,2)/xx(jx) + dv_EdR(jx,jz,jy,2) - dv_Edphi(jx,jz,jy,1)
                end do
            end do
        end do   
        
        pv_Ept_grid = (v_E_grid - v_E_old_grid)/dt
        
        v_E_old_grid = v_E_grid
    
    end subroutine calc_polarization_drift_var
    
    subroutine particle_orbit_polarization
        use distribution
        use var
        implicit none
        include 'mpif.h'

        
        do p=1,N
            do RK4=1,4
                select case(RK4)
                case(1)
                    Xold     = marker(p)%X
                    Vold     = marker(p)%v_para
                    Wold     = marker(p)%w
                    Xa       = marker(p)%X
                    Va       = marker(p)%v_para
                    Wa       = marker(p)%w
                    told     = time
                    timestep = 0.5*dt
                    coeff    = 1.0/6.0
                case(2)
                    timestep = 0.5*dt
                    coeff    = 1.0/3.0
                case(3)
                    timestep = dt
                    coeff    = 1.0/3.0
                case(4)
                    timestep = dt
                    coeff    = 1.0/6.0
                end select
                
                if(GYRO) then
                    call gyro_average_xyz(N_GY)
                else
                    
                    i = floor((marker(p)%X(1) - myxmin)/dxx) + 3
                    j = floor((marker(p)%X(2) - myymin)/dyy) + 3
                    k = floor((marker(p)%X(3) - myzmin)/dzz) + 3
                
                    if(i>=mx .or. i<1 .or. j>=my .or. j<1 .or. k>=mz .or. k<1) then
                        write(*,*)"zzzzzzzzzzzzzzzzzaaaaaaaaaaaaaa : ",i,j,k,marker(p)%X,marker(p)%id,nrank,time
                        write(*,*)"nnnnnnnnnnnnnnnnnbbbbbbbbbbbbbb : ",i,j,k,marker_initia(p)%X,marker_initia(p)%v_para,marker_initia(p)%v_perp,marker_initia(p)%id,nrank,time
                    end if
                
             
                    R   = xx(i)
                    phi = yy(j)
                    Z   = zz(k)
             
                    dR2 = xx(i+1)**2 - xx(i)**2
             
                    weight_line(1,1) = (marker(p)%X(1)**2 - R**2)/dR2
                    weight_line(2,1) = (marker(p)%X(2) - phi)/dyy
                    weight_line(3,1) = (marker(p)%X(3) - Z)/dzz
                
                    if(psi(i,k)<psmax .and. (psi(i+1,k)>psmax .or. psi(i,k+1)>psmax .or. psi(i+1,k+1)>psmax)) then
                        weight_line(1,1) = 0.0
                        weight_line(3,1) = 0.0
                    endif
             
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
                
                    if(psi(i,k)<psmax .and. (psi(i+1,k)>psmax .or. psi(i,k+1)>psmax .or. psi(i+1,k+1)>psmax)) then
                        weight_line(1,1) = 0.0
                        weight_line(2,1) = 0.0
                    endif
             
                    weight_line(1,2) = 1.0 - weight_line(1,1)
                    weight_line(2,2) = 1.0 - weight_line(2,1)

             
                    do ii=1,2
                        do jj=1,2
                            weight_square(ii,jj) = weight_line(1,3-ii)*weight_line(2,3-jj)
                        end do
                    end do                
                
                    B = 0
                    grad_B  = 0
                    curl_b  = 0
                    delta_E = 0
                    grad_RB_phi = 0
                    pBpt    = 0
                    do ii=1,2
                        do jj=1,2
                            do kk=1,2!note that B in CLT code is Bx(R,Z,phi)
                                B       = B + x(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
                                grad_B  = grad_B  + grad_B_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                                curl_b  = curl_b  + curl_b_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                                delta_E = delta_E + Ef(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                                pBpt    = pBpt    + xdif(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
                                grad_RB_phi(1) = grad_RB_phi(1) + (x(i-1+ii,k-1+kk,j-1+jj,7)+xx(i-1+ii)*xr(i-1+ii,k-1+kk,j-1+jj,7))*weight_cubic(ii,jj,kk)
                                grad_RB_phi(2) = grad_RB_phi(2) + (                                     xy(i-1+ii,k-1+kk,j-1+jj,7))*weight_cubic(ii,jj,kk)
                                grad_RB_phi(3) = grad_RB_phi(3) + (                          xx(i-1+ii)*xz(i-1+ii,k-1+kk,j-1+jj,7))*weight_cubic(ii,jj,kk)
                            end do
                        end do
                    end do
                    abs_B    = sqrt(dot_product(B,B))
                    b_unit   = B/abs_B
                    pabs_Bpt = dot_product(b_unit,pBpt)
                
    !polarization drift related variables                
                    if(VEB) then
                        curl_v_E   = 0
                        grad_v_E_2 = 0
                        pv_Ept     = 0
                        do ii=1,2
                            do jj=1,2
                                do kk=1,2!note that B in CLT code is Bx(R,Z,phi)
                                    curl_v_E   = curl_v_E   + curl_v_E_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                                    grad_v_E_2 = grad_v_E_2 + grad_v_E_2_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                            if(POL) pv_Ept     = pv_Ept     + pv_Ept_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                                end do
                            end do
                        end do
                        !calc pb_unitpt, which affect equation of dXdt and not affect equation of dv_||dt because of pb_unitpt * B = 0
                        pb_unitpt = (pBpt - pabs_Bpt*b_unit)/abs_B
                    end if

                
                    B_eq = 0
                    b_unit_eq = 0
                    grad_B_eq = 0
                    curl_b_eq = 0
                    psi_par   = 0
                    do ii=1,2
                        do kk=1,2
                            B_eq      = B_eq + xint(i-1+ii,k-1+kk,6:8)*weight_square(ii,kk)
                            b_unit_eq = b_unit_eq + b_unit_grid_eq(i-1+ii,k-1+kk,:)*weight_square(ii,kk)
                            grad_B_eq = grad_B_eq + grad_B_grid_eq(i-1+ii,k-1+kk,:)*weight_square(ii,kk)
                            curl_b_eq = curl_b_eq + curl_b_grid_eq(i-1+ii,k-1+kk,:)*weight_square(ii,kk)
                            psi_par   = psi_par  + psi(i-1+ii,k-1+kk)*weight_square(ii,kk)
                        end do
                    end do                
                    
                end if
                
               
                R          = marker(p)%X(1)
                v_parallel = marker(p)%v_para
                mu         = marker(p)%mu
                E          = 0.5*m*v_parallel**2 + mu*abs_B
                v          = sqrt(2*E/m)
                P_phi      = m*v_parallel*R*(B(2)/abs_B) - Zh*psi_par
                Lambda     = mu*B0/E
                
                if(BAS) then
                    v_banos    = mu/2/Zh*dot_product(b_unit,curl_b)
                    v_banos_eq = mu/2/Zh*dot_product(b_unit_eq,curl_b_eq)
                else
                    v_banos    = 0
                    v_banos_eq = 0
                end if
                
                B_star          = B + m_div_Zh*v_parallel*curl_b + m_div_Zh*curl_v_E
                B_star_parallel = dot_product(B_star,b_unit)
                E_star          = delta_E - m_div_Zh*(mu*grad_B/m + 0.5*grad_v_E_2 + v_parallel*pb_unitpt + pv_Ept)
        
                B_star_eq          = B_eq + m_div_Zh*v_parallel*curl_b_eq
                B_star_parallel_eq = dot_product(B_star_eq,b_unit_eq)
                E_star_eq          = -m_div_Zh*(mu*grad_B_eq/m)
 
                
                RHS_X    =  1.0/B_star_parallel*((v_parallel+v_banos)*B_star + cross_product(E_star,b_unit))
                RHS_V    =  Zh_div_m/B_star_parallel*dot_product(B_star,E_star)
                RHS_X_eq =  1.0/B_star_parallel_eq*((v_parallel+v_banos_eq)*B_star_eq + cross_product(E_star_eq,b_unit_eq))
                RHS_V_eq =  Zh_div_m/B_star_parallel_eq*dot_product(B_star_eq,E_star_eq)
    !delta-f part
                dXdt1          = RHS_X - RHS_X_eq 
                dv_paralleldt1 = RHS_V - RHS_V_eq
                
                
                v_d = 1.0/(Zh*B_star_parallel)*(m*v_parallel**2*curl_b + mu*cross_product(b_unit,grad_B))
        
    
                grad_P_phi(1) =   m*v_parallel/abs_B*grad_RB_phi(1) - m*v_parallel*R*B(2)/abs_B**2*grad_B(1) + Zh*B(3)*R
                grad_P_phi(2) =   m*v_parallel/abs_B*grad_RB_phi(2) - m*v_parallel*R*B(2)/abs_B**2*grad_B(2)
                grad_P_phi(3) =   m*v_parallel/abs_B*grad_RB_phi(3) - m*v_parallel*R*B(2)/abs_B**2*grad_B(3) - Zh*B(1)*R
                dP_phidv_parallel = m*R*B(2)/abs_B
        
        
                dP_phidt = dot_product(dXdt1,grad_P_phi) + dv_paralleldt1*dP_phidv_parallel
                dEdt     = m*v_parallel*RHS_V + mu*dot_product(RHS_X,grad_B) + mu*pabs_Bpt
!                dEdt     = Zh*dot_product(v_d + (v_parallel+v_banos)*B/B_star_parallel + m_div_Zh*v_parallel/B_star_parallel*curl_v_E,E_star + m_div_Zh*mu*grad_B) + mu*pabs_Bpt
                
                if(FGP) then
                    dLambdadt = -Lambda**2*dEdt
                end if
                
                if(v_parallel>0) then
                    sgn_v_parallel = +1
                else
                    sgn_v_parallel = -1
                end if
      
                if(ORBIT_AVERAGE_METHOD == 1) then
                    if((1-Lambda) > 0) then
                        bracket_psi = - P_phi/Zh + m_div_Zh*sgn_v_parallel*v*R0*sqrt(1-Lambda)
                    else
                        bracket_psi = - P_phi/Zh
                    end if
                else if(ORBIT_AVERAGE_METHOD == 2) then
                    bracket_psi = - P_phi/Zh
                end if
                
                    
    
        !Finite width of Gaussian distribution of Pitch angle  : f(Lambda) ~ exp[(Lambda-Lambda_0)^2/Delta_Lambda^2]
    
    
                if(DFR_type == 1) then
                    if(FGP) then
                        f0 = c_f*1.0/(v**3+v_c**3)*(1+erf((v_0-v)/Delta_v))*exp(-(bracket_psi/(c_1*Delta_psi)))*exp(-((Lambda-Lambda_0)/Delta_Lambda)**2)
                    else
                        f0 = c_f*1.0/(v**3+v_c**3)*(1+erf((v_0-v)/Delta_v))*exp(-(bracket_psi/(c_1*Delta_psi)))
                    end if
                else if(DFR_type == 2) then
                    bracket_psi_norm = (bracket_psi - psmin)/(psmax-psmin)
                    f0 = c_f*1.0/(v**3+v_c**3)*(1+erf((v_0-v)/Delta_v))*exp(-Delta_n/L_n*tanh((bracket_psi_norm - psi_n)/Delta_n))
                end if

                if(DFR_type == 1) then
                    df0dP_phi =  1.0/(Zh*c_1*Delta_psi)*f0
                else if(DFR_type == 2) then
                    df0dP_phi =  1.0/(Zh*L_n*Delta_psi)*sech((bracket_psi_norm - psi_n)/Delta_n)**2*f0
                end if
        
                if(ORBIT_AVERAGE_METHOD == 1) then
                    if((1-Lambda) > 0) then
                        if(DFR_type == 1) then
                            df0dE = -3/m*v/(v**3+ v_c**3)*f0 - exp(-((v_0-v)/Delta_v)**2)/(Delta_v*sqrt(PI*E*m/2)*(1+erf((v_0-v)/Delta_v)))*f0 - sgn_v_parallel*R0/(Zh*c_1*Delta_psi*sqrt(v*v-2*mu*B0/m))*f0
                        else if(DFR_type == 2) then
                            df0dE = -3/m*v/(v**3+ v_c**3)*f0 - exp(-((v_0-v)/Delta_v)**2)/(Delta_v*sqrt(PI*E*m/2)*(1+erf((v_0-v)/Delta_v)))*f0 - sgn_v_parallel*R0/(Zh*c_1*Delta_psi*sqrt(v*v-2*mu*B0/m))*sech((bracket_psi_norm - psi_n)/Delta_n)**2*f0
                        end if
                    else
                        df0dE = -3/m*v/(v**3+ v_c**3)*f0 - exp(-((v_0-v)/Delta_v)**2)/(Delta_v*sqrt(PI*E*m/2)*(1+erf((v_0-v)/Delta_v)))*f0
                    end if
                else if(ORBIT_AVERAGE_METHOD == 2) then
                    df0dE = -3/m*v/(v**3+ v_c**3)*f0 - exp(-((v_0-v)/Delta_v)**2)/(Delta_v*sqrt(PI*E*m/2)*(1+erf((v_0-v)/Delta_v)))*f0
                end if
                
                    
    
                if(FGP) then
                    df0dLambda =   -2*(Lambda-Lambda_0)/Delta_Lambda**2*f0
                end if
    
                if(FGP) then
                    df0dt = dP_phidt*df0dP_phi + dEdt*df0dE + dLambdadt*df0dLambda
                else
                    df0dt = dP_phidt*df0dP_phi + dEdt*df0dE
                end if
        
                RHS_W    = -1.0/marker(p)%g*df0dt
                
                if(abs(RHS_W) > 1E5) then
!                    write(*,*)"ahahahahaha",dP_phidt,df0dP_phi,dEdt,df0dE,dLambdadt,df0dLambda,bracket_psi,f0,v,c_f*1.0/(v**3+v_c*v_c*v_c),exp(-(bracket_psi/(c_1*Delta_psi))),Delta_psi,bracket_psi/(c_1*Delta_psi)
                    write(*,*)"ahahahahaha",f0,v,exp(-(bracket_psi/(c_1*Delta_psi))),bracket_psi,P_phi/Zh, m_div_Zh*sgn_v_parallel*v*R0*sqrt(1-Lambda),m*v_parallel*R*(B(2)/abs_B),Zh*psi_par
                    write(*,*)"gagagagga",i,j,k,marker(p)%X(1),marker(p)%X(2),marker(p)%X(3)

                    exit 
                endif
                
!                RHS_W    = m*v_parallel*RHS_V + mu*dot_product(RHS_X,grad_B)
        
        
                RHS_X(2) = RHS_X(2)/R
    
                marker(p)%X          = Xold + timestep*RHS_X
                marker(p)%v_para     = Vold + timestep*RHS_V
                marker(p)%w          = Wold + timestep*RHS_W
                Xa    = Xa + coeff*dt*RHS_X
                Va    = Va + coeff*dt*RHS_V
                Wa    = Wa + coeff*dt*RHS_W
                time  = told + timestep
                
                if(psi_par >= psmax_clt_k) then
                    marker(p)%OUT_OF_BOUNDARY = .true.
                    exit
                end if
                
            end do
            
            
            if(marker(p)%OUT_OF_BOUNDARY) then
                marker(p)%X          = Xold
                marker(p)%v_para     = Vold
                marker(p)%w          = Wold
                marker(p)%X(3)       =  - marker(p)%X(3)
            else
                marker(p)%X          = Xa
                marker(p)%v_para     = Va
                marker(p)%w          = Wa
            end if

!boundary condition
!		    do while(marker(p)%X(2) >= 2*PI)
!                marker(p)%X(2) = marker(p)%X(2) - 2*PI
!            end do
!            do while(marker(p)%X(2) < 0)
!                marker(p)%X(2) = marker(p)%X(2) + 2*PI
!            end do
            
            time = told
        end do
        
!        do p=1,N
!            if(marker(p)%id >N_initia*64) write(*,*)"333333hehehehehehehehe : ",marker(p)%id,time
!        end do
        
        if(RECYCLE_METHOD == 1) call recycle_particle
        if(RECYCLE_METHOD == 2) call remove_particle
!        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!        if(nrank == 0) write(*,*)"222222222222222222"
!        do p=1,N
!            if(marker(p)%id >N_initia*64) write(*,*)"hehehehehehehehehehehe : ",marker(p)%id,time
!        end do
        
        call update
!        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!        if(nrank == 0) write(*,*)"333333333333333333"
!        do p=1,N
!            if(marker(p)%id >N_initia*64) write(*,*)"121212121heheheheheheh : ",marker(p)%id,time,nrank
!        end do
        
        !output number of particle vs time
        call MPI_REDUCE(N,N_total,1,MPI_INTEGER,MPI_SUM, 0, MPI_COMM_WORLD, IERROR)
        if(nrank == 0 .and. (.not. COF)) write(*,*)"number of particle : ",N_total," at time = ",time
        if(N > N_initia+deltaN) write(*,*)"more than maxima of N in per process!!!!!   N: ",N," N_initia+deltaN :", N_initia+deltaN," nrank :",nrank
        

                    
    end subroutine particle_orbit_polarization
    
    
    subroutine calc_current_old
        use distribution
        use var
        implicit none
        include 'mpif.h'
        real*8, dimension(2) :: V_cell
        real*8               :: d1f2, d1fc
        real*8 :: f_jm1,f_jm2,f_j,f_jp1,f_jp2,x_jm1,x_j,x_jp1,coeff_a,coeff_b,coeff_c,coeff_d
        integer ix_first, ix_last,iz_first, iz_last,iy_first,iy_last
        integer jx,jy,jz

!  define statement functions
!  d1f2= d f / dx  with second-order accuracy central difference
        
        d1f2(f_jm1,f_j,f_jp1,x_jm1,x_j,x_jp1)= &
        ((x_jm1-x_j)/(x_jp1-x_j)*(f_jp1-f_j) &
        -(x_jp1-x_j)/(x_jm1-x_j)*(f_jm1-f_j))/(x_jm1-x_jp1)  
        
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(f_jm2,f_jm1,f_j,f_jp1,f_jp2,coeff_a,coeff_b,coeff_c,coeff_d)= &
       coeff_a*(f_jp1-f_j)+coeff_b*(f_j-f_jm1)+coeff_c*(f_jp2-f_j)+coeff_d*(f_j-f_jm2)
        
        
        ix_first=1
        ix_last=mx
        iz_first=1
        iz_last=mz
        iy_first=1
        iy_last=my
        
        J_h = 0
        M_h = 0
        do p=1,N
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
            
            if(psi(i,k)<psmax .and. (psi(i+1,k)>psmax .or. psi(i,k+1)>psmax .or. psi(i+1,k+1)>psmax)) then
                weight_line(1,1) = 0.0
                weight_line(3,1) = 0.0
            endif
             
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
            
            B = 0
            b_unit  = 0
            grad_B  = 0
            curl_b  = 0
            do ii=1,2
                do jj=1,2
                    do kk=1,2!note that B in CLT code is Bx(R,Z,phi)
                        B      = B + x(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
                        b_unit = b_unit + b_unit_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                        grad_B = grad_B + grad_B_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                        curl_b = curl_b + curl_b_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                    end do
                end do
            end do
            abs_B = sqrt(dot_product(B,B))
                
            R          = marker(p)%X(1)
            v_parallel = marker(p)%v_para
            mu         = marker(p)%mu
                
            !due to Jacobian ~ B_star_parallel, so denominator in v_d should be replaced by abs_B
            v_d        = 1.0/(Zh*abs_B)*(m*v_parallel**2*curl_b + mu*cross_product(b_unit,grad_B))
            
!polarization drift related variables                
            if(VEB) then
                curl_v_E = 0
                grad_v_E_2 = 0
                pv_Ept = 0
                do ii=1,2
                    do jj=1,2
                        do kk=1,2!note that B in CLT code is Bx(R,Z,phi)
                            curl_v_E   = curl_v_E   + curl_v_E_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                            grad_v_E_2 = grad_v_E_2 + grad_v_E_2_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                    if(POL) pv_Ept     = pv_Ept     + pv_Ept_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                        end do
                    end do
                end do
                
                !add v_EB induced by e*B flow in Lagrangian
                v_d = v_d + m_div_Zh/abs_B*(v_parallel*curl_v_E + 0.5*cross_product(b_unit,grad_v_E_2))
                
                !add polarization current
                if(POL) v_d = v_d + m_div_Zh/abs_B*cross_product(b_unit,pv_Ept)
                
            end if

            delta_f    = marker(p)%w*marker(p)%g

            
            
            dR2 = ((xx(i)+dxx+xx(i))/2)**2 - ((xx(i)+xx(i)-dxx)/2)**2
            V_cell(1) = dyy*dzz*dR2/2
            dR2 = ((xx(i+1)+dxx+xx(i+1))/2)**2 - ((xx(i+1)+xx(i+1)-dxx)/2)**2
            V_cell(2) = dyy*dzz*dR2/2
            
            do ii=1,2
                do jj=1,2
                    do kk=1,2!note that field in CLT code is Bx(R,Z,phi)
                       if(psi(i-1+ii,k-1+kk) <= (psmin_clt_k + 1.0*(psmax_clt_k - psmin_clt_k))) then
                            J_h(i-1+ii,k-1+kk,j-1+jj,:) = J_h(i-1+ii,k-1+kk,j-1+jj,:) + Zh*(v_parallel*b_unit + v_d)*delta_f/V_cell(ii)*weight_cubic(ii,jj,kk)
                            M_h(i-1+ii,k-1+kk,j-1+jj,:) = M_h(i-1+ii,k-1+kk,j-1+jj,:) + (-mu*b_unit)*delta_f/V_cell(ii)*weight_cubic(ii,jj,kk)
                       endif
                    end do
                end do
            end do
        end do
        
!update J_h and M_h on boundary of MPI sub-domian
        call update_vector_sub_domain_bounnary(J_h)
        call update_vector_sub_domain_bounnary(M_h)

!update J_h(:,:,2,:) and M_h(:,:,2,:) 
        call mpi_transfersm(M_h,3)
        call mpi_transfersm(J_h,3)

!magnezation current 
!heavy smooth!!!!!
!        call possion_solver_3D(M_h(:,:,:,1))
!        call possion_solver_3D(M_h(:,:,:,2))
!        call possion_solver_3D(M_h(:,:,:,3))

        
        do s=1,10
            call bndry3_ex(M_h, 0)
            call smooth26_old(M_h, 0.8d0)
        end do
        
        do s=1,3
            do jy=iy_first,iy_last
                do jz=iz_first+2,iz_last-2
                    do jx=ix_first+2,ix_last-2
                        if(FOUR_TH) then
                            dM_hdR(jx,jz,jy,s) = d1fc(M_h(jx-2,jz,jy,s),M_h(jx-1,jz,jy,s),M_h(jx,jz,jy,s),M_h(jx+1,jz,jy,s),M_h(jx+2,jz,jy,s),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
                            dM_hdZ(jx,jz,jy,s) = d1fc(M_h(jx,jz-2,jy,s),M_h(jx,jz-1,jy,s),M_h(jx,jz,jy,s),M_h(jx,jz+1,jy,s),M_h(jx,jz+2,jy,s),az1(jz),bz1(jz),cz1(jz),dz1(jz))
                        else
                            dM_hdR(jx,jz,jy,s) = d1f2(M_h(jx-1,jz,jy,s),M_h(jx,jz,jy,s),M_h(jx+1,jz,jy,s),xx(jx-1),xx(jx),xx(jx+1))
                            dM_hdZ(jx,jz,jy,s) = d1f2(M_h(jx,jz-1,jy,s),M_h(jx,jz,jy,s),M_h(jx,jz+1,jy,s),zz(jz-1),zz(jz),zz(jz+1))
                        end if
                  end do
              end do
            end do
        end do
        
        do s=1,3
            do jz=iz_first,iz_last
                do jx=ix_first,ix_last
                    do jy=iy_first+2,iy_last-2
                        if(FOUR_TH) then
                            dM_hdphi(jx,jz,jy,s) = d1fc(M_h(jx,jz,jy-2,s),M_h(jx,jz,jy-1,s),M_h(jx,jz,jy,s),M_h(jx,jz,jy+1,s),M_h(jx,jz,jy+2,s),ay1(jy),by1(jy),cy1(jy),dy1(jy))/xx(jx)
                        else
                            dM_hdphi(jx,jz,jy,s) = d1f2(M_h(jx,jz,jy-1,s),M_h(jx,jz,jy,s),M_h(jx,jz,jy+1,s),yy(jy-1),yy(jy),yy(jy+1))/xx(jx)
                        end if
                    end do
                end do
            end do
        end do
              
        J_M(:,:,:,1)   = dM_hdphi(:,:,:,3) - dM_hdZ(:,:,:,2)
        J_M(:,:,:,2)   = dM_hdZ(:,:,:,1)   - dM_hdR(:,:,:,3)
        do jz=iz_first,iz_last
            do jx=ix_first,ix_last
                do jy=iy_first,iy_last
                    J_M(jx,jz,jy,3)   = M_h(jx,jz,jy,2)/xx(jx) + dM_hdR(jx,jz,jy,2) - dM_hdphi(jx,jz,jy,1)
                end do
            end do
        end do
        
!energtic particle current is composed of guiding center current and magnezation current
        J_h = J_h + J_M
        
        if(ADI) then
            !substract adiabatic term : xi dot grad f0!!!!
            !calc plasma displayment 
            xi = xi + x(:,:,:,3:5)*dt
        
            do i=3,mx-2
                do j=3,mz-2
                    do k=3,my-2
                        J_h(i,j,k,1) = J_h(i,j,k,1) + dot_product(xi(i,j,k,:),grad_Jh0_R(i,j,k,:))
                        J_h(i,j,k,2) = J_h(i,j,k,2) + dot_product(xi(i,j,k,:),grad_Jh0_phi(i,j,k,:))
                        J_h(i,j,k,3) = J_h(i,j,k,3) + dot_product(xi(i,j,k,:),grad_Jh0_Z(i,j,k,:))
                    end do
                end do
            end do
        end if
        
        call possion_solver_3D(J_h(:,:,:,1),J_h_old(:,:,:,1))
        call possion_solver_3D(J_h(:,:,:,2),J_h_old(:,:,:,2))
        call possion_solver_3D(J_h(:,:,:,3),J_h_old(:,:,:,3))
        
        J_h_old = J_h
        
!3 times smooth!
        do s=1,3
            call bndry3_ex(J_h,0)
            call smooth26_old(J_h, coeff_smooth)
        end do
        
        
    end subroutine calc_current_old
    
    
    subroutine smooth26_old(data_in, coefficient_smooth)
        use var
        implicit none
        include 'mpif.h'
        
        real*8, dimension(mx,mz,my,3), intent(inout) :: data_in
        real*8, intent(in) :: coefficient_smooth
        
        real*8, dimension(3,3,3) :: weight
        real*8 :: distance, weight_total        
        real*8, dimension(mx,mz,my,3) :: data_smooth

        call mpi_transfersm(data_in,3)
        
        !points on boundary should not be smoothed!!
        data_smooth = data_in
        
        do i=3,mx-2
            weight_total = 0.0
            do ii=1,3
                do jj=1,3
                    do kk=1,3
                        if(ii /= 2 .or. jj /=2 .or. kk /= 2) then
                            distance = sqrt(abs(ii-2)*dxx**2 + abs(jj-2)*dzz**2 + abs(kk-2)*(xx(i)*dyy)**2)
                            weight(ii,jj,kk) = 1.0/distance
                            weight_total = weight_total + weight(ii,jj,kk)
                        end if
                    end do
                end do
            end do
            
            weight = weight/weight_total*coefficient_smooth
            weight(2,2,2) = (1.0 - coefficient_smooth)
            
            do j=3,mz-2
                do k=3,my-2
                    if(psi(i,j)<psmax .and. psi(i+1,j)<psmax .and. psi(i,j+1)<psmax .and. psi(i+1,j+1)<psmax .and. psi(i-1,j)<psmax .and. psi(i,j-1)<psmax .and. psi(i-1,j-1)<psmax .and. psi(i+1,j-1)<psmax .and. psi(i-1,j+1)<psmax) then
                        data_smooth(i,j,k,:) = 0.0
                        do ii=1,3
                            do jj=1,3
                                do kk=1,3
                                    data_smooth(i,j,k,:) = data_smooth(i,j,k,:) + weight(ii,jj,kk)*data_in(i+ii-2,j+jj-2,k+kk-2,:)
                                end do
                            end do
                        end do
                    end if
                end do
            end do
        end do
        
        data_in = data_smooth
        
        call mpi_transfersm(data_in,3)
        
    end subroutine smooth26_old
    
    subroutine calc_grad_cuurent_initia
        USE DECLARE
        implicit none
        real*8,dimension(mx,mz,my,3) :: dJh0dR, dJh0dphi, dJh0dZ
        real*8               :: d1f2, d1fc
        real*8 :: f_jm1,f_jm2,f_j,f_jp1,f_jp2,x_jm1,x_j,x_jp1,coeff_a,coeff_b,coeff_c,coeff_d
        include 'mpif.h'
      
!  define statement functions
!  d1f2= d f / dx  with second-order accuracy central difference
        
        d1f2(f_jm1,f_j,f_jp1,x_jm1,x_j,x_jp1)= &
        ((x_jm1-x_j)/(x_jp1-x_j)*(f_jp1-f_j) &
        -(x_jp1-x_j)/(x_jm1-x_j)*(f_jm1-f_j))/(x_jm1-x_jp1)  
        
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(f_jm2,f_jm1,f_j,f_jp1,f_jp2,coeff_a,coeff_b,coeff_c,coeff_d)= &
       coeff_a*(f_jp1-f_j)+coeff_b*(f_j-f_jm1)+coeff_c*(f_jp2-f_j)+coeff_d*(f_j-f_jm2)
        
        call mpi_transfersm(J_h0,3)
      
        do m=1,3
            do jy=iy_first,iy_last
                do jz=iz_first+2,iz_last-2
                    do jx=ix_first+2,ix_last-2
                        if(FOUR_TH) then
                            dJh0dR(jx,jz,jy,m) = d1fc(J_h0(jx-2,jz,jy,m),J_h0(jx-1,jz,jy,m),J_h0(jx,jz,jy,m),J_h0(jx+1,jz,jy,m),J_h0(jx+2,jz,jy,m),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
                            dJh0dZ(jx,jz,jy,m) = d1fc(J_h0(jx,jz-2,jy,m),J_h0(jx,jz-1,jy,m),J_h0(jx,jz,jy,m),J_h0(jx,jz+1,jy,m),J_h0(jx,jz+2,jy,m),az1(jz),bz1(jz),cz1(jz),dz1(jz))
                        else
                            dJh0dR(jx,jz,jy,m) = d1f2(J_h0(jx-1,jz,jy,m),J_h0(jx,jz,jy,m),J_h0(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))
                            dJh0dZ(jx,jz,jy,m) = d1f2(J_h0(jx,jz-1,jy,m),J_h0(jx,jz,jy,m),J_h0(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))
                        end if
                    end do
                end do
            end do
        end do
        
        do m=1,3
            do jz=iz_first,iz_last
                do jx=ix_first,ix_last
                    do jy=iy_first+2,iy_last-2
                        if(FOUR_TH) then
                            dJh0dphi(jx,jz,jy,m) = d1fc(J_h0(jx,jz,jy-2,m),J_h0(jx,jz,jy-1,m),J_h0(jx,jz,jy,m),J_h0(jx,jz,jy+1,m),J_h0(jx,jz,jy+2,m),ay1(jy),by1(jy),cy1(jy),dy1(jy))
                        else
                            dJh0dphi(jx,jz,jy,m) = d1f2(J_h0(jx,jz,jy-1,m),J_h0(jx,jz,jy,m),J_h0(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))
                        end if
                    end do
                end do
            end do
        end do
        
      
        grad_Jh0_R(:,:,:,1)   = dJh0dR(:,:,:,1)
        grad_Jh0_R(:,:,:,2)   = dJh0dphi(:,:,:,1)
        grad_Jh0_R(:,:,:,3)   = dJh0dZ(:,:,:,1)
      
        do jx=1,mx
            grad_Jh0_R(jx,:,:,2) = grad_Jh0_R(jx,:,:,2)/xx(jx)
        end do
      
        grad_Jh0_phi(:,:,:,1)   = dJh0dR(:,:,:,2)
        grad_Jh0_phi(:,:,:,2)   = dJh0dphi(:,:,:,2)
        grad_Jh0_phi(:,:,:,3)   = dJh0dZ(:,:,:,2)
      
        do jx=1,mx
            grad_Jh0_phi(jx,:,:,2) = grad_Jh0_phi(jx,:,:,2)/xx(jx)
        end do
      
        grad_Jh0_Z(:,:,:,1)   = dJh0dR(:,:,:,3)
        grad_Jh0_Z(:,:,:,2)   = dJh0dphi(:,:,:,3)
        grad_Jh0_Z(:,:,:,3)   = dJh0dZ(:,:,:,3)
      
        do jx=1,mx
            grad_Jh0_Z(jx,:,:,2) = grad_Jh0_Z(jx,:,:,2)/xx(jx)
        end do
      
        call mpi_transfersm(grad_Jh0_R,3)
        call mpi_transfersm(grad_Jh0_phi,3)
        call mpi_transfersm(grad_Jh0_Z,3)
                  
                  
        return
    end

    
!************************************************************************************************
!                                                                                               *
!          Input  : import scalar variable xi_psi(psi) from NOVA code                           *
!          Output : E_phi1(R,Z,phi)                                                             *
!          Method : superposition of all poloidal compmonet for xi_psi_m                        *
!                   E_phi1 = xi_psi*omega/R and xi_psi = sum xi_psi_m*cos(n*phi + m*theta)      *
!                                                                                               *
!     p.s. #1 here we omit omega in expression for E_phi1 because of a constant                 *
!          #2 igrid=1 should be set in NOVA because we use xi_psi uniformly in sqrt(psi)        *
!                                                                                               *
!************************************************************************************************    
    subroutine xi_from_NOVA(nosurf,m_tot)
        use DECLARE
        implicit none
        include 'mpif.h'
        integer status(mpi_status_size)
        
        integer, intent(in) :: nosurf,m_tot
        integer :: km,ix,r_index
        integer :: i,j,k
        real*8  :: drt_norm_psi,norm_psi,E_phi_NOVA(mxt,mzt,my),weight
        real*8, parameter :: small_coeff = 1.0e-6
        integer,allocatable :: mm(:)
        real*8,allocatable  :: zpsi(:,:)
        allocate(mm(m_tot))
        allocate(zpsi(nosurf,m_tot))
        
        open(unit=777,file='outzpsi_CLT')
        do km=1,m_tot
            read(777,703)mm(km),(zpsi(ix,km),ix=1,nosurf)
        end do
703     format(i4,<nosurf>e13.5)
        
        drt_norm_psi = 1.0/(nosurf - 1)
        E_phi_NOVA = 0;
        do i=1,mxt
            do j=1,mzt
                if(pst(i,j) < psmax) then
                    norm_psi = (pst(i,j) - psmin)/(psmax - psmin)
                    r_index  = floor(sqrt(norm_psi)/drt_norm_psi)
				    weight   = (sqrt(norm_psi) - r_index*drt_norm_psi)/drt_norm_psi
                    do m=1,m_tot
                        do k=1,my
                            E_phi_NOVA(i,j,k) = E_phi_NOVA(i,j,k) + ((1-weight)*zpsi(r_index,m) + weight*zpsi(r_index+1,m))*cos(n_filt*yy(k)+mm(m)*tht(i,j))
                        end do
                    end do
                end if
                
            end do
            E_phi_NOVA(i,:,:) = E_phi_NOVA(i,:,:)/xxt(i)
        end do
        
        E_phi_NOVA = small_coeff*E_phi_NOVA
        
        do jz=3,mz-2
            do jx=3,mx-2
                E_phi1(jx,jz,:) = E_phi_NOVA(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2,:)
            end do
        end do
        
        
        call mpi_transfersm(E_phi1,1)
        
        if(nrank==0) then
            open(unit=778,file='E_phi_NOVA.dat')
            do j=1,mzt
                write(778,704)(E_phi_NOVA(i,j,3),i=1,mxt)
            end do
704         format(<mxt>e13.5)
            
        end if


        
        return
        
    end subroutine xi_from_NOVA
    
!************************************************************************************************
!                          N points average on gyro ring                                        *
!************************************************************************************************    
    subroutine gyro_average(N_gyro)
        use var
        use distribution
        implicit none 
        include 'mpif.h'
        integer status(mpi_status_size)
        
        integer, intent(in) :: N_gyro
        real*8 :: rho_h,alpha_h !larmor_radius and gyro_angle 
        real*8 :: e1(3), e2(3), B_GC(3), b_unit_GC(3), rho(3), abs_B_GC, X_GC_vector(3), X_ring_vector(3)
        real*8 :: psi_par_ring
        real*8, allocatable :: X_ring(:,:)
        
        allocate(X_ring(N_gyro,3))
        
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
                
        if(psi(i,k)<psmax .and. (psi(i+1,k)>psmax .or. psi(i,k+1)>psmax .or. psi(i+1,k+1)>psmax)) then
            weight_line(1,1) = 0.0
            weight_line(3,1) = 0.0
        endif
             
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
             
                
        B_GC = 0

        do ii=1,2
            do jj=1,2
                do kk=1,2!note that B in CLT code is Bx(R,Z,phi)
                    B_GC       = B_GC + x(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
                end do
            end do
        end do
        abs_B_GC    = sqrt(dot_product(B_GC,B_GC))
        b_unit_GC   = B_GC/abs_B_GC
        
        e1(1) = - b_unit_GC(3)
        e1(2) =   0
        e1(3) =   b_unit_GC(1)
        
        e2 = cross_product(b_unit_GC, e1)
        
        !normalization e1 and e2
        e1 = e1/sqrt(dot_product(e1,e1))
        e2 = e2/sqrt(dot_product(e2,e2))

        rho_h = sqrt(2*m*marker(p)%mu/abs_B_GC)/Zh
        
        !get vector form of X_GC
        X_GC_vector = marker(p)%X
        X_GC_vector(2) = 0;

        
        if (time < 2*dt .and. nrank == 47) write(*,*)"*********gyro average process*********"
        if (time < 2*dt .and. nrank == 47) write(*,*)"rho_h : ",rho_h
        if (time < 2*dt .and. nrank == 47) write(*,*)marker(p)%X(1),marker(p)%X(2),marker(p)%X(3)  
        do n=1,N_gyro
            alpha_h = 2*pi/N_gyro*(n-1)
            rho = rho_h*cos(alpha_h)*e1 + rho_h*sin(alpha_h)*e2
            X_ring_vector = X_GC_vector + rho
            !X_ring_vector --> X_ring
            X_ring(n,1)  = sqrt(X_ring_vector(1)**2 + X_ring_vector(2)**2)
            X_ring(n,2)  = marker(p)%X(2) + atan2(X_ring_vector(2), X_ring_vector(1))
            X_ring(n,3)  = X_ring_vector(3)
            if (time < 2*dt .and. nrank == 47) write(*,*)X_ring(n,1),X_ring(n,2),X_ring(n,3) 
        end do
        if (time < 2*dt .and. nrank == 47) write(*,*)"*********gyro average process*********"
        
        !10 variables should be gyro-averaged
        B = 0
        grad_B  = 0
        curl_b  = 0
        delta_E = 0
        grad_RB_phi = 0
        pBpt    = 0
        
        B_eq = 0
        grad_B_eq = 0
        curl_b_eq = 0
        psi_par   = 0

        if(VEB) then
            curl_v_E   = 0
            grad_v_E_2 = 0
            pv_Ept     = 0
        end if
        
        do n=1,N_gyro
            
            i = floor((X_ring(n,1) - myxmin)/dxx) + 3
            j = floor((X_ring(n,2) - myymin)/dyy) + 3
            k = floor((X_ring(n,3) - myzmin)/dzz) + 3
                
            if(i>=mx .or. i<1 .or. j>=my .or. j<1 .or. k>=mz .or. k<1) then
                write(*,*)"zzzzzzzzzzzzzzzzzaaaaaaaaaaaaaa : ",i,j,k,marker(p)%X,marker(p)%id,nrank,time
                write(*,*)"nnnnnnnnnnnnnnnnnbbbbbbbbbbbbbb : ",i,j,k,marker_initia(p)%X,marker_initia(p)%v_para,marker_initia(p)%v_perp,marker_initia(p)%id,nrank,time
            end if
             
            R   = xx(i)
            phi = yy(j)
            Z   = zz(k)
             
            dR2 = xx(i+1)**2 - xx(i)**2
             
            weight_line(1,1) = (X_ring(n,1)**2 - R**2)/dR2
            weight_line(2,1) = (X_ring(n,2) - phi)/dyy
            weight_line(3,1) = (X_ring(n,3) - Z)/dzz
                
            if(psi(i,k)<psmax .and. (psi(i+1,k)>psmax .or. psi(i,k+1)>psmax .or. psi(i+1,k+1)>psmax)) then
                weight_line(1,1) = 0.0
                weight_line(3,1) = 0.0
            endif
             
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
            
            
            weight_line(1,1) = (X_ring(n,1) - R)/dxx
            weight_line(2,1) = (X_ring(n,3) - Z)/dzz
                
            if(psi(i,k)<psmax .and. (psi(i+1,k)>psmax .or. psi(i,k+1)>psmax .or. psi(i+1,k+1)>psmax)) then
                weight_line(1,1) = 0.0
                weight_line(2,1) = 0.0
            endif
             
            weight_line(1,2) = 1.0 - weight_line(1,1)
            weight_line(2,2) = 1.0 - weight_line(2,1)

             
            do ii=1,2
                do jj=1,2
                    weight_square(ii,jj) = weight_line(1,3-ii)*weight_line(2,3-jj)
                end do
            end do
            
            !if particle goes out of boundary, break!!!
            psi_par_ring = 0
            do ii=1,2
                do kk=1,2
                    psi_par_ring   = psi_par_ring  + psi(i-1+ii,k-1+kk)*weight_square(ii,kk)
                end do
            end do
            if(psi_par_ring >= psmax_clt_k) then
                marker(p)%OUT_OF_BOUNDARY = .true.
                exit
            else
                psi_par = psi_par + psi_par_ring
            end if
            
            do ii=1,2
                do kk=1,2
                    B_eq      = B_eq + xint(i-1+ii,k-1+kk,6:8)*weight_square(ii,kk)
                    grad_B_eq = grad_B_eq + grad_B_grid_eq(i-1+ii,k-1+kk,:)*weight_square(ii,kk)
                    curl_b_eq = curl_b_eq + curl_b_grid_eq(i-1+ii,k-1+kk,:)*weight_square(ii,kk)
                end do
            end do     

            do ii=1,2
                do jj=1,2
                    do kk=1,2!note that B in CLT code is Bx(R,Z,phi)
                        B       = B + x(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
                        grad_B  = grad_B  + grad_B_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                        curl_b  = curl_b  + curl_b_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                        delta_E = delta_E + Ef(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                        pBpt    = pBpt    + xdif(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
                        grad_RB_phi(1) = grad_RB_phi(1) + (x(i-1+ii,k-1+kk,j-1+jj,7)+xx(i-1+ii)*xr(i-1+ii,k-1+kk,j-1+jj,7))*weight_cubic(ii,jj,kk)
                        grad_RB_phi(2) = grad_RB_phi(2) + (                                     xy(i-1+ii,k-1+kk,j-1+jj,7))*weight_cubic(ii,jj,kk)
                        grad_RB_phi(3) = grad_RB_phi(3) + (                          xx(i-1+ii)*xz(i-1+ii,k-1+kk,j-1+jj,7))*weight_cubic(ii,jj,kk)
                    end do
                end do
            end do
            
            if(VEB) then
                do ii=1,2
                    do jj=1,2
                        do kk=1,2!note that B in CLT code is Bx(R,Z,phi)
                            curl_v_E   = curl_v_E   + curl_v_E_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                            grad_v_E_2 = grad_v_E_2 + grad_v_E_2_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                    if(POL) pv_Ept     = pv_Ept     + pv_Ept_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                        end do
                    end do
                end do
            end if
        end do
            
        
        B           = B/N_gyro
        grad_B      = grad_B/N_gyro
        curl_b      = curl_b/N_gyro
        delta_E     = delta_E/N_gyro
        grad_RB_phi = grad_RB_phi/N_gyro
        pBpt        = pBpt/N_gyro
        
        abs_B    = sqrt(dot_product(B,B))
        b_unit   = B/abs_B
        pabs_Bpt = dot_product(b_unit,pBpt)
        
        B_eq = B_eq/N_gyro
        grad_B_eq = grad_B_eq/N_gyro
        curl_b_eq = curl_b_eq/N_gyro
        psi_par   = psi_par/N_gyro
        b_unit_eq   = B_eq/sqrt(dot_product(B_eq,B_eq))
        
        if(VEB) pb_unitpt = (pBpt - pabs_Bpt*b_unit)/abs_B  !calc pb_unitpt, which affect equation of dXdt and not affect equation of dv_||dt because of pb_unitpt * B = 0
                              
        return
        
    end subroutine gyro_average

    
    
!************************************************************************************************
!                       N points average on gyro ring on (x,y,z) coordinate                     *
!************************************************************************************************    
    subroutine gyro_average_xyz(N_gyro)
        use var
        use distribution
        implicit none 
        include 'mpif.h'
        integer status(mpi_status_size)
        
        integer, intent(in) :: N_gyro
        real*8 :: rho_h,alpha_h !larmor_radius and gyro_angle 
        real*8 :: e1(3), e2(3), B_GC(3), b_unit_GC(3), b_unit_GC_xyz(3), rho(3), abs_B_GC, X_GC_xyz(3), X_ring_xyz(3)
        real*8 :: psi_par_ring
        real*8, allocatable :: X_ring(:,:)
        
        allocate(X_ring(N_gyro,3))
        
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
                
        if(psi(i,k)<psmax .and. (psi(i+1,k)>psmax .or. psi(i,k+1)>psmax .or. psi(i+1,k+1)>psmax)) then
            weight_line(1,1) = 0.0
            weight_line(3,1) = 0.0
        endif
             
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
             
                
        B_GC = 0

        do ii=1,2
            do jj=1,2
                do kk=1,2!note that B in CLT code is Bx(R,Z,phi)
                    B_GC       = B_GC + x(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
                end do
            end do
        end do
        abs_B_GC    = sqrt(dot_product(B_GC,B_GC))
        b_unit_GC   = B_GC/abs_B_GC

        !convert b_unit_GC(R,phi,Z) --> b_unit_GC_xyz(x,y,z)
        b_unit_GC_xyz(1) = b_unit_GC(1)*cos(marker(p)%X(2)) - b_unit_GC(2)*sin(marker(p)%X(2))
        b_unit_GC_xyz(2) = b_unit_GC(1)*sin(marker(p)%X(2)) + b_unit_GC(2)*cos(marker(p)%X(2))
        b_unit_GC_xyz(3) = b_unit_GC(3)        
        
        
        e1(1) = - b_unit_GC_xyz(3)
        e1(2) =   0
        e1(3) =   b_unit_GC_xyz(1)
        
        e2 = cross_product(b_unit_GC_xyz, e1)
        
        !normalization e1 and e2
        e1 = e1/sqrt(dot_product(e1,e1))
        e2 = e2/sqrt(dot_product(e2,e2))

        rho_h = sqrt(2*m*marker(p)%mu/abs_B_GC)/Zh
        
        !X_GC_xyz on (x,y,z) grid
        X_GC_xyz(1) = marker(p)%X(1)*cos(marker(p)%X(2))
        X_GC_xyz(2) = marker(p)%X(1)*sin(marker(p)%X(2))
        X_GC_xyz(3) = marker(p)%X(3)
        
        if (time < 2*dt .and. nrank == 47) write(*,*)"*********gyro average process*********"
        if (time < 2*dt .and. nrank == 47) write(*,*)"rho_h : ",rho_h
        if (time < 2*dt .and. nrank == 47) write(*,*)marker(p)%X(1),marker(p)%X(2),marker(p)%X(3)  
        do n=1,N_gyro
            alpha_h = 2*pi/N_gyro*(n-1)
            rho = rho_h*cos(alpha_h)*e1 + rho_h*sin(alpha_h)*e2
            X_ring_xyz = X_GC_xyz + rho
            !X_ring(x,y,z) --> X_ring(R,phi,Z)
            X_ring(n,1) = sqrt(X_ring_xyz(1)**2 + X_ring_xyz(2)**2)
            if(X_ring_xyz(2) >= 0) then
                X_ring(n,2) = atan2(X_ring_xyz(2),X_ring_xyz(1))
            else
                X_ring(n,2) = atan2(X_ring_xyz(2),X_ring_xyz(1)) + 2*PI
            end if
            X_ring(n,3) = X_ring_xyz(3)
            if (time < 2*dt .and. nrank == 47) write(*,*)X_ring(n,1),X_ring(n,2),X_ring(n,3) 
        end do
        if (time < 2*dt .and. nrank == 47) write(*,*)"*********gyro average process*********"
        
        !10 variables should be gyro-averaged
        B = 0
        grad_B  = 0
        curl_b  = 0
        delta_E = 0
        grad_RB_phi = 0
        pBpt    = 0
        
        B_eq = 0
        grad_B_eq = 0
        curl_b_eq = 0
        psi_par   = 0

        if(VEB) then
            curl_v_E   = 0
            grad_v_E_2 = 0
            pv_Ept     = 0
        end if
        
        do n=1,N_gyro
            
            i = floor((X_ring(n,1) - myxmin)/dxx) + 3
            j = floor((X_ring(n,2) - myymin)/dyy) + 3
            k = floor((X_ring(n,3) - myzmin)/dzz) + 3
                
            if(i>=mx .or. i<1 .or. j>=my .or. j<1 .or. k>=mz .or. k<1) then
                write(*,*)"zzzzzzzzzzzzzzzzzaaaaaaaaaaaaaa : ",i,j,k,marker(p)%X,marker(p)%id,nrank,time
                write(*,*)"nnnnnnnnnnnnnnnnnbbbbbbbbbbbbbb : ",i,j,k,marker_initia(p)%X,marker_initia(p)%v_para,marker_initia(p)%v_perp,marker_initia(p)%id,nrank,time
            end if
             
            R   = xx(i)
            phi = yy(j)
            Z   = zz(k)
             
            dR2 = xx(i+1)**2 - xx(i)**2
             
            weight_line(1,1) = (X_ring(n,1)**2 - R**2)/dR2
            weight_line(2,1) = (X_ring(n,2) - phi)/dyy
            weight_line(3,1) = (X_ring(n,3) - Z)/dzz
                
            if(psi(i,k)<psmax .and. (psi(i+1,k)>psmax .or. psi(i,k+1)>psmax .or. psi(i+1,k+1)>psmax)) then
                weight_line(1,1) = 0.0
                weight_line(3,1) = 0.0
            endif
             
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
            
            
            weight_line(1,1) = (X_ring(n,1) - R)/dxx
            weight_line(2,1) = (X_ring(n,3) - Z)/dzz
                
            if(psi(i,k)<psmax .and. (psi(i+1,k)>psmax .or. psi(i,k+1)>psmax .or. psi(i+1,k+1)>psmax)) then
                weight_line(1,1) = 0.0
                weight_line(2,1) = 0.0
            endif
             
            weight_line(1,2) = 1.0 - weight_line(1,1)
            weight_line(2,2) = 1.0 - weight_line(2,1)

             
            do ii=1,2
                do jj=1,2
                    weight_square(ii,jj) = weight_line(1,3-ii)*weight_line(2,3-jj)
                end do
            end do
            
            !if particle goes out of boundary, break!!!
            psi_par_ring = 0
            do ii=1,2
                do kk=1,2
                    psi_par_ring   = psi_par_ring  + psi(i-1+ii,k-1+kk)*weight_square(ii,kk)
                end do
            end do
            if(psi_par_ring >= psmax_clt_k) then
                marker(p)%OUT_OF_BOUNDARY = .true.
                exit
            else
                psi_par = psi_par + psi_par_ring
            end if
            
            do ii=1,2
                do kk=1,2
                    B_eq      = B_eq + xint(i-1+ii,k-1+kk,6:8)*weight_square(ii,kk)
                    grad_B_eq = grad_B_eq + grad_B_grid_eq(i-1+ii,k-1+kk,:)*weight_square(ii,kk)
                    curl_b_eq = curl_b_eq + curl_b_grid_eq(i-1+ii,k-1+kk,:)*weight_square(ii,kk)
                end do
            end do     

            do ii=1,2
                do jj=1,2
                    do kk=1,2!note that B in CLT code is Bx(R,Z,phi)
                        B       = B + x(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
                        grad_B  = grad_B  + grad_B_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                        curl_b  = curl_b  + curl_b_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                        delta_E = delta_E + Ef(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                        pBpt    = pBpt    + xdif(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
                        grad_RB_phi(1) = grad_RB_phi(1) + (x(i-1+ii,k-1+kk,j-1+jj,7)+xx(i-1+ii)*xr(i-1+ii,k-1+kk,j-1+jj,7))*weight_cubic(ii,jj,kk)
                        grad_RB_phi(2) = grad_RB_phi(2) + (                                     xy(i-1+ii,k-1+kk,j-1+jj,7))*weight_cubic(ii,jj,kk)
                        grad_RB_phi(3) = grad_RB_phi(3) + (                          xx(i-1+ii)*xz(i-1+ii,k-1+kk,j-1+jj,7))*weight_cubic(ii,jj,kk)
                    end do
                end do
            end do
            
            if(VEB) then
                do ii=1,2
                    do jj=1,2
                        do kk=1,2!note that B in CLT code is Bx(R,Z,phi)
                            curl_v_E   = curl_v_E   + curl_v_E_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                            grad_v_E_2 = grad_v_E_2 + grad_v_E_2_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                    if(POL) pv_Ept     = pv_Ept     + pv_Ept_grid(i-1+ii,k-1+kk,j-1+jj,:)*weight_cubic(ii,jj,kk)
                        end do
                    end do
                end do
            end if
        end do
            
        
        B           = B/N_gyro
        grad_B      = grad_B/N_gyro
        curl_b      = curl_b/N_gyro
        delta_E     = delta_E/N_gyro
        grad_RB_phi = grad_RB_phi/N_gyro
        pBpt        = pBpt/N_gyro
        
        abs_B    = sqrt(dot_product(B,B))
        b_unit   = B/abs_B
        pabs_Bpt = dot_product(b_unit,pBpt)
        
        B_eq = B_eq/N_gyro
        grad_B_eq = grad_B_eq/N_gyro
        curl_b_eq = curl_b_eq/N_gyro
        psi_par   = psi_par/N_gyro
        b_unit_eq   = B_eq/sqrt(dot_product(B_eq,B_eq))
        
        if(VEB) pb_unitpt = (pBpt - pabs_Bpt*b_unit)/abs_B  !calc pb_unitpt, which affect equation of dXdt and not affect equation of dv_||dt because of pb_unitpt * B = 0
                              
        return
        
    end subroutine gyro_average_xyz

    
    
    
    
!********************************************************************************************************************
!                                                                                                                   *
!     diagnostics for distribution function vs. P_phi and Energy E with fiexd pitch angle Lambda                    *
!                                                           (where Lambda = mu*B_0/E            = Lambda_0)         *
!                                                           (      E      = 0.5*m*(v_||)^2+mu*B = E_0	  )         *
!                                                                                                                   *
!    #file structrue : f_vs_P_phi_and_E_vs_Time_Lambda=XXX_windows_width=XXX.dat                                    *
!                         E E E                                                                                     *
!                       P X X X                                                                                     *
!                       P X X X                                                                                     *
!                       P X X X                                                                                     *
!                       p O O O                                                                                     *
!                       p O O O                                                                                     *
!                       P O O O                                                                                     *
!                                                                                                                   *
!            here    :  E : Energy  P : P_phi    X : t=t0    O : t=t0+dt(next time step)                            *
!                                                                                                                   *
!                                                                                                                   *
!*******************************************************************************************************************/
    
    subroutine diagnostics_output_f_vs_P_phi_and_E(Lambda_diag, windows_width, index)
        use var
        use distribution
        use diganostic
        use strings
        implicit none
        include 'mpif.h'
        
        integer :: STATUS(MPI_STATUS_SIZE)
        real*8, intent(in)  :: Lambda_diag, windows_width
        integer, intent(in) :: index
        real*8  :: my_P_phi_min, my_P_phi_max, my_E_min, my_E_max
        integer :: index_P_phi, index_E
        real*8  :: P_phi_floor, E_floor
        
        character*100 filename
        integer :: unit_number
        
      
        
!first step to get maxima and minimum of E and P_phi
        if(nstep == NSTEP_START) then
            
            my_P_phi_min =  1.0e9
            my_P_phi_max = -1.0e9
            my_E_min     =  1.0e9
            my_E_max     = -1.0e9
            
            do p=1,N
                if(DEB_type*marker(p)%v_para > 0 .or. DEB_type == 0 .or. DEP) then
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
                                B      = B + x(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
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

                    !update P_phi, E,Lambda
                    P_phi  = m*marker(p)%v_para*marker(p)%X(1)*(B(2)/abs_B) - Zh*psi_par
                    E      = 0.5*m*marker(p)%v_para**2 + marker(p)%mu*abs_B
                    Lambda = marker(p)%mu*B0/E
                
                    !<psi>-E 2D phase space
                    if(DEB) then
                        v          = sqrt(2*E/m)
                        if(marker(p)%v_para > 0) then
                            sgn_v_parallel = +1
                        else
                            sgn_v_parallel = -1
                        end if
      
                        if(ORBIT_AVERAGE_METHOD == 1) then
                            if((1-Lambda) > 0) then
                                bracket_psi = - P_phi/Zh + m_div_Zh*sgn_v_parallel*v*R0*sqrt(1-Lambda)
                            else
                                bracket_psi = - P_phi/Zh
                            end if
                        else if(ORBIT_AVERAGE_METHOD == 2) then
                            bracket_psi = - P_phi/Zh
                        end if
                        P_phi = bracket_psi
                    end if
                
        

                    if(abs(Lambda - Lambda_diag) < windows_width) then
                        my_P_phi_max = max(my_P_phi_max, P_phi)
                        my_P_phi_min = min(my_P_phi_min, P_phi)
                        my_E_max     = max(my_E_max, E)
                        my_E_min     = min(my_E_min, E)
                    end if
                    
                end if
            
            end do


        
            CALL MPI_ALLREDUCE(my_P_phi_max,P_phi_max(index),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
            CALL MPI_ALLREDUCE(my_P_phi_min,P_phi_min(index),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)
            
            CALL MPI_ALLREDUCE(my_E_max,E_max(index),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
            CALL MPI_ALLREDUCE(my_E_min,E_min(index),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)
            
            dP_phi(index) = (P_phi_max(index) - P_phi_min(index))/(GRID_P_PHI-1)
            
            dE(index) = (E_max(index) - E_min(index))/(GRID_E-1)
        
        end if 
        
        if(nstep == NSTEP_START) then
            myf_vs_P_phi_and_E(:,:,index)         = 0.0
            mydeltaf_vs_P_phi_and_E(:,:,index)    = 0.0
		    mynum_2d(:,:,index)                   = 0.0
            mynum_plot(index)                     = 0.0
        end if
    
        do p=1,N
            if(DEB_type*marker(p)%v_para > 0 .or. DEB_type == 0 .or. DEP) then
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
                            B      = B + x(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
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

                !update P_phi, E,Lambda
                P_phi  = m*marker(p)%v_para*marker(p)%X(1)*(B(2)/abs_B) - Zh*psi_par
                E      = 0.5*m*marker(p)%v_para**2 + marker(p)%mu*abs_B
                Lambda = marker(p)%mu*B0/E
            
                !<psi>-E 2D phase space
                if(DEB) then
                    v          = sqrt(2*E/m)
                    if(marker(p)%v_para > 0) then
                        sgn_v_parallel = +1
                    else
                        sgn_v_parallel = -1
                    end if
      
                    if(ORBIT_AVERAGE_METHOD == 1) then
                        if((1-Lambda) > 0) then
                            bracket_psi = - P_phi/Zh + m_div_Zh*sgn_v_parallel*v*R0*sqrt(1-Lambda)
                        else
                            bracket_psi = - P_phi/Zh
                        end if
                    else if(ORBIT_AVERAGE_METHOD == 2) then
                        bracket_psi = - P_phi/Zh
                    end if
                    P_phi = bracket_psi
                end if
            
                if(E > E_min(index) .and. E < E_max(index) .and. P_phi > P_phi_min(index) .and. P_phi < P_phi_max(index) .and. abs(Lambda - Lambda_diag) <= windows_width) then
                    index_P_phi        = floor((P_phi-P_phi_min(index))/dP_phi(index))
				    P_phi_floor        = P_phi_min(index) + index_P_phi*dP_phi(index)
				    weight_line(1,1)   = (P_phi - P_phi_floor)/dP_phi(index)
				    weight_line(1,2)   = 1 - weight_line(1,1)
                    index_E            = floor((E-E_min(index))/dE(index))
				    E_floor            = E_min(index) + index_E*dE(index)
                    weight_line(2,1)   = (E - E_floor)/dE(index)
				    weight_line(2,2)   = 1 - weight_line(2,1)
                
                    mynum_plot(index) = mynum_plot(index) + 1
                
                    do ii=1,2
                        do jj=1,2
                            weight_square(ii,jj) = weight_line(1,3-ii)*weight_line(2,3-jj)
                        end do
                    end do
                
                    do ii=1,2
                        do jj=1,2
                            mydeltaf_vs_P_phi_and_E(index_P_phi+ii,index_E+jj,index)   = mydeltaf_vs_P_phi_and_E(index_P_phi+ii,index_E+jj,index) + marker(p)%w*marker(p)%g*weight_square(ii,jj)
                            myf_vs_P_phi_and_E(index_P_phi+ii,index_E+jj,index)        = myf_vs_P_phi_and_E(index_P_phi+ii,index_E+jj,index)      +      marker(p)%f_over_g*weight_square(ii,jj)
                            mynum_2d(index_P_phi+ii,index_E+jj,index)                  = mynum_2d(index_P_phi+ii,index_E+jj,index)                +                     1.0*weight_square(ii,jj)
                        end do
                    end do  
                
                end if
            
            end if
            
        end do
        
        if(nstep_local==NSTEP_AVG .and. nstep .ne. NSTEP_START) then
            CALL MPI_REDUCE(myf_vs_P_phi_and_E(:,:,index),f_vs_P_phi_and_E,GRID_P_PHI*GRID_E,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
            CALL MPI_REDUCE(mynum_2d(:,:,index),num_2d,GRID_P_PHI*GRID_E,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
            CALL MPI_REDUCE(mydeltaf_vs_P_phi_and_E(:,:,index),deltaf_vs_P_phi_and_E,GRID_P_PHI*GRID_E,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
            
            CALL MPI_REDUCE(mynum_plot(index),num_plot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)

            f_vs_P_phi_and_E      = f_vs_P_phi_and_E/num_2d
            deltaf_vs_P_phi_and_E = deltaf_vs_P_phi_and_E/num_2d
            
            myf_vs_P_phi_and_E(:,:,index)         = 0.0
            mydeltaf_vs_P_phi_and_E(:,:,index)    = 0.0
		    mynum_2d(:,:,index)                   = 0.0
            mynum_plot(index)                     = 0.0
        end if
        
!output file        
        if(nrank == 0) then
            if(nstep == NSTEP_START) then
                if(DEP) filename = "DEP_parameter_Lambda="//trim(num2str(Lambda_diag))//".dat"
                if(DEB) filename = "DEB_parameter_Lambda="//trim(num2str(Lambda_diag))//".dat"
                unit_number = 1000 + floor(Lambda_diag*100)
                open(unit=unit_number,file=filename)
                
                if(DEP) filename = "DEP_Time_Lambda="//trim(num2str(Lambda_diag))//".dat"
                if(DEB) filename = "DEB_Time_Lambda="//trim(num2str(Lambda_diag))//".dat"
                open(unit=unit_number+1,file=filename,status='unknown',form='formatted',position='append')            
                
                if(DEP) filename = "DEP_f(E,P_phi,Lambda="//trim(num2str(Lambda_diag))//",t).dat"
                if(DEB) filename = "DEB_f(E,P_phi,Lambda="//trim(num2str(Lambda_diag))//",t).dat"
                open(unit=unit_number+2,file=filename,status='unknown',form='formatted',position='append')
                
                if(DEP) filename = "DEP_deltaf(E,P_phi,Lambda="//trim(num2str(Lambda_diag))//",t).dat"
                if(DEB) filename = "DEB_deltaf(E,P_phi,Lambda="//trim(num2str(Lambda_diag))//",t).dat"
                open(unit=unit_number+3,file=filename,status='unknown',form='formatted',position='append')
                
                !initial parameter output
                write(unit_number,'(e13.5)')E_min(index),E_max(index),P_phi_min(index),P_phi_max(index),Lambda_diag,windows_width
                write(unit_number,'(I13)')GRID_E,GRID_P_PHI
                
            else if(nstep_local==NSTEP_AVG) then
                unit_number = 1000 + floor(Lambda_diag*100)
                write(unit_number+1,'(e13.5,i13)')time,num_plot
                do i=1,GRID_P_PHI
                    write(unit_number+2,787)(f_vs_P_phi_and_E(i,j),j=1,GRID_E)
                    write(unit_number+3,787)(deltaf_vs_P_phi_and_E(i,j),j=1,GRID_E)
                end do
787             format(<GRID_E>e13.5)

            end if
            
        end if
        
        
    end subroutine diagnostics_output_f_vs_P_phi_and_E
    
    subroutine free_energy_vs_P_phi_and_E(Lambda_diag, windows_width, index)
        use var
        use distribution
        use diganostic
        use strings
        implicit none
        include 'mpif.h'
        
        integer :: STATUS(MPI_STATUS_SIZE)
        real*8, intent(in)  :: Lambda_diag, windows_width
        integer, intent(in) :: index
        real*8  :: my_P_phi_min, my_P_phi_max, my_E_min, my_E_max
        integer :: index_P_phi, index_E
        real*8  :: P_phi_floor, E_floor
        
        character*100 filename
        integer :: unit_number
        
      
        
!first step to get maxima and minimum of E and P_phi
            
        my_P_phi_min =  1.0e9
        my_P_phi_max = -1.0e9
        my_E_min     =  1.0e9
        my_E_max     = -1.0e9
            
        do p=1,N
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
                        B      = B + x(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
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

            !update P_phi, E,Lambda
            P_phi  = m*marker(p)%v_para*marker(p)%X(1)*(B(2)/abs_B) - Zh*psi_par
            E      = 0.5*m*marker(p)%v_para**2 + marker(p)%mu*abs_B
            Lambda = marker(p)%mu*B0/E

            if(abs(Lambda - Lambda_diag) < windows_width) then
                my_P_phi_max = max(my_P_phi_max, P_phi)
                my_P_phi_min = min(my_P_phi_min, P_phi)
                my_E_max     = max(my_E_max, E)
                my_E_min     = min(my_E_min, E)
            end if
            
        end do
        
        CALL MPI_ALLREDUCE(my_P_phi_max,P_phi_max(index),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
        CALL MPI_ALLREDUCE(my_P_phi_min,P_phi_min(index),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)
            
        CALL MPI_ALLREDUCE(my_E_max,E_max(index),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
        CALL MPI_ALLREDUCE(my_E_min,E_min(index),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)
            
        dP_phi(index) = (P_phi_max(index) - P_phi_min(index))/(GRID_P_PHI-1)
            
        dE(index) = (E_max(index) - E_min(index))/(GRID_E-1)
        
        mydf0dE_vs_P_phi_and_E(:,:,index)       = 0.0
        mydf0dP_phi_vs_P_phi_and_E(:,:,index)   = 0.0
		mynum_2d(:,:,index)                     = 0.0
    
        do p=1,N
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
                        B      = B + x(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
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

            !update P_phi, E,Lambda
            P_phi  = m*marker(p)%v_para*marker(p)%X(1)*(B(2)/abs_B) - Zh*psi_par
            E      = 0.5*m*marker(p)%v_para**2 + marker(p)%mu*abs_B
            Lambda = marker(p)%mu*B0/E
                               
            
            if(E > E_min(index) .and. E < E_max(index) .and. P_phi > P_phi_min(index) .and. P_phi < P_phi_max(index) .and. abs(Lambda - Lambda_diag) <= windows_width) then
                index_P_phi        = floor((P_phi-P_phi_min(index))/dP_phi(index))
				P_phi_floor        = P_phi_min(index) + index_P_phi*dP_phi(index)
				weight_line(1,1)   = (P_phi - P_phi_floor)/dP_phi(index)
				weight_line(1,2)   = 1 - weight_line(1,1)
                index_E            = floor((E-E_min(index))/dE(index))
				E_floor            = E_min(index) + index_E*dE(index)
                weight_line(2,1)   = (E - E_floor)/dE(index)
				weight_line(2,2)   = 1 - weight_line(2,1)
                
                do ii=1,2
                    do jj=1,2
                        weight_square(ii,jj) = weight_line(1,3-ii)*weight_line(2,3-jj)
                    end do
                end do

                !output free energy
                v          = sqrt(2*E/m)
                mu         = marker(p)%mu
                    
                if(marker(p)%v_para > 0) then
                    sgn_v_parallel = +1
                else
                    sgn_v_parallel = -1
                end if
                
                if(ORBIT_AVERAGE_METHOD == 1) then
                    if((1-Lambda) > 0) then
                        bracket_psi = - P_phi/Zh + m_div_Zh*sgn_v_parallel*v*R0*sqrt(1-Lambda)
                    else
                        bracket_psi = - P_phi/Zh
                    end if
                else if(ORBIT_AVERAGE_METHOD == 2) then
                    bracket_psi = - P_phi/Zh
                end if
    
    !Finite width of Gaussian distribution of Pitch angle  : f(Lambda) ~ exp[(Lambda-Lambda_0)^2/Delta_Lambda^2]
    
                if(DFR_type == 1) then
                    if(FGP) then
                        f0 = c_f*1.0/(v**3+v_c**3)*(1+erf((v_0-v)/Delta_v))*exp(-(bracket_psi/(c_1*Delta_psi)))*exp(-((Lambda-Lambda_0)/Delta_Lambda)**2)
                    else
                        f0 = c_f*1.0/(v**3+v_c**3)*(1+erf((v_0-v)/Delta_v))*exp(-(bracket_psi/(c_1*Delta_psi)))
                    end if
                else if(DFR_type == 2) then
                    bracket_psi_norm = (bracket_psi - psmin)/(psmax-psmin)
                    f0 = c_f*1.0/(v**3+v_c**3)*(1+erf((v_0-v)/Delta_v))*exp(-Delta_n/L_n*tanh((bracket_psi_norm - psi_n)/Delta_n))
                end if

                if(DFR_type == 1) then
                    df0dP_phi =  1.0/(Zh*c_1*Delta_psi)*f0
                else if(DFR_type == 2) then
                    df0dP_phi =  1.0/(Zh*L_n*Delta_psi)*sech((bracket_psi_norm - psi_n)/Delta_n)**2*f0
                end if
        
                if(ORBIT_AVERAGE_METHOD == 1) then
                    if((1-Lambda) > 0) then
                        if(DFR_type == 1) then
                            df0dE = -3/m*v/(v**3+ v_c**3)*f0 - exp(-((v_0-v)/Delta_v)**2)/(Delta_v*sqrt(PI*E*m/2)*(1+erf((v_0-v)/Delta_v)))*f0 - sgn_v_parallel*R0/(Zh*c_1*Delta_psi*sqrt(v*v-2*mu*B0/m))*f0
                        else if(DFR_type == 2) then
                            df0dE = -3/m*v/(v**3+ v_c**3)*f0 - exp(-((v_0-v)/Delta_v)**2)/(Delta_v*sqrt(PI*E*m/2)*(1+erf((v_0-v)/Delta_v)))*f0 - sgn_v_parallel*R0/(Zh*c_1*Delta_psi*sqrt(v*v-2*mu*B0/m))*sech((bracket_psi_norm - psi_n)/Delta_n)**2*f0
                        end if
                    else
                        df0dE = -3/m*v/(v**3+ v_c**3)*f0 - exp(-((v_0-v)/Delta_v)**2)/(Delta_v*sqrt(PI*E*m/2)*(1+erf((v_0-v)/Delta_v)))*f0
                    end if
                else if(ORBIT_AVERAGE_METHOD == 2) then
                    df0dE = -3/m*v/(v**3+ v_c**3)*f0 - exp(-((v_0-v)/Delta_v)**2)/(Delta_v*sqrt(PI*E*m/2)*(1+erf((v_0-v)/Delta_v)))*f0
                end if
                
                do ii=1,2
                    do jj=1,2
                        mydf0dE_vs_P_phi_and_E(index_P_phi+ii,index_E+jj,index)     = mydf0dE_vs_P_phi_and_E(index_P_phi+ii,index_E+jj,index)       +       df0dE*weight_square(ii,jj)
                        mydf0dP_phi_vs_P_phi_and_E(index_P_phi+ii,index_E+jj,index) = mydf0dP_phi_vs_P_phi_and_E(index_P_phi+ii,index_E+jj,index)   +   df0dP_phi*weight_square(ii,jj)
                        mynum_2d(index_P_phi+ii,index_E+jj,index)                   = mynum_2d(index_P_phi+ii,index_E+jj,index)                     +         1.0*weight_square(ii,jj)
                    end do
                end do  
                
            end if
            
        end do
        
        CALL MPI_REDUCE(mynum_2d(:,:,index),num_2d,GRID_P_PHI*GRID_E,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
            
        CALL MPI_REDUCE(mydf0dP_phi_vs_P_phi_and_E(:,:,index),df0dP_phi_vs_P_phi_and_E,GRID_P_PHI*GRID_E,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
        CALL MPI_REDUCE(mydf0dE_vs_P_phi_and_E(:,:,index),df0dE_vs_P_phi_and_E,GRID_P_PHI*GRID_E,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
                
        df0dP_phi_vs_P_phi_and_E  = df0dP_phi_vs_P_phi_and_E/num_2d
        df0dE_vs_P_phi_and_E      = df0dE_vs_P_phi_and_E/num_2d

        
!output file
        if(nrank == 0) then
            filename = "FEP_parameter_Lambda="//trim(num2str(Lambda_diag))//".dat"
            unit_number = 2000 + floor(Lambda_diag*100)
            open(unit=unit_number,file=filename)
                
            !initial parameter output
            write(unit_number,'(e13.5)')E_min(index),E_max(index),P_phi_min(index),P_phi_max(index),Lambda_diag,windows_width
            write(unit_number,'(I13)')GRID_E,GRID_P_PHI
                
            filename = "FEP_gamma(E,P_phi,Lambda="//trim(num2str(Lambda_diag))//").dat"
            open(unit=unit_number+1,file=filename,status='unknown',form='formatted')
                
            do i=1,GRID_P_PHI
                write(unit_number+1,797)(df0dP_phi_vs_P_phi_and_E(i,j),j=1,GRID_E)
            end do
            do i=1,GRID_P_PHI
                write(unit_number+1,797)(df0dE_vs_P_phi_and_E(i,j),j=1,GRID_E)
            end do
        
797         format(<GRID_E>e13.5)
            
        end if       
        
    end subroutine free_energy_vs_P_phi_and_E


!************************************************************************************************
!                                                                                               *
!           Calculate toroidal precession frequency : omega_phi   = Delta_phi/Delta_t           *
!                     poloidal transit    frequency : omega_theta = 2*pi/Delta_t                *
!                                                                                               *
!***********************************************************************************************/    
    subroutine calc_oribt_frequency(Lambda_initia, CO_GOING)
        use var
        use distribution
        use COF_var
        use strings
        implicit none
        include 'mpif.h'
        
        integer :: STATUS(MPI_STATUS_SIZE)
        real*8, intent(in)  :: Lambda_initia
        logical, intent(in) :: CO_GOING
        integer :: i_COF, j_COF        
        real*8  :: phi_old, phi_absolute
        real*8,dimension(2) :: data_diag
        real*8,allocatable :: data_diag_total(:,:)
        logical :: flag_diag
        logical,allocatable :: flag_diag_total(:)
        character*100 filename
        allocate(data_diag_total(2,nsize))
        allocate(flag_diag_total(nsize))
        
        
        
        if(nrank == 0) then
            if(CO_GOING) then
                filename = "oribt_frequency_Lambda="//trim(num2str(Lambda_initia))//"_co-going.dat"
            else
                filename = "oribt_frequency_Lambda="//trim(num2str(Lambda_initia))//"_counter-going.dat"
            end if
            
            open(unit=4000,file=filename,status='unknown',form='formatted',position='append')
            
            write(*,*)"Start calculating particle orbit frequency ......"
        end if
        
        call calc_grad_B
        call calc_polarization_drift_var
        
        do i_COF=1,NUM_E
            do j_COF=1,NUM_PPHI
                flag_diag = .false.
                R_initia    = R0 + (xmax-R0)*0.95*j_COF/(NUM_PPHI-1)
                if(Lambda_initia > 0.90) R_initia    = R0*Lambda_initia*1.05 + (xmax-R0*Lambda_initia*1.05)*0.95*j_COF/(NUM_PPHI-1)
                phi_initia  = 0.0
                Z_initia    = 0.0
                E_initia    = 0.0
                P_phi_initia = 0.0
                
                if(R_initia > myxmax .or. R_initia < myxmin .or. phi_initia > myymax .or. phi_initia < myymin .or. Z_initia > myzmax .or. Z_initia < myzmin) then
                    N = 0
                else
                    N = 1
                end if

                do p=1,N
                    E_initia    = 0.1 + (3.0-0.1)*i_COF/(NUM_E-1)
                    marker(p)%id = 1
                    marker(p)%X(1) = R_initia
                    marker(p)%X(2) = phi_initia
                    marker(p)%X(3) = Z_initia
                    
                
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
                                B      = B + x(i-1+ii,k-1+kk,j-1+jj,6:8)*weight_cubic(ii,jj,kk)
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
                    marker(p)%mu     = E_initia*Lambda_initia/B0
                    marker(p)%v_perp = sqrt(2*marker(p)%mu*abs_B/m)
                    !true for co-passing while false for counter-passing particle
                    if(CO_GOING) then
                        marker(p)%v_para = + sqrt(2*E_initia/m - marker(p)%v_perp**2)
                    else
                        marker(p)%v_para = - sqrt(2*E_initia/m - marker(p)%v_perp**2)
                    end if    
                    marker(p)%OUT_OF_BOUNDARY = .false.
                    write(*,*)"marker(p)%v_para : ",marker(p)%v_para," R_initia : ",R_initia
                    P_phi_initia  = m*marker(p)%v_para*marker(p)%X(1)*(B(2)/abs_B) - Zh*psi_par                        
                end do

                call MPI_ALLREDUCE(P_phi_initia, P_phi_initia, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERROR)
                call MPI_ALLREDUCE(E_initia, E_initia, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERROR)
    
                time = 0
                dt = 0.02
                FINISH_COF = .false.
                
                TIMESTEP_COF = 1
                Z_reverse_sign = .false.
                phi_absolute = phi_initia
                do while(.not. FINISH_COF)
                    do p=1,N
                        phi_old = marker(p)%X(2)
                        marker(p)%v_perp = phi_old
                    end do
                    
                    call particle_orbit_polarization
                    
                    do p=1,N
                        if(TIMESTEP_COF == 1) then
                            dZ_of_step1 = marker(p)%X(3) - Z_initia
                            marker(p)%f_over_g = dZ_of_step1
                            marker(p)%g        = marker(p)%v_perp                            
                        end if
                        

                        phi_old = marker(p)%v_perp
                        phi_absolute = marker(p)%g
                        phi_absolute = phi_absolute + marker(p)%X(2) - phi_old
                        if(abs(marker(p)%X(2)-phi_old)>abs(marker(p)%X(2)-2*PI-phi_old)) phi_absolute = phi_absolute - 2*PI
                        if(abs(marker(p)%X(2)-phi_old)>abs(marker(p)%X(2)+2*PI-phi_old)) phi_absolute = phi_absolute + 2*PI
                        marker(p)%g        = phi_absolute

                    end do
                    
                    call MPI_ALLREDUCE(N,N_total,1,MPI_INTEGER,MPI_SUM, MPI_COMM_WORLD, IERROR)
                    if(N_total == 0) then
                        omega_theta = 0
                        omega_phi   = 0
                        FINISH_COF = .true.
                    else
                        do p=1,N
                            if(marker(p)%id == 1) then
                                dZ_of_step1 = marker(p)%f_over_g
                                if(marker(p)%X(3)*dZ_of_step1 < 0) Z_reverse_sign = .true.
                                phi_absolute = marker(p)%g
                                if(abs(marker(p)%X(1)-R_initia) < eps_COF*a .and. abs(marker(p)%X(3)-Z_initia) < eps_COF*a .and. Z_reverse_sign) then
                                    omega_theta = 2*PI/time
                                    omega_phi   = (phi_absolute - phi_initia)/time
                                    FINISH_COF  = .true.
                                end if                
                            end if
                        end do                   
                    end if
                    
                    call MPI_ALLREDUCE(FINISH_COF, FINISH_COF, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, IERROR)
                    
                    TIMESTEP_COF = TIMESTEP_COF + 1
                    time = time + dt
                end do
                
                if(N>0) then
                    data_diag(1) = omega_theta
                    data_diag(2) = omega_phi
                    flag_diag    = .true.
                end if
                
                call MPI_GATHER(data_diag, 2, MPI_DOUBLE_PRECISION, data_diag_total(1:2,:), 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)
                call MPI_GATHER(flag_diag, 1, MPI_LOGICAL, flag_diag_total, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, IERROR)
                if(nrank == 0) then
                    write(*,*)"COF index :",i_COF,j_COF,"  time = ",time
                    if(N_total == 0) then
                        write(4000,'(4(e13.5))')E_initia,P_phi_initia,omega_theta,omega_phi
                    else
                        do s=1,nsize
                            if(flag_diag_total(s)) write(4000,'(4(e13.5))')E_initia,P_phi_initia,data_diag_total(:,s)
                        end do
                    end if
                end if
                
            end do
        end do
           

    end subroutine calc_oribt_frequency
    
    
!************************************************************************************************
!                                                                                               *
!                     output time evolution of E_phi vs psi(or r) at midplane                   *
!                                                                                               *
!***********************************************************************************************/ 
      subroutine output_E_phi_vs_psi
        USE DECLARE
        include 'mpif.h' 
        real*8, dimension(dgn_size_group) :: myE_phi_vs_psi,E_phi_vs_psi
      
        if(nrank == 0 .and. nstep == 0) then
            open(unit=166,file='E_phi_vs_psi.dat',position='append')
        end if

        myE_phi_vs_psi = 0
        E_phi_vs_psi   = 0

        !calc indexs of radial E_phi for diagnostic
        if(nstep == 0) then
            !do diagnostic for different radial location, span by q~(q_min*1.05, q_max*0.95)
            do jdgn=1,dgn_size_group
                qmode_group(jdgn) = qmin + (qmax-qmin)/(dgn_size_group-1)*(jdgn-1)
            end do
            
            do jdgn=1,dgn_size_group
                j = 1
                do while(q_NOVA(j) .lt. qmode_group(jdgn))
                    j=j+1
                end do
                if(j==1) j=j+2
                if(j==2) j=j+1

                call interp1d3l(xxst(3,j-2),xxst(3,j-1),xxst(3,j),xxst(3,j+1), &
                            q_NOVA(j-2),q_NOVA(j-1),q_NOVA(j),q_NOVA(j+1),qmode_group(jdgn),xxmode)
      
                jxtmode=floor((xxmode-xmin)/(xmax-xmin)*(mxt-1))+1
                jztmode=mzt/2
                nrkx_mode=(jxtmode-1)/mxm
                nrkz_mode=(jztmode-1)/mzm
                jxmode_group(jdgn)=jxtmode-nrkx_mode*mxm+2
                jzmode_group(jdgn)=jztmode-nrkz_mode*mzm+2
                nrank_mode_group(jdgn)=nrkz_mode*nprx+nrkx_mode
!                weight_group(jdgn)     = (xxmode-xxt(jxtmode))/(xxt(jxtmode+1) - xxt(jxtmode))
!                psi_group(jdgn)        = (1.0-weight_group(jdgn))*pst(jxtmode,jztmode) + weight_group(jdgn)*pst(jxtmode+1,jztmode)
                psi_group(jdgn)        = pst(jxtmode,jztmode)
            end do
            if(nrank == 0) write(166,7041)time, psi_group
            if(nrank == 0) write(166,7041)time, qmode_group
        end if
        
            
              
        do jdgn=1,dgn_size_group
            if(nrank==nrank_mode_group(jdgn)) then
!                myE_phi_vs_psi(jdgn) = (1.0-weight_group(jdgn))*Ef(jxmode_group(jdgn),jzmode_group(jdgn),1,2) + weight_group(jdgn)*Ef(jxmode_group(jdgn)+1,jzmode_group(jdgn),1,2)
                myE_phi_vs_psi(jdgn) = Ef(jxmode_group(jdgn),jzmode_group(jdgn),1,2)
            end if
        end do
        call MPI_REDUCE(myE_phi_vs_psi, E_phi_vs_psi, dgn_size_group, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, IERROR)
        if(nrank == 0) write(166,7041)time, E_phi_vs_psi

7041    format(<dgn_size_group+1>e13.5)
        return
      end
