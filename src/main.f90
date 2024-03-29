program main
  use modules, only : rk, init_meshgeom, calc_meshgeom, init_Riemann, max_vel, boundary_walls, boundary_grad, boundary_movewalls
  use io, only      : output_field, output_monitor, output_monitor_x
  use godunov, only : flux_godunov
  use vanleer, only : flux_vanleer
  use TVD, only     : reconst
  implicit none
  
  !Parameters
  character(*),parameter   :: input_file='params.nml', monitor_file='_mon.dat', output_file='_sol.dat'
  character(len=30)        :: start_file
  integer                  :: scheme, grani, tvdim=0
  real(kind=rk)            :: L, adi_k, c_P, CFL, dt, t_stop
  real(kind=rk)            :: U_L, U_R, rho_L, rho_R, p_L, p_R
  real(kind=rk)            :: freql=0, hl=0, freqr=0, hr=0, phaser=0
  integer                  :: nx, nt, L_size, imon, mon_tstep
  !Work arrays & variables
  real(kind=rk),allocatable:: x_cent(:), omega(:), u(:), rho(:), p(:), T(:)
  real(kind=rk),allocatable:: x_f(:), sigma_f(:), u_f(:)
  real(kind=rk),allocatable:: w(:,:), w_n(:,:), F(:,:), RHS(:,:)  
  real(kind=rk)            :: dx
  !Local variables
  real(kind=rk)            :: C_v, R_m, total_t=0
  real(kind=rk)            :: p_tilde(4), u_tilde(4), rho_tilde(4)
  real(kind=rk)            :: u_left=0, u_right=0 
  real(kind=rk),parameter  :: pi=4*atan(1.0_rk)
  character(len=30)        :: time_str
  !Service variables
  integer                  :: i, k, iu
  !Procedure pointers
  procedure(flux_godunov),pointer   :: flux=>flux_godunov
  procedure(boundary_walls),pointer :: boundary=>boundary_walls
  
 
  !/////////////////////////////////////////////////////////// 
  !Read input data
  namelist /params/ start_file, mon_tstep, scheme, grani, tvdim 
  namelist /params_zad/ L, adi_k, c_p,                &
                    nx, nt, CFL, t_stop, L_size,      &
                    U_L, U_R, rho_L, rho_R,           &
                    p_L, p_R, freql, hl, freqr, hr, phaser                    
  open(newunit=iu, file=input_file)
  read(iu, nml=params)
  close(iu)
  print*, 'Reading:', (trim(start_file) // ".nml")
  open(newunit=iu, file=(trim(start_file) // ".nml"))
  read(iu, nml=params_zad)
  close(iu)
  
  !Set scheme and boundary conditions
  if (scheme==2) flux=>flux_vanleer
  if (grani==2)  boundary=>boundary_grad
  if (grani==3)  boundary=>boundary_movewalls
  
  imon = L_size/2
  R_m  = C_p * (1-1.0_rk/adi_k)
  C_v  = C_p - R_m
  
  !Prepare work arrays, set initial conditions
  allocate(x_cent(0:nx+1),omega(0:nx+1),u(0:nx+1),rho(0:nx+1),p(0:nx+1),T(0:nx+1))
  allocate(x_f(1:nx+1),sigma_f(1:nx+1),u_f(1:nx+1))
  allocate(w(3,0:nx+1),w_n(3,0:nx+1),F(3,0:nx+1),RHS(3,0:nx+1))  
  
  call init_meshgeom(L,nx,dx,x_cent,x_f,sigma_f,omega)
  u = 0 ; rho = 0 ; p = 0 ; T = 0
  w = 0 ; w_n = 0 ; F = 0 ; u_f = 0
  call init_Riemann(u,rho,p,T,                          &
                    L_size,U_L,U_R,rho_L,rho_R,p_L,p_R, &
                    adi_k,C_p)
  call boundary(u,rho,p,T,adi_k,C_p,u_left,u_right)
  
  !Open monitor file
  open(newunit=iu,file=(trim(start_file) // monitor_file))
  
  !Solve equation
  time: do k = 1, nt
    if (total_t==t_stop) then
        print*, 'Saving'
        print*, 'Enter new end time:'
        write(time_str,'(es10.2)') total_t 
        !Output solution
        call output_field(x_cent,u,rho,p,T,trim(start_file) // trim(adjustl(time_str)) // output_file)        
        read(*,*) t_stop
        if (t_stop<total_t) exit
    end if 
    !Calculate time step
    dt = CFL * dx/(max_vel(u,rho,p,adi_k)+abs(2*pi*freql*hl))
    if (k<=5) dt = 0.2*dt
    total_t = total_t + dt
    if (total_t>t_stop) then
        dt = dt - (total_t-t_stop)
        total_t = t_stop
    end if
    print*, 'Step=', k, 'dt=', dt, 'Time=', total_t, 'Max velocity=', max_vel(u,rho,p,adi_k)
    !Conservative variables
    w(1,:)   = rho
    w(2,:)   = rho * u
    w(3,:)   = rho * (C_v*T + u**2/2)
    F(1,:)   = rho * u
    F(2,:)   = rho * u**2 + p
    F(3,:)   = rho * u * (C_p*T + u**2/2)
    !Added sources
    RHS(1,:) = 0
    RHS(2,:) = p*(sigma_f(1:nx)-sigma_f(2:nx+1))
    RHS(3,:) = 0
    
    !Moving mesh
    if (grani==3) then
        u_left = 2*pi*freql*hl*sin(2*pi*freql*(total_t-dt))
        u_right = 2*pi*freqr*hr*sin(2*pi*freqr*(total_t-dt)+phaser*pi/180)
        
        !Calculate new volumes and motion in each cell
        call calc_meshgeom(x_f,sigma_f,u_f,omega,x_cent,L,u_left,u_right,dt)
        call boundary(u,rho,p,T,adi_k,C_p,u_left,u_right)
    end if
    
    !Solve for conservative variables
    w_n(:,1) = w(:,1) - dt/omega(1) * &
        (flux(adi_k,C_p,p(1),p(2),rho(1),rho(2),u(1),u(2),u_f(2))*sigma_f(2)   &
       - flux(adi_k,C_p,p(0),p(1),rho(0),rho(1),u(0),u(1),u_f(1))*sigma_f(1))  &
       - dt/omega(1) * RHS(:,1)
    !TVD/non-TVD scheme far from boundaries
    if (tvdim==1) then
        do i=2, nx-1
            call reconst(p(i-2:i+2),rho(i-2:i+2),u(i-2:i+2),p_tilde, u_tilde, rho_tilde)
            w_n(:,i) = w(:,i) - dt/omega(i) * &
                (flux(adi_k,C_p,p_tilde(3),p_tilde(4),rho_tilde(3),rho_tilde(4),u_tilde(3),u_tilde(4),u_f(i+1))*sigma_f(i+1)   &
                - flux(adi_k,C_p,p_tilde(1),p_tilde(2),rho_tilde(1),rho_tilde(2),u_tilde(1),u_tilde(2),u_f(i))*sigma_f(i))    &
                - dt/omega(i) * RHS(:,i)
        end do
    else
        do i=2, nx-1
            w_n(:,i) = w(:,i) - dt/omega(i) * &
                (flux(adi_k,C_p,p(i),p(i+1),rho(i),rho(i+1),u(i),u(i+1),u_f(i+1))*sigma_f(i+1)   &
            - flux(adi_k,C_p,p(i-1),p(i),rho(i-1),rho(i),u(i-1),u(i),u_f(i))*sigma_f(i))    &
            - dt/omega(i) * RHS(:,i)
        end do
    end if
    w_n(:,nx) = w(:,nx) - dt/omega(nx) * &
        (flux(adi_k,C_p,p(nx),p(nx+1),rho(nx),rho(nx+1),u(nx),u(nx+1),u_f(nx+1))*sigma_f(nx+1)   &
       - flux(adi_k,C_p,p(nx-1),p(nx),rho(nx-1),rho(nx),u(nx-1),u(nx),u_f(nx))*sigma_f(nx))    &
       - dt/omega(nx) * RHS(:,nx)
       
    !Physical variables
    rho = w_n(1,:)
    u   = w_n(2,:)/rho
    T   = (w_n(3,:)/rho-u**2/2)/C_v
    p   = R_m * T * rho
    call boundary(u,rho,p,T,adi_k,C_p,u_left,u_right)
    !Monitor points
    if (mod(k,mon_tstep) == 0) then
        if (any(isnan(u))) stop 'Error:NaN'
        call output_monitor_x(iu,0.5_rk,total_t,u,rho,p,T,x_cent)
    end if
  end do time
  
  !Close monitor file
  close(iu)

end program main