program main
  use modules, only : rk, init_meshgeom, init_Riemann, boundary_walls 
  use io, only      : output_field
  use godunov, only : flux_godunov
  implicit none
  
  !Parameters
  character(*),parameter   :: input_file='params.nml'
  real(kind=rk)            :: L, adi_k, c_P, dt
  real(kind=rk)            :: U_L, U_R, rho_L, rho_R, p_L, p_R
  integer                  :: nx, nt, L_size
  !Work arrays & variables
  real(kind=rk),allocatable:: x_cent(:), omega(:), u(:), rho(:), p(:), T(:)
  real(kind=rk),allocatable:: x_f(:), sigma_f(:)
  real(kind=rk),allocatable:: w(:,:), w_n(:,:), F(:,:), RHS(:,:)  
  real(kind=rk)            :: dx
  !Local variables
  real(kind=rk)            :: C_v, R_m
  !Service variables
  integer                  :: i, k, iu
 
  !/////////////////////////////////////////////////////////// 
  !Read input data
  namelist /params/ L, adi_k, c_p,          &
                    nx, nt, dt, L_size,     &
                    U_L, U_R, rho_L, rho_R, &
                    p_L, p_R                    
  open(newunit=iu, file=input_file)
  read(iu, nml=params)
  close(iu)  
  R_m = C_p * (1-1.0_rk/adi_k)
  C_v = C_p - R_m
  
  !Prepare work arrays, set initial conditions
  allocate(x_cent(0:nx+1),omega(0:nx+1),u(0:nx+1),rho(0:nx+1),p(0:nx+1),T(0:nx+1))
  allocate(x_f(1:nx+1),sigma_f(1:nx+1))
  allocate(w(3,0:nx+1),w_n(3,0:nx+1),F(3,0:nx+1),RHS(3,0:nx+1))  
  
  call init_meshgeom(L,nx,dx,x_cent,x_f,sigma_f,omega)
  u = 0 ; rho = 0 ; p = 0 ; T = 0
  w = 0 ; w_n = 0 ; F = 0
  call init_Riemann(u,rho,p,T,                          &
                    L_size,U_L,U_R,rho_L,rho_R,p_L,p_R, &
                    adi_k,C_p)
  call boundary_walls(u,rho,p,T,adi_k,C_p)

  !Output initial field
  call output_field(x_cent,u,rho,p,T,'init.dat')  
  
  !Solve equation
  time: do k = 1, nt
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
    !Solve for conservative variables
    do i=1, nx
        w_n(:,i) = w(:,i) - dt/omega(i) * &
            (flux_godunov(adi_k,C_p,p(i),p(i+1),rho(i),rho(i+1),u(i),u(i+1),0.0_rk)*sigma_f(i+1)   &
           - flux_godunov(adi_k,C_p,p(i-1),p(i),rho(i-1),rho(i),u(i-1),u(i),0.0_rk)*sigma_f(i))    &
           - dt/omega(i) * RHS(:,i)
    end do
    !Physical variables
    rho = w_n(1,:)
    u   = w_n(2,:)/rho
    T   = (w_n(3,:)/rho-u**2/2)/C_v
    p   = R_m * T * rho
    call boundary_walls(u,rho,p,T,adi_k,C_p)    
  end do time
  
  !Output solution
  call output_field(x_cent,u,rho,p,T,'sol.dat')
  

end program main