module modules
  implicit none

  integer,parameter           :: rk=8
  
  contains
  
  
  subroutine init_Riemann(u,rho,p,T,                            &
                          L_size,U_L,U_R,rho_L,rho_R,p_L,p_R,   &
                          adi_k,C_p)
  real(kind=rk),intent(in)    :: u_L, u_R, rho_L, rho_R, p_L, p_R
  real(kind=rk),intent(in)    :: adi_k, C_p
  integer,intent(in)          :: L_size
  real(kind=rk),intent(inout) :: u(0:), rho(0:), p(0:), T(0:)
  integer                     :: nx
  real(kind=rk)               :: R_m
  
  nx = size(u)-2
  R_m = C_p*(1-1.0_rk/adi_k)
  
  u(0:L_size) = u_L ; rho(0:L_size) = rho_L ; p(0:L_size) = p_L
  u(L_size+1:nx+1)= u_R ; rho(L_size+1:nx+1)= rho_R ; p(L_size+1:nx+1)= p_R
  
  T = p/(rho*R_m)
  
  end subroutine init_Riemann
  
  
  subroutine boundary_walls(u,rho,p,T,adi_k,C_p)
  real(kind=rk),intent(inout) :: u(0:), rho(0:), p(0:), T(0:)
  real(kind=rk),intent(in)    :: adi_k, C_p
  integer                     :: nx
  real(kind=rk)               :: R_m
  
  nx  = size(u)-2
  R_m = C_p*(1-1.0_rk/adi_k)
  
  !По скорости - нулевая скорость на грани
  u(0)      = -u(1)
  u(nx+1)   = -u(nx)
  !По давлению - нулевой градиент
  p(0)      = p(1)
  p(nx+1)   = p(nx)
  !По плотности - нулевой градиент
  rho(0)    = rho(1)
  rho(nx+1) = rho(nx)
  !По температуре - расчёт из плотности, давления
  T         = p/(rho*R_m)
  
  end subroutine boundary_walls
  
  
  elemental function cross_area(x,L) result(res)
  real(kind=rk),intent(in) :: x, L
  real(kind=rk)            :: res
  real(kind=rk),parameter  :: alpha = 12, pi=atan(1.0_rk)
  
  res = 2*x*tan(alpha*pi/180)+L/10
  
  end function cross_area


end module modules