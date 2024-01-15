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

  function max_vel(u,rho,p,adi_k)
  real(kind=rk)               :: max_vel
  real(kind=rk),intent(in)    :: adi_k
  real(kind=rk),intent(in)    :: u(0:), rho(0:), p(0:)
  integer                     :: nx
  
  nx = size(u)-2
  max_vel = maxval(abs(u(1:nx))+sqrt(adi_k*p(1:nx)/rho(1:nx)))
  
  end function max_vel  
  
  subroutine boundary_walls(u,rho,p,T,adi_k,C_p,u_l)
  real(kind=rk),intent(inout) :: u(0:), rho(0:), p(0:), T(0:)
  real(kind=rk),intent(in)    :: adi_k, C_p
  real(kind=rk),intent(in)    :: u_l
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
  
  subroutine boundary_moveleft(u,rho,p,T,adi_k,C_p,u_l)
  real(kind=rk),intent(inout) :: u(0:), rho(0:), p(0:), T(0:)
  real(kind=rk),intent(in)    :: adi_k, C_p
  real(kind=rk),intent(in)    :: u_l
  integer                     :: nx
  real(kind=rk)               :: R_m
  
  nx  = size(u)-2
  R_m = C_p*(1-1.0_rk/adi_k)
  
  !По скорости - скорость поршня на грани
  u(0)      = 2*u_l-u(1)
  u(nx+1)   = -u(nx)
  !По давлению - нулевой градиент
  p(0)      = p(1)
  p(nx+1)   = p(nx)
  !По плотности - нулевой градиент
  rho(0)    = rho(1)
  rho(nx+1) = rho(nx)
  !По температуре - расчёт из плотности, давления
  T         = p/(rho*R_m)
  
  end subroutine boundary_moveleft
  
  subroutine boundary_grad(u,rho,p,T,adi_k,C_p,u_l)
  real(kind=rk),intent(inout) :: u(0:), rho(0:), p(0:), T(0:)
  real(kind=rk),intent(in)    :: adi_k, C_p
  real(kind=rk),intent(in)    :: u_l
  integer                     :: nx
  real(kind=rk)               :: R_m
  
  nx  = size(u)-2
  R_m = C_p*(1-1.0_rk/adi_k)
  
  !По скорости - нулевой градиент на грани
  u(0)      = u(1)
  u(nx+1)   = u(nx)
  !По давлению - нулевой градиент
  p(0)      = p(1)
  p(nx+1)   = p(nx)
  !По плотности - нулевой градиент
  rho(0)    = rho(1)
  rho(nx+1) = rho(nx)
  !По температуре - расчёт из плотности, давления
  T         = p/(rho*R_m)
  
  end subroutine boundary_grad
  
  
  subroutine init_meshgeom(L,nx,dx,x_cent,x_f,sigma_f,omega)
  real(kind=rk),intent(in)    :: L
  integer,intent(in)          :: nx
  real(kind=rk),intent(out)   :: dx, x_cent(0:), x_f(:), sigma_f(:), omega(0:)
  integer                     :: i
  
  dx            = L/nx  
  x_cent        = [(dx*(i+0.5_rk), i=-1,nx)]
  x_f           = [(dx*i, i=0,nx)]
  sigma_f       = cross_area(x_f,L)
  !Объём как площадь трапеции
  omega(1:nx)   = (sigma_f(1:nx)+sigma_f(2:nx+1))/2 * dx
  omega(0)      = 0 ; omega(nx+1)       = 0
  
  end subroutine init_meshgeom
  
  
  subroutine calc_meshgeom(x_f,sigma_f,u_f,omega,x_cent,L,u_l,u_r,dt)
  real(kind=rk),intent(in)    :: sigma_f(:), L, u_l, u_r, dt
  real(kind=rk),intent(inout) :: x_f(:), u_f(:), omega(0:), x_cent(0:)
  integer                     :: i, nx
  
  nx = size(sigma_f)-1
  !Расчёт скорости каждой каждой грани
  u_f(1:nx+1) = u_l + (u_r-u_l)/(x_f(nx+1)-x_f(1))*(x_f(1:nx+1)-x_f(1))
  !Перемещение граней
  x_f(1:nx+1) = x_f(:) + dt*u_f(1:nx+1)
  !Объём как площадь трапеции
  omega(1:nx)   = (sigma_f(1:nx)+sigma_f(2:nx+1))/2 * (x_f(2:nx+1)-x_f(1:nx))
  omega(0)      = 0 ; omega(nx+1)       = 0
  !Центр ячейки
  x_cent(1:nx) = (x_f(2:nx+1)+x_f(1:nx))/2
  x_cent(0)    = x_cent(1)  - 2*(x_cent(1)-x_f(1))
  x_cent(nx+1) = x_cent(nx) + 2*(x_f(nx+1)-x_cent(nx))
  
  end subroutine calc_meshgeom
  
  
  elemental function cross_area(x,L) result(res)
  real(kind=rk),intent(in) :: x, L
  real(kind=rk)            :: res
!  real(kind=rk),parameter  :: alpha = 0.0_rk, pi=4*atan(1.0_rk)
  real(kind=rk),parameter  :: xa = 0.49, xb = 0.8
  real(kind=rk),parameter  :: wallr=1, a=-33.5671847202175, b=64.9525024336209, c=-39.4750092309758, d=8.19680440401464
  
!  res = 2*x*tan(alpha*pi/180)+L
  res = L/2
!  if (x>xb) res = wallr
!  if (x>xa.and.x<xb) res = (a*x**3 + b*x**2 + c*x + d)
  res = res*2
  end function cross_area
  
  
end module modules