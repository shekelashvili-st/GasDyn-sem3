module TVD
  use modules, only: rk
  implicit none
  
  contains
  
  subroutine reconst(p, rho, u, p_tilde, u_tilde, rho_tilde)
  real(kind=rk), intent(in) :: p(1:5), rho(1:5), u(1:5)
  real(kind=rk)             :: r_p(1:3), r_u(1:3), r_rho(1:3)
  real(kind=rk)             :: small = 1e-7
  real(kind=rk)             :: p_R, p_L, p_LL, p_RR, &
                               u_R, u_L, u_LL, u_RR, &
                               rho_R, rho_L, rho_LL, rho_RR
  real(kind=rk), intent(out):: p_tilde(1:4), u_tilde(1:4), rho_tilde(1:4)
  !Индексы переменных: 1 - (i-2), 2 - (i-1), 3 - (i), 4 - (i+1), 4 - (i+2)
  !Индексы r - (i-1), i, (i+1)
  r_p(1:3)   = (p(3:5) - p(2:4))/(p(2:4)-p(1:3)+small)
  r_u(1:3)   = (u(3:5) - u(2:4))/(u(2:4)-u(1:3)+small)
  r_rho(1:3) = (rho(3:5) - rho(2:4))/(rho(2:4)-rho(1:3)+small)
  
  !Центральная ячейка
  p_R   = p(3) - 0.5_rk * psi(r_p(2)) * (p(3)-p(2))
  p_L   = p(3) + 0.5_rk * psi(r_p(2)) * (p(3)-p(2))
  u_R   = u(3) - 0.5_rk * psi(r_u(2)) * (u(3)-u(2))
  u_L   = u(3) + 0.5_rk * psi(r_u(2)) * (u(3)-u(2))
  rho_R = rho(3) - 0.5_rk * psi(r_rho(2)) * (rho(3)-rho(2))
  rho_L = rho(3) + 0.5_rk * psi(r_rho(2)) * (rho(3)-rho(2))  
  
  !Ячейка слева, поток слева 
  P_LL   = p(2) + 0.5_rk * psi(r_p(1)) * (p(2)-p(1))
  u_LL   = u(2) + 0.5_rk * psi(r_u(1)) * (u(2)-u(1))
  rho_LL = rho(2) + 0.5_rk * psi(r_rho(1)) * (rho(2)-rho(1))
  
  !Ячейка справа, поток справа
  p_RR   = p(4) - 0.5_rk * psi(r_p(3)) * (p(4)-p(3))
  u_RR   = u(4) - 0.5_rk * psi(r_u(3)) * (u(4)-u(3))
  rho_RR = rho(4) - 0.5_rk * psi(r_rho(3)) * (rho(4)-rho(3))
  
  !Вывод
  p_tilde(1:4)   = [p_LL, p_R, p_L, p_RR]
  u_tilde(1:4)   = [u_LL, u_R, u_L, u_RR]
  rho_tilde(1:4) = [rho_LL, rho_R, rho_L, rho_RR]
  end subroutine reconst
  
  function psi(r)
  real(kind=rk), intent(in):: r
  real(kind=rk)            :: psi, psir
  
  psi = 0
  !van Albada
  if (r>=0) then
    psir = 2.0_rk/(1+r)
    psi = (r*r+r)/(r*r+1)
    psi = min(psi,psir)
  end if
  
  end function psi
  
end module TVD