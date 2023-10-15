program test_godunov
  use modules, only : rk
  use godunov, only : flux_godunov
  implicit none
  
  real(kind=rk)     :: adi_k=1.4, C_p=1005, p_L=100000, p_R=4000 
  real(kind=rk)     :: rho_L=1.2, rho_R=0.047, U_L=0, U_R=0
  real(kind=rk)     :: flux(3), error, eps=1
  
  flux = flux_godunov(adi_k,C_p,P_L,P_R,rho_L,rho_R,U_L,U_R,0.0_rk)
  
  error = maxval(abs( &
               flux - (/137.26732542817240,66979.595688964619,33363586.467703633/)))
  print*, error
  if (error >= eps ) call exit(1) 
  print*, flux
  
  
end program
