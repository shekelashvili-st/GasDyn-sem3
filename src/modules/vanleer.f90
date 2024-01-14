module vanleer
  use modules, only: rk
  implicit none
  
  contains
  
  
  function flux_vanleer(AK,CP,P1,P2,R1,R2,U1,U2,UDOT) result(flux)
  real(kind=rk), intent(in):: AK,CP,P1,P2,R1,R2,U1,U2,UDOT
  real(kind=rk)            :: flux(3), w(3)
  real(kind=rk)            :: fl_plus(3), fl_minus(3) !Pf, Uf, Rf
  real(kind=rk)            :: RM, CV, C1, C2, M1, M2
  
  
  RM=CP*(1.0-1.0/AK)
  CV=CP-RM
  !Calculate speed of sound and M
  C1 = SQRT(AK*P1/R1); M1 = U1/C1
  C2 = SQRT(AK*P2/R2); M2 = U2/C2
  
  !Calculate F+ for physical variables in cell (i)
  if ((u1+c1)<=0) then
    fl_plus= [0.0_rk, 0.0_rk, 0.0_rk]
  else if ((u1-c1)>=0) then
    fl_plus = [R1*U1,                   &
                R1*U1*U1 + P1,  &
                R1*U1*(U1*U1/2+CP*P1/(R1*RM))]       
  else if (M1<=1) then
    fl_plus = [1.0_rk,                   &
                2*C1/AK*((AK-1)/2*M1+1),  &
                2*C1*C1/(AK*AK-1)*((AK-1)/2*M1+1)**2]
    fl_plus = fl_plus * R1/4 * C1 * (M1+1)**2
  end if
          
  !Calcualte F- for physical variables in cell(i+1)
  if ((u2-c2)>=0) then
    fl_minus= [0.0_rk, 0.0_rk, 0.0_rk]
  else if ((u2+c2)<=0) then
    fl_minus = [R2*U2,                   &
                R2*U2*U2 + P2,  &
                R2*U2*(U2*U2/2+CP*P2/(R2*RM))]     
  else if (M2<=1) then
    fl_minus = [1.0_rk,                  &
                2*C2/AK*((AK-1)/2*M2-1),  &
                2*C2*C2/(AK*AK-1)*((AK-1)/2*M2-1)**2]
    fl_minus = -fl_minus * R2/4 * C2 * (1-M2)**2
  end if
  
  !Calculate F on cell face (i+1/2)
  w(1)   = (R1 + R2)/2
  w(2)   = (R1*U1 + R2*U2)/2
  w(3)   = (R1*(CV*P1/(R1*RM) + U1**2/2) + R2*(CV*P2/(R2*RM) + U2**2/2))/2
  flux   = fl_plus + fl_minus - UDOT * w
 
  end function flux_vanleer
  
end module vanleer