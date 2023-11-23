module vanleer
  use modules, only: rk
  implicit none
  
  contains
  
  
  function flux_vanleer(AK,CP,P1,P2,R1,R2,U1,U2,UDOT) result(flux)
  real(kind=rk), intent(in):: AK,CP,P1,P2,R1,R2,U1,U2,UDOT
  real(kind=rk)            :: flux(3)
  real(kind=rk)            :: fl_plus(3), fl_minus(3) !Pf, Uf, Rf
  real(kind=rk)            :: C1, C2, M1, M2
  
  C1 = SQRT(AK*P1/R1); M1 = U1/C1
  C2 = SQRT(AK*P2/R2); M2 = U2/C2
  
  !Calculate F+ for physical variables in cell (i)
  fl_plus = [1.0_rk,                   &
             2*C1/AK*((AK-1)/2*M1+1),  &
             2*C1*C1/(AK*AK-1)*((AK-1)/2*M1+1)**2]
  fl_plus = fl_plus * R1/4 * C1 * (M1+1)**2
             
  !Calcualte F- for physical variables in cell(i+1)
  fl_minus = [1.0_rk,                  &
             2*C2/AK*((AK-1)/2*M2-1),  &
             2*C2*C2/(AK*AK-1)*((AK-1)/2*M2-1)**2]
  fl_minus = -fl_minus * R2/4 * C2 * (M2-1)**2  
   
  !Calculate F on cell face (i+1/2)
  flux = fl_plus + fl_minus
 
  end function flux_vanleer
  
end module vanleer