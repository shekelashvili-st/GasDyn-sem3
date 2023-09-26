module modules
  implicit none

  integer,parameter		      :: rk=8
  
  contains
  
  
  subroutine init_Riemann(u,L_size,u_L,u_R)
  real(kind=rk),intent(in)    :: u_L, u_R
  integer,intent(in)	      :: L_size
  real(kind=rk),intent(inout) :: u(0:)
  
  u(0:L_size) = u_L
  u(L_size+1:size(u)-1)= u_R
  
  end subroutine init_Riemann
  
  
  elemental function flux_cent(left,right,a) result(flux)
  real(kind=rk),intent(in)  :: left,right,a
  real(kind=rk)	            :: flux
  
  flux = a*(right+left)/2
  
  end function flux_cent
  
  elemental function flux_upwind(left,right,a) result(flux)
  real(kind=rk),intent(in)  :: left,right,a
  real(kind=rk)	            :: flux
  
  if (a>=0) then
	flux = a*left
  else
    flux = a*right
  end if
  
  end function flux_upwind


end module modules