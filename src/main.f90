program main
  use modules, only : rk, init_Riemann, flux_cent, flux_upwind
  use io, only      : output_field
  use godunov, only : flux_godunov
  implicit none
  
  !Parameters
  character(*),parameter   :: input_file='params.nml'
  real(kind=rk)			   :: L, CFL, a, U_L, U_R
  integer	   			   :: nx, nt, L_size
  !Work arrays & variables
  real(kind=rk),allocatable:: x_cent(:), u(:), u_n(:)
  real(kind=rk)			   :: dx, dt
  !Local variables
  integer				   :: lala
  !Service variables
  integer	   			   :: i, j, k, iu
 
  !/////////////////////////////////////////////////////////// 
  !Read input data
  namelist /params/ L, CFL, a, nx, nt, U_L, U_R, L_size
  open(newunit=iu, file=input_file)
  read(iu, nml=params)
  close(iu)  
  
  !Prepare work arrays, set initial conditions
  allocate(x_cent(0:nx+1),u(0:nx+1),u_n(0:nx+1))
  dx = L/nx
  dt = CFL * dx / a
  x_cent = [(dx*(i+0.5_rk), i=-1,nx)]
  u = 0 ; u_n = 0
  call init_Riemann(u,L_size,U_L,U_R)

  !Output initial field
  call output_field(x_cent,u,'init.dat')  
  
  !Solve equation
  time: do k = 1, nt
	do i=1, nx
		u_n(i) = u(i) - dt/dx*(flux_upwind(u(i),u(i+1),a) - flux_upwind(u(i-1),a*u(i),a))
	end do
	u(1:nx) = u_n(1:nx)
  end do time
  
  !Output solution
  call output_field(x_cent,u,'sol.dat')
  
  

end program main