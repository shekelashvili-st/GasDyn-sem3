module io
  use modules, only: rk
  implicit none
  
  contains
  
  subroutine output_field(x_cent,u,rho,p,T,file_name)
  real(kind=rk), intent(in) :: x_cent(0:), u(0:), rho(0:), p(0:), T(0:)
  character(*), intent(in)  :: file_name
  integer                   :: iu, i
  
  open(newunit=iu, file=file_name)
  do i=1, size(x_cent)-2
    write(iu,'(5(es23.16,x))') x_cent(i), u(i), rho(i), p(i), T(i)
  end do
  close(iu)
 
  end subroutine output_field
  
  
  subroutine output_monitor(iu,i_x,Time,u,rho,p,T)
  integer,intent(in)       :: iu, i_x
  real(kind=rk),intent(in) :: Time
  real(kind=rk),intent(in) :: u(0:), rho(0:), p(0:), T(0:)
  
  write(iu,'(5(es23.16,x))') Time, u(i_x), rho(i_x), p(i_x), T(i_x)
  
  end subroutine output_monitor
  
end module io