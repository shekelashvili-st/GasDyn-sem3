module io
  use modules, only: rk
  implicit none
  
  contains
  
  subroutine output_field(x_cent,u,file_name)
  real(kind=rk), intent(in) :: x_cent(0:), u(0:)
  character(*), intent(in)  :: file_name
  integer			        :: iu, i
  
  open(newunit=iu, file=file_name)
  do i=0, size(x_cent)-1
	write(iu,'(2(es23.16,x))') x_cent(i), u(i)
  end do
  close(iu)
 
  end subroutine output_field
  
end module io