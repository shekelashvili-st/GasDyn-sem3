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
  
  subroutine output_monitor_x(iu,xm,Time,u,rho,p,T,x_cent)
  integer,intent(in)       :: iu
  real(kind=rk),intent(in) :: Time, xm
  real(kind=rk),intent(in) :: u(0:), rho(0:), p(0:), T(0:), x_cent(0:)
  integer                  :: ipos, nx
  real(kind=rk)            :: um, rhom, pm, Tm
  nx = size(u) - 2
  
  ipos = locate(x_cent(1:nx),xm)
  um = interp1d(x_cent(1:nx),u(1:nx),xm,ipos)
  rhom = interp1d(x_cent(1:nx),rho(1:nx),xm,ipos)
  pm = interp1d(x_cent(1:nx),p(1:nx),xm,ipos)
  Tm = interp1d(x_cent(1:nx),T(1:nx),xm,ipos)
  
  write(iu,'(5(es23.16,x))') Time, um, rhom, pm, Tm
  
  end subroutine output_monitor_x
  
  integer function locate(xx,x)
  ! Locate a value in a sorted array
  real(kind=rk), dimension(:), intent(in) :: xx
  real(kind=rk), intent(in) :: x
  integer :: n,jl,jm,ju
  logical :: ascnd
  n=size(xx)
  ascnd = (xx(n) >= xx(1))
  jl=0
  ju=n+1
  do
     if (ju-jl <= 1) exit
     jm=(ju+jl)/2
     if (ascnd .eqv. (x >= xx(jm))) then
        jl=jm
     else
        ju=jm
     end if
  end do

  if (x == xx(1)) then
     locate = 1
  else if (x == xx(n)) then
     locate = n-1
  else if(ascnd.and. (x > xx(n) .or. x < xx(1))) then
     locate = -1
  else if(.not.ascnd.and. (x < xx(n) .or. x > xx(1))) then
     locate = -1
  else
     locate = jl
  end if

  end function locate
  
  real(kind=rk) function interp1d_lin(x1,y1,x2,y2,xval) result(yval)
    real(kind=rk),intent(in) :: x1,y1,x2,y2,xval
    real(kind=rk) :: frac
    frac = ( xval - x1 ) / ( x2 - x1 )
    yval = y1 + frac * ( y2 - y1 )
  end function interp1d_lin
  
  real(kind=rk) function interp1d(x,y,xval,ipos) result(yval)
    integer :: n
    ! the size of the array

    real(kind=rk),dimension(:),intent(in) :: x,y
    ! the x and y arrays

    real(kind=rk),intent(in) :: xval
    ! the value at which to interpolate y
    integer, intent(in)      :: ipos
    ! position of x value in x array 

    n = size(x)


    if(ipos == -1) then
        write(0,'("ERROR: Interpolation out of bounds : ",ES11.4," in [",ES11.4,":",ES11.4,"]")') xval,x(1),x(n)
    end if

    if( ipos < n .and. ipos > 0) then
       yval = interp1d_lin(x(ipos), y(ipos), x(ipos+1), y(ipos+1), xval)
    else if(ipos == n) then
       yval = y(n)
    else if(ipos == 0) then
       yval = y(1)
    else
       write(0,'("ERROR: Unexpected value of ipos : ",I0)') ipos
    end if
  end function interp1d 
  
  
end module io