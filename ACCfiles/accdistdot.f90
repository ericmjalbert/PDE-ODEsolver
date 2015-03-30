
real function accdistdot(n,x,ix,y,iy)
!----------------------------
implicit none
 integer :: n,ix,iy
 real, dimension(1+(n-1)*ix) :: x
 real, dimension(1+(n-1)*iy) :: y
 integer :: i 
 real :: ddt

ddt=0.
!$acc kernels loop private(ddt) reduction(+: ddt)
do i=1,n
 ddt=ddt+x(i)*y(i)
enddo
accdistdot=ddt
end function
