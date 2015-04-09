real function dotProd(n,u,v)
!----------------------------
implicit none
 integer, intent(in) :: n
 real, dimension(n), intent(in) :: u,v
 integer :: i 

 real :: sol

sol=0.
!$acc kernels loop private(sol) reduction(+: sol)
do i=1,n
 sol=sol+u(i)*v(i)
enddo
dotProd=sol
end function 
