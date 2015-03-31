subroutine dotProd(n,u,v,sol)
implicit none
    integer, intent(in) :: n
    real, dimension(n), intent(in) :: u,v
    integer :: i

    real, intent(out) :: sol

    sol = 0.
    !$omp parallel do reduction(+:sol)
    do i=1,n
        sol = sol + u(i) * v(i)
    enddo
    !$omp end parallel do

end subroutine
