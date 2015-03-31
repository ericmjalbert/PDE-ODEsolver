subroutine amuxd (n,x,y,diag,idiag,ioff) 
!-----------------------------------------------------------------------
!        A times a vector in Diagonal storage format (DIA) 
!        f90/f95 version of the sparskit f77 subroutine
!----------------------------------------------------------------------- 
! multiplies a matrix by a vector when the original matrix is stored 
! in the diagonal storage format.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! ndiag  = integer. The first dimension of array adiag as declared in
!         the calling program.
!         (obscolete (=n always)
! idiag  = integer. The number of diagonals in the matrix.
! diag   = real array containing the diagonals stored of A.
! idiag  = number of diagonals in matrix.
! diag   = real array of size (ndiag x idiag) containing the diagonals
!          
! ioff   = integer array of length idiag, containing the offsets of the
!   	   diagonals of the matrix:
!          diag(i,k) contains the element a(i,i+ioff(k)) of the matrix.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=A*x
!
!-----------------------------------------------------------------------
implicit none
  integer, intent(in)::  n, idiag
  integer, intent(in),dimension(idiag) :: ioff
  real, dimension(n), intent(in) :: x    
  real, dimension(n,idiag), intent(in) :: diag
  real, dimension(n), intent(out) :: y
  integer :: j, io, i1, i2, i       

      !$acc kernels loop present(y)
      do i=1,n
        y(i)=0.  
      enddo
      !$acc end kernels

      do j=1, idiag
         io = ioff(j)
         i1 = max0(1,1-io)
         i2 = min0(n,n-io)

         !$acc parallel loop copyout(y) present(x,diag) firstprivate(io)
         do i=i1,i2
           y(i) = y(i)+diag(i,j)*x(i+io)
         enddo
         !$acc end parallel loop  
      enddo
    
end subroutine amuxd
