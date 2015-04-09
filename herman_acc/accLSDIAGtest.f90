program test_linear_solvers_in_diagonal_format
!----------------------------------------------------------------------------
! hje, for M6020
!----------------------------------------------------------------------------
! a test bed for solvers of linear systems Ax=b 
! where A is stored in sparse diagonal format
!----------------------------------------------------------------------------
! in the main program
! (1) a test problem is set up:  subroutine genDIAG
! (2) problem is solved: subroutine solveLSDIAG
! (3) report on results is given
!----------------------------------------------------------------------------
! the linear solver is supplied externally, regardless of the method used
! but they all should use the same interface to exchange data with the
! calling program; this is explicitly defined in an interface
!----------------------------------------------------------------------------
! parameters
!----------------------------------------------------------------------------
implicit none
!----------------------------------------------------
! this is where the variable declaration of he actual
! program begins
!----------------------------------------------------
 integer :: n,ndiag
 parameter (n=1000000, ndiag=5)
 real,dimension(n,ndiag) ::A 
 integer,dimension(ndiag) :: ioff
 real,dimension(n) :: rhs,sol,res,y,x,defect,truesol
 integer :: nit,stopcrit,i
 real :: err1,err2

 integer :: tstart,tfinish,clock_rate
 real :: elapsed_time


 ! initialise: set tolerances, max no iterations
 ! and initial guess
 !----------------------------------------------
 err1=1e-12;  err2=1e-12;  nit=n*100;  x=0.


 !(1) setup a test case: prescribe solution sol
 ! and generate a matrix A and rhs, such that A*sol=rhs
 !-----------------------------------------------------
 call genDIAG(n,ndiag,ioff,A,rhs,truesol)

 !(2) call jacobi method
 !-----------------------------------------------------------
 
 sol=0.
 call system_clock(tstart)
 call solve_linear_system(n,ndiag,ioff,A,rhs,sol) 
 call system_clock(tfinish, clock_rate)
 elapsed_time = float(tfinish-tstart) / float(clock_rate)


 !(3) report results: 
 !-------------------------------
 defect=truesol-sol ! compute defect
 write(*,'(A40,5(E14.7,X))')'max defekt, max/min sol/rhs: ',maxval(abs(defect)),maxval(sol),minval(sol),maxval(rhs),minval(rhs)
 write(*,'(A40,E14.7)') 'CPU time in seconds: ',elapsed_time
end program




subroutine genDIAG(n,ndiag,ioff,A,rhs,x)
!---------------------------------------
! author: hje for M6020
!---------------------------------------
! generates a matrix A in diagonal format
! and determines the rhs to go with known
! solution x
!----------------------------------------
! input
!   n:     problem size
!   ndiag: number of off-diagonals
! output
!   ioff:  diagonal offsets
!   A:     matrix values as n x ndiag array
!   rhs:   rhs for matrix A and solution x
!   x:     solution for Ax=b
!-------------------------------------------
implicit none
 integer, intent(in) :: n,ndiag
 integer,intent(out),dimension(ndiag) :: ioff
 real, dimension(n),intent(out) :: rhs
 real, dimension(n,ndiag),intent(out) :: a
 real,dimension(n),intent(out) :: x

 integer :: i,j,sw,nd
 real :: rn

  write(*,'(A30,I9,A10,I3,A10)') 'genDIAG creates fixed matrix with ',n,' rows and ',ndiag,' diagonals'

  ioff=(/ -1000,-1,0,1,1000 /)
  a(:,1) =-1.; 
  a(:,2) =-2.; 
  a(:,3) =6.01;
  a(:,4) = a(:,2)
  a(:,5) = a(:,1)
  


  ! define solution vector x
  x=1.

  ! compute rhs as rhs=A*x
  write(*,'(A31)') 'genDIAG creates right hand side'
  call amuxd(n,x,rhs,A,ndiag,ioff)
  write(*,'(A15)') 'genDIAG is done'


end subroutine



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
! diag   = real array containing the diagonals stored of A.
! idiag  = number of diagonals in matrix.
! diag   = real array of size (ndiag x idiag) containing the diagonals
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

      y=0.

      do j=1, idiag
         io = ioff(j)
         i1 = max0(1,1-io)
         i2 = min0(n,n-io)
         do i=i1,i2
           y(i) = y(i)+diag(i,j)*x(i+io)
         enddo
      enddo  
end subroutine amuxd
