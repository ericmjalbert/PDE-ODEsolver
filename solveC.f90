!==============================================================================
!   Solves the solution, C
!-----------------------------------------------------------------------------
!   Uses the trapizoidal rule for an ODE and solves for C.
!   Different gSelects give different values for b and c. 
!   gSelect = 1 --> g = nu*C/(k+C)
!   gSelect = 2 --> g = nu*C/k
!   gSelect = 3 --> g = nu*C/(k*M+C)
!
!   For when C_t = (r+s)*g(C)*f(C,M)*M - s*g(C)*M
!   gSelect = 4
!
!   gSelect = 5 --> C_t = -kC
!   gSelect = 6 --> C_t = -C/(k+C)*M/nu
!
!==============================================================================
subroutine solveC(M,Mnew,Csol,Cnew,row,col,n,k,y,m0,c0,gam,fSelect,&
                  tDel,nu,gSelect)
    implicit none 
    integer,intent(in) :: row,col,n,gSelect,fSelect
    real,intent(in) :: k,tDel,nu,y,m0,c0,gam
    real,dimension(n),intent(in) :: M,Mnew,Csol
    real,dimension(n),intent(out) :: Cnew

    integer :: i,j    ! grid index
    real :: b,c       ! quadratic equation terms
    real :: aux 
    
    aux = tDel*0.5*nu
    
    if (gSelect == 1) then
      !#omp parallel do private(b,c) shared(k,Csol,aux,Mnew,M,Cnew)
      do i = 1,n 
        b = k - Csol(i) + aux*(Mnew(i) + Csol(i)*M(i)/(k+Csol(i)) )
        c = -k*Csol(i) + aux*k*Csol(i)*M(i)/(k + Csol(i))
        Cnew(i) = 0.5 * (-b + SQRT(b*b - 4 * c))
      enddo
      !#omp end parallel do
  end if    
end subroutine solveC


