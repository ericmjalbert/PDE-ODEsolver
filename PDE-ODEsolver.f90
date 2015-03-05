!==============================================================================
!   PDE - ODE Coupled Solver
!-----------------------------------------------------------------------------
!     This fortran code solves the PDE - ODE coupled system that describes 
!   the growth of Clostridium Thermocellum and its consumption of the carbon
!   substrait.
!   
!   The system is based off the following system:
!     M_t = nabla (D(M) nabla M ) + f(C,M) 
!     C_t = - g(C,M)
!   where M is the biomass of C. Thermocellum and C is the concentration of
!   Carbon. D(M) is the diffusion coeffient for the biomass movement. f(C,M) 
!   is the growth and death term for the biomass. g(C,M) is the consumption
!   of carbon substrait.
!-----------------------------------------------------------------------------
!     The method of solving is by trapzidral rule for C, and by finite
!   difference for M 
!==============================================================================

program cThermoPDEODE
!    use omp_lib
    implicit none
    INTERFACE
        subroutine solveLSDIAG(n,ndiag,ioff,a,sol,rhs,nit,err1,err2,stopcrit)
        !---------------------------------------------------------------------
        ! input:  n       problem size
        !         ndiag:  number of diagonals
        !         ioff:   offsets (distance of sub diagonals to main diagonal)
        !         Mnew:      matrix values
        !         sol:    initial guess for iteration ( overwritten by result)
        !         rhs:    the righ hand side of the linear system
        !         nit:    max num of iter to be carried out (overwritten)
        !         err1:   tol for 1st stopping criterion (will be overwritten)
        !         err2:   tol for 2nd stopping criterion (will be overwritten)
        ! output: sol:    solution of Ax=b
        !         nit:    number of iterations taken
        !         err1:   computed value for 1st stopping criterion
        !         err2:   computed value for 2nd stopping criterion
        !         stopcrit: value that tells which stopping criti activated
        !---------------------------------------------------------------------
        implicit none
        integer, intent(in)               :: n,ndiag
        real,dimension(n,ndiag),intent(in):: a
        integer, dimension(ndiag),intent(in)::ioff
        real,dimension(n),intent(in)      :: rhs
        real,dimension(n),intent(inout)   :: sol
        real,intent(inout)                :: err1,err2
        integer,intent(inout)             :: nit
        integer,intent(out)               :: stopcrit
        end subroutine
    END INTERFACE
  
    !=======================
    ! Variable Descriptions
    !=======================
    
    ! filename Variables
    character(100) :: filename

    ! Problem Parameters
    real :: length,lambda,depth,height,m0,c0,gam
    real :: k,delta,nu,mu,y
    integer :: alpha,beta
    integer :: fSelect, dSelect, gSelect, MinitialCond

    ! Numerical Method Parameters
    integer :: pSize,row,col,n,ndiag,nit,nOuts
    real :: tEnd,xLen,yLen,tDel,yDel,xDel,e1,e2,esoln

    ! Solution variables
    real,dimension(:),allocatable  :: C, M
    
    ! Reporting variables
    real :: avgIters, maxIters
    real :: avgNit, maxNit
    real :: startTime, endTime
    
    write(*,*) "Enter parameter file name: ex. 'parameter.txt' "
    write(*,'(A)', advance="no") "    "
    read(*,*) filename

    write(*,*) "Setting Parameters that shouldn't change"
    ndiag = 5
    write(*,'(A, I5)') "     ndiag = ", ndiag

    write(*,*) "Getting problem size from file"
    call getProbSize(pSize, filename, len(filename)) 
    write(*,*) "    pSize = ", pSize
    row = pSize
    col = 4
    n = pSize * col
    write(*,*) "    row   = ", row
    write(*,*) "    col   = ", col        
    write(*,*) "    n     = ", n

    write(*,*) "Opening Parameter file"
    call paramSet(length, lambda, depth, height, m0, c0, alpha, beta, gam, k, &
                  mu, y, nu, delta, nOuts, tEnd, xLen, yLen, tDel, e1, e2, & 
                  esoln, fSelect, dSelect, gSelect, MinitialCond, filename, &
                  len(filename))
    yDel = yLen/real(row)
    xDel = yDel
    nit = n * 100
    write(*,*) "Parameters set:"
    write(*,*) "    length      = ", length
    write(*,*) "    lambda      = ", lambda
    write(*,*) "    depth       = ", depth
    write(*,*) "    height      = ", height
    write(*,*) "    m0          = ", m0
    write(*,*) "    c0          = ", c0
    write(*,*) "    alpha       = ", alpha
    write(*,*) "    beta        = ", beta
    write(*,*) "    gam         = ", gam
    write(*,*) "    y           = ", y
    write(*,*) "    mu          = ", mu
    write(*,*) "    nu          = ", nu
    write(*,*) "    delta       = ", delta
    write(*,*) "    nOuts       = ", nOuts
    write(*,*) "    tEnd        = ", tEnd
    write(*,*) "    xLen        = ", xLen
    write(*,*) "    yLen        = ", yLen
    write(*,*) "    tDel        = ", tDel
    write(*,*) "    e1          = ", e1
    write(*,*) "    e2          = ", e2
    write(*,*) "    nit         = ", nit
    write(*,*) "    xDel        = ", xDel
    write(*,*) "    yDel        = ", yDel
    write(*,*) "    fSelect     = ", fSelect
    write(*,*) "    dSelect     = ", dSelect
    write(*,*) "    gSelect     = ", gSelect
    write(*,*) "    MinitialCond= ", MinitialCond
    
    write(*,*) "Allocating the size of arrays"
    allocate(C(n),M(n))
    write(*,'(A,I6,A)') "    C and M are now dimension(", n, ") arrays"

    write(*,*) "Setting Initial Conditions"
    call setInitialConditions(C,M,row,col,n,depth,height,yDel,MinitialCond)
        write(*,*) "    C = 1"
        if (MinitialCond == 1) then
            write(*,*) "    M has non-homogenous Initial Conditions"
        else if (MinitialCond == 2) then
            write(*,*) "    M has homogenous initial conditions"
            write(*,*) "        M = ", height
        endif
        

    write(*,*) "Starting Solver"
    call cpu_time(startTime)
    call solveOrder2(tEnd,nOuts,tDel,n,row,col,M,C,yLen,xDel,yDel,&
       c0,m0,y,nu,gam,alpha,beta,k,delta,ndiag,e1,e2,nit,esoln,dSelect,&
       fSelect,gSelect,avgIters,maxIters,avgNit,maxNit)
    call cpu_time(endTime)
    
!    call solveOrder1(tEnd,nOuts,tDel,n,row,col,M,C,xLen,yLen,xDel,yDel&
!       ,c0,m0,nu,gam,alpha,beta,k,delta,ndiag,e1,e2,nit,dSelect,fSelect)
  
    ! Report Statistics
    write(*,*) "-------------------------------------------------------------"
    write(*,*) "Statsitcs:"
    write(*,*) "-------------------------------------------------------------"
    write(*,*) "Time to compute = ", endTime-startTime
    write(*,*) ""
    write(*,*) "Avg Iters for iterating betn. soln. =", avgIters
    write(*,*) "Max Iters for iterating betn. soln. =", maxIters
    write(*,*) "Avg Iters for linear solver =", avgNit
    write(*,*) "Max Iters for linear solver =", maxNit  
    write(*,*) "Writing statistics to file"
    call reportStats(avgIters,maxIters,avgNit,maxNit,endTime-startTime)
    
    write(*,*) "Execution completed"
  
end program


!==============================================================================
!   Sets the problem size and number of diagonal elements
!-----------------------------------------------------------------------------
!     This must be done first since the problem size is used in the variable 
!   declaration of the arrays.
!==============================================================================
subroutine getProbSize(pSize, filename, nameLen) 
    implicit none
    integer,intent(out) :: pSize
    integer,intent(in)  :: nameLen
    character(nameLen),intent(in) :: filename
    character :: dum   ! Dummy variable
    write(*,*) filename
    open(UNIT = 19, FILE = filename, STATUS = "old", ACTION = "read")
    read(19,*); read(19,*); read(19,*)  ! Skip first 3 lines
    read(19,*) dum, dum, pSize
    close(19)
end subroutine getprobSize


!==============================================================================
!   Sets the all the parameter values
!-----------------------------------------------------------------------------
!     Opens the parameter.txt file, which holds all the parameter values and 
!   function uses, and sets the variables accordingly.
!     Takes in all the parameters on entry and outputs them with the appropriate
!   value.
!==============================================================================
subroutine paramSet(length, lambda, depth, height, m0, c0, alpha, beta, gam, &
                    k, mu,y, nu, delta, nOuts, tEnd, xLen, yLen, tDel, e1, e2,&
                    esoln, fSelect, dSelect, gSelect, MinitialCond, filename, &
                    nameLen)
    implicit none
    real, intent(out) :: lambda, depth, height, m0, c0, gam, k, delta, nu, mu
    real, intent(out) :: tEnd, xLen, yLen, tDel, e1, e2, length, esoln, y
    integer, intent(out) :: alpha, beta, nOuts
    integer, intent(out) :: fSelect, dSelect, gSelect, MinitialCond
    integer, intent(in) :: nameLen
    character(nameLen), intent(in) :: filename

    real :: kBar, nuBar, delBar
    character :: dum, dum2      ! Dummy variable
    
    open(UNIT = 19, FILE = filename, STATUS = "old", ACTION = "read")
    read(19,*); read(19,*); read(19,*); read(19,*)  ! Skip first 4 lines

    read(19,*) dum, dum2, length
    read(19,*) dum, dum2, lambda
    read(19,*) dum, dum2, depth
    read(19,*) dum, dum2, height
    read(19,*) dum, dum2, m0
    read(19,*) dum, dum2, c0
    read(19,*) dum, dum2, alpha
    read(19,*) dum, dum2, beta
    read(19,*) dum, dum2, gam
    read(19,*) dum, dum2, mu
    read(19,*) dum, dum2, y
    read(19,*) dum, dum2, kBar
    read(19,*)                          ! Skip nuBar 
    read(19,*) dum, dum2, delBar
    read(19,*); read(19,*); read(19,*)  ! Skip k, delta, and nu
    read(19,*) dum, dum2, nOuts
    read(19,*) dum, dum2, tEnd
    read(19,*) dum, dum2, yLen
    read(19,*)                          ! Skip xLen
    read(19,*) dum, dum2, tDel
    read(19,*); read(19,*)              ! Skip yDel and xDel
    read(19,*) dum, dum2, e1
    read(19,*) dum, dum2, e2
    read(19,*) dum, dum2, esoln
    read(19,*)                          ! Skip nit
    read(19,*); read(19,*); read(19,*)  ! Skip 3 lines
    read(19,*) dum, dum2, fSelect
    read(19,*) dum, dum2, dSelect
    read(19,*) dum, dum2, gSelect
    read(19,*) dum, dum2, MinitialCond
    
    nuBar = mu*m0*1.5873    ! 1/0.63 = 1.5873
    k = kBar/c0
    delta = delBar/(mu*length*length)
    nu = nuBar/(mu*c0)
    xLen = yLen*lambda

    close(19)
    
end subroutine paramSet


!==============================================================================
!   Sets the initial conditions for the system
!-----------------------------------------------------------------------------
!   For C, we have homogenous initial conditions, trivial to do
!   For M, the innoculation point has a smooth curve to avoid sharp numerical 
!     artifacts in the system. This is done with a polynomial f(x) = a*x^8+b, 
!     this function is computed based on the depth and height parameters, 
!     calculating b = height and a = b/(depth)^8
!==============================================================================
subroutine setInitialConditions(C,M,row,col,n,depth,height,yDel,MinitialCond)
    implicit none
    integer,intent(in) :: row,col,n,MinitialCond
    real,intent(in) :: depth, height, yDel
    real,dimension(n),intent(out) :: C,M

    integer :: i,j
    real :: f,a       ! function for IC curve
    
    C = 1.; j = 0

    if (MinitialCond == 1) then
        M = 0
        a = -height/(depth)**4
        f = a*j*yDel + height
        do while(f > 0)
            do i = 1,col
                M(i+j*col) = f
            enddo
            j = j+1
            f = a*(j*yDel)**4 + height
        enddo
    else if (MinitialCond == 2) then
        M = height
    else if (MinitialCond == 3) then
        j = row ! have multiple innoculation points
    endif


end subroutine setInitialConditions


!==============================================================================
!   Function f(C,M)
!----------------------------------------------------------------------
!   fSelect == 7 should be used with gSelect == 2
!   fSelect == 8 should be used with gSelect == 3
!   otherwise gSelect == 1
!==============================================================================
subroutine fFunc(M,C,f,k,y,nu,m0,c0,gam,fSelect)
    implicit none    
    integer, intent(in) :: fSelect
    real, intent(in) :: M,C,k,m0,c0,gam,nu,y
    real, intent(out) :: f
    real :: eps
    eps = 0.00000001
    if (fSelect == 1) f = y*C/ (k + C) * (1 - (M*m0/(C*c0+eps))**gam)
    if (fSelect == 2) f = 1.
    if (fSelect == 3) f = (1. - M)
    if (fSelect == 4) f = (1. - C)*(C/(M+eps))
    if (fSelect == 5) f = y*nu*C/(k+C)*(1. - (M/(C+eps))**gam)
    if (fSelect == 6) f = y*nu*C/(k+C)*(1. - (M*(k+C)/(nu*C+eps))**gam)
    if (fSelect == 7) f = y*nu/k*C*(1. - (M*k/(nu*C+eps))**gam)
    if (fSelect == 8) f = y*nu*C/(k*M+C)*(1. -(M*(k*M+C)/(nu*C+eps))**gam)
    if (fSelect == 9) f = y*nu*(C**gam - M**gam)* C**(1-gam)/(k+C)
    if (fSelect == 10) f = y*nu*((C/(k+C))**gam - M**gam)*(C/(k+C))**(1-gam)
    if (fSelect == 11) f = y*nu*(C**gam - M**gam)*C/(k+C)/k**gam
    if (fSelect == 12) f = C/(k+C)-0.1*nu 
    if (fSelect == 13) f = C/(k+C)
    if (fSelect == 14) f = C/(k*M+C+eps)-0.1*nu
    if (fSelect == 15) f = C/(k*M+C+eps)
end subroutine fFunc


!==============================================================================
!   Function d(M)
!==============================================================================
subroutine dFunc(M,d,delta,alpha,beta,dSelect)
    implicit none    
    integer,intent(in) :: alpha, beta, dSelect
    real, intent(in) :: M, delta
    real, intent(out) :: d
    if (dSelect == 1) d = delta
    if (dSelect == 2) d = delta * M** alpha
    if (dSelect == 3) d = delta * M**alpha / ((1 - M)**beta)
end subroutine dFunc


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
    integer :: g      ! current grid point
    real :: b,c       ! quadratic equation terms
    real :: f         ! f(C,M) value at gridpoint
    real :: aux 
    
    real :: r,s
    r = 1
    s = 0.5

    aux = tDel*0.5*nu
    
    do i = 1,row
        do j = 1,col
            g = j + (i - 1) * col
            
            if (gSelect == 1) then
                b = k - Csol(g) + aux*(Mnew(g) + Csol(g)*M(g)/(k+Csol(g)) )
                c = -k*Csol(g) + aux*k*Csol(g)*M(g)/(k + Csol(g))
                Cnew(g) = 0.5 * (-b + SQRT(b*b - 4 * c))
            else if (gSelect == 2) then
                Cnew(g) = (k - aux*M(g))*Csol(g) / (aux*Mnew(g) + k)
            else if (gSelect == 3) then
                b = k*Mnew(g) - Csol(g) + aux*(Mnew(g) &
                    + Csol(g)*M(g)/(k*Mnew(g)+Csol(g)) )
                c = aux*k*Mnew(g)*Csol(g)*M(g)/(k*Mnew(g) + Csol(g)) &
                    - k*Mnew(g)*Csol(g)
                Cnew(g) = 0.5 * (-b + SQRT(b*b - 4 * c))
            else if (gSelect == 4) then
                call fFunc(Mnew(g),Csol(g),f,k,y,nu,m0,c0,gam,fSelect)
                b = k - Csol(g) - tDel*Mnew(g)*((s-r)*f - s)
                c = -k*Csol(g)
                Cnew(g) = 0.5 * (-b + SQRT(b*b - 4*c))
            else if (gSelect == 5) then
                Cnew(g) = Csol(g)*(1-aux)/(1+aux)
            endif

            if (Cnew(g) .le. 0) then
                Cnew(g) = 0
            endif
        enddo
    enddo

end subroutine solveC


!==============================================================================
!   Generates matrix M in diagonal format for the current timestep
!==============================================================================
subroutine GenMatrixM(M,C,MatrixM,Mioff,Mrhs,row,col,n,ndiag,delta,nu,alpha, &
                      beta,k,m0,c0,gam,tDel,xDel,yDel,dSelect,fSelect,yConst)
    implicit none
    integer,intent(in) :: row,col,n,ndiag,alpha,beta,dSelect,fSelect
    real,intent(in) :: delta,k,m0,c0,gam,xDel,yDel,tDel,nu,yConst
    real,dimension(n),intent(in) :: M,C

    real,dimension(n,ndiag),intent(out) :: MatrixM
    real,dimension(n),intent(out) :: Mrhs
    integer,dimension(ndiag),intent(out) :: Mioff

    real :: xCof,yCof
    real :: f
    real,dimension(n) :: diff
    integer :: x,y,g

    real :: tDela
    tDela = tDel      ! Used for testing purposes

    xCof = 1/(xDel*xDel)
    yCof = 1/(yDel*yDel)

    Mioff = (/ -col, -1, 0, 1, col/)
    MatrixM(:,:) = 0

    ! Compute all the diffusion coefficients
    do y = 1,n
        call dFunc(M(y), diff(y), delta, alpha, beta, dSelect)
    enddo
  
    do y = 1,row
        write(*,*) "y = ", y
        do x = 1,col
            g = x + (y - 1) * col
      
            if (y .NE. 1) then
                MatrixM(g,1) = MatrixM(g,1) - yCof*0.5*(diff(g-col)+diff(g))
                MatrixM(g,3) = MatrixM(g,3) + yCof*0.5*(diff(g-col)+diff(g))
            else
                MatrixM(g,5) = MatrixM(g,5) - yCof*0.5*(diff(g+col)+diff(g))
                MatrixM(g,3) = MatrixM(g,3) + yCof*0.5*(diff(g+col)+diff(g))
            endif
              
            if (x .NE. 1) then
                MatrixM(g,2) = MatrixM(g,2) - xCof*0.5*(diff(g-1)+diff(g))
                MatrixM(g,3) = MatrixM(g,3) + xCof*0.5*(diff(g-1)+diff(g))
            else
                MatrixM(g,4) = MatrixM(g,4) - xCof*0.5*(diff(g+1)+diff(g))
                MatrixM(g,3) = MatrixM(g,3) + xCof*0.5*(diff(g+1)+diff(g))
            endif

            if (x .NE. col) then
                MatrixM(g,4) = MatrixM(g,4) - xCof*0.5*(diff(g+1)+diff(g))
                MatrixM(g,3) = MatrixM(g,3) + xCof*0.5*(diff(g+1)+diff(g))
            else
                MatrixM(g,2) = MatrixM(g,2) - xCof*0.5*(diff(g-1)+diff(g))
                MatrixM(g,3) = MatrixM(g,3) + xCof*0.5*(diff(g-1)+diff(g))
            endif
              
            if (y .NE. row) then
                MatrixM(g,5) = MatrixM(g,5) - yCof*0.5*(diff(g+col)+diff(g))
                MatrixM(g,3) = MatrixM(g,3) + yCof*0.5*(diff(g+col)+diff(g))
            else
                MatrixM(g,1) = MatrixM(g,1) - yCof*0.5*(diff(g-col)+diff(g))
                MatrixM(g,3) = MatrixM(g,3) + yCof*0.5*(diff(g-col)+diff(g))
            endif
            
            call fFunc(M(g),C(g),f,k,yConst,nu,m0,c0,gam,fSelect)
            MatrixM(g,3) = MatrixM(g,3) - f + (1/tDel)
            Mrhs(g) = M(g)/tDel
write(*,"(I5,F10.6,F10.6,F10.6,F10.6,F10.6,F10.6,F10.6)") g,MatrixM(g,1),MatrixM(g,2),MatrixM(g,3),MatrixM(g,4),MatrixM(g,5),Mrhs(g)
        enddo
    enddo
    write(*,*) "[Debug 496] Printing Matrix"
    call exit(1)

end subroutine GenMatrixM


!==============================================================================
!   Solve Order 2
!-----------------------------------------------------------------------------
!   1. Solves for M_{i+1} using C_i and M_i
!   2. Solves for C_{i+1} using C_i, M_i, and M_{i+1}
!   ... repeat until convergence
!==============================================================================
subroutine solveOrder2(tEnd,nOuts,tDel,n,row,col,M,C,yLen,xDel,yDel,&
    c0,m0,y,nu,gam,alpha,beta,k,delta,ndiag,e1,e2,nit,eSoln,dSelect,fSelect,&
    gSelect,avgIters,maxIters,avgNit,maxNit)
  implicit none
  integer,intent(in) :: nOuts,n,row,col,alpha,beta,ndiag
  integer,intent(in) :: dSelect,fSelect,gSelect
  real,intent(in) :: tEnd,tDel,yLen,xDel,yDel,c0,m0,gam,k,delta,nu,y
  real,intent(in) :: eSoln
  real,intent(inout) :: e1,e2
  integer,intent(inout) :: nit
  real,dimension(n),intent(out) :: M,C
  real,intent(out) :: avgIters,maxIters,avgNit,maxNit
  
  integer :: counter        ! Counts the num. of iter. in system solver
  integer :: endLoop        ! endLoop = 1 -> solving loop can stop
  integer :: filter         ! Controls frequency of outputs written
  
  real :: totalMassM        ! Total Biomass over region
  real :: totalMassC        ! Total Substrait over region
  
  real,dimension(n) :: Mnew
  real,dimension(n) :: Cnew
  real,dimension(n,ndiag) :: MatrixM
  real,dimension(n) :: Mrhs
  real,dimension(n) :: Cprev, Mprev
  integer,dimension(ndiag) :: Mioff 
  integer :: stopcritria 
  integer :: stat 
  real :: diffC, diffM 
  integer :: countIters 
  real :: peak, height, intfac  ! Track Interface and wave peak
  
  filter = int(tEnd/(nOuts*tDel)) +1
  endLoop = 0
  counter = 0
  countIters = 0
  avgIters = 1 
  avgNit = 1 
  
  open(UNIT = 124, IOSTAT = stat, FILE = "total.dat", STATUS = "old")
  if (stat .EQ. 0) close(124, STATUS = "delete")
  open(UNIT = 120, FILE = "total.dat", POSITION = "append", ACTION = "write")
  
  open(UNIT = 124, IOSTAT = stat, FILE = "peakInfo.dat", STATUS = "old")
  if (stat .EQ. 0) close(124, STATUS = "delete")
  open(UNIT = 121, FILE = "peakInfo.dat", POSITION = "append", ACTION = "write")
    
  write(*,*) "   time    avgIter  maxIter      avgNit  maxNit        avgM        avgC"

  do while(counter * tDel <= tEnd)
    ! Output every 100 times more then nOuts for the peak info
    if (MOD(counter, int(filter/100+1)) == 0) then
      call calcMass(M,totalMassM,n,row,col)
      call calcMass(C,totalMassC,n,row,col)
      call calcPeakInterface(M, row, col, peak, height, intfac, yLen)
      write(120,*) counter*tDel, totalMassM, totalMassC
      write (121,*) tDel*counter, peak, height, intfac
    endif 
  
    ! Write to file / report Total Mass
    if (MOD(counter, filter) == 0) then
      call printToFile2D(n,row,col,M,C,yLen)
      if (counter == 0) then
        write(*,'(F8.2,F12.2,I8,F12.2,I8,F12.6,F12.6)') 0.0, 0.0, &
          int(maxIters), 0.0, int(maxNit), totalMassM, totalMassC
      else
        write(*,'(F8.2,F12.2,I8,F12.2,I8,F12.6,F12.6)') tDel*counter, &
          real(avgIters/counter), int(maxIters), real(avgNit/avgIters), &
          int(maxNit),totalMassM, totalMassC
      endif
    endif

    ! Output every 100 times more then nOuts for the peak info
    if (MOD(counter, int(filter/100+1)) == 0) then
      call calcPeakInterface(M, row, col, peak, height, intfac, yLen)
      write (121,*) tDel*counter, peak, height, intfac
    endif
 
    diffC = 1
    diffM = 1
    countIters = 0
    Cprev = C
    Mprev = M

    do while(diffC + diffM > eSoln)
        nit = 100*n
        e1 = 1.e-12
        e2 = e1
        
        ! Solve M
        call GenMatrixM(Mprev,C,MatrixM,Mioff,Mrhs,row,col,n,ndiag,delta,nu,&
                    alpha,beta,k,m0,c0,gam,tDel,xDel,yDel,dSelect,fSelect,y)
        call solveLSDIAG(n,ndiag,Mioff,MatrixM,Mnew,Mrhs,nit,e1,e2,stopcritria)
  
        ! Solve C
        call solveC(Mprev,Mnew,Cprev,Cnew,row,col,n,k,y,m0,c0,gam,fSelect,&
                  tDel,nu,gSelect)

        if(nit > maxNit) maxNit = nit
        avgNit = avgNit + nit

        call calcDiff(diffC, C, Cnew, row, col)
        call calcDiff(diffM, M, Mnew, row, col)

        C = Cnew
        M = Mnew
        countIters = countIters+1
        write(*,*) countIters, diffC, diffM, C(12),M(12)
    enddo
    if(countIters > maxIters) maxIters = countIters
    avgIters = avgIters + countIters
    
    counter = counter + 1
  enddo
  avgNit = avgNit/(avgIters) ! avgIters right now is the total
  avgIters = avgIters/counter

  close(120)
  close(121)
  
end subroutine solveOrder2


!==============================================================================
!   Calculate the Total Mass
!-----------------------------------------------------------------------------
!   Uses Riemann Sums kind of concept to check if the program runs correctly
!   since the total mass of the region can be computed and checked.
!==============================================================================
subroutine calcMass(X,totalMass,n,row,col)
    implicit none
    integer,intent(in) :: n,row,col
    real,dimension(n),intent(in) :: X
  
    real,intent(out) :: totalMass
  
    integer :: i,j,g
  
    totalMass = 0
  
    do i = 1,row
        do j = 1,col
            g = j +(i - 1) * col
            totalMass = totalMass + X(g)
        enddo
    enddo
  
    totalMass = totalMass / n
  
end subroutine calcMass


!==============================================================================
!   Prints out the solution in a format accepted by gnuplot, used for graphing
!-----------------------------------------------------------------------------
!   Runs through the grid, row-by-row. The MOD and filter act to reduce the 
!     number of grid points written, useful when comparing different grid 
!     sizes.
!-----------------------------------------------------------------------------
!   Input:
!     n      =  The problem size
!     row    =  The number of rows
!     col    =  The number of columns
!     M      =  The solution vector for biomass
!     C      =  The solution for substrait
!     xLen   =  The length of each x-grid slot
!     yLen   =  The length of each y-grid slot
!   Output:
!     none
!==============================================================================
subroutine printToFile(n,row,col,M,C,xLen,yLen)
    implicit none

    integer,intent(in) :: n, row, col
    real,intent(in) ::xLen, yLen
    real,dimension(n),intent(in) :: M,C

    integer :: p
    integer :: i,j
    integer :: stat
  
    !-------------------------------------------
    ! Deletes the old output file if it exist
    !-------------------------------------------
    open(UNIT = 123, IOSTAT = stat, FILE = "output.dat", STATUS = "old")
    if (stat .EQ. 0) close(123, STATUS = "delete")

    open(UNIT = 11, FILE = "output.dat", POSITION = "append", ACTION = "write")

    do, i=1,row
        do, j=1,col
            p = (j + (i-1)*col)
            write(11,*) real(j-1)/real(col-1)*xLen, &
                        real(i-1)/real(row-1)*yLen, M(p), C(p)
        enddo
        write(11,*) ' '
    enddo
    write(11,*) 
    
      
end subroutine printToFile


!==============================================================================
!   Prints out the solution in a 2D format, used for graphing the Travelling
!     Wave Example.
!-----------------------------------------------------------------------------
!   Runs through the grid, row-by-row. The MOD and filter act to reduce the 
!     number of grid points written, useful when comparing different grid 
!     sizes.
!   Unique here is that the average of the x-axis is taken so that the system 
!     can be reduced to just y. Also written are the max and min for each y 
!     value; this is used for showing that the system can be reduced.
!-----------------------------------------------------------------------------
!   Input:
!     n      =  The problem size
!     row    =  The number of rows
!     col    =  The number of columns
!     M      =  The solution vector for biomass
!     C      =  The solution for substrait
!     xLen   =  The length of each x-grid slot
!     yLen   =  The length of each y-grid slot
!   Output:
!     none
!==============================================================================
subroutine printToFile2D(n,row,col,M,C,yLen)
    implicit none

    integer,intent(in) :: n, row, col
    real,intent(in) :: yLen
    real,dimension(n),intent(in) :: M,C

    integer :: p
    integer :: i,j
    integer :: stat
    real :: averageM, averageC
    real :: maxM, minM, maxC, minC
    real :: y
  
    !-------------------------------------------
    ! Deletes the old output file if it exist
    !-------------------------------------------
    open(UNIT = 124, IOSTAT = stat, FILE = "2D_output.dat", STATUS = "old")
    if (stat .EQ. 0) close(124, STATUS = "delete")

    open(UNIT = 12, FILE = "2D_output.dat", POSITION = "append", ACTION = "write")
    
    do, i=0,row
        averageM = 0
        averageC = 0
        maxM = 0
        maxC = 0
        minM = 1
        minC = 1
        do, j=0,col
            if (i == 0 .and. j == 0) then
                p = 1
            else if (i == 0) then
                p = j
            else if (j == 0) then
                p = 1 + (i-1)*col
            else
                p = j + (i-1)*col
            endif
            averageM = averageM + M(p)
            averageC = averageC + C(p)
            if(M(p) .ge. maxM) then 
                maxM = M(p)
            endif
            if(M(p) .le. minM) then 
                minM = M(p)
            endif
            if(C(p) .ge. maxC) then
                maxC = C(p)
            endif
            if(C(p) .le. minC) then
                minC = C(p)
            endif
        enddo
        averageM = averageM/(col+1)
        averageC = averageC/(col+1)
        y = real(i)/real(row)*yLen
        write(12,'(f20.12,f20.12,f20.12)') y,averageM,averageC
!write(12,'(f14.10,f14.10,f14.10,f14.10,f14.10,f14.10,f14.10)') y, averageM, averageC, minM, maxM, minC, maxC
  enddo
  write(12,*) 
      
      
end subroutine printToFile2D


!==============================================================================
!   Calculates the difference between each grid point for the solutions at
!       different iterations
!-----------------------------------------------------------------------------
!   Input:
!     row    =  The number of rows
!     col    =  The number of columns
!     C      =  The previous solution for substrait
!     Cnew   =  The current solution for substrait
!   Output:
!     diff   =  The average difference between the two C's
!==============================================================================
subroutine calcDiff(diff, C, Cnew, row, col)
    integer, intent(in) :: row, col
    real, dimension(row*col), intent(in) :: C, Cnew
    real, intent(out) :: diff

    integer :: i

    diff = 0
    do i=1,row*col
        diff = diff + abs(C(i) - Cnew(i))
    enddo
    diff = diff/real(row*col)
    
end subroutine calcDiff


!==============================================================================
!   Calculates the peak and interface info at a single timestep
!-----------------------------------------------------------------------------
!   Input:
!     row    =  The number of rows
!     col    =  The number of columns
!     M      =  The solution for biomass, array size (n)
!   Output:
!     peak   =  Peak location
!     height =  Peak height
!     intfac =  Interface location
!==============================================================================
subroutine calcPeakInterface(M, row, col, peak, height, intfac, yLen)
    implicit none
    integer, intent(in):: row,col
    real, intent(in) :: yLen
    real, dimension(row*col), intent(in) :: M
    real, intent(out) :: peak, height, intfac

    real :: hei
    real :: y
    integer :: i,j,p

    hei = 0
    do, i=1,row
      y = real(i-1)/real(row-1)*yLen
      do, j=1,col
          p = (j + (i-1)*col)
          if (M(p) >= hei) then 
              peak = y
              hei = M(p)
          endif
          if (M(p) > 0.1) then
              intfac = y
          endif
      enddo
    enddo
    height = hei

end subroutine 


!==============================================================================
!   Writes a bunch of statistics to file
!-----------------------------------------------------------------------------
!   Input:
!       avgIters = average number of iterations from between solutions
!       maxIters = maximum number of iterations from between solutions
!       avgNit   = average number of iterations from linear solver
!       maxNit   = maximum number of iterations from linear solver
!       time     = time to complete solveOrder
!   Output: 
!       write everything to the file.
!==============================================================================
subroutine reportStats(avgIters,maxIters,avgNit,maxNit,time)
  implicit none
  real,intent(in)::avgIters,maxIters,avgNit,maxNit,time
  integer :: stat

  open(UNIT = 125, IOSTAT = stat, FILE = "statReport.dat", STATUS = "old")
  if (stat .EQ. 0) close(125, STATUS = "delete")
  open(UNIT = 128, FILE = "statReport.dat", POSITION = "append", ACTION = "write")
  
  write(128,*) "Statsitcs:"
  write(128,*) "-------------------------------------------------------------"
  write(128,*) "Time to compute = ", time
  write(128,*) ""
  write(128,*) "Avg Iters for iterating betn. soln. =", avgIters
  write(128,*) "Max Iters for iterating betn. soln. =", maxIters
  write(128,*) "Avg Iters for linear solver =", avgNit
  write(128,*) "Max Iters for linear solver =", maxNit  
  
  close(128)
  
end subroutine



subroutine amuxd (n,x,y,diag,idiag,ioff) 
!-----------------------------------------------------------------------
!        Mnew times a vector in Diagonal storage format (DIA) 
!        f90/f95 version of the sparskit f77 subroutine
!----------------------------------------------------------------------- 
! multiplies a matrix by a vector when the original matrix is stored 
! in the diagonal storage format.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n     = row dimension of Mnew
! x     = real array of length equal to the column dimension of
!         the Mnew matrix.
! ndiag  = integer. The first dimension of array adiag as declared in
!         the calling program.
!         (obscolete (=n always)
! idiag  = integer. The number of diagonals in the matrix.
! diag   = real array containing the diagonals stored of Mnew.
! idiag  = number of diagonals in matrix.
! diag   = real array of size (ndiag x idiag) containing the diagonals
!          
! ioff   = integer array of length idiag, containing the offsets of the
!   	   diagonals of the matrix:
!          diag(i,k) contains the element a(i,i+ioff(k)) of the matrix.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=Mnew*x
!
!-----------------------------------------------------------------------
implicit none
    integer, intent(in)::  n, idiag
    integer, intent(in),dimension(idiag) :: ioff
    real, dimension(n), intent(in) :: x    
    real, dimension(n,idiag), intent(in) :: diag
    real, dimension(n), intent(out) :: y
    integer :: j, io, i1, i2, i       

    !$omp parallel shared(y,diag,x,n) private(j,io,i1,i2)

    !!$omp workshare
    !$omp do
    do i=1,n
        y(i)=0.  
    enddo
    !$omp enddo
    !!$omp end workshare    

    do j=1, idiag
        io = ioff(j)
        i1 = max0(1,1-io)
        i2 = min0(n,n-io)
        !$omp do
        do i=i1,i2
            y(i) = y(i)+diag(i,j)*x(i+io)
        enddo
        !$omp end do
    enddo  
    !$omp end parallel
    
end subroutine amuxd





