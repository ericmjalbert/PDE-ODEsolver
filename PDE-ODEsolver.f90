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
    real :: depth,height
    real :: kappa,delta,nu,gama
    integer :: alpha,beta
    integer :: fSelect, dSelect, gSelect, MinitialCond

    ! Numerical Method Parameters
    integer :: pSize,row,col,n,ndiag,nit,nOuts
    real :: tEnd,tDel,xDel,e1,e2,eSoln
    integer :: true2D

    ! Solution variables
    real,dimension(:),allocatable  :: C, M
    
    ! Reporting variables
    real :: avgIters, maxIters
    real :: avgNit, maxNit
    integer :: startTime, endTime, clock_rate, clock_max
    real :: elapsedTime
    
!    write(*,*) "Enter parameter file name: ex. 'parameter.txt' "
!    write(*,'(A)', advance="no") "    "
    read(*,*) filename

!   write(*,*) "Setting Parameters that shouldn't change"
    ndiag = 5
!    write(*,'(A, I5)') "     ndiag = ", ndiag

!    write(*,*) "Getting problem size from file"
    call getProbSize(pSize, true2D, filename, len(filename)) 
    if (true2D == 1) then 
!      write(*,*) "    Problem is 2D"
      col = 4
    else if (true2D == 0) then 
!      write(*,*) "    Problem is 3D"
      col = pSize
    end if
!    write(*,*) "    pSize = ", pSize
    row = pSize
    n = pSize * col
!    write(*,*) "    row   = ", row
!    write(*,*) "    col   = ", col        
    write(*,*) "    n     = ", n

!    write(*,*) "Opening Parameter file"
    call paramSet(depth, height,alpha, beta, kappa, &
                  gama, nu, delta, nOuts, tEnd, tDel, e1, e2, & 
                  eSoln, fSelect, dSelect, gSelect, MinitialCond, filename, &
                  len(filename))
    xDel = 1/real(row)
    nit = n * 100
!    write(*,*) "Parameters set:"
!    write(*,*) "    depth       = ", depth
!    write(*,*) "    height      = ", height
!    write(*,*) "    alpha       = ", alpha
!    write(*,*) "    beta        = ", beta
    write(*,*) "    kappa       = ", kappa
    write(*,*) "    gama        = ", gama
    write(*,*) "    nu          = ", nu
    write(*,*) "    delta       = ", delta
    write(*,*) "    nOuts       = ", nOuts
    write(*,*) "    tEnd        = ", tEnd
    write(*,*) "    tDel        = ", tDel
    write(*,*) "    e1          = ", e1
    write(*,*) "    e2          = ", e2
!    write(*,*) "    nit         = ", nit
!    write(*,*) "    xDel        = ", xDel
!    write(*,*) "    fSelect     = ", fSelect
!    write(*,*) "    dSelect     = ", dSelect
!    write(*,*) "    gSelect     = ", gSelect
!    write(*,*) "    MinitialCond= ", MinitialCond
!   write(*,*) "Allocating the size of arrays"
    allocate(C(n),M(n))
!    write(*,'(A,I12,A)') "    C and M are now dimension(", n, ") arrays"

    call setInitialConditions(C,M,row,col,n,depth,height,xDel,MinitialCond)
        
!    write(*,*) "Starting Solver"
    call system_clock(COUNT_RATE=clock_rate, COUNT_MAX=clock_max)
    call system_clock(startTime)
    call solveOrder2(tEnd,nOuts,tDel,n,row,col,M,C,xDel,&
       gama,nu,alpha,beta,kappa,delta,ndiag,e1,e2,nit,eSoln,dSelect,&
       fSelect,gSelect,avgIters,maxIters,avgNit,maxNit,true2D)
    call system_clock(endTime)
   
    ! Convert times into seconds
    elapsedTime = real(endTime - startTime)/real(clock_rate)


    ! Report Statistics
    write(*,*) "-------------------------------------------------------------"
    write(*,*) "Statsitcs:"
    write(*,*) "-------------------------------------------------------------"
    write(*,*) "Time to compute = ", elapsedTime
    write(*,*) ""
    write(*,*) "Avg Iters for iterating betn. soln. =", avgIters
    write(*,*) "Max Iters for iterating betn. soln. =", maxIters
    write(*,*) "Avg Iters for linear solver =", avgNit
    write(*,*) "Max Iters for linear solver =", maxNit  
    write(*,*) "Writing statistics to file"
    call reportStats(avgIters,maxIters,avgNit,maxNit,elapsedTime)
    
    write(*,*) "Execution completed"
  
end program


!==============================================================================
!   Sets the problem size and number of diagonal elements
!-----------------------------------------------------------------------------
!     This must be done first since the problem size is used in the variable 
!   declaration of the arrays.
!==============================================================================
subroutine getProbSize(pSize, true2D, filename, nameLen) 
    implicit none
    integer,intent(out) :: pSize, true2D
    integer,intent(in)  :: nameLen
    character(nameLen),intent(in) :: filename
    character :: dum   ! Dummy variable
    write(*,*) filename
    open(UNIT = 19, FILE = filename, STATUS = "old", ACTION = "read")
    read(19,*); read(19,*); read(19,*)  ! Skip first 3 lines
    read(19,*) dum, dum, pSize
    read(19,*) dum, dum, true2D
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
subroutine paramSet(depth, height, alpha, beta,  &
                    kappa, gama, nu, delta, nOuts, tEnd, tDel, e1, e2,&
                    eSoln, fSelect, dSelect, gSelect, MinitialCond, filename, &
                    nameLen)
    implicit none
    real, intent(out) :: depth, height, kappa, delta, nu 
    real, intent(out) :: tEnd, tDel, e1, e2, eSoln, gama
    integer, intent(out) :: alpha, beta, nOuts
    integer, intent(out) :: fSelect, dSelect, gSelect, MinitialCond
    integer, intent(in) :: nameLen
    character(nameLen), intent(in) :: filename

    character :: dum, dum2      ! Dummy variable
    
    open(UNIT = 19, FILE = filename, STATUS = "old", ACTION = "read")
    ! Skip first 5 lines
    read(19,*); read(19,*); read(19,*); read(19,*); read(19,*)  

    read(19,*) dum, dum2, depth
    read(19,*) dum, dum2, height
    read(19,*) dum, dum2, alpha
    read(19,*) dum, dum2, beta
    read(19,*) dum, dum2, nu
    read(19,*) dum, dum2, kappa 
    read(19,*) dum, dum2, gama 
    read(19,*) dum, dum2, delta
    read(19,*) dum, dum2, nOuts
    read(19,*) dum, dum2, tEnd
    read(19,*) dum, dum2, tDel
    read(19,*)                          ! Skip xDel 
    read(19,*) dum, dum2, e1
    read(19,*) dum, dum2, e2
    read(19,*) dum, dum2, eSoln
    read(19,*)                          ! Skip nit
    read(19,*); read(19,*); read(19,*)  ! Skip 3 lines
    read(19,*) dum, dum2, fSelect
    read(19,*) dum, dum2, dSelect
    read(19,*) dum, dum2, gSelect
    read(19,*) dum, dum2, MinitialCond

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
subroutine setInitialConditions(C,M,row,col,n,depth,height,xDel,MinitialCond)
    implicit none
    integer,intent(in) :: row,col,n,MinitialCond
    real,intent(in) :: depth, height, xDel
    real,dimension(n),intent(out) :: C,M

    integer :: i,j,jstop,x,y
    real :: f,a       ! function for IC curve
    
    C = 1.; j = 0

    if (MinitialCond == 1) then
        M = 0
        a = -height/(depth)**4
        f = a*(j*xDel)**4 + height
        jstop = INT(depth/xDel)+1
        !$omp parallel do private(f) shared(height,a) 
        do j = 0, jstop
            do i = 1,col
                M(i+j*col) = f
            enddo
            f = a*(j*xDel)**4 + height
        enddo
        !$omp end parallel do
    else if (MinitialCond == 2) then
        M = height

    ! 3D initial conditions
    else if (MinitialCond == 3) then
        M = 0
        a = -height/(depth)**2
        f = a*(j*xDel)**2 + height
        jstop = INT(depth/xDel)+1
        !$omp parallel do private(f,x,y) shared(height,a) 
        do i = 1, n
          x = MOD(i, col)
          y = i / row
          f = a*((x*xDel-0.5)*(x*xDel-0.5) + (y*xDel-0.5)*(y*xDel-0.5)) + height
          M(i) = f
          if (M(i) .LE. 0) M(i) = 0
        enddo
        !$omp end parallel do

    endif
end subroutine setInitialConditions


!==============================================================================
!   Function f(C,M)
!==============================================================================
subroutine fFunc(M,C,f,kappa,nu,fSelect)
    implicit none    
    integer, intent(in) :: fSelect
    real, intent(in) :: M,C,kappa,nu
    real, intent(out) :: f
    real :: eps
    eps = 0.0000000001
    if (fSelect == 1) f = C/(kappa+C)-nu 
    if (fSelect == 2) f = C/(kappa+C)
    if (fSelect == 3) f = C/(kappa*M+C+eps)-nu
    if (fSelect == 4) f = C/(kappa*M+C+eps)
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
!   gSelect = 1 --> g = gama*C/(kappa+C)
!==============================================================================
subroutine solveC(M,Mnew,Csol,Cnew,n,kappa,gama,tDel,gSelect)
    implicit none 
    integer,intent(in) :: n,gSelect
    real,intent(in) :: kappa,tDel,gama
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

    aux = tDel*0.5*gama
    
    if (gSelect == 1) then
     !$omp parallel do private(b,c) shared(kappa,Csol,aux,Mnew,M)
      do i = 1,n 
        b = kappa - Csol(i) + aux*(Mnew(i) + Csol(i)*M(i)/(kappa+Csol(i)) )
        c = -kappa*Csol(i) + aux*kappa*Csol(i)*M(i)/(kappa + Csol(i))
        Cnew(i) = 0.5 * (-b + SQRT(b*b - 4 * c))
      enddo
     !$omp end parallel do
    end if    
      
    end subroutine solveC


!====================================================================
!   Generates matrix M in diagonal format for the current timestep
!====================================================================
subroutine GenMatrixM(M,C,MatrixM,Mioff,Mrhs,row,col,n,ndiag,delta,nu,&
                      alpha,beta,kappa,tDel,xDel,dSelect,fSelect)
    implicit none
    integer,intent(in) :: row,col,n,ndiag,alpha,beta,dSelect,fSelect
    real,intent(in) :: delta,kappa,xDel,tDel,nu
    real,dimension(n),intent(in) :: M,C

    real,dimension(n,ndiag),intent(out) :: MatrixM
    real,dimension(n),intent(out) :: Mrhs
    integer,dimension(ndiag),intent(out) :: Mioff

    real :: xCof
    real :: f
    real,dimension(n) :: diff
    !integer :: x,y,g
    integer :: i

    real :: tDela
    tDela = tDel      ! Used for testing purposes

    xCof = 1/(xDel*xDel)

    Mioff = (/ -col, -1, 0, 1, col/)
    MatrixM(:,:) = 0

    ! Compute all the diffusion coefficients
    !$acc kernels
    !$omp parallel do shared(M,diff,delta,alpha,beta,dSelect)
    do i = 1,n
        call dFunc(M(i), diff(i), delta, alpha, beta, dSelect)
    enddo
    !$omp end parallel do
    !$acc end kernels

    !$acc kernels
    !$omp parallel do private(f)
    do i = 1,n
        if (i .LE. col) then
            MatrixM(i,5) = MatrixM(i,5) - xCof*0.5*(diff(i+col)+diff(i))
            MatrixM(i,3) = MatrixM(i,3) + xCof*0.5*(diff(i+col)+diff(i))
        else
            MatrixM(i,1) = MatrixM(i,1) - xCof*0.5*(diff(i-col)+diff(i))
            MatrixM(i,3) = MatrixM(i,3) + xCof*0.5*(diff(i-col)+diff(i))         
        endif
          
        if (MOD(i,col) == 1) then
            MatrixM(i,4) = MatrixM(i,4) - xCof*0.5*(diff(i+1)+diff(i))
            MatrixM(i,3) = MatrixM(i,3) + xCof*0.5*(diff(i+1)+diff(i))
        else
            MatrixM(i,2) = MatrixM(i,2) - xCof*0.5*(diff(i-1)+diff(i))
            MatrixM(i,3) = MatrixM(i,3) + xCof*0.5*(diff(i-1)+diff(i))         
        endif

        if (MOD(i,col) == 0) then
            MatrixM(i,2) = MatrixM(i,2) - xCof*0.5*(diff(i-1)+diff(i))
            MatrixM(i,3) = MatrixM(i,3) + xCof*0.5*(diff(i-1)+diff(i))
        else
            MatrixM(i,4) = MatrixM(i,4) - xCof*0.5*(diff(i+1)+diff(i))
            MatrixM(i,3) = MatrixM(i,3) + xCof*0.5*(diff(i+1)+diff(i))         
        endif
          
        if  (i .GE. n-col) then
            MatrixM(i,1) = MatrixM(i,1) - xCof*0.5*(diff(i-col)+diff(i))
            MatrixM(i,3) = MatrixM(i,3) + xCof*0.5*(diff(i-col)+diff(i))
        else
            MatrixM(i,5) = MatrixM(i,5) - xCof*0.5*(diff(i+col)+diff(i))
            MatrixM(i,3) = MatrixM(i,3) + xCof*0.5*(diff(i+col)+diff(i))         
        endif
        
        call fFunc(M(i),C(i),f,kappa,nu,fSelect)
        MatrixM(i,3) = MatrixM(i,3) - f + (1/tDel)
        Mrhs(i) = M(i)/tDel
    enddo
    !$omp end parallel do
    !$acc end kernels

end subroutine GenMatrixM


!==============================================================================
!   Solve Order 2
!-----------------------------------------------------------------------------
!   1. Solves for M_{i+1} using C_i and M_i
!   2. Solves for C_{i+1} using C_i, M_i, and M_{i+1}
!   ... repeat until convergence
!==============================================================================
subroutine solveOrder2(tEnd,nOuts,tDel,n,row,col,M,C,xDel,&
    gama,nu,alpha,beta,kappa,delta,ndiag,e1,e2,nit,eSoln,dSelect,fSelect,&
    gSelect,avgIters,maxIters,avgNit,maxNit,true2D)
  implicit none
  integer,intent(in) :: nOuts,n,row,col,alpha,beta,ndiag
  integer,intent(in) :: dSelect,fSelect,gSelect,true2D
  real,intent(in) :: tEnd,tDel,xDel,kappa,delta,nu,gama
  real,intent(in) :: eSoln
  real,intent(inout) :: e1,e2
  integer,intent(inout) :: nit
  real,dimension(n),intent(out) :: M,C
  real,intent(out) :: avgIters,maxIters,avgNit,maxNit
  
  integer :: counter        ! Counts the num. of iter. in system solver
  integer :: endLoop        ! endLoop = 1 -> solving loop can stop
  integer :: filter         ! Controls frequency of outputs written
  
  real :: stored_e1, stored_e2
 
  real :: totalMassM        ! Total Biomass over region
  real :: totalMassC        ! Total Substrait over region
  
  integer :: startTime, endTime, clock_rate, clock_max
  real :: mTime, cTime 

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
  
  stored_e1 = e1
  stored_e2 = e2 
 
  filter = int(tEnd/(nOuts*tDel)) +1
  endLoop = 0
  counter = 0
  countIters = 0
  avgIters = 1 
  avgNit = 1 
mTime = 0  
cTime = 0
!  open(UNIT = 124, IOSTAT = stat, FILE = "total.dat", STATUS = "old")
!  if (stat .EQ. 0) close(124, STATUS = "delete")
!  open(UNIT = 120, FILE = "total.dat", POSITION = "append", ACTION = "write")
!  
!  if (true2D == 1) then
!    open(UNIT = 124, IOSTAT = stat, FILE = "peakInfo.dat", STATUS = "old")
!    if (stat .EQ. 0) close(124, STATUS = "delete")
!    open(UNIT = 121, FILE = "peakInfo.dat", POSITION = "append", ACTION = "write")
!  end if    

  write(*,*) "   time    avgIter  maxIter      avgNit  maxNit        avgM        avgC"

  do while(counter * tDel <= tEnd)
    ! Output every 100 times more then nOuts for the peak info
    if (MOD(counter, int(filter/100+1)) == 0) then
      call calcMass(M,totalMassM,n,row,col)
      call calcMass(C,totalMassC,n,row,col)
     ! write(120,*) counter*tDel, totalMassM, totalMassC
     ! if (true2D == 1) then
     !   call calcPeakInterface(M, row, col, peak, height, intfac)
     !   write (121,*) tDel*counter, peak, height, intfac
     ! end if
    endif 
  
    ! Write to file / report Total Mass
    if (MOD(counter, filter) == 0) then
!      if (true2D == 1) then
!        call printToFile2D(n,row,col,M,C)
!      else
!        call printToFile(n,row,col,M,C)
!      end if

      if (counter == 0) then
        write(*,'(F8.2,F12.2,I8,F12.2,I8,F12.6,F12.6)') 0.0, 0.0, &
          int(maxIters), 0.0, int(maxNit), totalMassM, totalMassC
      else
        write(*,'(F8.2,F12.2,I8,F12.2,I8,F12.6,F12.6)') tDel*counter, &
          real(avgIters/counter), int(maxIters), real(avgNit/avgIters), &
          int(maxNit),totalMassM, totalMassC
      endif
    endif

    diffC = 1
    diffM = 1
    countIters = 0
    Cprev = C
    Mprev = M

    do while(diffC + diffM > eSoln)
        nit = 100*n
        e1 = stored_e1 
        e2 = stored_e2

  !$acc data       
         
        ! Solve M
  call system_clock(COUNT_RATE=clock_rate, COUNT_MAX=clock_max)
  call system_clock(startTime)
        call GenMatrixM(Mprev,C,MatrixM,Mioff,Mrhs,row,col,n,ndiag,delta,nu,&
                    alpha,beta,kappa,tDel,xDel,dSelect,fSelect)
        call solveLSDIAG(n,ndiag,Mioff,MatrixM,Mnew,Mrhs,nit,e1,e2,stopcritria)
  call system_clock(endTime)
  mTime = mTime + real(endTime - startTime)/real(clock_rate) 
  !$acc end data
        ! Solve C

  call system_clock(COUNT_RATE=clock_rate, COUNT_MAX=clock_max)
  call system_clock(startTime)
  call solveC(Mprev,Mnew,Cprev,Cnew,n,kappa,gama,tDel,gSelect)
  call system_clock(endTime)
  cTime = cTime + real(endTime - startTime)/real(clock_rate)        

        if(nit > maxNit) maxNit = nit
        avgNit = avgNit + nit

        call calcDiff(diffC, C, Cnew, row, col)
        call calcDiff(diffM, M, Mnew, row, col)

        C = Cnew
        M = Mnew
        countIters = countIters+1

        if (countIters > 100) then
          write(*,*) "[!] Over 100 iterations in one timestep"
          write(*,*) "[!] Solutions not converging. Exit!"
          exit
        end if
!        write(*,*) countIters, diffC, diffM, C(12),M(12)
    enddo
    if(countIters > maxIters) maxIters = countIters
    avgIters = avgIters + countIters
    
    counter = counter + 1
  enddo
  avgNit = avgNit/(avgIters) ! avgIters right now is the total
  avgIters = avgIters/counter

  ! Print final solution
!  if (true2D == 1) then
!    call printToFile2D(n,row,col,M,C)
!  else
!    call printToFile(n,row,col,M,C)
!  end if
!
!  close(120)
!  close(121)
  write(*,*) "Solve M: ", mTime
  write(*,*) "Solve C: ", cTime
  
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

    !$omp parallel do reduction(+:totalMass) 
    do i=1,n
      totalMass = totalMass + X(i)
    enddo
    !$omp end parallel do 
    
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
!   Output:
!     none
!==============================================================================
subroutine printToFile(n,row,col,M,C)
    implicit none

    integer,intent(in) :: n, row, col
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
            write(11,*) real(j-1)/real(col-1), &
                        real(i-1)/real(row-1), M(p), C(p)
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
!   Output:
!     none
!==============================================================================
subroutine printToFile2D(n,row,col,M,C)
    implicit none

    integer,intent(in) :: n, row, col
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
        y = real(i)/real(row)
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
    !$omp parallel do reduction(+:diff) 
    do i=1,row*col
        diff = diff + abs(C(i) - Cnew(i))
    enddo
    !$omp end parallel do
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
subroutine calcPeakInterface(M, row, col, peak, height, intfac)
    implicit none
    integer, intent(in):: row,col
    real, dimension(row*col), intent(in) :: M
    real, intent(out) :: peak, height, intfac

    real :: hei
    real :: y
    integer :: i,j,p

    hei = 0
    do, i=1,row
      y = real(i-1)/real(row-1)
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
  real,intent(in)::avgIters,maxIters,avgNit,maxNit
  real,intent(in)::time
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



