program qs_biofilms
!----------------------------------------------------------------------------
implicit none
!----------------------------------------------------------------------------
! parameters
!----------------------------------------------------------------------------
  real ::  eps,dt,dx,t,tout,to,tmax !,L
  integer :: n,alfa,beta, i,j,it,ntstep,dtok,m,nm,p,oi,jspec,nspec, simit,nr,nsimit
  integer, dimension(5) :: ioff
  integer, dimension (2:5) :: pp
  parameter(nspec=4,eps=1e-3, beta=4,alfa=4, n=200,m=n/2, ntstep=20000000)
  parameter(nm=n*m, ioff=(/0,1,-1,m,-m /))
  parameter(nr=99)
  real, dimension(nr) ::pr
  real, dimension(1:n,1:m,nspec) :: u,R1,R2
  real, dimension(n*m,5)     :: a
  real, dimension (n,m,nspec)    :: d
  real, dimension(1:n*m)     :: dconc,af
  real, dimension(1:n*m,nspec)   :: uu
  real, dimension(1:n,1:m)   :: tf
  real :: occupancy,maxoccu,u1tot,u2tot,bulkconc,tstart,grmax,u4tot
  integer :: nfilm,maxheight,height,abbr,outsw1,outsw2
  character(len=3) :: suf
  character(len=6) :: filename
  parameter(tmax=14, tout=50, maxoccu=0.75,maxheight=1000,nsimit=1)

  call random_seed()

! nsimit simulations are carried out in one run
do  simit=1,nsimit
    write(*,*) "+++++++++++++++++++++++++++++++++++++++++++++++++++++"
    write(*,*) " new simulation ",simit
    write(*,*) "+++++++++++++++++++++++++++++++++++++++++++++++++++++"

    ! some intialisation
    u=0.; a=0.0; af=0.; d=0.; uu=0;  oi=0; occupancy=0.; abbr=0
    call prefix(simit,filename)
    call parameters(pr,nr,n,m,simit,nsimit)
    call initial_data_II(n,m,2,pr(14),pr(15),u,simit)
    dx=pr(1)/(1.*n)

    open(13,file=filename//'tme', action='write')
    open(15,file=filename//'bf',action='write')
    do i=1,nr
       if (pr(i)>10e-20) write(13,*) '#', i,pr(i)
    enddo 
    write(*,'(I6,1X,E9.4,1X,E9.4,1X,E9.4,1X,E9.4,1X,E9.4,1X,E9.4,1X,E9.4,1X,E9.4)')  &
        0, 0.,0.,maxval(u(:,:,1)),maxval(u(:,:,2)),minval(u(:,:,3)),maxval(u(:,:,4)),occupancy, &
        maxval(u(:,:,1)+u(:,:,2))

    write(13,'(I6,1X,E9.4,1X,E9.4,1X,E9.4,1X,E9.4,1X,E9.4,1X,E9.4)') &
        0, 0.,0.,maxval(u(:,:,1)),maxval(u(:,:,2)),minval(u(:,:,3)),maxval(u(:,:,4))

     ! time step initialisation
   it=0;  t=0.; to=tout; occupancy=0.; tstart=-1.; height=0
   !-----------------------------------------------------------------
   ! time loop
   write(*,*) 'start time loop'
   do while (abbr==0)
         
      if (maxval(u(:,:,4))<1.) then 
         outsw1=0
      else
         outsw1=1
      endif
         
      if (u4tot<1.) then
         outsw2=0
      else 
         outsw2=1
      endif

      if (tstart<0. .and. occupancy>=pr(20)) tstart=t
      call volfrac(n,m,nspec,u,tf)
      if (maxval(tf)>1) then
          write(*,*) "biomass too large"
          stop
      endif
      it=it+1
      ! determine local diffusion coefficients

      do jspec=1,nspec,2 
         if (jspec==1 .or. jspec==2) then
           ! biomass fractions
           d(1:n,1:m,jspec)= pr(60)*(tf)**beta   * (1.-tf)**(-alfa)
        else
        ! dissolved substrates (linear interpolation)
          d(1:n,1:m,jspec)=pr(13+jspec)*(1.- (1.-pr(15+jspec))*tf)
        endif 
      enddo ! loop over-species
      d=d/(dx*dx)

      ! determine reaction rates
      call reactions(n,m,nspec,u,nr,pr,R1,R2)

       ! estimate time-step
      dt=100./maxval(d(:,:,1)+d(:,:,2))
      grmax=maxval(u(:,:,3), tf>1e-6)
      grmax=1./(pr(4)*grmax/(grmax+pr(3)))
      dt=min(0.99*grmax,dt)
      dt=min(dt,0.002)
      dtok=1.

      do while (dtok>0.)  ! carry out time step
        t=t+dt
        ! loop over all species and substrates
        do jspec=1,nspec
         !matrix A setup
             A=0.0; af=0.
             do i=1,n
             do j=1,m
               dconc((i-1)*m+j)=d(i,j,jspec) 
             enddo
             enddo

             do i=1,n
             do j=1,m
               p= (i-1)*m + j
               pp(2:5)= p + ioff(2:5)
               IF (i==1) pp(5)  = p
               IF (j==1) pp(3)  = p
               IF (i==n) pp(4) = p
               IF (j==m) pp(2) = p

               !diffusion 
               a(p,1)   = 1.+0.5*dt*(4.*dconc(p)+SUM(dconc(pp(2:5))) )
               a(p,2:5) = -0.5*dt*(dconc(pp(2:5))+dconc(p))
          
               ! reaction
               a(p,1)=a(p,1)-dt*R1(i,j,jspec)
               af(p)=u(i,j,jspec)+dt*R2(i,j,jspec)
              enddo
              enddo 
      
              ! boundary condition corrections
              select case(jspec)

                ! in/out
                case(1:4) 
                   DO j=1,m
                     p=(n-1)*m+j
                     a(j,1)= a(j,1)+a(j,5)              
                     a(j,5) = 0.
                     a(p,1)=a(p,1)+a(p,4)
                     a(p,4) = 0.
                   ENDDO
                !case(3) !dirichlet
                !     Do j=1,m
                !      bulkconc=pr(11+jspec)
                !      p=(n-1)*m+j
                !      af(p)=af(p)+dt*dconc(p)*bulkconc
                !      af(j)=af(j)+dt*dconc(j)*bulkconc
                !      a(j,5)=0. 
                !      a(p,4)=0.
                !     enddo
               end select
                 

                !top bottom 
                select case (jspec) 
                  case(1:2)                   
                   Do i=1,n
                     p=(i-1)*m+1
                     a(p,1)=a(p,1)+a(p,3)
                     a(p,3)=0.

                     p=i*m
                     a(p,1)=a(p,1)+a(p,3)
                     a(p,2)=0.            
                    enddo

                   case(3:4) !dirichlet bottom, neumann top
                     Do i=1,n
                      bulkconc=pr(11+jspec)
                      p=(i-1)*m+1
                      a(p,1)=a(p,1)+a(p,3)
                      a(p,3)=0.
 
                      p=i*m
                      af(p)=af(p)+dt*dconc(p)*bulkconc
                      a(p,2)=0.                     
                     enddo



                   case(-1) !dirichlet , neumann bottom
                     Do i=1,n
                      bulkconc=pr(11+jspec)
                      p=(i-1)*m+1
                      af(p)=af(p)+dt*dconc(p)*bulkconc
                      a(p,3)=0.
 
                      p=i*m
                      a(p,1)=a(p,1)+a(p,3)
                      a(p,2)=0.                        
                     enddo

                 endselect
 

      

              ! solve system

              call solve_linear_system(n*m,5,ioff,A,af,uu(:,jspec))  ! solve A * uu = af
              uu(:,jspec)=abs(uu(:,jspec))

          enddo !loop over all species and substrate

      

         if (maxval(uu(:,1)+uu(:,2))>=1) then
             t=t-dt
             dt=0.5*dt
             write(*,*) 'new time step:',dtok,dt         
             if (dt<1e-13) then
               write(*,*) "time step too small" !"u>0.9999"
               stop
             endif  
          else ! evaluation and output
             dtok=-1.
             ! transfer back to naive variables
             do jspec=1,nspec
                where (abs(uu(:,jspec))<1e-30) uu(:,jspec)=0.
                do i=1,m
                do j=1,n
                   u(j,i,jspec)=uu((j-1)*m+i,jspec)
                enddo
                enddo
             enddo


         

             ! occupancy: volume fraction occupied by biofilm

             nfilm=0;  u1tot=0; u2tot=0.; u4tot=0.
             do i=1,n
             do j=1,m  
               if (u(i,j,1)+u(i,j,2)>1e-6) nfilm=nfilm+1
               if (u(i,j,1)+u(i,j,2)>1e-6 .and. j>height) height=j
               u1tot=u1tot+u(i,j,1)
               u2tot=u2tot+u(i,j,2) 
               u4tot=u4tot+u(i,j,4)   
             enddo
             enddo
             occupancy=1.*nfilm/(1.*n*m)
             u1tot=u1tot/nm
             u2tot=u2tot/nm
             u4tot=u4tot/nm
             
             ! output
             write(13,'(I6,1X,E9.4,1X,E9.4,1X,E9.4,1X,E9.4,1X,E9.4,1X,E9.4)') &
               it, t,dt,maxval(u(:,:,1)),maxval(u(:,:,2)),minval(u(:,:,3)),maxval(u(:,:,4))
             write(15,'(E9.4,1X,E9.4,1X,E9.4,1X,E9.4,1X,E9.4)')  t,occupancy,u1tot,u2tot,u4tot
             write(*,'(I6,1X,E9.4,1X,E9.4,1X,E9.4,1X,E9.4,1X,E9.4,1X,E9.4,1X,E9.4,1X,E9.4)') &
               it, t,dt,maxval(u(:,:,1)),maxval(u(:,:,2)),minval(u(:,:,3)),maxval(u(:,:,4)) &
               ,occupancy,maxval(u(:,:,1)+u(:,:,2))
             endif

             !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             ! output: changed to minval(u(:,:,:,3)) since substrate c is consumed
             !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       enddo !dtok loop

   

       if (t>tmax)              abbr=abbr-1
       if (occupancy>maxoccu)   abbr=abbr-10
       if (height>maxheight)    abbr=abbr-100
       if (minval(u(:,:,4))>20.)  abbr=abbr-1000
       if (abs(u1tot)<(1e-3)*abs(u2tot))  abbr=abbr-10000   

       if (abbr<0) write(*,*) 'STOP  ',abbr

       if (maxval(u(:,:,4))>=1. .and. outsw1==0) outsw1=2
       
       if (u4tot >=1. .and. outsw2==0) outsw2=2

       !  write results into file
       to=to+dt
       if (to>tout .or. abbr<0 .or. outsw1==2 .or. outsw2==2) then 
         oi=oi+1
         call suffix(oi,suf)
         open(11,file=filename//suf  ,action='write')
         write(11,*) ' # t=',t
         do j=1,n
         do i=1,m
           write(11,'(4(E10.4,1X))') u(j,i,1),u(j,i,2),u(j,i,3),u(j,i,4)!,uf(j,i),vf(j,i),pf(j)
         enddo 
         write(11,*) !<- include for gnuplot visualisation
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! note: I non-dimensionalised the substrate on output for better visulisation
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         enddo
         to=0.
         close(11)
       endif

     
    enddo !time-loop
    close(13)
enddo !simit

end program



subroutine suffix(n,suf)
!--------------------------------------------
! defines the suffix of the output file names
!--------------------------------------------
implicit none
  integer, intent(in) :: n
  integer:: hv
  character(len=3),intent(out)  :: suf
  character :: s1,s2,s3
  
  suf='000'

  if (n<0) then
    write(*,*) 'negative counter'
    stop
  else  
    if (n>999) then
      write(*,*) 'warning:  n>999 --- set suf=XXX'
      suf='XXX'  
    else
      hv=modulo(n,10) 
      s1=achar(48+hv)
      hv=(n-hv)/10
      hv=modulo(hv,10)
      s2=achar(48+hv)

      hv=(n-modulo(n,100))/100 
      s3=achar(48+hv)

      suf=s3//S2//S1
    endif
  endif

end subroutine suffix




subroutine parameters(pr,nr,n,m,simit,nsimit)
!------------------------------------------------
! this defines a number of reaction parameters,
! such as reaction rates, diffusion constants etc
!------------------------------------------------
implicit none
  integer, intent(in) :: nr,n,m,simit,nsimit
  real, dimension(nr), intent(inout) :: pr
  real :: rn 

pr=0.

!--------------------------------------------------
! biofilm growth parameters from BM1 / ES, JTB:2008
! limiting dissolved substrate:  carbon
!
!--------------------------------------------------
pr(1)=2. *2.*.5e-3    !  L=0.0005  channel length [m]

pr(2)=1.0 	  !kp1
pr(12)=10000. !max cell density, M_infty
pr(3)=4.0	  !kp2   monod half saturation contant [g/m3]
pr(4)=6.0	  !kp3   max growth rate [1/d]
pr(5)=0.4	  !kp4   decay rate [1/5]
pr(6)=1.0	  !kp5   QS regulation rate
pr(7)=0.!(2.3e-10)*24.*1e12/pr(12)   ! alfa
pr(8)=0.!(2.3e-9)*24.*1e12/pr(12)	  !beta
pr(9)=10.	  !tau
pr(10)=2.5	  !n
pr(11)=0.005545*24.	  !gamma
pr(13)=0.63	  !yield coefficient 

! bc + diff
pr(14:15)=(/ 20.0, 0.0 /)     !  cbulk, tbulk
pr(16:17)=(/ 1e-4, 7.7605008e-5 /)  ! diffusion coefficients c, AHL
pr(18:19)=(/ 0.8, 0.5/)        ! diff coeff ratios water:biofilm

! flow
pr(25)= 24.*3.6e-3  ! nu
pr(26)= 1000.*1000.   ! rho
pr(27)= 5e-3  ! Re

! hydrodynamics
rn = pr(1)*m/n                            ! H=L*M/n
pr(40)= pr(25)*pr(27)/rn                  ! Uinfty
pr(41)=2./3.*pr(40)*rn                    ! Q
pr(42)= -8.*pr(40)*pr(25)/(rn*rn)*pr(26)  ! dpdxinfty

!auxiliaries
pr(50) = 1.                         ! new kappa3   - growth rate
pr(51) = pr(3)/pr(14)               ! new kappa2   - Monod
pr(52) = pr(5)/pr(4)                ! new kappa4   - decay rate
pr(53) = pr(6)*pr(9)**pr(10)/pr(4)  ! new kappa5   - upregulation rate
pr(54) = pr(12)/(pr(14)*pr(13))     ! new kappa1   - consumtion rate 
pr(55) = pr(11)/pr(4)               ! new gamma    - AHL decay rate
pr(56) = pr(7)*pr(12)/(pr(4)*pr(9)) ! new alfa     - AHL prod rate
pr(57) = pr(8)*pr(12)/(pr(4)*pr(9)) ! new beta     - AHL prod rate
pr(58) = 0.!=1.5  ! concentration boundary layer thichness, relative to system length

pr(80)=0.8*pr(50)
pr(81)=0.5*pr(51)
pr(82)=pr(52)
pr(84)=pr(54)

pr(16:17)=pr(16:17)/(pr(1)*pr(1)*pr(4))
pr(60)= (1e-12)/(pr(1)*pr(1)*pr(4))   ! new epsilon!!!!!
pr(1)=1.
pr(14)=1.
end subroutine parameters


subroutine volfrac(n,m,nspec,u,tf)
!-------------------------------------------
!  computes the volume fraction occupied in
!  every grid cell
!-------------------------------------------
implicit none
 integer, intent(in) :: n,m,nspec
 real, dimension(n,m,nspec), intent(in) ::u
 real, dimension(n,m), intent(out) :: tf

 tf=u(:,:,1)+u(:,:,2)
end subroutine


subroutine initial_data_II(n,m,sw,cbulk,tbulk,u,simit)
!-------------------------------------------------------------------------------------
!   inoculation sites are randomly chosen
!   biomass density in these sites is also chosen randomly
!-------------------------------------------------------------------------------------
implicit none
  integer, intent(in) :: n,m,sw,simit
  real, dimension(n,m,4), intent(out) :: u
  real, intent(in) :: cbulk, tbulk
  integer :: i,j,k,r
  real :: rn,r1,r2,rr,dx,posx,posy,u1,u2
  
  write(*,*) 'initial data entry'


select case(sw)
case(0)
  u=0.; k=0  
  dx=1./n
  r1=12.

  do k=1,3
  select case(k)
   case(1)  
     posx=1
     posy=1
     r1=0.05*n
     u1=0.35
     u2=0.
   case(2)
     posx=n/2.
     posy=1
     r1=0.02*n
     u1=0.25
     !u2=0.175
   case(3)
     posx=n
     posy=1.
     r1=0.05*n
     u1=0.15
     !u2=0.35
  end select

  do i=1,n
  do j=1,m
    rr=sqrt( ((i-0.5)-posx)**2 + ((j-0.5)-posy)**2 )
    if (rr<=r1) then
       u(i,j,1)=u1
!       u(i,j,2)=u2
    endif
  enddo
  enddo
  enddo


  u(:,:,3)=cbulk
!  u(:,:,4)=tbulk
  

case(1) ! read from file
 open(39,file='init.dat',action='read')
 read(39,*)
        do j=1,n
         do i=1,m
           read(39,'(4(E10.4,1X))') u(j,i,1),u(j,i,2),u(j,i,3),u(j,i,4)
         enddo 
          read(39,*) !<- include for gnuplot visualisation
         enddo

case(2)
  u=0.
  do i=1,n
     posx=i*1./n
     posy= m*(0.5*(1.-cos(posx*2.*6.2830)))**2 * posx *0.25
     j=1
     do while(j<posy)
       u(i,j,1)=0.3
       j=j+1
     enddo
  enddo
  u(:,:,2)=0.
  u(:,:,3)=cbulk
  u(:,:,4)=tbulk
endselect  

end subroutine initial_data_II



subroutine reactions(n,m,nspec,u,nr,pr,R1,R2)
!--------------------------------------------
! reaction terms 
!--------------------------------------------
implicit none 
   integer, intent(in) :: n,m,nspec,nr
   real, dimension(nr), intent(in) :: pr
   real, dimension(n,m,nspec), intent(in) :: u
   real, dimension(n,m,nspec), intent(out) :: R1, R2
   integer :: i,j

   R1=0.
   R2=0.   

   ! M_0: u(i,j,1)
   ! M_1: u(i,j,2)
   ! C: u(i,j,3)
   ! AHL: u(i,j,4)	
   !$omp parallel do	
   do i=1,n
   do j=1,m
     if (u(i,j,1)>1e-8) then
       R1(i,j,1)=pr(50)*u(i,j,3)/(pr(51)+u(i,j,3)) - pr(52) !- pr(53)*u(i,j,4)**pr(10)
       !R2(i,j,1)=pr(53)*u(i,j,2)
     endif
!     if (u(i,j,2)>1e-8) then
!        R1(i,j,2)=pr(80)*u(i,j,3)/(pr(81)+u(i,j,3))-pr(82)  !-pr(53)
!     endif
     !if (u(i,j,1)>1e-8) then 
     !   R2(i,j,2)=pr(53)*u(i,j,1)*u(i,j,4)**pr(10)
     !endif
     R1(i,j,3)=-pr(54)*(u(i,j,1))/(pr(51) + u(i,j,3))  !-pr(84)*(u(i,j,2))/(pr(81) + u(i,j,3))

     !R1(i,j,4)=-pr(55)
     !if (u(i,j,1)>1e-8) R2(i,j,4)=pr(56)*u(i,j,1)
     !if (u(i,j,2)>1e-8) R2(i,j,4)=R2(i,j,4)+pr(57)*u(i,j,2)
   enddo
   enddo 
end subroutine
