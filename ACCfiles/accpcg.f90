subroutine solve_linear_system(n,idiag,ioff,A,rhs,sol)
!-----------------------------------------------------
! solves      A*sol = rhs
! using PCG, follows Templates
!-------------------------------------
implicit none
  integer, intent(in) :: n,idiag 
  real, dimension(n,idiag), intent(in) :: A
  integer, dimension(idiag), intent(in) :: ioff
  real, dimension(n), intent(in) :: rhs
  real, dimension(n), intent(inout) :: sol
  real,dimension(n):: p,r,z,dsol,q,hvv,rm,diff
  real :: rho,rhoold,beta,alfa,tol1,tol2,tol3,err1,err2,err3,err4
  real :: dd2,hvr,diffmin,diffmax,diffavg
  integer :: abbr,it,i,maxit
  real :: accdistdot
  external accdistdot

!  write(*,*) 'start'

  beta=0.; rhoold=0.

  tol1=1e-12
  tol2=1e-6
  tol3=1e-12
  maxit=100*n
  rho=1.
  alfa=1.
  abbr=0
  err1=2.*tol1
  it=0

  
  !$acc data copyin(rhs,A) copy(sol) create(dsol,r,p,q,rm,z,hvv,diff) 
  call accamuxd(n,sol,r,A,idiag,ioff)


  !$acc kernels present(r,p,rm,a,rhs)

  err2=0.
  !$acc loop reduction(max:err2)
  do i=1,n
    r(i)=rhs(i)-r(i)
    p(i)=0.
    rm(i)=1./a(i,1)
    err2=max(err2,abs(r(i)))
  enddo
  !$acc end kernels

  if (err2<tol2) then
    abbr=-100
  else
    abbr=0
  endif
  !write(*,*) 'while loop starts now'
  do while(abbr>=0)
    it=it+1
    rhoold=rho

    !solve Mz=r
    !$acc parallel loop present(z,r,rm)
    do i=1,n
     z(i)=r(i)*rm(i)
    enddo


    rho=accdistdot(n,z,1,r,1)

    if (it==1)then
     !$acc parallel loop present(z,p) 
     do i=1,n
       p(i)=z(i)
     enddo
    else
      beta=(rho/rhoold)
      !$acc parallel loop present(z,p) firstprivate(beta)
      do i=1,n
        p(i)=z(i)+beta*p(i)
      enddo
     endif
     call accamuxd(n,p,q,A,idiag,ioff)

    dd2=accdistdot(n,p,1,q,1)

    alfa=rho/dd2
   
    err1=0.; err2=0.; err3=0.
    !$acc parallel loop present(dsol,sol,p,r,q) firstprivate(alfa) reduction(max:err1,err2) reduction(+:err3)
    do i=1,n
      dsol(i)=alfa*p(i)
      sol(i)=sol(i)+dsol(i)
      r(i)=r(i)-alfa*q(i)
      err1=max(err1,abs(dsol(i)))
      err2=max(err2,abs(r(i)))
      err3=err3+dsol(i)*dsol(i)
    enddo
    err3=sqrt(err3)/n


    if (err1<tol1) abbr=abbr-10
    if (err2<tol2) abbr=abbr-100
    if (err3<tol3) abbr=abbr-1000
    if (it==maxit) abbr=abbr-1
!     write(*,*) 'dsol',maxval(dsol),minval(dsol)
!     write(*,*) 'p   ',maxval(p),minval(p)
!     write(*,*) 's   ',maxval(s),minval(s)
!     write(*,*) 'r   ',maxval(r),minval(r)
!     write(*,*) 'v   ',maxval(v),minval(v)
!     write(*,*) 't   ',maxval(t),minval(t)
!     write(*,*) alfa,omega,beta


!   call accamuxd(n,sol,hvv,A,idiag,ioff)
!
!    err4=0.
!    !$acc parallel loop present(RHS,hvv) reduction(max:err4) 
!    do i=1,n
!      err4=max(err4,abs(hvv(i)-rhs(i)))
!    enddo
!   write(*,'(I5,4(X,E14.7))') it,err1,err2,err3,maxval(sol)

  enddo

!  call accamuxd(n,sol,hvv,A,idiag,ioff)
!
!  diffmin=1e12; diffmax=-1e12; diffavg=0.
!  !$acc parallel loop present(hvv,rhs,diff) reduction(min:diffmin) reduction(max:diffmax) reduction(+:diffavg)
!  do i=1,n
!    diff(i)=rhs(i)-hvv(i)
!    diffmin=min(diffmin,diff(i))
!    diffmax=max(diffmax,diff(i))
!    diffavg=diffavg+abs(diff(i))
!  enddo  


  !$acc end data


!! write(*,'(I8,X,E14.7,1X,I8,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)') it,err1,abbr,err2,err3,maxval(sol),minval(sol)

! ! write(*,*) '+---------------------------------------------'
!  diffavg=diffavg/n
! ! write(*,'(I5,X,6(E14.7,X))') it,diffmin,diffavg,diffmax,sum(abs(rhs))/n,minval(abs(rhs))

end subroutine solve_linear_system



