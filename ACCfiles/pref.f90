subroutine prefix(n,pre)
!--------------------------------------------
! defines the prefix of the output file names
!--------------------------------------------
implicit none
  integer, intent(in) :: n
  integer:: hv
  character(len=6),intent(out)  :: pre
  character :: s1,s2,s3
  
  pre='AM'
  if (n<0) then
    write(*,*) 'negative counter'
    stop
  else  
    if (n>999) then
      write(*,*) 'warning:  n>999 --- set suf=XXX'
      pre='XXX'  
    else
      hv=modulo(n,10) 
      s1=achar(48+hv)
      hv=(n-hv)/10
      hv=modulo(hv,10)
      s2=achar(48+hv)
      hv=(n-modulo(n,100))/100 
      s3=achar(48+hv)
      pre=trim(pre)//s3//s2//s1
    endif
  endif

  pre=trim(pre)//'.'
  pre=trim(pre)
  !write(*,*) n,pre,s3,s2,s1
end subroutine prefix


