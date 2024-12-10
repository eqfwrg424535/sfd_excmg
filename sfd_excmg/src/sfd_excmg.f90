!  sfd_excmg.f90 
!
!  FUNCTIONS:
!  sfd_excmg - Entry point of console application.

!****************************************************************************
!
!  PROGRAM: sfd_excmg
!
!  PURPOSE: Large-scale MT forward modeling
!! version 1.0
! author: Jinxuan Wang (jinxuanwang@cug.edu.cn)
! realease date: 08/12/24
! Notes:    
    ! the original code calculates the MT responses frequency by frequency
    ! here the ex array contains the solutions of all frequencies, which means
    ! the post process is called only once
!****************************************************************************

    program sfd_excmg
    ! sfd_excmg
    
     use global_parameter
     use file_operation
     use model_selection
     use fwd3d
     
    implicit none

    ! Variables
    
    real*8,allocatable :: a(:),b(:),c(:),sigma(:,:),freq(:),recv(:,:),topo(:)
    complex(8),allocatable :: ex(:,:)
    integer :: out_type,nair,level, slv,abu(6),n
    real*8 eps
    character(len = 30) name,str
 
    n = iargc()
    call GETARG(1,name)
  
!    name = 'dtm1'
!    read_input(model_name,a,b,c,nair,level,eps,solver,sigma,freq,recv,out_type,abu,topo)
    call read_input(name,a,b,c,nair,level,eps,slv,sigma,freq,recv,out_type,abu)
! use the command line arguments instead    
    if (n > 1) then
      call getarg(2,str)
      read(str,*) level
      call getarg(3,str)
      read(str,*) eps
    endif 
    
    call fwd_solver(name,a,b,c,nair,level,slv,eps,sigma,freq,ex)
  
    call post_process(name,a,b,c,nair,abu,freq,recv,ex,out_type)
    
   !call output_data(model_name,freq,recv,out_type,arho)
    
    end program sfd_excmg

