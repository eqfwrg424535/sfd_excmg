module file_operation
   ! module mainly for reading and writing files
     use global_parameter
    implicit none
    
    contains
    
    subroutine read_input(model_name,a,b,c,nair,level,eps,solver,sigma,freq,recv,out_type,abu,topo)     
    implicit none
     real*8,allocatable::a(:),b(:),c(:),sigma(:,:),freq(:),recv(:,:)
     real*8,allocatable,optional ::topo(:)
     integer ::nair,nf,nrec, level,out_type,topo_indic,i,j,k,nx,ny,nz,abu(6),conduc_type,solver
     real*8 ::eps,ext(2),sig_tmp
     integer ncell1,ncell2
     character(len=30):: model_name
    character(len=200):: string
    
    open(20, file = trim(adjustl(model_name))//'.fwd',status='old')
    ! some description
    read(20,*) string
    read(20,*) level,eps,solver,out_type
    if(level.lt.3)then
      write(*,*) 'Error: the grid level should be larger than 3!' 
      stop
    endif
    if(eps.ge.1e-8)  write(*,*) 'Warning: Too big tolerance, the SFD results might be wrong!' 
    if(solver.lt.0.or.solver.ge.5)then
      write(*,*) 'Error: Unsupported solver type!'
      stop
    endif
    if(out_type.lt.0.or.out_type.ge.3)then
      write(*,*) 'Error: Unsupported ouput type!'
      stop
    endif
    read(20,*) nx,ny,nz,nair,topo_indic
    if(nx/2**level.lt.2)then
      write(*,*) 'Error: the grid level is too big, please check'
      stop
    endif
    ncell1 = nx*ny*nair
    ncell2 = nx*ny*(nz+nair)
    allocate(a(nx),b(ny),c(nz+nair))
    read(20,*) a         
    read(20,*) b
    read(20,*) c(nair+1:nair+nz)
    write(*,*) 'size of the domain(excluding air)',sum(a),sum(b),sum(c(nair+1:nair+nz))
    read(20,*) ext(1),ext(2)
    read(20,*) abu ! the boundary of uniform cells, also the boundary for topography
   !  do i = 1,nair
   !     c(i) = c(nair+1)*ext(2)**(nair+1-i)
   !  enddo
    c(1:nair) = c(nz+nair:nz+1:-1)  
   write(*,*) 'air thickness',sum(c(1:nair))
   read(20,*) conduc_type
    allocate(sigma(nx*ny*(nz+nair),6))
    do k = 1,nz
      do j = 1,ny
         do i = 1,nx
             if (conduc_type.eq.1)then
               read(20,*) sig_tmp
               sigma((nair+k-1)*nx*ny+(j-1)*nx+i,1:3) = sig_tmp
               sigma((nair+k-1)*nx*ny+(j-1)*nx+i,4:6) = 0
             elseif(conduc_type.eq.2)then
               read(20,*) sigma((nair+k-1)*nx*ny+(j-1)*nx+i,1:3)
               sigma((nair+k-1)*nx*ny+(j-1)*nx+i,4:6) = 0
            else
              read(20,*) sigma((nair+k-1)*nx*ny+(j-1)*nx+i,:) !sigma_xx sigma_yy sigma_zz sigma_xy sigma_xz sigma_yz
            endif
         enddo
      enddo
    enddo
    sigma(1:nair*nx*ny,1:3) = sigma_air
   write(*,*) 'sigma range(without air)',minval(sigma(ncell1+1:ncell2,1)),maxval(sigma(ncell1+1:ncell2,1))
   ! sigma(1:nair*nx*ny,4) = sigma_air
   ! sigma(1:nair*nx*ny,6) = sigma_air
    close(20)
    ! read in the topo data and modify the sigma array
    if(topo_indic.ne.0)then
        call read_topo(model_name,nx,ny,c,nair,sigma,topo)
    endif
    ! read in the frequency and receiver location    
    open(21, file = trim(adjustl(model_name))//'.dat',status='old')
    read(21,*) nf
    allocate(freq(nf))
    do i = 1,nf
       read(21,*) freq(i)
    enddo
    read(21,*) nrec
    allocate(recv(nrec,3))
    do i = 1,nrec
       read(21,*) recv(i,:)
    enddo    
    close(21)
       
    end subroutine
    
    subroutine read_topo(model_name,nx,ny, c,nair,sigma,topo)
      implicit none
      real*8,allocatable::sigma(:,:),topo(:),c(:),ccz(:)
      integer:: n1,n2,n3,n4,nair,i,j,k,nx,ny
      character(len=30):: model_name
      real*8:: cstep
      
      cstep = c(nair+1)
     ! allocate(ccz(size(c)+1))
     ! ccz(1) = -sum(c(1:nair))
     ! do i = 1,nz
     !    ccz(i+1) =  -sum(c(1:nair))+sum(c(1:i))
     ! enddo
      open(22, file = trim(adjustl(model_name))//'.topo',status='old')
      read(22,*) n1,n2  ! the number of cells covered by the topo(n1,n2 are the begining and ending index of a&b
      ! n3,n4 are the begining and ending index of c)
      allocate(topo((n2-n1)**2))
      ! the topo array is useful in the post process, for the location of certain receivers
        do j = n1,n2
          do i = n1,n2
             read(22,*) topo((j-n1)*(n2-n1+1)+i+1-n1)
             n3 = topo((j-n1)*(n2-n1+1)+i+1-n1)/cstep
             do k = 1,n3  
                 sigma((nair-n3+k-1)*nx*ny+(j-1)*nx+i,1)= sigma_hhs 
                 sigma((nair-n3+k-1)*nx*ny+(j-1)*nx+i,4)= sigma_hhs
                 sigma((nair-n3+k-1)*nx*ny+(j-1)*nx+i,6)= sigma_hhs
             enddo
           enddo
         enddo 
      close(22)
    
    end subroutine
    
    
    subroutine output_data(model_name,freq,recv,out_type,arho,tipper)
      implicit none
      real*8 :: freq(:),recv(:,:)
      integer ::out_type,i,j,k,nf,nrec
      real*8,allocatable :: arho(:,:),tipper(:,:)
      character(len=30) model_name
      
      nf = size(freq)
      nrec =size(recv,1)
      open(19,file =trim(adjustl(model_name))//'.res',status='replace',action='write')
          do k = 1,nf
              do i= 1,nrec
                 write(19,145) recv(i,:),freq(k),j,arho((k-1)*nrec+i,:)   
               enddo
            enddo
      close(19)
      
      if(output_tipper)then
          open(19,file = trim(adjustl(model_name))//'.res',status='old',position ='append')
            do k = 1,nf
              do i= 1,nrec
                 write(19,145) recv(i,:),freq(k),j,tipper((k-1)*nrec+i,:)   
               enddo
            enddo
          close(19)
      endif
145 format(4(e12.5,2x),i1,2x,4(e12.5,2x))      
    end subroutine
end module
