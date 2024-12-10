module fwd3d
    ! main module
    use global_parameter
    use solvers
    use model_selection
    use excmg
  !  use mkl_pardiso
  !  use divergence_corr
    use file_operation
   implicit none
   
  ! integer,parameter :: solver=1,src_type =0
  ! logical,parameter :: hybrid = .false., reg=.false., sma =.true.
  ! logical,parameter :: exc=.true.
  ! integer,parameter :: n_water =32
  ! real(8),dimension(5),parameter:: cor_ts=(/50,0,0,1,1/)*1.d0
contains
  subroutine mesh_coarsen(a,b,c,cont,k,ncont,a1,b1,c1)
  !! another way of generating nested grids using coasening
  !! volume-averaged conductivity is used for simplicity
  !! k is the time of coasening, at least 1
  !! for simplicity, only 8->1 merging is considered 
    implicit none
    real(8),intent(in)::a(:),b(:),c(:),cont(:,:)
    integer nx,ny,nz,k,i,j,ne,m,asf,r,s,t,ind
   !  integer,allocatable:: abr1(:)
    real(8),allocatable,intent(out)::a1(:),b1(:),c1(:),ncont(:,:) 
    real(8) vbig,vsmall,tmp(6)

    nx =size(a); ny=size(b); nz= size(c)
    asf = 2**k 
    if (mod(nx,asf).ne.0.or.mod(ny,asf).ne.0.or.mod(nz,asf).ne.0.)then
      print *,'error in mesh_coasen, please check the input mesh'
      stop
    endif
    nx = nx/asf
    ny = ny/asf
    nz = nz/asf
    allocate(a1(nx),b1(ny),c1(nz))
   do i= 1,Nx
       a1(i) = sum(a((i-1)*asf+1:i*asf))
    enddo
    do i= 1,Ny
       b1(i) = sum(b((i-1)*asf+1:i*asf))
    enddo
    do i= 1,Nz
       c1(i) = sum(c((i-1)*asf+1:i*asf))
    enddo

    ne = nx*ny*nz
    print *,'current grid size',nx,ny,nz
   ! print *,size(cont)
    
    allocate(ncont(ne,6))
    do m = 1,nz
      do j =1,ny
        do i = 1,nx
           vbig = a1(i)*b1(j)*c1(m) 
          ! print *,vbig 
           tmp = 0
           do r = 1,asf
             do s = 1,asf
                do t =1,asf
                   ind =((m-1)*asf+r)*nx*ny*asf**2+((j-1)*asf+s)*nx*asf+(i-1)*asf+t
                   vsmall =a((i-1)*asf+t)*b((j-1)*asf+s)*c((m-1)*asf+r) 
                   tmp = tmp+cont(ind,:)*vsmall
                enddo 
             enddo
          enddo
          ncont((m-1)*nx*ny+(j-1)*nx+i,:) = tmp/vbig
         enddo
       enddo  
     enddo
   
    ! allocate(abr1(size(abr)))
    ! do i = 1,size(abr)
       !   if(mod(i,2).eq.1) then
       !     abr1(i) =(abr(i)-1)/asf+1
       !   else
    !        abr1(i) = abr(i)/asf
       !   endif
    !  enddo
 
 end subroutine
  
subroutine mesh_refine(hc,a,b,c,Nx,Ny,Nz,abr,abu,k,N_air,a1,b1,c1)
  implicit none
    real(8),intent(in)::a(:),b(:),c(:)
    integer,intent(inout)::Nx,Ny,Nz,N_air,abr(:),abu(:)
    integer,intent(in)::k,hc
    integer::i,asd
    real(8),allocatable,intent(out)::a1(:),b1(:),c1(:)

    asd = hc**k
    Nx = Nx*hc
    Ny = Ny*hc
    Nz = Nz*hc
    allocate(a1(Nx))
    allocate(b1(Ny))
    allocate(c1(Nz))

    do i = 1,Nx/asd
        a1(asd*(i-1)+1:asd*i) = a(i)/asd
    enddo
    do i = 1,Ny/asd
        b1(asd*(i-1)+1:asd*i) = b(i)/asd
    enddo
    do i = 1,Nz/asd
        c1(asd*(i-1)+1:asd*i) = c(i)/asd
    enddo
      do i = 1,size(abr)
        if(mod(i,2).eq.1)then
           abr(i) = abr(i)*hc-(hc-1)
        else
          abr(i) = abr(i)*hc
        endif
      enddo
      do i = 1,size(abu)
        if(mod(i,2).eq.1)then
           abu(i) = abu(i)*hc-(hc-1)
        else
          abu(i) = abu(i)*hc
        endif
      enddo
    N_air = N_air*hc

 end subroutine

   subroutine avgsig(dt,np,a,b,c,sigma,sigv)
  !!function for calculating average conductivity at nodes
  !!follow the steps put forward by Li et al(2022,Geophysics)
    implicit none
    integer dt,nx,ny,nz,ne,nl,np,xnl,ynl
    real*8 a(:),b(:),c(:),h(6),sig1(8,3),sig2(6),ve(6)
    real*8,allocatable:: sigma(:,:),a1(:),b1(:),c1(:),sigv(:)
    integer i,j,k,m,ind,eid,eid2,eid3,eid4,ids(8),idl(6)

          nx = size(a); ny= size(b); nz =size(c)
          ne = nx*ny*nz; np = (nx+1)*(ny+1)*(nz+1)
          nl = nx*(ny+1)*(nz+1)+(nx+1)*ny*(nz+1)+(nx+1)*(ny+1)*nz
          xnl = nx*(ny+1)*(nz+1)
          ynl = (nx+1)*ny*(nz+1)
         ! allocate(avgsig(np))
           if(dt.eq.1)then
             allocate(sigv(np))
           else
             allocate(sigv(nl))
           endif
          ! the domaied is manually enlarged to ensure all nodes are connected
          ! with 6 edges, leading to a standard calculation process
          allocate(a1(nx+2),b1(ny+2),c1(nz+2))
          a1(2:nx+1) = a; a1(1) = a(1); a1(nx+2) = a(nx) 
          b1(2:ny+1) = b; b1(1) = b(1); b1(ny+2) = b(ny)
          c1(2:nz+1) = c; c1(1) = c(1); c1(nz+2) = c(nz)
          do k = 1,nz+1
            do j = 1,ny+1
               do i = 1,nx+1
                  eid = (k-2)*nx*ny+(j-2)*nx+i-1
                  h = (/a1(i),a1(i+1),b1(j),b1(j+1),c1(k),c1(k+1)/)
                  ids = eid+(/0,1,nx,nx+1,nx*ny,nx*ny+1,nx*ny+nx,nx*ny+nx+1/)
                  ve = (/(a1(i)+a1(i+1))*(b1(j)+b1(j+1))*c1(k),(a1(i)+a1(i+1))*(b1(j)+b1(j+1))*c1(k+1),&
                  (c1(k)+c1(k+1))*(b1(j)+b1(j+1))*a1(i),(c1(k)+c1(k+1))*(b1(j)+b1(j+1))*a1(i+1),&
                 (a1(i)+a1(i+1))*(c1(k)+c1(k+1))*b1(j),(a1(i)+a1(i+1))*(c1(k)+c1(k+1))*b1(j+1)/)/4
                  if(eid.gt.0)then
                    if(k.eq.nz+1) ids(5:8) = ids(1:4)
                    if(j.eq.ny+1)then
                        ids(3:4) = ids(1:2); ids(7:8)=ids(5:6)
                    endif
                    if(i.eq.nx+1) ids(2:8:2) = ids(1:7:2)
                  else
                    if(k.eq.1) ids(1:4) = ids(5:8) 
                    if(j.eq.1)then
                       ids(1:2) = ids(3:4); ids(5:6)=ids(7:8)
                    endif
                    if(i.eq.1) ids(1:7:2) = ids(2:8:2)            
                  endif 
                  sig1 = sigma(ids,1:3)
                if(dt.eq.2)then
                  !volume average over 6 edges
                  sig2(1) =(sig1(1,3)*h(1)*h(3)+sig1(2,3)*h(2)*h(3)+sig1(3,3)*h(1)*h(4)+&
                  sig1(4,3)*h(2)*h(4))/((h(1)+h(2))*(h(3)+h(4)))
                  sig2(2) =(sig1(5,3)*h(1)*h(3)+sig1(6,3)*h(2)*h(3)+sig1(7,3)*h(1)*h(4)+&
                  sig1(8,3)*h(2)*h(4))/((h(1)+h(2))*(h(3)+h(4)))
                  sig2(3) =(sig1(1,1)*h(3)*h(5)+sig1(3,1)*h(4)*h(5)+sig1(5,1)*h(3)*h(6)+&
                  sig1(7,1)*h(4)*h(6))/((h(3)+h(4))*(h(5)+h(6))) 
                  sig2(4) =(sig1(2,1)*h(3)*h(5)+sig1(4,1)*h(4)*h(5)+sig1(6,1)*h(3)*h(6)+&
                  sig1(8,1)*h(4)*h(6))/((h(3)+h(4))*(h(5)+h(6)))
                  sig2(5) =(sig1(1,2)*h(1)*h(5)+sig1(2,2)*h(2)*h(5)+sig1(5,2)*h(1)*h(6)+&
                  sig1(6,2)*h(2)*h(6))/((h(1)+h(2))*(h(5)+h(6)))
                  sig2(6) =(sig1(3,2)*h(1)*h(5)+sig1(4,2)*h(2)*h(5)+sig1(7,2)*h(1)*h(6)+&
                  sig1(8,2)*h(2)*h(6))/((h(1)+h(2))*(h(5)+h(6)))  
!                  !average at the nodes
                 idl =(/xnl+ynl+(k-1)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i,xnl+ynl+k*(nx+1)*(ny+1)+(j-1)*(nx+1)+i,&
                    (k-1)*nx*(ny+1)+(j-1)*nx+i-1, (k-1)*nx*(ny+1)+(j-1)*nx+i, &
                    xnl+(k-1)*ny*(nx+1)+(j-2)*(nx+1)+i, xnl+(k-1)*ny*(nx+1)+(j-1)*(nx+1)+i/) 
                 do m = 1,6
                    if(idl(m)<=0) idl(m)= 1
                    if(idl(m)>nl) idl(m) = nl
                 enddo 
                 if(k.ne.1.and.k.ne.nz+1.and.j.ne.1.and.j.ne.ny+1.and.i.ne.1.and.i.ne.nx+1) sigv(idl) =sig2
              else
                  do m  =1,3
                     sig2(m) =(sig1(1,m)*h(1)*h(3)*h(5)+sig1(2,m)*h(2)*h(3)*h(5)+&
                    sig1(3,m)*h(1)*h(4)*h(5)+sig1(4,m)*h(2)*h(4)*h(5)+sig1(5,m)*h(1)*h(3)*h(6)+&
                    sig1(6,m)*h(2)*h(3)*h(6)+sig1(7,m)*h(1)*h(4)*h(6)+&
                    sig1(8,m)*h(2)*h(4)*h(6))/((h(1)+h(2))*(h(3)+h(4))*(h(5)+h(6)))
                   enddo

                 ind = (k-1)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i
                if(i.ne.1.and.i.ne.nx+1.and.j.ne.1.and.j.ne.ny+1.and.k.ne.1.and.k.ne.nz+1)then
                    sigv(ind) = sum(sig2(1:3))/3
                else
                   sigv(ind) = 1.d0
                endif
             
               endif
             enddo
            enddo
          enddo

          deallocate(a1,b1,c1)

        end subroutine


 subroutine fd_assemb(ff,a,b,c,sigma,v,vr,vc) 
  implicit none
  real*8,allocatable :: a(:),b(:),c(:),sigma(:,:),aax(:),bby(:),ccz(:) 
  complex(8),allocatable ::v(:),p1(:),p2(:)
  complex(8) coef
  integer,allocatable ::vc(:),vr(:)
  real*8 ff,w,hx,hy,hz,sigmas,ve,am,sigmav(4)
  integer sid,nair,axn,byn,czn,h,nl,xnl,ynl,znl,nx,ny,nz,np,i,j,k, &
      ind2,l,m,n,ie,kms
  integer(8) nnz,ind

      w = 2*pi*ff
      coef = cmplx(0,1,8)*w*mu
    
      nx = size(a); ny =size(b); nz =size(c)
      xnl = nx*(ny+1)*(nz+1)
      ynl = ny*(nz+1)*(nx+1)
      nl = xnl+ynl+nz*(nx+1)*(ny+1)
      nnz = int(nl,8)*13
      allocate(v(nnz),vr(nnz),vc(nnz))    
      v = 0
      !!Ex
      do k = 2,nz
         do j = 2,ny
            do i = 1,nx
               ie = (k-2)*nx*ny+(j-2)*nx+i
	        ve = a(i)*(b(j-1)+b(j))*(c(k-1)+c(k))/4
           ! print *,a(i),b(j-1),c(k-1),b(j),c(k)
               sigmas = (sigma(ie,1)*b(j-1)*c(k-1)+sigma(ie+nx,1)*b(j)*c(k-1)+&
	             sigma(ie+nx*ny,1)*b(j-1)*c(k)+sigma(ie+nx*ny+nx,1)*b(j)*c(k))/ &
			    ((b(j-1)+b(j))*(c(k-1)+c(k)))
               sigmav = (/(sigma(ie,4)*b(j-1)*c(k-1)+sigma(ie+nx*ny,4)*b(j-1)*c(k))/2,&
                         (sigma(ie+nx,4)*b(j)*c(k-1)+sigma(ie+nx*ny+nx,4)*b(j)*c(k))/2,&
                         (sigma(ie,5)*b(j-1)*c(k-1)+sigma(ie+nx,5)*b(j)*c(k-1))/2,&
               (sigma(ie+nx*ny,5)*b(j-1)*c(k)+sigma(ie+nx+nx*ny,5)*b(j)*c(k))/2/)/((b(j-1)+b(j))*(c(k-1)+c(k)))
                ind = int(k-1,8)*nx*(ny+1)+int(j-1,8)*nx+int(i,8) 
                vr((ind-1)*13+1:ind*13) = ind            
                v((ind-1)*13+1) = -2/(c(k-1)*(c(k-1)+c(k)))
                vc((ind-1)*13+1) = ind -nx*(ny+1)
                v((ind-1)*13+2) = -2/(b(j-1)*(b(j-1)+b(j)))
                vc((ind-1)*13+2) = ind -nx
                v((ind-1)*13+4) = -2/(b(j)*(b(j-1)+b(j)))
                vc((ind-1)*13+4) = ind +nx
                v((ind-1)*13+5) = -2/(c(k)*(c(k-1)+c(k)))
                vc((ind-1)*13+5) = ind +nx*(ny+1)
                v((ind-1)*13+3) = -sum(v((ind-1)*13+1:(ind-1)*13+5))-coef*sigmas
                vc((ind-1)*13+3) = ind
                v((ind-1)*13+6) = 2/(a(i)*(b(j-1)+b(j)))
                vc((ind-1)*13+6) = xnl+(k-1)*(nx+1)*ny+(j-2)*(nx+1)+i
                v((ind-1)*13+8) = -v((ind-1)*13+6); v((ind-1)*13+7) =-v((ind-1)*13+6); 
                v((ind-1)*13+9) = v((ind-1)*13+6)
                v((ind-1)*13+6:(ind-1)*13+7) = v((ind-1)*13+6:(ind-1)*13+7)-coef*sigmav(1)
                v((ind-1)*13+8:(ind-1)*13+9) = v((ind-1)*13+8:(ind-1)*13+9)-coef*sigmav(2)
                vc((ind-1)*13+7) = xnl+(k-1)*(nx+1)*ny+(j-2)*(nx+1)+i+1
                vc((ind-1)*13+8) = xnl+(k-1)*(nx+1)*ny+(j-1)*(nx+1)+i
                vc((ind-1)*13+9) = xnl+(k-1)*(nx+1)*ny+(j-1)*(nx+1)+i+1
                v((ind-1)*13+10) = 2/(a(i)*(c(k-1)+c(k)))
                vc((ind-1)*13+10) = xnl+ynl+(k-2)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i
                v((ind-1)*13+12) = -v((ind-1)*13+10); v((ind-1)*13+11) =-v((ind-1)*13+10); 
                v((ind-1)*13+13) = v((ind-1)*13+10)
                vc((ind-1)*13+11) = xnl+ynl+(k-2)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i+1
                vc((ind-1)*13+12) = xnl+ynl+(k-1)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i
                vc((ind-1)*13+13) = xnl+ynl+(k-1)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i+1
                v((ind-1)*13+10:(ind-1)*13+11) =v((ind-1)*13+10:(ind-1)*13+11)-coef*sigmav(3)
                v((ind-1)*13+12:(ind-1)*13+13) =v((ind-1)*13+12:(ind-1)*13+13)-coef*sigmav(4)

               ! v((ind-1)*13+1:ind*13) = v((ind-1)*13+1:ind*13)*ve
            enddo
        enddo
      enddo
    !!Ey
      do k = 2,nz
        do j = 1,ny
           do i = 2,nx
            ie = (k-2)*nx*ny+(j-1)*nx+i-1
            ve = b(j)*(a(i-1)+a(i))*(c(k-1)+c(k))/4
            sigmas = (sigma(ie,2)*a(i-1)*c(k-1)+sigma(ie+1,2)*a(i)*c(k-1)+&
                  sigma(ie+nx*ny,2)*a(i-1)*c(k)+sigma(ie+nx*ny+1,2)*a(i)*c(k))/ &
                 ((a(i-1)+a(i))*(c(k-1)+c(k)))
            sigmav = (/(sigma(ie,4)*a(i-1)*c(k-1)+sigma(ie+nx*ny,4)*a(i-1)*c(k))/2,&
                       (sigma(ie+1,4)*a(i)*c(k-1)+sigma(ie+nx*ny+1,4)*a(i)*c(k))/2 ,&
                       (sigma(ie,6)*a(i-1)*c(k-1)+sigma(ie+1,6)*a(i)*c(k-1))/2 ,&
               (sigma(ie+nx*ny,6)*a(i-1)*c(k)+sigma(ie+nx*ny+1,6)*a(i)*c(k))/2/)/((a(i-1)+a(i))*(c(k-1)+c(k)))

             ind =int(xnl,8)+int(k-1,8)*(nx+1)*ny+int(j-1,8)*(nx+1)+int(i,8) 
                 vr((ind-1)*13+1:ind*13) = ind            
                 v((ind-1)*13+1) = 2/(b(j)*(a(i-1)+a(i)))
                 vc((ind-1)*13+1) = (k-1)*nx*(ny+1)+(j-1)*nx+i-1
                 v((ind-1)*13+2) = -v((ind-1)*13+1); v((ind-1)*13+3) =-v((ind-1)*13+1); 
                 v((ind-1)*13+4) = v((ind-1)*13+1)
                 vc((ind-1)*13+2) = (k-1)*nx*(ny+1)+(j-1)*nx+i
                 vc((ind-1)*13+3) = (k-1)*nx*(ny+1)+(j-1)*nx+i-1+nx
                 vc((ind-1)*13+4) = (k-1)*nx*(ny+1)+(j-1)*nx+i+nx
                v((ind-1)*13+1:(ind-1)*13+3:2) =v((ind-1)*13+1:(ind-1)*13+3:2)-coef*sigmav(1)
                v((ind-1)*13+2:(ind-1)*13+4:2)=v((ind-1)*13+2:(ind-1)*13+4:2)-coef*sigmav(2)
                 v((ind-1)*13+5) = -2/(c(k-1)*(c(k-1)+c(k)))
                 vc((ind-1)*13+5) = ind -(nx+1)*ny
                 v((ind-1)*13+6) = -2/(a(i-1)*(a(i-1)+a(i)))
                 vc((ind-1)*13+6) = ind -1
                 v((ind-1)*13+8) = -2/(a(i)*(a(i-1)+a(i)))
                 vc((ind-1)*13+8) = ind +1
                 v((ind-1)*13+9) = -2/(c(k)*(c(k-1)+c(k)))
                 vc((ind-1)*13+9) = ind +(nx+1)*ny
                 v((ind-1)*13+7) = -sum(v((ind-1)*13+5:(ind-1)*13+9))-coef*sigmas
                 vc((ind-1)*13+7) = ind

                 v((ind-1)*13+10) = 2/(b(j)*(c(k-1)+c(k)))
                 vc((ind-1)*13+10) = xnl+ynl+(k-2)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i
                 v((ind-1)*13+12) = -v((ind-1)*13+10); v((ind-1)*13+11) =-v((ind-1)*13+10); 
                 v((ind-1)*13+13) = v((ind-1)*13+10)
                 vc((ind-1)*13+11) = xnl+ynl+(k-2)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i+nx+1
                 vc((ind-1)*13+12) = xnl+ynl+(k-1)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i
                 vc((ind-1)*13+13) = xnl+ynl+(k-1)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i+nx+1
                  v((ind-1)*13+10:(ind-1)*13+11)=v((ind-1)*13+10:(ind-1)*13+11)-coef*sigmav(3)
                v((ind-1)*13+12:(ind-1)*13+13)=v((ind-1)*13+12:(ind-1)*13+13)-coef*sigmav(4)
               !  v((ind-1)*13+1:ind*13) = v((ind-1)*13+1:ind*13)*ve
           enddo
       enddo
     enddo
    !!Ez  
     do k = 1,nz
        do j = 2,ny
           do i = 2,nx
            ie = (k-1)*nx*ny+(j-2)*nx+i-1
            ve = c(k)*(b(j-1)+b(j))*(a(i-1)+a(i))/4
            sigmas = (sigma(ie,3)*b(j-1)*a(i-1)+sigma(ie+nx,3)*b(j)*a(i-1)+&
                  sigma(ie+1,3)*b(j-1)*a(i)+sigma(ie+nx+1,3)*b(j)*a(i))/ &
                  ((b(j-1)+b(j))*(a(i-1)+a(i)))              
            sigmav =(/(sigma(ie,5)*b(j-1)*a(i-1)+sigma(ie+nx,5)*b(j)*a(i-1))/2,&
                      (sigma(ie+1,5)*b(j-1)*a(i)+sigma(ie+nx+1,5)*b(j)*a(i))/2,&
                      (sigma(ie,6)*b(j-1)*a(i-1)+sigma(ie+1,6)*b(j-1)*a(i))/2,&
           (sigma(ie+nx,6)*b(j)*a(i-1)+sigma(ie+nx+1,6)*b(j)*a(i))/2/)/((b(j-1)+b(j))*(a(i-1)+a(i)))   
               ind = int(xnl,8)+int(ynl,8)+int(k-1,8)*(nx+1)*(ny+1)+int(j-1,8)*(nx+1)+int(i,8) 
                  vr((ind-1)*13+1:ind*13) = ind            
                  v((ind-1)*13+1) = 2/(c(k)*(a(i-1)+a(i)))
                  vc((ind-1)*13+1) = (k-1)*nx*(ny+1)+(j-1)*nx+i-1
                  v((ind-1)*13+2) = -v((ind-1)*13+1); v((ind-1)*13+3) =-v((ind-1)*13+1); 
                  v((ind-1)*13+4) = v((ind-1)*13+1)
                  vc((ind-1)*13+2) = (k-1)*nx*(ny+1)+(j-1)*nx+i
                  vc((ind-1)*13+3) = (k-1)*nx*(ny+1)+(j-1)*nx+i-1+nx*(ny+1)
                  vc((ind-1)*13+4) = (k-1)*nx*(ny+1)+(j-1)*nx+i+nx*(ny+1)
                  v((ind-1)*13+1:(ind-1)*13+3:2)=v((ind-1)*13+1:(ind-1)*13+3:2)-coef*sigmav(1)
                v((ind-1)*13+2:(ind-1)*13+4:2)=v((ind-1)*13+2:(ind-1)*13+4:2)-coef*sigmav(2)

                  v((ind-1)*13+5) = 2/(c(k)*(b(j-1)+b(j)))
                  vc((ind-1)*13+5) = xnl+(k-1)*(nx+1)*ny+(j-2)*(nx+1)+i
                  v((ind-1)*13+6) = -v((ind-1)*13+5); v((ind-1)*13+7) =-v((ind-1)*13+5); 
                  v((ind-1)*13+8) = v((ind-1)*13+5)
                  vc((ind-1)*13+6) = xnl+(k-1)*(nx+1)*ny+(j-1)*(nx+1)+i
                  vc((ind-1)*13+7) = xnl+k*(nx+1)*ny+(j-2)*(nx+1)+i
                  vc((ind-1)*13+8) = xnl+k*(nx+1)*ny+(j-1)*(nx+1)+i
                 v((ind-1)*13+5:(ind-1)*13+7:2)=v((ind-1)*13+5:(ind-1)*13+7:2)-coef*sigmav(3)
                v((ind-1)*13+6:(ind-1)*13+8:2)=v((ind-1)*13+6:(ind-1)*13+8:2)-coef*sigmav(4)

                  v((ind-1)*13+9) = -2/(b(j-1)*(b(j-1)+b(j)))
                  vc((ind-1)*13+9) = ind -(nx+1)
                  v((ind-1)*13+10) = -2/(a(i-1)*(a(i-1)+a(i)))
                  vc((ind-1)*13+10) = ind -1
                  v((ind-1)*13+12) = -2/(a(i)*(a(i-1)+a(i)))
                  vc((ind-1)*13+12) = ind +1
                  v((ind-1)*13+13) = -2/(b(j)*(b(j-1)+b(j)))
                  vc((ind-1)*13+13) = ind +(nx+1)
                  v((ind-1)*13+11) = -sum(v((ind-1)*13+9:(ind-1)*13+13))-coef*sigmas
                  vc((ind-1)*13+11) = ind
               !   v((ind-1)*13+1:ind*13) = v((ind-1)*13+1:ind*13)*ve
           enddo
       enddo
     enddo

    ! print *,v(1194170)
 end subroutine


 subroutine dcfast(ff,ax,by,cz,siga,v,vr,vc,etype,v2,vc2,ex,kmat,ia,ja,p1,p2)
    implicit none
    integer nx,ny,nz,nl,xnl,ynl,znl,np,abu(6),ids(12)
    complex(8),allocatable :: v(:),v2(:),ex(:),kmat(:),p1(:),p2(:)
    integer,allocatable::me(:,:),vr(:),vc(:),vc2(:),ja(:),intsec(:)
    integer(8),allocatable::ia(:)
    integer(8) k,nnz,ind
    logical,allocatable:: bori(:),etype(:)
    complex(8) tmp,tmp2,rhs(3),rhs2(3),rhs3(3)
    real*8 ff,ax(:),by(:),cz(:),coor(12,3),w
    real*8,allocatable::siga(:,:), sigv(:)
    integer i,j,m,n,l,ind1,ind2,ind3,il1

    nx= size(ax); ny =size(by); nz =size(cz)
    np = (nx+1)*(ny+1)*(nz+1) 
    xnl = nx*(ny+1)*(nz+1)  
    ynl =(nx+1)*ny*(nz+1)
    znl = (nx+1)*(ny+1)*nz	
    nl = xnl+ynl+znl 
    w =2*pi*ff
   ! print *,xnl,xnl+ynl,nl,(nx+1)*(ny+1)
     allocate(me(nl,12))
    do l = 1,nz
       do m = 1,ny
          do n =1,nx
              il1 =  (l-1)*ny*nx+(m-1)*nx+n
             ind1 = (l-1)*nx*(ny+1)+(m-1)*nx+n
             ind2 = xnl+(l-1)*(nx+1)*ny+(m-1)*(nx+1)+n
             ind3 = xnl+ynl+(l-1)*(nx+1)*(ny+1)+(m-1)*(nx+1)+n
             me(il1,:) = (/ind1,ind1+nx,ind1+nx*(ny+1),ind1+nx*(ny+1)+nx,&
               ind2,ind2+(nx+1)*ny,ind2+1,ind2+(nx+1)*ny+1, &
               ind3,ind3+1,ind3+nx+1,ind3+nx+2/)
        enddo
      enddo
    enddo

     call avgsig(2,np,ax,by,cz,siga,sigv)
     allocate(bori(nl))
    bori = .false.
    !! index of boundary dofs, all 6 faces
    do k = 1,nz+1 
      do j = 1,ny+1
	 do i = 1,nx
	     ind = (k-1)*nx*(ny+1)+(j-1)*nx+i
           if(k.eq.1.or.k.eq.nz+1.or.j.eq.1.or.j.eq.ny+1)then
              bori(ind) = .true.
           endif
	  enddo
	enddo
    enddo  
    do k = 1,nz+1 
	  do j = 1,ny
	    do i = 1,nx+1
	       ind = xnl+(k-1)*(nx+1)*ny+(j-1)*(nx+1)+i
           if(k.eq.1.or.k.eq.nz+1.or.i.eq.1.or.i.eq.nx+1)then
              bori(ind) = .true.
           endif
	  enddo
	enddo
     enddo  
     do k = 1,nz
	do j = 1,ny+1
	  do i = 1,nx+1
	   ind = xnl+ynl+(k-1)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i
           if(i.eq.1.or.i.eq.nx+1.or.j.eq.1.or.j.eq.ny+1)then
              bori(ind) = .true.
           endif
	  enddo
       enddo
     enddo  
   ! deal with Dirichlet boundaries
    allocate(p1(nl),p2(nl))
    p1 = 0; p2 = 0
    do i = 1,ny+1
        p1((i-1)*nx+1:i*nx) = (1.d0,0.d0)
    enddo
     do k = 1,nz
        p1(k*nx*(ny+1)+1:k*nx*(ny+1)+nx) = ex(k+1)
        p1((k+1)*nx*(ny+1)-nx+1:(k+1)*nx*(ny+1)) = ex(k+1)
     enddo
      p1(xnl-nx*(ny+1)+1:xnl) = ex(nz+1)
  
  
   ! write(44,*) p1
   ! stop

     do i = 1,ny
        p2(xnl+(i-1)*(nx+1)+1:xnl+i*(nx+1)) = (1.d0,0.d0)
     enddo
     do k = 1,nz
       do i = 1,ny
         p2(xnl+k*(nx+1)*ny+(i-1)*(nx+1)+1) = ex(nz+1+k+1)
         p2(xnl+k*(nx+1)*ny+i*(nx+1)) = ex(nz+1+k+1)
       enddo
     enddo
     p2(xnl+ynl-(nx+1)*ny+1:xnl+ynl) = ex(2*nz+2)

     allocate(ia(nl+1))
     ia(1) = 1
     k = 0
     do i = 1,nl
        tmp = 0; tmp2 =0
        if(.not.bori(i))then
         do j= 1,13
          ind = int(i-1,8)*13+j
         ! if(vc(ind).gt.nl) print *,i,j
          if(vc(ind).gt.0)then
            if(.not.bori(vc(ind)).and.v(ind).ne.0)then
                k = k +1
                v(k) = v(ind)
                vc(k) = vc(ind)
            else	   
                tmp =tmp+v(ind)*p1(vc(ind))
                tmp2 =tmp2+v(ind)*p2(vc(ind))
            endif
           endif
        enddo
        p1(i) = p1(i)-tmp
        p2(i) = p2(i)-tmp2
    else
        k =k +1
        v(k) = 1.d0
        vc(k) = i
    endif
    ia(i+1) = k+1
   enddo

   nnz = k
   print *,'number of nonzeros',nnz
   allocate(kmat(nnz),ja(nnz))
    kmat = v(1:nnz)
    ja = vc(1:nnz)
!   open(11,file='fdfe_mat.dat') 
!   do i =1,nl
!    do j =ia(i),ia(i+1)-1
!     write(11,605) i,ja(j),real(kmat(j)),imag(kmat(j))
!     if(isnan(abs(kmat(j))))then
!         print *,i,j,ja(j)
!         stop
!      endif
!    enddo 
!   enddo
!605 format(2(i8,2x),2(e9.2,2x))
!   close(11)
!   stop
   deallocate(etype,bori)
   deallocate(v,vc,vr)
 end subroutine
 
                         ! (name,a,b,c,nair,abu,freq,recv,ex,out_type,topo)
 subroutine post_process(modelname,a,b,c,nair,abr,freq,recv,u,out_type,topo)
    !!calculating rho and phi in double-mode
    implicit none

    real(8),intent(in)::a(:),b(:),c(:),freq(:),recv(:,:)
    integer,intent(in)::nair,abr(:),out_type
    real(8),allocatable,optional:: topo(:)
     complex(8),allocatable,intent(in) :: u(:,:)
    integer :: nx,ny,nz,i,j,k,k1,m,n,halfi,unx,n1,n2,nl,idx,idx2,idx3,nrec,nf
    complex(8) ::exz,ezx,eyz,ezy,exy,eyx,tmp,eht(8),zz(2),tx,ty
    complex(8),allocatable::x1(:),x2(:),ex(:,:,:),ey(:,:,:),hx(:,:,:),hy(:,:,:),zxy(:,:),&
       zyx(:,:),arho(:,:),arho2(:,:),ez(:,:,:),hz(:,:,:)
    real(8),allocatable::resxy(:,:),phasexy(:,:),resyx(:,:),phaseyx(:,:),&
       results(:,:),aax(:),bby(:),ccz(:),result2(:,:)
    real(8) ff,omega
    character(len=30) modelname
    character(len=12) str1,str2
         !print *,x1(2703)
         k = nair
         m = k
        Nx=size(a);Ny=size(b);Nz=size(c);
        unx = nx-1
        nrec = size(recv,1);nf = size(freq)
       ! m = NAIR*(NX-1)*(NY-1)+(NAIR-1)*(NX*(NY+1)+(NX+1)*NY)
       ! n1 = (nx+1)*(ny+1)+(nx+1)*ny+(ny+1)*nx
       ! n2 = (nx+1)*ny+(ny+1)*nx
        n1 = nx*(ny+1)*(nz+1)
        n2 = ny*(nz+1)*(nx+1)
        nl = n1+n2+nz*(nx+1)*(ny+1)
        halfi = int(ny/2)
        allocate(ex(nx,ny,2),ey(nx,ny,2),ez(nx,ny,1))
        allocate(hx(nx,ny,2),hy(nx,ny,2),hz(nx,ny,2))
        allocate(zxy(nx,ny),zyx(nx,ny),resxy(nx,ny),resyx(nx,ny))
        allocate(phasexy(nx,ny),phaseyx(nx,ny),results(nrec*nf,4),result2(nrec*nf,4))
        allocate(arho(nrec,5),arho2(nrec,5))
        allocate(x1(nl),x2(nl))
        !           do i = 1,nrec
         allocate(aax(nx),bby(ny),ccz(nz))
         do i = 1,nx
           aax(i) = -sum(a)/2+sum(a(1:i))-a(i)/2
          enddo
         do i = 1,ny
           bby(i) = -sum(b)/2+sum(b(1:i))-b(i)/2
         enddo
         do i = 1,nz
            ccz(i) = -sum(c(1:nair))+sum(c(1:i))-c(i)/2
         enddo
       ! print *,topography 
      do k = 1,size(freq)
          ff = freq(k)
          x1 = u(1:nl,k)
          x2 = u(nl+1:2*nl,k)
        do j = 1,ny
           do i = 1,nx
            !! topo modification
            if(present(topo))then
              if(i.ge.abr(1).and.i.le.abr(2).and.j.ge.abr(3).and.j.le.abr(4))then
                do k1 = abr(5),abr(6)
                  ! if(ccz(k1).ge.-topo_func(aax(i),bby(j)))then 
                    if(ccz(k1).ge.-topo((j-abr(3))*(abr(4)-abr(3)+1)+i+1-abr(1)))then
                       m = k1 
                       exit 
                   endif 
                enddo
              else
                m = nair
              endif
            endif
          !  print *,m
            idx = m*nx*(ny+1)+(j-1)*nx+i
            idx2 = n1+m*ny*(nx+1)+(j-1)*(nx+1)+i
            idx3 = n1+n2+m*(nx+1)*(ny+1)+(j-1)*(nx+1)+i
            ex(i,j,1) = (x1(idx)+x1(idx+nx))/2.d0
            ex(i,j,2) = (x2(idx)+x2(idx+nx))/2.d0
            ey(i,j,1) = (x1(idx2)+x1(idx2+1))/2.d0
            ey(i,j,2) = (x2(idx2)+x2(idx2+1))/2.d0
            ez(i,j,1) = (x1(idx3+1)+x1(idx3)+x1(idx3+nx+1)+x1(idx3+nx+2))/4
   
            exz = ((x1(idx+nx*(ny+1))+x1(idx+nx*(ny+1)+nx))/2.d0-ex(i,j,1))/c(nair+1)
            ezx = ((x1(idx3+1)+x1(idx3+nx+2))/2.d0-(x1(idx3+nx+1)+x1(idx3))/2.d0)/a(i)
            eyz = ((x1(idx2+ny*(nx+1))+x1(idx2+ny*(nx+1)+1))/2.d0-ey(i,j,1))/c(nair+1)
            ezy = (-(x1(idx3)+x1(idx3+1))/2.d0+(x1(idx3+nx+1)+x1(idx3+nx+2))/2.d0)/b(j)
            exy = (x1(idx+nx)-x1(idx))/b(j)
            eyx = (x1(idx2+1)-x1(idx2))/a(i) 
            hy(i,j,1) = (exz-ezx)/(cmplx(0,1.d0,8)*2.d0*pi*ff*mu)
           ! if(j==ny/2.and.i==2) print *,exz,ezx
            hx(i,j,1) = (ezy-eyz)/(cmplx(0,1.d0,8)*2.d0*pi*ff*mu)
            hz(i,j,1) = (exy-eyx)/(cmplx(0,1.d0,8)*2.d0*pi*ff*mu)
          ! if(j.eq.45.and.i.ge.28.and.i.le.38) print *,hy(i,j,1)
         enddo
      enddo
     call arho_interp2(1,2,a,b,ex,ey,hx,hy,hz,recv,arho)
   !  results(:,1:3) = abs(arho(:,1:3))
   !  results(:,4:6) = -atan(imag(arho(:,1:3))/real(arho(:,1:3)))*180.d0/pi
    !  results(:,1:3) = real(arho(:,1:3))  
    !  results(:,4:6) = imag(arho(:,1:3)) 
  !   open(9,file='casca_ey_yx_ssor.dat')
  !   do j = 1,ny
  !     do i =1,nx
  !        write(9,1236) real(ey(i,j,2)), imag(ey(i,j,2))
  !     enddo
  !   enddo
  !   close(9)

      do j = 1,ny
        do i = 1,nx
        !! topo modification
            if(present(topo))then
              if(i.ge.abr(1).and.i.le.abr(2).and.j.ge.abr(3).and.j.le.abr(4))then
                do k1 = abr(5),abr(6)
                   if(ccz(k1).ge.-topo((j-abr(3))*(abr(4)-abr(3)+1)+i+1-abr(1)))then
                       m = k1
                       exit
                   endif
                enddo
              else
                m = nair
              endif
            endif

            idx = m*nx*(ny+1)+(j-1)*nx+i
            idx2 = n1+m*ny*(nx+1)+(j-1)*(nx+1)+i
            idx3 = n1+n2+m*(nx+1)*(ny+1)+(j-1)*(nx+1)+i

            exz =((x2(idx+nx*(ny+1))+x2(idx+nx*(ny+1)+nx))/2.d0-ex(i,j,2))/c(nair+1)
            ezx =((x2(idx3+1)+x2(idx3+nx+2))/2.d0-(x2(idx3+nx+1)+x2(idx3))/2.d0)/a(i)
            eyz =((x2(idx2+ny*(nx+1))+x2(idx2+ny*(nx+1)+1))/2.d0-ey(i,j,2))/c(nair+1)
            ezy =(-(x2(idx3)+x2(idx3+1))/2.d0+(x2(idx3+nx+1)+x2(idx3+nx+2))/2.d0)/b(j)
            exy = (x2(idx+nx)-x2(idx))/b(j)
            eyx = (x2(idx2+1)-x2(idx2))/a(i) 
            hy(i,j,2) = (exz-ezx)/(cmplx(0,1.d0,8)*2.d0*pi*ff*mu)
            hx(i,j,2) = (ezy-eyz)/(cmplx(0,1.d0,8)*2.d0*pi*ff*mu)
            hz(i,j,2) = (exy-eyx)/(cmplx(0,1.d0,8)*2.d0*pi*ff*mu)
           enddo
       enddo

       call arho_interp2(2,2,a,b,ex,ey,hx,hy,hz,recv,arho2)
  
       do i = 1,nrec
        tmp =arho(i,3)*arho2(i,4)-arho2(i,3)*arho(i,4)
        zz(1) = (arho2(i,1)*arho(i,3)-arho(i,1)*arho2(i,3))/tmp
        zz(2) = (arho(i,2)*arho2(i,4)-arho2(i,2)*arho(i,4))/tmp
     !  zxy(i,j) = x1(idx)/hy(i,j,1)
       !  zyx(i,j) = ey(i,j,2)/hx(i,j,2)
       if(out_type.eq.2)then 
        results((k-1)*nrec+i,1) = abs(zz(1)**2*cmplx(0,1,8))/(2.d0*pi*ff*mu)
         results((k-1)*nrec+i,2) = abs(zz(2)**2*cmplx(0,1,8))/(2.d0*pi*ff*mu)
         results((k-1)*nrec+i,3) = -atan(imag(zz(1))/real(zz(1)))*180.d0/pi
         results((k-1)*nrec+i,4) = -atan(imag(zz(2))/real(zz(2)))*180.d0/pi
        if(i.eq.nrec/2+1) print *,results(i,:)
       else
         results((k-1)*nrec+i,1:4) = (/real(zz(1)),imag(zz(1)),real(zz(2)),imag(zz(2))/)  
       endif
       ! tippers
       tx = (arho2(i,4)*arho(i,5)-arho(i,4)*arho2(i,5))/tmp
       ty = (arho(i,3)*arho2(i,5)-arho2(i,3)*arho(i,5))/tmp
       result2((k-1)*nrec+i,1:4) = (/real(tx),imag(tx),real(ty),imag(ty)/)
       enddo
      enddo

    call output_data(modelname,freq,recv,out_type,results,result2)  
      
    deallocate(phasexy,phaseyx,arho)
    deallocate(arho2)
    deallocate(zxy,zyx,resxy,resyx)
    deallocate(ex,ey,ez,hx,hy,hz)
    deallocate(aax,bby,ccz)
    end subroutine
      
    subroutine arho_interp2(mode,rd,ax,by,ex,ey,hx,hy,hz,xdet,arho)
        implicit none
        complex(8),allocatable::ex(:,:,:),ey(:,:,:),hx(:,:,:),hy(:,:,:),hz(:,:,:),&
             arho(:,:)
        real*8,allocatable ::aax(:),bby(:)
        real*8 ypos,coef1(3),coef2(3),ax(:),by(:),xdet(:,:)
        complex(8) ehtmp(5)
        integer mode,rd,i,j,k,m,n,nx,ny,jrec,nrec,irec
  
         nx =size(ax); ny =size(by);nrec =size(xdet,1)
         allocate(aax(nx),bby(ny))
         aax(1) = (ax(1)-sum(ax))/2; bby(1) = (by(1)-sum(by))/2;
         do i = 2,nx
           aax(i)=-sum(ax)/2+sum(ax(1:i-1))+ax(i)/2
         enddo
         do i = 2,ny
           bby(i)=-sum(by)/2+sum(by(1:i-1))+by(i)/2
         enddo
       
  
         do i = 1,nrec
           ehtmp = 0
           do j = 1,nx-1
            if(xdet(i,1).ge.aax(j).and.xdet(i,1).le.aax(j+1))then
                 irec = j
                 exit
            endif
           enddo
          do j =1,ny-1
            if(xdet(i,2).ge.bby(j).and.xdet(i,2).le.bby(j+1))then
              jrec = j
              exit
            endif
          enddo
          !bi-quadratic interp
           coef1 = coef(xdet(i,1),aax(irec-1:irec+1))
           coef2 = coef(xdet(i,2),bby(jrec-1:jrec+1))
           do m = 1,3
             do n = 1,3
               if(rd.eq.1)then !along y-dir(default)
                ehtmp(1) = ehtmp(1)+coef1(m)*coef2(n)*ex(m+irec-2,n+jrec-2,mode)
                ehtmp(2) = ehtmp(2)+coef1(m)*coef2(n)*ey(m+irec-2,n+jrec-2,mode)
                ehtmp(3) = ehtmp(3)+coef1(m)*coef2(n)*hx(m+irec-2,n+jrec-2,mode)
                ehtmp(4) = ehtmp(4)+coef1(m)*coef2(n)*hy(m+irec-2,n+jrec-2,mode)
                ehtmp(5) = ehtmp(5)+coef1(m)*coef2(n)*hz(m+irec-2,n+jrec-2,mode)
              else !along x-dir
               ehtmp(1) = ehtmp(1)+coef1(m)*coef2(n)*ex(m+jrec-2,n+irec-2,mode)
                ehtmp(2) = ehtmp(2)+coef1(m)*coef2(n)*ey(m+jrec-2,n+irec-2,mode)
                ehtmp(3) = ehtmp(3)+coef1(m)*coef2(n)*hx(m+jrec-2,n+irec-2,mode)
                ehtmp(4) = ehtmp(4)+coef1(m)*coef2(n)*hy(m+jrec-2,n+irec-2,mode)
                ehtmp(5) = ehtmp(5)+coef1(m)*coef2(n)*hz(m+irec-2,n+jrec-2,mode)
              endif
             enddo
           enddo
            arho(i,:) =ehtmp
         !  arho(i,1:7:2) = real(ehtmp)
         !  arho(i,2:8:2) = aimag(ehtmp)
         enddo
  
         deallocate(aax,bby)
       contains
         function coef(x,x0)
           !! quadratic interpolation using Lagragian basis
           real*8 x,x0(3),coef(3)
           !real*8 val(3),l2
  
           coef = (/(x-x0(2))*(x-x0(3))/(x0(1)-x0(2))/(x0(1)-x0(3)),&
                    (x-x0(1))*(x-x0(3))/(x0(2)-x0(1))/(x0(2)-x0(3)),&
                    (x-x0(1))*(x-x0(2))/(x0(3)-x0(1))/(x0(3)-x0(2))/)
         !  l2 = dot_product(coef,val)
  
          end function
      end subroutine

  subroutine sig_interp(patm,a,b,c,n_air,conv,sigma)
  !! subroutine for interpolation sigma from caorse grid to fine grid
   real*8 patm(:),a(:),b(:),c(:)
   real*8,allocatable:: conv(:,:),sigma(:,:)
   integer n_air,nx,ny,nz,ne, nx0, ny0,nz0
   integer i,j,k,ind,m,ie(8),ide(8)
   integer,allocatable::ife(:)
   real*8 imat1(64,8),imat2(512,8),imat3(4096,8)

   nx0 =int(size(conv,2)**(1.d0/3))
   ny0 = nx0; nz0 = nx0
   nx = size(a); ny =size(b); nz =size(c)
   allocate(sigma(nx*ny*nz,6))
   sigma = 0
   select case(nx/nx0)
   case(2)
     imat1 = gen_imat(2) 
   case(4)
     imat2 = gen_imat(4) 
   case(8)
     imat3 = gen_imat(8)
   end select

    allocate(ife((nx/nx0*2)**3))
     do k = 1,nz/nz0*2
        do j = 1,ny/ny0*2
          do i =1,nx/nx0*2
            ife((k-1)*nx*ny/nx0/ny0*4+(j-1)*nx/nx0*2+i) =(k-1)*nx*ny+(j-1)*nx+i
          enddo
        enddo
     enddo

   !! air 
   sigma(1:n_air*nx*ny,1:3) = 1./patm(1)
   ie = (/1,2,1+nx0,2+nx0, 1+nx0*ny0,2+nx0*ny0,1+nx0*ny0+nx0,2+nx0*ny0+nx0/) 

   do k = 1,nz0/2
      do j = 1,ny0/2
        do i = 1,nx0/2
           ide = ie+(k-1)*2*nx0*ny0+(j-1)*2*nx0+(i-1)*2
           ind = n_air*nx*ny+(k-1)*nx*ny*nz/nz0*2+(j-1)*nx*ny/ny0*2+(i-1)*nx/nx0*2
         do m = 1,3
          select case(nx/nx0)
          case(2)
           sigma(ind+ife,m) = matmul(imat1,conv(1,ide))
          case(4)
           sigma(ind+ife,m) = matmul(imat2,conv(1,ide)) 
          case(8)
           sigma(ind+ife,m) = matmul(imat3,conv(1,ide))
          end select
         enddo
        enddo
      enddo
   enddo
 
  deallocate(ife)

  contains
    function gen_imat(as)
      integer as,i,j,k,l,m,n,ind,ind2
      real*8 gen_imat((2*as)**3,8),coor(16) 

         select case(as)       
          case(2)
            coor(1:4) =(/-3:3:2/)/4.d0  
          case(4)
            coor(1:8) =(/-7:7:2/)/8.d0 
          case(8)
            coor = (/-15:15:2/)/16.d0
         end select

       do l = 1,2
         do m = 1,2
           do n =1,2
              ind =(l-1)*4+(m-1)*2+n             
              do k = 1,as*2
                do j = 1,as*2
                  do i = 1,as*2
                     ind2 =(k-1)*as**2*4+(j-1)*as*2+i
                     gen_imat(ind2,ind) =l1(n,coor(i))*l1(m,coor(j))*l1(l,coor(k)) 
                  enddo
                enddo
              enddo
            enddo
          enddo
       enddo 
    end function

    function l1(isd,x)
       real*8 x,l1
       integer isd

       if(isd.eq.1)then
         l1 = x+0.5d0
       else
         l1 = -(x-0.5d0) 
       endif
    end function
  end subroutine

   subroutine doublemode_nor(modelname,ff,a,b,c,N_air,conv,mid,slv,tol,u1,u2)
        implicit none
           real(8),intent(in)::ff
           !real,intent(out):: tsol
           integer,intent(in)::N_air,slv,mid
           complex(8),allocatable::xx(:),e_dc(:),hx(:)
           complex(8),allocatable,intent(out) ::U1(:),u2(:)
           real(8),allocatable ::a(:),b(:),c(:),sigma(:,:),rho(:), arho(:,:),sigma_s(:,:),&
            conv(:,:),cont1(:,:),rmat(:)
           complex(8),allocatable::V(:),v2(:),k_mat(:),P1(:),P2(:),EX(:),x0(:)
           integer(4),allocatable::VC(:),vc2(:),VR(:),vr2(:),ndc(:),dc(:),ia(:),rr(:),rc(:)
           integer,allocatable::ia_n(:),Ja(:),xyz(:,:),idx(:),idx2(:),idx3(:),idx4(:)
           integer(8),allocatable::Ian(:)
           logical,allocatable::etype(:)
           real(8) ::rho_bj,tol,tol1d,contmp(3,3),patm(2)
           complex(8) imp(4)
           integer::maxit,i,j,k,ind,N1,N2,N3,N4,NE,NL,np,Nx,Ny,Nz,nab,nnz,ndcn,idd(7),idz(3)
           integer(8)::N_less
           integer::t1,t2,t3,t4,cr
           real :: iter1(3),iter2(3),tsol1,tsol2
            character(len=30) modelname
           character(len=2) :: string
   
           !xdet = (/-40:40:2/)*1.d3
           !rho_bj = patm(2)
          ! idz = (/1,4,6/)
           Nx=size(a);Ny=size(b);Nz=size(c)
           NE = Nx*Ny*Nz
           np = (nx+1)*(ny+1)*(nz+1)
           n1 = nx*(ny+1)*(nz+1)
           tol1d = 1.d-16
   
           ! tol = tol*10**(-n_grid)
           nl = NX*(NY+1)*(NZ+1)+(NX+1)*NY*(NZ+1)+(NX+1)*(NY+1)*NZ
           allocate(U1(nl),u2(nl))
           u1 = 0; u2 = 0
          ! call modelrho(PaTM,Nx,Ny,Nz,N_air,abr,rho)
         ! call load_sigma(PaTM,a,b,c,N_air,abr,conv,sigma)
         !  allocate(sigma_s(ne,6))
         !  if(csem)then
         !    sigma_s(:,1:3) = sigma(:,1:3)-1.d0/patm(2)
         !  else
         !    sigma_s =0
         !  endif
          print *,'sigma loaded'
   
          allocate(cont1(nz,6))
           cont1 = 0
           do i=1,nz
             cont1(i,1:3) = conv((i-1)*nx*ny+1,1:3)
           enddo
            !do i = 1,nab
            ! !! only upper part of sigma is needed since it's symmetric
            ! if(abr((i-1)*6+1).eq.1.and.abr((i-1)*6+3).eq.1)then
            !   contmp = conten_cal(conv(i,:))
            !   do j = abr((i-1)*6+5),abr(i*6)
            !     cont1(j,:) =(/contmp(1,1),contmp(2,2),contmp(3,3),&
            !     contmp(1,2), contmp(1,3),contmp(2,3)/)
            !   !  print *,cont1(j,:)
            !   !  stop
            !   enddo
            ! endif
            !enddo
          
        ! write(58,*) cont1
       !  allocate(ex(2*nz+2))
       if(bgf_type.eq.1)then
          call mt1d_anl(ff,cont1,c,e_dc,hx,imp)
          deallocate(hx)
       else
           patm = (/cont1(1,1),cont1(nz,1)/) ! assume the last layer is HHS
           call mt1dfem_n(ff,c,nz,n_air,patm,cont1,e_dc)
       endif
           !   write(59,*) e_dc
       deallocate(cont1)
   
      ! stop
         write(*,*) 'background model loaded'
           call system_clock(count_rate=cr)
           call system_clock(t1)
           call fd_assemb(ff,a,b,c,conv,V,vr,vc)
           allocate(v2(size(v)),vc2(size(v)),etype(nl))
           call dcfast(ff,a,b,c,conv,v,vr,vc,etype,v2,vc2,e_dc,k_mat,ian,ja,p1,p2)              
           call system_clock(t2)
           allocate(xx(size(p1)),ia(size(p1)+1))
           print *,'number of edges and inteerior edges',nl,size(ndc)
           ia = ian
           n_less = ian(size(p1)+1)-1
           nnz= n_less
           xx = 0
        !   do j = 1,ny+1
        !      xx((j-1)*nx+1:j*nx) = e_dc(1)
        !   enddo
        !   xx(1:nx*(ny+1)) = e_dc(1)
        !   do k = 1,nz
        !     do j = 1,ny+1
        !      ind = k*nx*(ny+1)+(j-1)*nx
        !      xx(ind+1:ind+nx) = e_dc(k+1)  ! using background field as initial solution
        !     enddo
        !   enddo
          print *,slv,tol
           maxit = 5000
           idd =(/1,nl,nl,np,75,ny,nz/)
           nnz = n_less
           call  system_clock(t1)
           if(mid.eq.1.or.mid.eq.2)then
            call pardiso_spd(k_mat,ia,ja,p1,1,xx,size(p1),nnz)
           else
            select case(slv)
             case(1)
               call bicg_stab_P(modelname,ff,'xy',a,b,c,conv,conv,e_dc,n_less,k_mat,ian,ja,p1, xx,tol,size(p1),idd,maxit,iter1)
             case(2)
               call ptfqmr(modelname,ff,'xy',n_less,a,b,c,conv,conv,e_dc,idd,K_mat, Ian, Ja, p1, xx,tol, nl,maxit, iter1)  
             case(3)
               call pgpbicg(modelname,ff,'xy',n_less, nl, a,b,c,conv,conv,e_dc,idd, k_mat, ian, ja, xx, p1, tol,maxit,iter1)  
             case(4)
               call pbicg(modelname, ff,'xy',idd,a,b,c,conv,conv,e_dc,k_mat,ian,ja,p1,xx,tol,nl,maxit,iter1(1))
            ! case(5) 
             !  call zagmg(nl,k_mat,ja,ia,p1,xx,10,6,100,maxit,tol)
             end select
           endif
           call system_clock(t2)
           write(*,"(/,' Time of slove(xy mode)  = ',F9.2,'s.')"), (t2 -t1)/real(cr)
           write(*,"('iter=',f8.1)"), iter1(1)
           tsol1 = (t2 -t1)/real(cr)
           U1 = xx
          ! open(94,file='fdfe_sol.dat')
          ! do i = 1,nl
          !   write(94,'(2(e12.5,2x))') real(u1(i)),imag(u1(i))
          ! enddo
          ! close(94) 
          ! stop
          maxit = 5000
          xx = 0
          !!yx mode
           call  system_clock(t1)
          if(mid.eq.1.or.mid.eq.2)then
            call pardiso_spd(k_mat,ia,ja,p2,1,xx,size(p1),nnz)
           else
            select case(slv)
             case(1)
               call bicg_stab_P(modelname,ff,'yx',a,b,c,conv,conv,e_dc,n_less,k_mat,ian,ja,p2, xx,tol,size(p1),idd,maxit,iter2)
             case(2)
               call ptfqmr(modelname,ff,'yx',n_less,a,b,c,conv,conv,e_dc,idd,K_mat, Ian, Ja, p2, xx,tol, nl,maxit, iter2)  
             case(3)
               call pgpbicg(modelname,ff,'yx',n_less, nl,a,b,c,conv,conv,e_dc,idd, k_mat, ian, ja, xx, p2, tol,maxit,iter2)  
             case(4)
               call pbicg(modelname, ff,'yx',idd,a,b,c,conv,conv,e_dc,k_mat,ian,ja,p2,xx,tol,nl,maxit,iter2(1))
            ! case(5)
            !   call zagmg(nl,k_mat,ja,ia,p2,xx,10,6,100,maxit,tol)
             end select
           endif
           call system_clock(t2)
           write(*,"(/,' Time of slove(yx mode)  = ',F9.2,'s.')"), (t2 -t1)/real(cr)
           write(*,"('iter=',f8.1)"), iter2(1)
           tsol2 = (t2 -t1)/real(cr)
           U2 =  xx
          ! write(108,*) u2
          !  call resistivity_cal(ff,a,b,c,n_air,1,u1,arho)
         !  stop 
         !  call post_process(ff,a,b,c,n_air,abr,u1,u2,xdet,arho)
         !  arho(1,1) = tsol1; arho(1,2) = tsol2
          ! print *,arho()
          ! stop
          ! deallocate(sigma_s,sigma)
           deallocate(xx)
           deallocate(k_mat,ia,ian,ja)
           deallocate(p1,p2)
           if(allocated(v2)) deallocate(v2,vc2)      
           if(allocated(vr2)) deallocate(vr2)

    end subroutine

!   subroutine singlemode_nor(modelname,ff,PATM,a,b,c,N_air,abr,abu,conv,grid_id,xdet,tol,maxit,u1,res)
!
!    implicit none
!    real(8),intent(in)::ff,patm(:),xdet(:)
!    !real,intent(out):: tsol
!    integer,intent(in)::N_air,grid_id,abr(:),abu(:)
!    complex(8),allocatable::xx(:),e_dc(:),xin(:),xx2(:),xx3(:)
!    complex(8),allocatable,intent(out) ::U1(:)
!    real(8),allocatable ::a(:),b(:),c(:),conv(:,:),sigma(:,:), arho(:,:),sigma_s(:,:), &
!          rho(:), sig1d(:),cont1(:,:)
!    !real(8),allocatable,optional:: mcpu(:,:)
!    real(8),allocatable,optional ::res(:,:)
!    complex(8),allocatable::V(:),v2(:),k_mat(:),P1(:),P2(:),EX(:),x0(:),hx(:)
!    integer(4),allocatable::VC(:),VR(:),vc2(:),vr2(:),ndc(:),dc(:),db(:)
!    integer,allocatable::ia(:),Ja(:),xyz(:,:),idx(:),idx2(:),idx3(:),idx4(:)
!    logical,allocatable::etype(:)
!    integer(8),allocatable::Ian(:)    
!    real(8) ::w,rho_bj,tol,tol1d,rec_loc(101),layh(3),layh2(5),patm2(3),&
!      patm3(4),contmp(3,3)
!    complex(8) extmp,hytmp,imp(4)
!    !parameter(cor_ts=(/0,0,950,1,1/)*1.d0)
!    integer::i,j,k,ind,N1,N2,N3,N4,NE,NL,Nx,Ny,Nz,nnz,tmp1(12),tmp2(9),d54(54),&
!           ne_air,np,ne2,ne3,nab,maxit
!    integer(8)::N_less
!    integer::t1,t2,t3,t4,cr,idd(7),nrec
!    real :: iter1(3),tsol
!    character(len=20) :: modelname
!    character(len=2) :: string
!
!   ! layh = (/-1.d5,0.d0,2.d2,1.7d3,1.8d3/) 
!   !  layh = (/-1.d5,0.d0,6.d2,8.5d2,3.15d3/) 
!    layh = (/-1.d4,0.d0,1.d3/)
!    rec_loc = (/-5000:5000:100/)*2.d0
!   ! patm2(1:4) = 1.d0/patm;patm2(5) =1.d0/patm(3)
!    layh2 = (/-1.d4,0.d0,1.d3,1.4d3,1.5d3/)
!    patm3 = (/1.d-9,1.d-1,1.d-2,2.d-3/) 
!    patm2 = 1.d0/patm(1:3)
!    idd = 0
!    w = 2*pi*ff
!    rho_bj = patm(2)
!    Nx=size(a);Ny=size(b);Nz=size(c);nab=size(abr)/6
!    NE = Nx*Ny*nz
!    ne_air=nx*ny*(abr(5)-1); ne2 =nx*ny*(abr(6)); ne3 = nx*ny*abr(12)
!    tol1d = 1.d-16
!    nl = NX*(NY+1)*(NZ+1)+(NX+1)*NY*(NZ+1)+(NX+1)*(NY+1)*NZ
!    np = (nx+1)*(ny+1)*(nz+1)
!    n1 = nx*(ny+1); n2 = nx
!     !write(32,*) e_dc
!       ! tol = tol*10**(-n_grid)
!        allocate(U1(nl),sigma_s(ne,6),rho(ne))
!        !call modelrho(PaTM,Nx,Ny,Nz,N_air,abr,rho)
!    
!      !    call sig_interp(patm,a,b,c,n_air,conv,sigma)
!      ! else
!        call load_sigma(PaTM,a,b,c,N_air,abr,conv,sigma) 
!      
!         rho = 1.d0/sigma(:,1)
!        !write(*,*) 'resistivity model loaded'
!        sigma_s = 0 
!
!     
!        !  allocate(sig1d(nz))
!        !  do i = 1,nz
!        !     sig1d(i) = sigma((i-1)*nx*ny+1,1)
!        !  enddo
!        !  call mt1dfem(ff,n_air,patm,nz+1,c,abr,TOL1d,MAXit,E_DC)          
!        allocate(cont1(nz,6))
!        cont1 = 0
!        do i=1,nz
!           if(i.le.n_air)then
!             cont1(i,1:3) = 1.d0/patm(1)
!           else
!             cont1(i,1:3) = 1.d0/patm(2)
!           endif
!        enddo
!    
!         do i = 1,nab
!          !! only upper part of sigma is needed since it's symmetric
!          if(abr((i-1)*6+1).eq.1.and.abr((i-1)*6+3).eq.1)then
!          contmp = conten_cal(conv(i,:))
!          do j = abr((i-1)*6+5),abr(i*6)
!            cont1(j,:) =(/contmp(1,1),contmp(2,2),contmp(3,3),&
!            contmp(1,2), contmp(1,3),contmp(2,3)/)
!          enddo
!          endif
!         enddo
!       
!          call mt1d_anl(ff,cont1,c,e_dc,hx,imp) 
!         deallocate(cont1)
!       !  write(63,*) e_dc
!   
!       ! (ff,N_air,PaTM,np,c,abr,tol,itmax,EX)
!  
!      
!    !  stop
!       call system_clock(count_rate=cr)
!       ! call system_clock(t1)
!        call fd_assemb(ff,a,b,c,conv,V,vr,vc)
!        allocate(v2(size(v)),vc2(size(v)),etype(nl))
!        call dcfast(ff,a,b,c,conv,v,vr,vc,etype,v2,vc2,e_dc,k_mat,ian,ja,p1,p2)
!    
!      !print *,vnorm(p1,nl) 
!      !if(grid_id.eq.2) stop
!       n_less = ian(nl+1)-1
!       ! allocate(xx(size(p1)),ia_n(size(p1)+1))
!        allocate(xx(nl))
!        xx =0
!       ! n1 = (nx+1)*(ny+1)+nx*(ny+1)+(nx+1)*ny
!       ! n2 = 2*nx+1 
!  
!     !  write(73,*) p1        
!!        print *,'dof', size(p1)
!      ! maxit = 1000
!       idd = (/1,nl,nl,np,33,ny,nz/)
!       allocate(ia(nl+1))
!       ia = ian
!       call  system_clock(t1)
!       if(nl.ge.1000000)then
!         call bicg_stab_P(modelname,ff,'yx',a,b,c,sigma,sigma_s,e_dc,n_less,k_mat,ian,ja,p2, xx,tol,size(p1),idd,maxit,iter1)
!        !  call pbicg(ff,idd,2,a,b,c,sigma,sigma_s,e_dc,K_mat, Ian, Ja, P1, xx,&
!        ! tol,size(p1), maxit, iter1(1))
!        !  call BiCGSTAB_pre_64(K_mat, Ian, Ja, p1, xx, tol, size(p1), maxit) 
!       !   call zagmg(size(p1),k_mat,ja,ia,p1,xx,10,6,50,maxit,tol)
!      ! call ptfqmr(n_less,a,b,c,sigma,sigma_s,e_dc,idd,K_mat, Ian, Ja, p2,&
!      !   xx,tol, nl,maxit, iter1)
!       ! call pidrs(n_less,idd, k_mat, ian, ja, p1, xx, tol, nl,10, maxit, iter1)
!        !  call p_qmr(n_less,k_mat,ian,ja,nl,idd, xx,p1,tol,maxit,iter1)
!       else
!        call pardiso_spd(k_mat,ia,ja,p2,1,xx,size(p1),nnz)
!       endif
!       deallocate(ia)
!
!        call system_clock(t2)
!        write(*,"(/,' Time of slove = ',F9.2,'s.')"), (t2 -t1)/real(cr)
!        write(*,"('iter=',i8)"), iter1(1)
!        tsol = (t2 -t1)/real(cr)
!      !  mcpu(1,1) = tsol; !mcpu(1,2) = iter1
!        u1 = xx
!       !! calculate some norm      
!    
!        do i = 1,nl
!           if(isnan(abs(u1(i))))then
!              print *,'nan in x',i
!           endif
!        enddo
!       ! stop
!       ! write(89,*) imag(xx)
!
!!        call post_process(ff,a,b,c,abr(6),abr,xx,xx,xdet,arho)
!
!        deallocate(e_dc,sigma_s,rho)
!      !  if(allocated(arho)) deallocate(arho)
!        deallocate(xx)
!        deallocate(k_mat,ian,ja)
!        deallocate(p1,p2)
!        if(allocated(xx2)) deallocate(xx2,xx3)
!       if(allocated(ndc)) deallocate(ndc)
!        if(allocated(v2)) deallocate(v2,vc2)
!        if(allocated(vr2)) deallocate(vr2)
!   
!    end subroutine

  
    subroutine u1u2(modelname,ff,a,b,c,n_air,ng,sigf,slv,tol,u01,u02,u11,u12)
      implicit none
      integer,intent(in) :: ng
      real(8),intent(in) :: ff
      integer :: nx,ny,nz,n1,n2,n_air,slv
      integer,allocatable:: abr1(:)
      real(8),allocatable ::a(:),b(:),c(:),a1(:),b1(:),c1(:),conv(:,:),sig(:,:),sigf(:,:)
      complex(8),allocatable,intent(out) :: u01(:),u02(:),u11(:),u12(:)
      real(8) tol
      character(len=30) modelname

      nx = size(a); ny = size(b); nz = size(c)
     ! write(*,*) '-------On level 1-------'
     !if(.not.uniform_refine)then
     ! call load_sigma(PaTM,a,b,c,N_air,abr,conv,sigf)
      call mesh_coarsen (a,b,c,sigf,ng-1,sig,a1,b1,c1)
   !   print *,abr1
      call doublemode_nor(modelname,ff,a1,b1,c1,N_air/2**(ng-1),sig,1,slv,tol,u01,u02)
    ! write(171,*) u0  
      deallocate(a1,b1,c1,sig) 
      call mesh_coarsen (a,b,c,sigf,ng-2,sig,a1,b1,c1)
      ! print *,size(a1),size(b1),size(c1),n_air
     ! write(*,*) '-------On level 2-------'
      call doublemode_nor(modelname,ff,a1,b1,c1,N_air/2**(ng-2),sig,2,slv,tol,u11,u12)
       deallocate(a1,b1,c1,sig) 
    !  write(183,*) u1
    ! print *,u0(1035),u1(8009)

   end subroutine

    subroutine doublemode_acc(modelname,ff,a,b,c,N_air,N_grid,conv,slv,tol,u21,u22)

        implicit none
        real(8),intent(in)::ff
        integer::N_air,maxit,N_grid,idd(7),slv
       ! real(8),allocatable :: arho(:,:)
        complex(8),allocatable::V(:),v2(:),k_mat(:),P1(:),p2(:),EX(:),xx(:),e_dc(:),hx(:)
        integer(4),allocatable::VC(:),VR(:),vc2(:),vr2(:),abr1(:)
        logical,allocatable::etype(:)
        complex(8),allocatable::u01(:),U02(:),U11(:),u12(:),u21(:),u22(:),xx2(:),uexp(:),utmp(:)
        integer,allocatable::Ja(:),xyz(:,:),idx(:),idx2(:),NDC(:),dc(:),db(:),ia(:)
        integer(8),allocatable::Ian(:)
        real(8)::tol,tol1d,rho_bj,w,res(4),layh(5),layh2(3),patm(2),patm2(3),contmp(3,3)
        !parameter(cor_ts=(/0,0,950,1,1/)*1.d0)
        integer::i,j,k,ind,ind2,ind3,N1,N2,NE,ne_air,ne2,ne3,NL,np,N_tol,Nx,Ny,Nz,mit,xnl
        integer(8)::N_less
        integer::nnz,l,t1,t2,cr,itmax1d,err,nab,nrec,idz(3)
        real :: iter1(3),iter2(3),t_sol1,t_sol2,wu !! (iter|cpu on each level)
        real(8),allocatable::a(:),b(:),c(:),rho(:),a1(:),b1(:),c1(:),arho(:,:),&
           conv(:,:),sigma(:,:),sigma_s(:,:),sigf(:,:),assf(:,:),sig1d(:),cont1(:,:)
        complex(8) extmp,hytmp,imp(4)
        character(len=30) modelname
       
         layh2 = (/-1.d4,0.d0,1.d3/) 
         layh = (/-1.d5,0.d0,1.d3,1.4d3,1.5d3/)
         idz = (/1,4,6/)
       !  patm2(1:4) = 1.d0/patm; patm2(5) = 1.d0/patm(3)
       !  patm2 = 1.d0/patm(1:3)
          idd = 0
         t_sol1 = 0
         Nx=size(a);Ny=size(b);Nz=size(c);
         nl = NX*(NY+1)*(NZ+1)+(NX+1)*NY*(NZ+1)+(NX+1)*(NY+1)*NZ
         NE = Nx*Ny*Nz
         tol1d= 1.d-10
       !  w = 2*pi*ff
         ALLOCATE(assf(n_grid,7))
         !  call  mesh_refine(a,b,c,Nx,Ny,Nz,abr,1,N_air,a1,b1,c1)
            !allocate(a1(nx),b1(ny),c1(nz))
          !  allocate(db(nx*(ny+1)))
          !  do i = 1,ny+1
          !     db((i-1)*nx+1:i*nx) = (/1:nx/)+(2*nx+1)*(i-1)
          !  enddo
            !a1 = a; b1 = b; c1 = c
        !   call load_sigma(PaTM,a,b,c,N_air,abr,conv,sigma)
            write(*,*) '------------------SFD with EXCMG(doublemode)-----------------'
            write(*,*) 'initial mesh size:  ',Nx,Ny,Nz 
             write(*,*) 'number of unknowns(edges):  ',nl
            write(*,"(  'frequency:',e9.2,'Hz')"),ff
            write(*,*) 'solver: ',slv
            write(*,*) 'preconditioner:  ',pre_type
             write(*,*) 'grid level:  ',n_grid
            write(*,*) 'tolerance:  ',tol
            write(*,*) 'tolerance and interval for DC:  ',corr_tol,corr_intv
         !   print *,n_grid
            tol = tol*10.d0**(-n_grid+1)
            !!loop of CMG
           
            call system_clock(count_rate=cr)
            call system_clock(t1)
            !(modelname,ff,a,b,c,n_air,ng,sigf,tol,u01,u02,u11,u12)
            call u1u2(modelname,ff,a,b,c,n_air,n_grid,conv,slv,tol,u01,u02,u11,u12)         
            !deallocate(a1,b1,c1)
            call system_clock(t2)
          !  deallocate(sigma) 
          !  nx = nx*3; ny = ny*3; nz = nz*3; 
            tol = tol*10
            t_sol1 =(t2 - t1)/real(cr)/2
            t_sol2 =(t2 - t1)/real(cr)/2
            allocate(sigma(size(conv,1),size(conv,2)))
            sigma = conv

         do i = 3,N_grid
              
            tol = tol*10
         !  if(.not.uniform_refine)then
           if(i.ne.n_grid)then 
             call mesh_coarsen (a,b,c,sigma,n_grid-i,sigf,a1,b1,c1)   
           else
             allocate(sigf(size(conv,1),size(conv,2)))
             allocate(a1(size(a)),b1(size(b)),c1(size(c)))
             a1 =a; b1 =b;c1 =c
             sigf = sigma
           endif
         !    call mesh_refine(2,a,b,c,Nx,Ny,Nz,abr,abu,i-1,N_air,a1,b1,c1)
         !    allocate(abr1(SIZE(abr)))
         !    abr1 = abr
         !  endif
          !  print *,abr
            nx = size(a1); ny = size(b1); nz =size(c1); !nab = size(abr)/6
            nl = NX*(NY+1)*(NZ+1)+(NX+1)*NY*(NZ+1)+(NX+1)*(NY+1)*NZ
            xnl = NX*(NY+1)*(NZ+1)
             ne =nx*ny*nz; 
            print *,'on level',i,nx,ny,nz,size(sigf) 
           !ne_air = nx*ny*(abr(5)-1);
            !ne2 =nx*ny*abr(6); !ne3 = nx*ny*abr(18)
            np =(nx+1)*(ny+1)*(nz+1)
          !  print *,nx
          !  if(allocated(u2).or.allocated(db)) print *,'yes'
            allocate(U21(nl),U22(NL))
            u21 = 0; u22 = 0 
          !  call linear_interp(Nx,Ny,nz,n_air,nl,u1,u2)
          !  call quadra_interp(Nx,Ny,nz,n_air,nl,u1,u2)
          !  call quad_interp_n(Nx,Ny,nz,n_air,nl,u1,u2)
          !  print *,size(u0),size(u1)
          ! if(anlcase)then
          !   u2 = extint_v(nx,ny,nz,nl,u0,u1,3,2)
          ! else
             u21 = extint_v(nx,ny,nz,nl,u01,u11,2,2)
             u22 = extint_v(nx,ny,nz,nl,u02,u12,2,2)
          ! endif
          ! write(18,*) u2
          ! write(19,*) u1
          ! stop   
           !   allocate(e_dc(nz+1),sig1d(nz))
           !   sig1d = sigma(1:ne:nx*ny,1)    
              ! print *,sig1d  
           !   call mt1dfem(ff,n_air,patm,nz+1,c1,abr,TOL1d,MAXit,E_DC)
           !   deallocate(sig1d)
            allocate(cont1(nz,6))
            cont1 = 0
           do k = 1,nz
             cont1(k,1:3) = sigf((k-1)*nx*ny+1,1:3)
           enddo
          ! print *,c1
          !if(i.eq.n_grid) write(77,*) cont1
          if(bgf_type.eq.1)then
           call mt1d_anl(ff,cont1,c1,e_dc,hx,imp)
          ! if(i.eq.n_grid) write(78,*) e_dc
           deallocate(hx)
          else
           patm = (/cont1(1,1),cont1(nz,1)/) ! assume the last layer is HHS
           call mt1dfem_n(ff,c1,nz,n_air/2**(n_grid-i),patm,cont1,e_dc)
          endif
           deallocate(cont1)

           ! write(*,*) 'resistivity model loaded'     
           call fd_assemb(ff,a1,b1,c1,sigf,V,vr,vc)
           allocate(v2(size(v)),vc2(size(v)),etype(nl))
           call dcfast(ff,a1,b1,c1,sigf,v,vr,vc,etype,v2,vc2,e_dc,k_mat,ian,ja,p1,p2)
           ! debug para
           do k = 1,nl
              if(isnan(abs(p1(k))))then                
                 print *,'error: NaN value found in the RHS(p1)',k
                 stop
              endif     
           enddo
           allocate(xx(size(p1)),ia(size(p1)+1))
           xx = u21; ia = ian
         !  write(92,*) xx
          ! xx = 0
           !if(i.eq.n_grid-1) xx= 0
          ! do k = 1,nz+1
          !    xx(xnl+(k-1)*(nx+1)*ny:xnl+k*(nx+1)*ny) = e_dc(k)
          ! enddo
           n_less = ian(nl+1)-1; nnz = n_less 
           idd =(/1,np,nl,np,33,ny,nz/)
           maxit = 1000
          ! mit = 10*8**(n_grid-i)  
          ! print *,mit
           call  system_clock(t1)
          ! call pbicg(ff,idd,1,a1,b1,c1,sigma,sigma_s,e_dc, K_mat, Ian, Ja, P1, xx,&
          !  tol,size(p1), maxit, iter1(1))
           select case(slv)
          case(1) !default
           call  Bicg_stab_p(modelname,ff,'xy',a1,b1,c1,sigf,sigf,e_dc,n_less, K_mat, Ian, Ja, p1, xx, tol, size(p1),idd, maxit, iter1)
          case(2)
            call ptfqmr(modelname,ff,'xy',n_less,a1,b1,c1,sigf,sigf,e_dc,idd,K_mat, Ian, Ja, p1, xx,tol, nl,maxit, iter1)  
          case(3)
            call pgpbicg(modelname,ff,'xy',n_less, nl, a1,b1,c1,sigf,sigf,e_dc,idd, k_mat, ian, ja, xx, p1, tol,maxit,iter1)  
          case(4)
            call pbicg(modelname, ff,'xy',idd,a1,b1,c1,sigf,sigf,e_dc,k_mat,ian,ja,p1,xx,tol,nl,maxit,iter1(1))
          end select
          call system_clock(t2)
          t_sol1 = t_sol1 + (t2 - t1)/real(cr)
          ! write(*,*) 'dof of A',size(p1)
           write(*,"(/,' Time of slove (xy)  = ',F8.2,'s.')"), t_sol1
           write(*,"(/,'iter=',f8.1)"), iter1(1)
           u21 = xx
         ! stop 
           
           xx = u22
          call  system_clock(t1)
             select case(slv)
           case(1) !default
           call  Bicg_stab_p(modelname,ff,'yx',a1,b1,c1,sigf,sigf,e_dc,n_less, K_mat, Ian, Ja, p2, xx, tol, size(p1),idd, maxit, iter2)
          case(2)
            call ptfqmr(modelname,ff,'yx',n_less,a1,b1,c1,sigf,sigf,e_dc,idd,K_mat, Ian, Ja, p2, xx,tol, nl,maxit, iter2)  
          case(3)
            call pgpbicg(modelname,ff,'yx',n_less, nl, a1,b1,c1,sigf,sigf,e_dc,idd, k_mat, ian, ja, xx, p2, tol,maxit,iter2)  
          case(4)
            call pbicg(modelname, ff,'yx',idd,a1,b1,c1,sigf,sigf,e_dc,k_mat,ian,ja,p2,xx,tol,nl,maxit,iter2(1))
           end select
          call system_clock(t2)
          ! write(*,*) 'dof of A',size(p1)
          t_sol2 = t_sol2 + (t2 - t1)/real(cr)
           write(*,"(/,' Time of slove (yx)  = ',F8.2,'s.')"), t_sol2
           write(*,"(/,'iter=',f8.1)"), iter2(1)
          ! m_sol(i,1) = iter1(1); m_sol(i,2) = t_sol1;
          ! stop
          ! write(*,*) 'l-inf of the intial guess ',maxval(abs(u2-xx))
          !wu = wu+iter1(1)*27.d0**(i-n_grid)
          ! allocate(db(nx*(ny+1)))
!          if(anlcase)then
!            res(1) = maxval(abs(xx-e_dc))
!            res(2) = vnorm(xx-e_dc,nl)/sqrt(nl*1.d0)
!           !  write(*,*) 'l-inf and l-2 norm of the FD solution',res(1:2)
!             assf(i,1) = res(1)
!             assf(i,3) = res(2)
!          !   assf(i,5) = maxval(abs(u2-xx))
!             assf(i,5) = vnorm(xx-u2,nl)/sqrt(nl*1.d0)
!            if(i.gt.2)then
!              assf(i,7) = assf(i,5)/res(2)
!              assf(i,2:6:2) = log(assf(i-1,1:5:2)/assf(i,1:5:2))/log(3.d0)
!             endif
!           !! extrapolate for higher order
!            allocate(uexp(nl),utmp(nl))
!           !  call cubic_interp(nx,ny,nz,nl,u1,utmp,2)
!             call qc_interp(nx,ny,nz,nl,u1,utmp,2)
!             uexp = xx*9/8-utmp/8
!            print *,"|u_l-u|_inf", maxval(abs(uexp-e_dc))
!            print *,"|u_l-u|_2", vnorm(uexp-e_dc,nl)/sqrt(nl*1.d0)
!            deallocate(uexp,utmp)
!           !! accuracy validation
!            if(i.eq.n_grid)then
!              write(*,*) '|uh-u|_inf    order   |uh-u|_2    order   |uh-wh|_inf &
!                order     rh'
!               do l = 1,n_grid
!                 write(*,989) assf(l,:)
!               enddo
!989     format(3(e12.5,2x,f5.2,2x),e10.3)
!             endif
!            endif

            u22 = xx
          ! if(i.eq.n_grid)  print *,ff,'total CPU(s)',t_sol1
            !  print *,'WU',wu
            !  return
           
            !u2 = 0
          !  do j = 1,ny+1
          !     db((j-1)*nx+1:j*nx) = (/1:nx/)+(2*nx+1)*(j-1)
          !  enddo
          
          ! u2(db) = (1.d0,0.d0)
          !!call resistivity_cal()
   
          if(i<N_grid)then
           ! if(i.eq.n_grid-1) call post_process(ff,a1,b1,c1,n_air/2,abr1,xx,xx,rec_loc,arho)
              deallocate(u01,u02)
              allocate(u01(size(u11)),u02(size(u12)))
              u01 = u11
              u02 = u12
              deallocate(u11,u12)
              allocate(U11(nl),u12(nl)) 
               U11 = U21
               u12 = u22
              deallocate(U21,u22)
            !  deallocate(V,VR,VC)
              deallocate(k_mat,ia,ja,ian)
              deallocate(p1,p2,xx,e_dc)
              deallocate(a1,b1,c1)
             ! deallocate(sigma_s)
             if(allocated(v2)) deallocate(v2,vc2)
             if(allocated(vr2)) deallocate(vr2)
             if(allocated(sigf)) deallocate(sigf)
           !  if(uniform_refine) deallocate(sigma)
           ! else
            !nrec = size(rec_loc)
            ! allocate(xx2(size(rec_loc)))
            !!  xx2 = 0
            !  if(csem) xx = xx+e_dc
            !  call pri_dipole(ff,2,a1,b1,c1,layh,n_air,cor_ts,rec_loc,1./patm,sigma,xx2)
            !!  call resistivity_cal(ff,a1,b1,c1,n_air,1,xx,arho)
            !!  res(1:2) = arho(nx/2,1:2) 
            !!  xx = xx+xx2
            !  call post_process(ff,a1,b1,c1,abr(6),abr,xx,xx,rec_loc,arho)
            !!  write(12,*) arho
            !!  write(14,*) xx2  
            !! arho(:,1) =arho(:,1)+real(xx2)
            !! arho(:,2) =arho(:,2)+imag(xx2)
            ! open(71,file='mcsem_model1_0.5hz_inl.dat')
            ! do j = 1,nrec
            !    write(71,'(12(e9.2,2x))') arho(j,1:6),abs(xx2(j)),abs(xx2(j+nrec)),&
            !     abs(xx2(j+2*nrec)),-atan(imag(xx2(j))/real(xx2(j)))*180/pi, &
            !     -atan(imag(xx2(j+nrec))/real(xx2(j+nrec)))*180/pi, &
            !     -atan(imag(xx2(j+2*nrec))/real(xx2(j+2*nrec)))*180/pi
            ! enddo
            ! close(71)
            ! deallocate(xx2)
!                extmp =cmplx(arho(size(rec_loc)/2+1,1),arho(size(rec_loc)/2+1,2),8)
!                hytmp = cmplx(arho(size(rec_loc)/2+1,7),arho(size(rec_loc)/2+1,8),8)
!                res(1) = abs(extmp/hytmp)**2/w/mu
!                res(2) =atan(imag(extmp/hytmp)/real(extmp/hytmp))*180/pi
!                res(3) = arho(size(rec_loc)/2+1,1)
!                res(4) = arho(size(rec_loc)/2+1,2)
            !    print *,res
          endif    
         enddo
   
        deallocate(e_dc)
        deallocate(u01,u02,u11,U12)   
        deallocate(k_mat,ia,ja,ian)
        deallocate(p1,p2,xx)
        deallocate(a1,b1,c1)
        deallocate(sigma,sigf)
    end subroutine 

    subroutine fwd_solver(modelname,a,b,c,nair,level,slv,eps,sigma,freq,ex)
      implicit none
      real*8,allocatable ::a(:),b(:),c(:),sigma(:,:),freq(:),recv(:,:)
     ! real*8,allocatable,optional :: topo(:)
      integer nair, level, nx,ny,nz, nf,nrec,nl,i,j,k,slv
      complex(8),allocatable :: ex(:,:),u1(:),u2(:)
      real*8 eps, ff
      character(len=30) modelname
     
      nx = size(a); ny = size(b); nz = size(c); nf =size(freq)
      nl = nx*(ny+1)*(nz+1)+(nx+1)*ny*(nz+1)+(nx+1)*(ny+1)*nz
       
      allocate(ex(nl*2,nf))
      do i = 1,nf
          ff = freq(i)
          print *,'Current frequenncy(Hz):',ff,i,'of',nf
          ! select excmg or not
          if (fwd_acc.eq.0)then   
             call doublemode_nor(modelname,ff,a,b,c,Nair,sigma,0,slv,eps,u1,u2)
             ex(1:nl,i) =  u1
             ex(nl+1:nl*2,i) =  u2
          else
              call doublemode_acc(modelname,ff,a,b,c,Nair,level,sigma,slv,eps,u1,u2)
              ex(1:nl,i) =  u1
              ex(nl+1:nl*2,i) =  u2
          endif
         deallocate(u1,u2)
      enddo    
        
    end subroutine
    
 end module

!program fdfe_main
!  use global
!
!  implicit none
!        real(8) ::f(26),aa(12)
!    ! real(8) ::ax(35),by(32),cz(37),f(26),aa(12)
!    integer :: N_air
!    integer maxit, i, j,k, m,N_g, NX, NY, NZ,nl,nrec
!    real(8) :: R_tol,tol,ff,skd,extcoef(4),mgs,res(4)
!    ! real(8) ::PaTM(3) = (/1.d12,100.d0,0.5d0/)
!    real*16 c
!    complex(8),allocatable::xx(:),u1(:),u2(:),xx2(:)
!    real(8),allocatable ::arhowe(:,:),a1(:),b1(:),c1(:),cpum(:,:),cpusts(:,:),&
!            ax(:),by(:),cz(:),xdet(:),conv(:,:),xr(:),xi(:),patm(:),arho(:,:),&
!            sigma(:,:),sigma_new(:,:)
!    real(8)  xmin, ymin, ztop,one
!    CHARACTER(len=32) :: arg
!    Character(len=9) :: string
!    real t_sol1,t_sol2
!    real,allocatable ::m_sol(:,:),accr(:,:),cpu_nor(:,:),cpus(:),iters(:)
!    integer,allocatable:: itersts(:,:),abr_in(:),abr(:),abu(:),abr1(:)
!!       character(len=13) CTemp
!    extcoef =(/1.d0,1.2d0,1.5d0,1.85d0/)
!    do i = 1,25
!        f(i) = -3.0+(i-1)*0.25
!    enddo
!    f = 10.d0**f
!    ! do i=1,21
!    !    f(i) = 2.0**(i-1)
!    ! enddo
!  ! c = mptest(3)
!   if(.not.anlcase)then
!     call model_set( ax, by, cz,abr_in,abu, PaTM,xdet,conv)
!  !! layered model for validating accuracy(Chen et al,2016,CJG)
!!    allocate(ax(52),by(52),cz(48),xdet(21))
!!    allocate(patm(3),conv(2,6),abr_in(12),abu(6))
!!    patm = (/1.d9, 5.d2, 1.d2/)
!!    ax(25:28) = 50.d0 
!!    do i = 1,24
!!      ax(28+i) = 50.d0*1.12**i
!!    enddo
!!    ax(1:24) = ax(52:29:-1) 
!!    by  =ax
!!    cz(11:18) = (/50,75,125,250,400,400,400,300/)*1.d0
!!    cz(19:23) = 300.d0; cz(22) =480.d0; cz(24:30) = 300.d0; cz(27)=200.d0
!!    do i = 1,18
!!       cz(30+i) = 300.d0*1.15**i
!!    enddo
!!    do i =1,10
!!       cz(11-i) = 50.d0*1.35**i
!!    enddo
!!    abr_in =(/1,52,1,52,19,23, 1,52,1,52,24,30/)
!!    conv=0
!!    !! layered model from Bai etal(2022,COGSC), anisotropic
!!     conv(1,:) = (/2.d2,1.d3,2.d2,pi/18,pi/6,pi/9/)
!!     conv(2,:) = (/1.d3,2.d2,1.d3,pi/18,pi/9,pi/6/)
!!    abu = (/19,34,19,34,6,26/)
!!    xdet = (/-10:10:1/)*1.d2
!    !! math model with analytical solution
!   else 
!     ! allocate(ax(8),by(8),cz(8),xdet(21))
!      allocate(ax(3),by(3),cz(3),xdet(21))
!      allocate(patm(3),conv(1,6),abr_in(6),abu(6))
!       ax = 1.d0/3; by=ax; cz= ax
!       patm = (/1.d0, 1.d0, 1.d2/)
!       abr_in = (/1,3,1,3,1,3/)
!      conv =0
!      conv(1,1:3) = (/1.d0,1.d0,1.d0/)
!      xdet = (/0:20:1/)/20.d0
!  endif
!
!  select case(model)
!  case(0)
!     n_air = 3
!  case(1)
!      if(topography==0) then
!        N_air = 1
!      else
!        N_air = 12
!      endif
!  case(2,4,16)
!       N_air = 4
!  case(3)
!       N_air = 9
!  case(5)
!       n_air = 5
!  case(6)
!       n_air = 2
!  case(7)
!       n_air = 12
!  case(10)
!       n_air = 256
!  case(11)
!      n_air = 20 
!  case(14)
!      n_air = 30
!  case(15)
!      n_air = 128
!end select
!
! nx = size(ax); ny = size(by); nz = size(cz); nrec =size(xdet)
! nl = nx*(ny+1)*(nz+1)+(nx+1)*ny*(nz+1)+(nx+1)*(ny+1)*nz
!
! maxit = 2000
! r_tol = 1.d-10
! n_g = 6
!! ff =1.d0/1325
!! ff = 1.d0/112.47
! ff  = 1.d-2
!  !   n_air = 16
! allocate(arhowe(25,6))
!
! abr = abr_in
! print *,'physical size of computational domain(x,z)',sum(ax),sum(cz)
!! print *,topography
! !stop
!! abu = 0
!! stop 
!! do j = 1,4
!!   ff = 10.d0**(j-4)
!!   call singlemode_nor(ff,PATM,ax,by,cz,N_air,abr,abu,conv,1,xdet,r_tol,maxit,u1,arho) 
!!  call doublemode_nor(ff,patm,ax,by,cz,N_air,abr,abu,conv,xdet,R_tol,maxit,U1,U2,arho)
!!   deallocate(u1)
!!   if(allocated(arho)) deallocate(arho)
!! enddo
!! stop
! if(.not.exc)then
! do i = 1,n_g-1
!   call mesh_refine(2,ax,by,cz,Nx,Ny,Nz,abr,abu,i,N_air,a1,b1,c1)
!  if(i.eq.n_g-1)then 
!!     call load_sigma(PaTM,ax,by,cz,N_air,abr,conv,sigma)
! !    do j = 1,4
! !       ff = 10.d0**((j-1)/2.d0-2)  
!    !    ff = 2.d0**(j-1)
!  !      ff = 10.d0**(j-4)
!  !     ff = f(j)
!      ! skd = -(log10(ff)+3)*0.4/6+1.4    
!       print *,'f(Hz)',ff
!!       do k = 1,18 
!!          cz(42+k) = 1.d2*skd**k
!!       enddo 
!!        call doublemode_nor(ff,patm,a1,b1,c1,N_air,abr,abu,conv,xdet,R_tol,maxit,U1,u2,arho)
!       call singlemode_nor(ff,PATM,a1,b1,c1,N_air,abr,abu,conv,1,xdet,r_tol,maxit,u1,arho)
! !      call singlemode_acc(ff,patm,a1,b1,c1,N_air,N_g+1,abr,abu,conv,xdet,r_tol,maxit,res) 
!  !    arhowe(j,1:4) = arho(nrec/2+1,:)
!   !    arhowe(j,5) =arho(1,1); arhowe(j,6) = arho(1,2) 
!      ! print *,arho
!  !   stop
!       deallocate(u1)
!      if(allocated(arho)) deallocate(arho)
!  !   enddo
!
!   endif
!   deallocate(a1,b1,c1)
!  enddo
! else
!!   call mesh_refine(2,ax,by,cz,Nx,Ny,Nz,abr,abu,1,N_air,a1,b1,c1)
!!   call load_sigma(PaTM,ax,by,cz,N_air,abr,conv,sigma)
!!   call mesh_coarsen(ax,by,cz,abr,sigma,1,sigma_new,a1,b1,c1,abr1)
!!   n_air = 2 
!   print *,maxit
!   call singlemode_acc(ff,patm,ax,by,cz,N_air,N_g,abr,abu,conv,xdet,r_tol,maxit,res)
! endif 
! !stop
!
!! open(125,file='fd_anislayer_21f.dat',position='append') 
!  open(125,file='fd_dtm1.0_17f.dat')
!  do j = 1,17
!   write(125,1029) arhowe(j,:)
! enddo
!1029 format(6(e12.5,2x))
! close(125)
!
!end program






