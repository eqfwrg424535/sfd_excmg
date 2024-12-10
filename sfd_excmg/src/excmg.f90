module excmg
 use global_parameter
 implicit none 
  !real(8), parameter::pi=3.1592653589793238462643d0
   real(8), parameter::order = 2.0d0

  contains
  subroutine nexcmg(nx,ny,nz,nair,u1,u2,u3,u4)
  !! four layer excmg for aphi3dfe
  implicit none
  complex(8),allocatable :: u1(:),u2(:),u3(:),u4(:),tmp(:),tmp2(:),tmp3(:),&
         tmp4(:)
  integer nx,ny,nz,nair,i,j,k,np1,np2,np3,np4,n1,n2,n3,npm(6)
  !integer,allocatable :: iix(:),iixy(:)
 
   n1 =size(u1); n2 =size(u2); n3 =size(u3)
   np1 = (nx/8+1)*(ny/8+1)*(nz/8+1) 
   np2 = (nx/4+1)*(ny/4+1)*(nz/4+1) 
   np3 = (nx/2+1)*(ny/2+1)*(nz/2+1)
   np4 = (nx+1)*(ny+1)*(nz+1)
   allocate(tmp(np1),tmp2(np2),tmp3(np3),tmp4(np4))
   do i = 1,3
      npm = (/(i-1)*np1+1,i*np1,(i-1)*np2+1,i*np2,(i-1)*np3+1,i*np3/)
      tmp = u1(npm(1):npm(2)); tmp2 =u2(npm(3):npm(4))
      tmp3 = u3(npm(5):npm(6)); !tmp4 = u4((i-1)*np4+1:i*np4)
      tmp4 = extrap_new(np4,nx,ny,nz,tmp,tmp2,tmp3)
      call quart_interp(nx,ny,nz,tmp4)
      u4((i-1)*np4+1:i*np4) = tmp4
   enddo
   deallocate(tmp,tmp2,tmp3,tmp4)

   npm = (/3*np1+1,n1,3*np2+1,n2,3*np3+1,n3/)
   allocate(tmp(n1-3*np1),tmp2(n2-3*np2),tmp3(n3-3*np3),tmp4(size(u4)-3*np4))
   tmp = u1(npm(1):npm(2)); tmp2 =u2(npm(3):npm(4))
   tmp3 = u3(npm(5):npm(6));
   tmp4 =extrap_new(size(u4)-3*np4,nx,ny,nz-nair,tmp,tmp2,tmp3)
   call quart_interp(nx,ny,nz-nair,tmp4)
   u4(3*np4+1:size(u4)) = tmp4 
   deallocate(tmp,tmp2,tmp3,tmp4)
   end subroutine

  

   function extrap_new(n,nx,ny,nz,u1,u2,u3)
   integer n,nx,ny,nz,nair,idx(8,4),idx2(27,2),idx3(125,2)
   complex(8),allocatable :: u1(:),u2(:),u3(:)
   complex(8) extrap_new(n),tmp0(8),tmp1(27),tmp2(125),inv(3)
   real*8 beta(3),l1(27,8),l2(125,27)
   parameter(beta=(/20,-21,1/)/64.d0)
   integer i,j,k,m,p,q
   
   !! two interpolation matrices
   l1 = cal_co8() 
   l2 = cal_co27()
   !! extrapolation cell by cell
   do k = 1,nz/8
      do j = 1,ny/8
         do i = 1,nx/8 
          !! find those node indexes on each grid
            do m = 1,2
               do p =1,2
                 do q = 1,2
                    idx((m-1)*4+(p-1)*2+q,1) =(k+m-2)*(nx/8+1)*(ny/8+1)+&
                                           (j+p-2)*(nx/8+1)+i+q-1
                    idx((m-1)*4+(p-1)*2+q,2) = 2*(k+m-2)*(nx/4+1)*(ny/4+1)+&
                                           2*(j+p-2)*(nx/4+1)+2*(i+q-2)+1
                    idx((m-1)*4+(p-1)*2+q,3) =4*(k+m-2)*(nx/2+1)*(ny/2+1)+&
                                           4*(j+p-2)*(nx/2+1)+4*(i+q-2)+1
                    idx((m-1)*4+(p-1)*2+q,4) =8*(k+m-2)*(nx+1)*(ny+1)+&
                                           8*(j+p-2)*(nx+1)+8*(i+q-2)+1
                 enddo
              enddo
            enddo
        
            do m = 1,5
               do p =1,5
                  do q = 1,5
                    idx3((m-1)*25+(p-1)*5+q,1) =(4*k+m-5)*(nx/2+1)*(ny/2+1)+&
                                           (4*j+p-5)*(nx/2+1)+4*i+q-4
                    idx3((m-1)*25+(p-1)*5+q,2)=(8*(k-1)+2*(m-1))*(nx+1)*(ny+1)+&
                                           (8*(j-1)+2*(p-1))*(nx+1)+8*(i-1)+2*(q-1)+1                                   
                  enddo
                enddo
            enddo
            
            do m = 1,8
               inv = (/u3(idx(m,3)),u2(idx(m,2)),u1(idx(m,1))/)
                tmp0(m) = dot_product(beta,inv) 
            enddo
             tmp1 = matmul(l1,tmp0)
             tmp2 = matmul(l2,tmp1)
             extrap_new(idx3(:,2)) = u3(idx3(:,1))+tmp2
            enddo
          enddo
       enddo
                        
  end function

  subroutine quart_interp(nx,ny,nz,u)
  !! quartic Lagragian interpolation    
   complex(8) ::u(:) 
   integer nx,ny,nz,idx1(125),idx2(729) 
   real*8 l4(729,125)
   integer i,j,k,m,p,q,ik
   
   l4 = cal_co125()
   do k = 1,nz/8
      do j = 1,ny/8
         do i = 1,nx/8
                   ! ik = 0
            do m = 1,5
              do p = 1,5
                do q =1,5
                  idx1((m-1)*25+(p-1)*5+q) =(8*(k-1)+2*(m-1))*(nx+1)*(ny+1)+&
                           (8*(j-1)+2*(p-1))*(nx+1)+8*(i-1)+2*(q-1)+1
                enddo
              enddo
            enddo
            do m = 1,9
              do p = 1,9
                do q =1,9
                  idx2((m-1)*81+(p-1)*9+q) =(8*(k-1)+(m-1))*(nx+1)*(ny+1)+&
                           (8*(j-1)+(p-1))*(nx+1)+8*(i-1)+q
                 enddo
               enddo
            enddo
            u(idx2) = matmul(l4,u(idx1))
          enddo
        enddo
     enddo
         
    end subroutine

   function cal_co125()
    !coor(729,125)——插值矩阵的计算
    implicit none
    real*8::cal_co125(729,125)

    integer i,j,k,l,m,n,row,col
    real*8 coeff1D(9,5)

    coeff1D = coeff4_1D()
    do i=1,9           ! z-direction
        do j=1,9       ! y-direction
            do k=1,9   ! x-direction
                row = (i-1)*81 + (j-1)*9 + k
                do l=1,5         ! z-direction
                    do m=1,5     ! y-direction
                        do n=1,5 ! x-direction
                            col = (l-1)*25 + (m-1)*5 + n
                            cal_co125(row,col) =coeff1D(i,l)*coeff1D(j,m)*coeff1D(k,n)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
  
    end function

    function coeff4_1D()
    !!coor(9,5)——插值矩阵的计算
    implicit none
    real*8::coeff4_1D(9,5)

    integer i
    real*8 x(5), y(9)                           

    forall(i=1:5) x(i) = -1.0d0 +  0.5*(i-1)   
    forall(i=1:9) y(i) = -1.0d0 + 0.25*(i-1)  

    do i=1,9
        coeff4_1D(i,:) = coeff5(y(i))
    enddo
    contains
    
        function coeff5(y)
        real*8::coeff5(5)
        real*8::y

        integer j,k
        coeff5 = 1.0d0
        do j=1,5
            do k=1,5
                if (j /= k) then
                    coeff5(j) = coeff5(j)*(y-x(k))/(x(j)-x(k))
                endif
            enddo
        enddo
        end function

    end function


    function extrap(U1,U2,Nx,Ny,Nz,ND)
    integer::Nx,Ny,Nz,ND
    complex(8)::U1(ND)
    complex(8)::U2(ND)
    complex(8)::extrap(ND)

    integer j,k,dis,Nxy,dis1,dis2,dis3,dis4
    integer,allocatable::IIx(:),index1(:),IIxy(:)

    real*8 c1, c2, c3, c4
    c1 = 2**order
    c2 = 2**(order+1)
    c3 = 2**(order+2)
    c4 = 2**(order+3)

    !coarse mesh nodes*********************************************************************
    Nxy = (Nx/4+1)*(Ny/4+1)
    allocate(index1((Nx/4+1)*(Ny/4+1)*(Nz/4+1)),IIx(Nx/4+1),IIxy(Nxy))
    IIx=(/1:Nx+1:4/)
    do j=1,Ny/4+1
        index1((j-1)*(Nx/4+1)+(/1:Nx/4+1/))=(j-1)*(Nx+1)*4+IIx
    enddo
    IIxy = index1(1:Nxy)
    do k=2,Nz/4+1
        index1((k-1)*Nxy+(/1:Nxy/))=(k-1)*(Nx+1)*(Ny+1)*4+IIxy
    enddo
    extrap(index1)=U2(index1)+(U2(index1)-U1(index1))/c1
    deallocate(index1,IIx,IIxy)
    !middle points at x-dir******************************************************************
    Nxy=(Nx/4)*(Ny/4+1)
    allocate(index1((Nx/4)*(Ny/4+1)*(Nz/4+1)),IIx(Nx/4),IIxy(Nxy))
    IIx=(/1:Nx:4/)+2
    do j=1,Ny/4+1
        index1((j-1)*(Nx/4)+(/1:Nx/4/))=(j-1)*(Nx+1)*4+IIx
    enddo
    IIxy=index1(1:Nxy)
    do k=2,Nz/4+1
        index1((k-1)*Nxy+(/1:Nxy/))=(k-1)*(Nx+1)*(Ny+1)*4+IIxy
    enddo
    extrap(index1)=U2(index1)+(U2(index1-2)-U1(index1-2)+U2(index1+2)-U1(index1+2))/c2
    deallocate(index1,IIx,IIxy)
    !middle points at y-dir******************************************************************
    Nxy=(Ny/4)*(Nx/4+1)
    allocate(index1((Nx/4+1)*(Ny/4)*(Nz/4+1)),IIx(Nx/4+1),IIxy(Nxy))
    IIx=(/1:Nx+1:4/)
    do j=1,Ny/4
        index1((j-1)*(Nx/4+1)+(/1:Nx/4+1/))=(j-0.5)*(Nx+1)*4+IIx
    enddo
    IIxy=index1(1:Nxy)
 do k=2,Nz/4+1
        index1((k-1)*Nxy+(/1:Nxy/))=(k-1)*(Nx+1)*(Ny+1)*4+IIxy
    enddo
    dis=2*(Nx+1)
    extrap(index1)=U2(index1)+(U2(index1-dis)-U1(index1-dis)+U2(index1+dis)-U1(index1+dis))/c2
    deallocate(index1,IIx,IIxy)
    !middle points at z-dir******************************************************************
    Nxy=(Ny/4+1)*(Nx/4+1)
    allocate(index1((Nx/4+1)*(Ny/4+1)*(Nz/4)),IIx(Nx/4+1),IIxy(Nxy))
    IIx=(/1:Nx+1:4/)+2*(Nx+1)*(Ny+1)
    do j=1,Ny/4+1
        index1((j-1)*(Nx/4+1)+(/1:Nx/4+1/))=(j-1)*(Nx+1)*4+IIx
    enddo
    IIxy=index1(1:Nxy)
    do k=2,Nz/4
        index1((k-1)*Nxy+(/1:Nxy/))=(k-1)*(Nx+1)*(Ny+1)*4+IIxy
    enddo
    dis=2*(Nx+1)*(Ny+1)
    extrap(index1)=U2(index1)+(U2(index1-dis)-U1(index1-dis)+U2(index1+dis)-U1(index1+dis))/c2
    deallocate(index1,IIx,IIxy)

    !center points at xoy planes******************************************************************
    Nxy=(Ny/4)*(Nx/4)
    allocate(index1((Nx/4)*(Ny/4)*(Nz/4+1)),IIx(Nx/4),IIxy(Nxy))
    IIx=(/1:Nx:4/)+2+2*(Nx+1)
    do j=1,Ny/4
        index1((j-1)*(Nx/4)+(/1:Nx/4/))=(j-1)*(Nx+1)*4+IIx
    enddo
    IIxy=index1(1:Nxy)
    do k=2,Nz/4+1
        index1((k-1)*Nxy+(/1:Nxy/))=(k-1)*(Nx+1)*(Ny+1)*4+IIxy
    enddo
    dis1=2+2*(Nx+1);    dis2=2*(Nx+1)-2

    extrap(index1) = U2(index1)+(U2(index1-dis1)-U1(index1-dis1)+U2(index1+dis1)-U1(index1+dis1))/c3 &
        & + (U2(index1+dis2)-U1(index1+dis2)+U2(index1-dis2)-U1(index1-dis2))/c3
    deallocate(index1,IIx,IIxy)
    !center points at xoz planes******************************************************************
    Nxy=(Ny/4+1)*(Nx/4)
    allocate(index1((Nx/4)*(Ny/4+1)*(Nz/4)),IIx(Nx/4),IIxy(Nxy))
    IIx=(/1:Nx:4/)+2+2*(Nx+1)*(Ny+1)
    do j=1,Ny/4+1
        index1((j-1)*(Nx/4)+(/1:Nx/4/))=(j-1)*(Nx+1)*4+IIx
    enddo
    IIxy=index1(1:Nxy)
    do k=2,Nz/4
        index1((k-1)*Nxy+(/1:Nxy/))=(k-1)*(Nx+1)*(Ny+1)*4+IIxy
    enddo
    dis1=2+2*(Nx+1)*(Ny+1);    dis2=2*(Nx+1)*(Ny+1)-2
    extrap(index1) = U2(index1)+(U2(index1-dis1)-U1(index1-dis1)+U2(index1+dis1)-U1(index1+dis1))/c3 &
        & + (U2(index1+dis2)-U1(index1+dis2)+U2(index1-dis2)-U1(index1-dis2))/c3
    deallocate(index1,IIx,IIxy)
    !center points at xoy planes******************************************************************
    Nxy=(Nx/4+1)*(Ny/4)
    allocate(index1((Nx/4+1)*(Ny/4)*(Nz/4)),IIx(Nx/4+1),IIxy(Nxy))
    IIx=(/1:Nx+1:4/)+2*(Nx+1)+2*(Nx+1)*(Ny+1)
    do j=1,Ny/4
        index1((j-1)*(Nx/4+1)+(/1:Nx/4+1/))=(j-1)*(Nx+1)*4+IIx
    enddo
    IIxy=index1(1:Nxy)
    do k=2,Nz/4
        index1((k-1)*Nxy+(/1:Nxy/))=(k-1)*(Nx+1)*(Ny+1)*4+IIxy
    enddo
    dis1=2*(Nx+1)*(Ny+1)+2*(Nx+1);    dis2=2*(Nx+1)*(Ny+1)-2*(Nx+1)
    extrap(index1) = U2(index1)+(U2(index1-dis1)-U1(index1-dis1)+U2(index1+dis1)-U1(index1+dis1))/c3 &
        & + (U2(index1+dis2)-U1(index1+dis2)+U2(index1-dis2)-U1(index1-dis2))/c3
    deallocate(index1,IIx,IIxy)
    !element center points****************************************************************
    Nxy=(Nx/4)*(Ny/4)
    allocate(index1((Nx/4)*(Ny/4)*(Nz/4)),IIx(Nx/4),IIxy(Nxy))
    IIx=(/1:Nx:4/)+2+2*(Nx+1)+2*(Nx+1)*(Ny+1)
    do j=1,Ny/4
        index1((j-1)*(Nx/4)+(/1:Nx/4/))=(j-1)*(Nx+1)*4+IIx
    enddo
    IIxy=index1(1:Nxy)
    do k=2,Nz/4
        index1((k-1)*Nxy+(/1:Nxy/))=(k-1)*(Nx+1)*(Ny+1)*4+IIxy
    enddo
    dis1=2*(Nx+1)*(Ny+1)+2*(Nx+1)+2;    dis2=2*(Nx+1)*(Ny+1)-2*(Nx+1)+2
    dis3=2*(Nx+1)*(Ny+1)-2*(Nx+1)-2;    dis4=2*(Nx+1)*(Ny+1)+2*(Nx+1)-2
    extrap(index1) = U2(index1)+(U2(index1-dis1)-U1(index1-dis1)+U2(index1+dis1)-U1(index1+dis1))/c4 &
        & + (U2(index1+dis2)-U1(index1+dis2)+U2(index1-dis2)-U1(index1-dis2))/c4 &
        & + (U2(index1+dis3)-U1(index1+dis3)+U2(index1-dis3)-U1(index1-dis3))/c4 &
        & + (U2(index1+dis4)-U1(index1+dis4)+U2(index1-dis4)-U1(index1-dis4))/c4
    deallocate(index1,IIx,IIxy)

    end function

    subroutine quad_interp27(U3,Nx,Ny,Nz,ND)
    integer::Nx,Ny,Nz,ND
    complex(8)::U3(ND)

    integer i,j,k,jN,N1,N2,Nx4,Ny4,Nz4
    integer D27(27),D125(125),temp5(5),temp25(25),I27(27),I125(125)
    real*8 L0(125,27)

    L0 = cal_co27()

    N1=(Nx+1);      N2=(Nx+1)*(Ny+1)
    Nx4=Nx/4;       Ny4=Ny/4;       Nz4=Nz/4
    D27=(/1, 3, 5, 2*N1+1, 2*N1+3, 2*N1+5, 4*N1+1, 4*N1+3, 4*N1+5, &
        &(/1, 3, 5, 2*N1+1, 2*N1+3, 2*N1+5, 4*N1+1, 4*N1+3, 4*N1+5/)+2*N2, &
        &(/1, 3, 5, 2*N1+1, 2*N1+3, 2*N1+5, 4*N1+1, 4*N1+3, 4*N1+5/)+4*N2 /)
    temp5=(/1:5/)
    temp25=(/temp5,  temp5+N1,  temp5+2*N1,  temp5+3*N1,  temp5+4*N1/)
    D125 =(/temp25, temp25+N2, temp25+2*N2, temp25+3*N2, temp25+4*N2/)


    do k=1,Nz4
        do j=1,Ny4
            do i=1,Nx4
                jN=4*(k-1)*(Nx+1)*(Ny+1)+4*(j-1)*(Nx+1)+4*(i-1)
                I27=D27+jN
                I125=D125+jN
                U3(I125)=matmul(L0,U3(I27))
            enddo
        enddo
    enddo

    end subroutine

    function cal_co27()

    double precision::cal_co27(125,27)

    integer i, j, k, ijk,ind
    double precision coor(125,3)

    do k=1,5
        do j=1,5
            do i=1,5
                ijk=(k-1)*25+(j-1)*5+i
                coor(ijk,:)=(/(i-3)*0.5,(j-3)*0.5,(k-3)*0.5/)
            enddo
        enddo
    enddo
    do ijk=1,125
        do k=1,3
            do j=1,3
                 do i=1,3
                    ind = (k-1)*9+(j-1)*3+i
                    cal_co27(ijk,ind) = l2(i,coor(ijk,1))*l2(j,coor(ijk,2))*l2(k,coor(ijk,3))
                enddo
            enddo
        enddo
    enddo

    end function

    function extint_v(nx,ny,nz,nl,u0,u1,opt,codsch)
     !! prolongation for VFEM  
      integer :: nx,ny,nz,nl,opt,codsch
      complex(8) :: u0(*),u1(*),extint_v(nl)
      complex(8),allocatable:: utmp(:),w0(:),w1(:)
      integer i,j,k,nl_new,n1,nl2,nl3,xnl,ynl,xnl2,xnl3,ynl2,ynl3
      real(8) co1,co2,co3


      n1 = ((nx/2+1)*ny/2+(ny/2+1)*nx/2)*(nz/2+1)+(nx/2+1)*(ny/2+1)*nz/2
     !n1= ((nx/3+1)*ny/3+(ny/3+1)*nx/3)*(nz/3+1)+(nx/3+1)*(ny/3+1)*nz/3 
      nl2 = nx/4*(ny/4+1)*(nz/4+1)+(nx/4+1)*ny/4*(nz/4+1)+(nx/4+1)*(ny/4+1)*nz/4
      xnl = nx*(ny+1)*(nz+1);ynl =ny*(nx+1)*(nz+1)
      xnl2 = (ny/2+1)*nx/2*(nz/2+1); ynl2 = (nx/2+1)*ny/2*(nz/2+1)
      xnl3 = (ny/4+1)*nx/4*(nz/4+1); ynl3 = (nx/4+1)*ny/4*(nz/4+1)
      co1 = 1.d0+1.d0/2.d0**order
      co2 = 1.d0/2.d0**order
      allocate(utmp(n1))
      utmp = 0
   
    !  w0 = cmplx(real(u0(1:nl2))/2,imag(u0(1:nl2))/4) 
     !  w0 = u0(1:nl2)
     !  w1 = u1(1:n1)
     ! print *,w0(1),w1(1)
     ! stop
     ! w0(1:xnl3) = u0(1:xnl3)*nx/4;w1(1:xnl2) = u1(1:xnl2)*nx/2
     !  w0(xnl3+1:xnl3+ynl3) = u0(xnl3+1:xnl3+ynl3)*ny/4
     !   w1(xnl2+1:xnl2+ynl2) = u1(xnl2+1:xnl2+ynl2)*ny/2
     !   w0(xnl3+ynl3+1:nl2) = u0(xnl3+ynl3+1:nl2)*nz/4
     !   w1(xnl2+ynl2+1:n1) = u1(xnl2+ynl2+1:n1)*nz/2
    if(opt.eq.1)then
     ! call lin_interp_e(Nx/2,Ny/2,Nz/2,n1,w0,utmp,codsch)
       call quad_interp_e(nx/2,ny/2,nz/2,n1,u0,utmp,codsch)
    elseif(opt.eq.2)then 
     ! call cubic_interp(nx/2,ny/2,nz/2,n1,u0,utmp,codsch)
      call lin_interp_e(Nx/2,Ny/2,Nz/2,n1,u0,utmp,codsch) 
    elseif(opt.eq.0)then
      call lin_interp_e(Nx/2,Ny/2,Nz/2,n1,u0,utmp,codsch)
    else
      call qc_interp(nx/3,ny/3,nz/3,n1,u0,utmp,codsch)
    endif
       do i = 1,n1
          utmp(i) = co1*u1(i)-co2*utmp(i)
         ! if(abs(utmp(i)).le.1.d-50.and.utmp(i).ne.0) print *,i,w1(i),utmp(i) 
       enddo
     if(opt.eq.1)then
       call quad_interp_e(nx,ny,nz,nl,utmp,extint_v,codsch)
    !  call cubic_interp(nx,ny,nz,nl,utmp,extint_v,codsch)
     elseif(opt.eq.2)then
       call cubic_interp(nx,ny,nz,nl,utmp,extint_v,codsch)
     elseif(opt.eq.0)then
     !  call lin_interp_e(Nx,Ny,Nz,nl,utmp,extint_v,codsch)
      call quad_interp_e(nx,ny,nz,nl,utmp,extint_v,codsch)
     else
       call qc_interp(nx,ny,nz,nl,utmp,extint_v,codsch)
     endif
    ! print *,'prolong completed'
     !   extint_v(1:xnl) = extint_v(1:xnl)/nx
     !   extint_v(xnl+1:xnl+ynl) = extint_v(xnl+1:xnl+ynl)/ny
     !   extint_v(xnl+ynl+1:nl) = extint_v(xnl+ynl+1:nl)/nz
    ! extint_v = cmplx(real(extint_v)/2,imag(extint_v)/4)
    ! extint_v =  extint_v/2
     do i = 1,nl
         if(isnan(abs(extint_v(i)))) print *,i
     enddo
      deallocate(utmp) 
    end function

    subroutine qc_interp(nx,ny,nz,nl,u1,u2,cod_sch)
     !! mixed quadratic-cubic interpolation for trichotomous refinement
        implicit none
       integer nx,ny,nz,nl,cod_sch
       complex(8)::u1(*),u2(nl)
       integer i,j,k,ind,ind2,nd,m1,n2,n3,n4,i900(900),i48(48),d16(16),&
            d100(100),d900(900,3),d10(10),d48(48,3),xnl1,ynl1,&
            xnl2,ynl2,ind3,ind4  !d900 and d48 are excursion matrixes
       real(8) L1(900,48)

       l1 = cal_co48()

        nd =(nx+1)*(ny+1)*(nz+1)
       if(cod_sch.eq.1)then
         m1 = (nx/3+1)*ny/3+(ny/3+1)*nx/3+(nx/3+1)*(ny/3+1)
         n2 = (nx/3+1)*ny/3+(ny/3+1)*nx/3
         n3 = (nx+1)*ny+(ny+1)*nx+(nx+1)*(ny+1)
         n4 = (nx+1)*ny+(ny+1)*nx
       else
         m1 = (ny/3+1)*nx/3*(nz/3+1) !xnl_coarse
         n2 = (nx/3+1)*ny/3*(nz/3+1) !ynl_coarse
         n3 = (ny+1)*nx*(nz+1)  !xnl_fine
         n4 = (nx+1)*(nz+1)*ny  !ynl_fine
       endif
      !   print *,m1,n2,n3,n4
       if(cod_sch.eq.1)then
        do i = 1,4
           d16((i-1)*4+1:i*4) = (i-1)*m1+(/0,nx*2/3+1,4*nx/3+2,2*nx+3/)
        enddo
        do i = 1,3
           d48((i-1)*16+1:i*16,1) = d16+i-1
        enddo
        do i = 1,4
           d16((i-1)*4+1:i*4) = (i-1)*m1+(/0,1,2,3/)+nx/3
        enddo
        do i = 1,3
           d48((i-1)*16+1:i*16,2) = d16+(i-1)*(nx*2/3+1)
        enddo
        do i = 1,4
           d16((i-1)*4+1:i*4) = (i-1)*(nx/3+1)+(/0,1,2,3/)
        enddo
        do i = 1,3
           d48((i-1)*16+1:i*16,3) = d16+(i-1)*m1+n2
        enddo

        do i = 1,10
           d10(i) = (i-1)*(2*nx+1)
        enddo
        do i = 1,10
           d100((i-1)*10+1:i*10) = (i-1)*n3 +d10
        enddo
        do i = 1,9
           d900((i-1)*100+1:i*100,1) = d100+i-1
        enddo
        do i = 1,10
           d10(i) = i-1+nx
        enddo
        do i =1,10
           d100((i-1)*10+1:i*10) = (i-1)*n3 +d10
        enddo
        do i =1,9
           d900((i-1)*100+1:i*100,2) = d100+(i-1)*(2*nx+1)
        enddo
        d10 = (/0:9:1/)
        do i =1,10
           d100((i-1)*10+1:i*10) = (i-1)*(nx+1) +d10
        enddo
        do i =1,9
           d900((i-1)*100+1:i*100,3) = d100+(i-1)*n3+n4
        enddo
      else
        do i = 1,4
           d16((i-1)*4+1:i*4) = (i-1)*(ny/3+1)*nx/3+(/0,nx/3,2*nx/3,nx/)
        enddo
        do i = 1,3
           d48((i-1)*16+1:i*16,1) = d16+i-1
        enddo
        do i = 1,4
           d16((i-1)*4+1:i*4) = m1+(i-1)*(nx/3+1)*ny/3+(/0,1,2,3/)
        enddo
        do i = 1,3
           d48((i-1)*16+1:i*16,2) = d16+(i-1)*(nx/3+1)
        enddo
        do i = 1,4
           d16((i-1)*4+1:i*4) = m1+n2+(i-1)*(nx/3+1)+(/0,1,2,3/)
        enddo
        do i = 1,3
           d48((i-1)*16+1:i*16,3) = d16+(i-1)*(nx/3+1)*(ny/3+1)
        enddo

        do i = 1,10
           d10(i) = (i-1)*nx
        enddo
        do i = 1,10
           d100((i-1)*10+1:i*10) = (i-1)*(ny+1)*nx +d10
        enddo
        do i = 1,9
           d900((i-1)*100+1:i*100,1) = d100+i-1
        enddo
        do i = 1,10
           d10(i) = i-1+n3
        enddo
        do i =1,10
           d100((i-1)*10+1:i*10) = (i-1)*(nx+1)*ny +d10
        enddo
        do i =1,9
           d900((i-1)*100+1:i*100,2) = d100+(i-1)*(nx+1)
        enddo
        d10 = (/0:9:1/)
        do i =1,10
           d100((i-1)*10+1:i*10) = (i-1)*(nx+1) +d10+n3+n4
        enddo
        do i =1,9
           d900((i-1)*100+1:i*100,3) = d100+(i-1)*(nx+1)*(ny+1)
        enddo
      endif

       do k = 1,nz/9
          do j = 1,ny/9
            do i = 1,nx/9
              !! ex
              if(cod_sch.eq.1)then
               ind = (k-1)*3*m1+(j-1)*(2*nx+3)+3*i-2
               ind2 = (k-1)*9*n3+(j-1)*(18*nx+9)+9*i-8
              else
                ind = (k-1)*3*nx/3*(ny/3+1)+(j-1)*nx+3*i-2
               ind2 = (k-1)*9*nx*(ny+1)+(j-1)*9*nx+9*i-8
              endif
               i48 = d48(:,1)+ind
               i900 = d900(:,1)+ind2
               u2(i900) = matmul(l1,u1(i48))
              !! ey
              if(cod_sch.eq.2)then
                 ind = (k-1)*3*ny/3*(nx/3+1)+(j-1)*(nx+3)+3*i-2
                 ind2 = (k-1)*9*ny*(nx+1)+(j-1)*(9*nx+9)+9*i-8
               endif
               i48 = d48(:,2)+ind
               i900 = d900(:,2)+ind2
               u2(i900) = matmul(l1,u1(i48))
              !! ez
              if(cod_sch.eq.1)then
                ind3 =  (k-1)*3*m1+(j-1)*(nx+3)+3*i-2
               ind4 =  (k-1)*9*n3+(j-1)*(9*nx+9)+9*i-8
              else
               ind3 =  (k-1)*3*(nx/3+1)*(ny/3+1)+(j-1)*(nx+3)+3*i-2
               ind4 =  (k-1)*9*(nx+1)*(ny+1)+(j-1)*(9*nx+9)+9*i-8
             endif
               i48 = d48(:,3)+ind3
               i900 = d900(:,3)+ind4
               u2(i900) = matmul(l1,u1(i48))
            enddo
          enddo
        enddo

      contains 
         function cal_co48()
          integer i,j,k,ind,ijk
          real(8) cx1(4),cx2(3),coor(900,3), cal_co48(900,48)

         cal_co48 = 0
         cx1 =(/-1.d0,-1.d0/3,1.d0/3,1.d0/)
         cx2 =(/-2.d0/3,0.d0,2.d0/3/)
           do k = 1,9
               do j= 1,10
                   do i = 1,10
                       ijk = (k-1)*100+(j-1)*10+i
                       coor(ijk,1) = -1.d0+(i-1)*2.d0/9
                       coor(ijk,2) = -1.d0+(j-1)*2.d0/9
                       coor(ijk,3)  = -8.d0/9+(k-1)*2.d0/9
                   enddo
                enddo
           enddo

         do ijk = 1,900
            do k  = 1,3
               do j = 1,4
                  do i = 1,4
                     ind = (k-1)*16+(j-1)*4+i
                     cal_co48(ijk,ind)=lbf(3,cx1,i,coor(ijk,1))*lbf(3,cx1,j,coor(ijk,2))*&
                                  lbf(2,cx2,k,coor(ijk,3))
                   enddo
                enddo
            enddo
        enddo

         end function
     end subroutine


    function extint_n(nx,ny,nz,nl,np,u0,u1,opt,cod_sch)
    !! prolongation for NEFEM
      implicit none
      integer :: nx,ny,nz,nl,np,opt,cod_sch
      complex(8) :: u0(*),u1(*),extint_n(nl+np)
      complex(8),allocatable:: utmp(:),w01(:),w11(:),w02(:),w12(:)
      integer i,j,k,nl_new,n1,n2,nl1,nl2,nl3,np1,np2,nxy1,nxy2
      integer,allocatable :: iiy(:),iixy(:),index1(:) 
     real(8) co1,co2,co3

      nl1 = ((nx/2+1)*ny/2+(ny/2+1)*nx/2)*(nz/2+1)+(nx/2+1)*(ny/2+1)*nz/2
      np1 = (nx/2+1)*(ny/2+1)*(nz/2+1)
      nl2 = nx/4*(ny/4+1)*(nz/4+1)+(nx/4+1)*ny/4*(nz/4+1)+(nx/4+1)*(ny/4+1)*nz/4
      np2 = (nx/4+1)*(ny/4+1)*(nz/4+1) 
      nxy1 = (nx/2+1)*(ny/2+1); nxy2=(nx/4+1)*(ny/4+1)
      n1= nl1+np1; n2 = nl2+np2
     co1 = 1.d0+1.d0/2.d0**order
      co2 = 1.d0/2.d0**order

      allocate(utmp(n1),w01(nl2),w11(nl1),w02(np2),w12(np1))
      utmp = 0
      w01 = u0(1:nl2); w02 = u0(nl2+1:n2)
      w11 = u1(1:nl1); w12 = u1(nl1+1:n1)

     !! interpolation-1(only linear on coarsest grid)
       if(opt.eq.1)then
         call lin_interp_e(Nx/2,Ny/2,Nz/2,nl1,w01,utmp(1:nl1),cod_sch)
         call lin_interp_n(nx/2,ny/2,nz/2,np1,w02,utmp(nl1+1:n1))
       else
         call quad_interp_e(nx/2,ny/2,nz/2,nl1,w01,utmp(1:nl1),cod_sch)
         allocate(index1(np2))
         allocate(IIy(Nx/4+1),IIxy(Nxy2))
         ! index1=0;index2=0;
         IIy=(/1:Nx/2+1:2/) !! different from matlab
         do i=1,Ny/4+1
          index1((i-1)*(Nx/4+1)+(/1:Nx/4+1/))=(i-1)*(Nx/2+1)*2+IIy
         enddo
         IIxy=index1(1:Nxy2)
         do k=2,Nz/4+1
           index1((k-1)*Nxy2+(/1:Nxy2/))=(k-1)*Nxy1*2+IIxy
         enddo
         index1 = index1+nl1
         deallocate(IIy,IIxy)
         utmp(index1) =w02 
         deallocate(index1)
         call quad_interp27(utmp(nl1+1:n1),Nx/2,Ny/2,Nz/2,np1)
       endif
     !! global extrapolation
        do i = 1,nl1
          utmp(i) = co1*w11(i)-co2*utmp(i)
        enddo
        do i = 1,np1
           utmp(nl1+i) =co1*w12(i)-co2*utmp(nl1+i)
           if(isnan(real(utmp(nl1+i))).or.isnan(imag(utmp(nl1+i))))then
                print *,i,w12(i),utmp(nl1+i)
                stop
           endif
        enddo
      !!interpolation-2
       call quad_interp_e(nx,ny,nz,nl1,utmp(1:nl1),extint_n(1:nl),cod_sch)
         allocate(index1(np1))
         allocate(IIy(Nx/2+1),IIxy(Nxy1))
         ! index1=0;index2=0;
         IIy=(/1:Nx+1:2/) !! different from matlab
         do i=1,Ny/2+1
          index1((i-1)*(Nx/2+1)+(/1:Nx/2+1/))=(i-1)*(Nx+1)*2+IIy
         enddo
         IIxy=index1(1:Nxy1)
         do k=2,Nz/2+1
           index1((k-1)*Nxy1+(/1:Nxy1/))=(k-1)*(Nx+1)*(Ny+1)*2+IIxy
         enddo
         index1 = index1+nl
         deallocate(IIy,IIxy)
         extint_n(index1) = utmp(nl1+1:n1)
         deallocate(index1)
     
       call quad_interp27(extint_n(nl+1:nl+np),Nx,Ny,Nz,np)
        
       deallocate(utmp,w01,w02,w11,w12)
    end function

    subroutine quad_interp_e(nx,ny,nz,nl,u1,u2,opt)
        implicit none
        integer,intent(in) :: nx,ny,nz,nl,opt
        complex(8):: u1(*),u2(nl)
        integer i,j,k,ind,ind2,nd,N1,N2,N3,N4,n5,n6,ind3,ind4, d9(9),d27(27),&
                d5(5),d25(25),d100(100),i27(27),g27(27),i100(100),xnl1,ynl1,&
                xnl2,ynl2,i18(18)
         real(8) L1(100,27),L2(100,27),l3(100,18)

        ! L1 = cal_co100(1)
        ! L2 = cal_co100(2)
        l3 = cal_co18()
        nd  = (nx+1)*(ny+1)*(nz+1)
        xnl1 = nx/2*(ny/2+1)*(nz/2+1); ynl1 = (nx/2+1)*ny/2*(nz/2+1);
        xnl2 = nx*(ny+1)*(nz+1); ynl2 = (nx+1)*ny*(nz+1)
        !allocate(utmp1(nd),utmp2(nd),utmp3(nd))
        !utmp1 = 0; utmp2= 0; utmp3 =0
       if(opt.eq.1)then
          N1 = (nx/2+1)*ny/2+(ny/2+1)*nx/2+(nx/2+1)*(ny/2+1)
         n2 = (nx/2+1)*ny/2+(ny/2+1)*nx/2
         n3 = (nx+1)*ny+(ny+1)*nx+(nx+1)*(ny+1)
         n4 = (nx+1)*ny+(ny+1)*nx
       else
         n1 = (ny/2+1)*nx/2; n2 = nx*(ny+1)
         n3 = (nx/2+1)*ny/2; n4 = ny*(nx+1)
         n5 = (nx/2+1)*(ny/2+1); n6 = (nx+1)*(ny+1)
       endif
         !! x-edges
       if(opt.eq.1)then
         d9 = (/0,nx+1, 2*nx+2, n1,n1+nx+1,n1+2*nx+2, 2*n1, &
                 2*n1+nx+1, 2*n1+2*nx+2/)
       else
         d9 = (/0,nx/2,nx, n1,n1+nx/2, n1+nx, 2*n1,2*n1+nx/2,2*n1+nx/)
       endif
       do i= 1,5
         if(opt.eq.1)then
            d25((i-1)*5+1:5*i) = (i-1)*n3+(/0,2*nx+1,4*nx+2,6*nx+3,8*nx+4/)
         else
            d25((i-1)*5+1:5*i) = (i-1)*n2+(/0,nx,2*nx,3*nx,4*nx/)
         endif
        enddo

        do i= 1,4
           d100((i-1)*25+1:25*i) = d25+i-1
        enddo
        do k = 1,Nz/4
            do j = 1,Ny/4
                do i = 1,Nx/4
                    if(opt.eq.1)then
                      ind = (k-1)*4*n3+(j-1)*(8*nx+4)+4*i-3
                       ind2 = (k-1)*2*n1+(j-1)*(2*nx+2)+2*i-1
                    else
                      ind = (k-1)*4*n2+(j-1)*(4*nx)+4*i-3
                       ind2 = (k-1)*2*n1+(j-1)*nx+2*i-1
                    endif
                    i100 = d100+ind
                    i18 = ind2+(/d9,d9+1/)
                    u2(I100) = matmul(L3,u1(i18))
!                    if(i.eq.1) then
!                       i27 = ind2+(/d9,d9+1,d9+2/)
!                        u2(I100) = matmul(L1,u1(i27))
!                     elseif(i.eq.nx/4)then
!                        i27 = ind2+(/d9-1,d9,d9+1/)
!                        u2(I100) = matmul(L2,u1(i27))
!                     else
!                        i27 = ind2+(/d9,d9+1,d9+2/)
!                        g27 = ind2+(/d9-1,d9,d9+1/)
!                        u2(i100) =(matmul(L1,u1(i27))+matmul(L2,u1(g27)))/2
!                    endif
                enddo
             enddo
        enddo
      !! y-edges
        if(opt.eq.1)then
          d9 = (/nx/2,nx/2+1,nx/2+2,n1+nx/2,n1+nx/2+1,n1+nx/2+2, 2*n1+nx/2,&
                  2*n1+nx/2+1,2*n1+nx/2+2/)
        else
          d9 = xnl1+(/0,1,2,n3,n3+1,n3+2, 2*n3,2*n3+1,2*n3+2/)
        endif
        do i = 1,5
          if(opt.eq.1)then
            d25((i-1)*5+1:5*i) = (i-1)*n3+(/nx,nx+1,nx+2,nx+3,nx+4/)
          else
            d25((i-1)*5+1:5*i) = xnl2+(i-1)*n4+(/0:4/)
          endif
        enddo
        do i =1,4
          if(opt.eq.1)then
            d100((i-1)*25+1:25*i) = d25+(i-1)*(2*nx+1) !(2*nx+1)
          else
           d100((i-1)*25+1:25*i) = d25+(i-1)*(nx+1)
          endif
        enddo
        do k  = 1,nz/4
            do j = 1,ny/4
                do i = 1,nx/4
                 if(opt.eq.1)then
                     ind = (k-1)*4*n3+(j-1)*(8*nx+4)+4*i-3
                     ind2 = (k-1)*2*n1+(j-1)*(2*nx+2)+2*i-1
                     i100 =d100+ind
                      i18 = ind2+(/d9,d9+nx+1/)
                      u2(i100) = matmul(l3,u1(i18))
                  else
                     ind = (k-1)*4*n4+(j-1)*(4*nx+4)+4*i-3
                    ind2 = (k-1)*2*n3+(j-1)*(nx+2)+2*i-1
                     i100 = d100+ind 
                      i18 = ind2+(/d9,d9+nx/2+1/)
                      u2(i100) = matmul(l3,u1(i18))
                   endif                  
!                   i100 = d100+ind
!                   if(j.eq.1) then
!                     if(opt.eq.1)then
!                      i27 = ind2+(/d9,d9+nx+1,d9+2*nx+2/)
!                     else
!                      i27 = ind2+(/d9,d9+nx/2+1,d9+nx+2/)
!                     endif
!                      u2(I100) = matmul(L1,u1(i27))
!                   elseif(j.eq.ny/4)then
!                     if(opt.eq.1)then
!                      i27 =ind2+(/d9-nx-1,d9,d9+nx+1/)
!                     else
!                      i27 =ind2+(/d9-nx/2-1,d9,d9+nx/2+1/)
!                     endif
!                      u2(I100) = matmul(L2,u1(i27))
!                   else
!                     if(opt.eq.1)then
!                      i27 = ind2+(/d9,d9+nx+1,d9+2*nx+2/)
!                      g27 = ind2+(/d9-nx-1,d9,d9+nx+1/)
!                     else 
!                      i27 = ind2+(/d9,d9+nx/2+1,d9+nx+2/)
!                      g27 = ind2+(/d9-nx/2-1,d9,d9+nx/2+1/)
!                     endif
!                      u2(i100) =(matmul(L1,u1(i27))+matmul(L2,u1(g27)))/2
!                   endif
                enddo
            enddo
        enddo

       !! z- edges
       if(opt.eq.1)then
         d9 = (/0,1,2, nx/2+1,nx/2+2,nx/2+3, nx+2,nx+3,nx+4/)
       else
         d9 = xnl1+ynl1+(/0,1,2, nx/2+1,nx/2+2,nx/2+3, nx+2,nx+3,nx+4/)
       endif
        do i = 1,5
           d25((i-1)*5+1:5*i) = (i-1)*(nx+1)+(/0:4/)
        enddo
         do i = 1,4
          if(opt.eq.1)then
           d100((i-1)*25+1:25*i) = d25+(i-1)*n3
         else
           d100((i-1)*25+1:25*i) = d25+(i-1)*n6+xnl2+ynl2
         endif
        enddo
        do k  = 1,nz/4
            do j = 1,ny/4
                do i = 1,nx/4
                   if(opt.eq.1)then
                     ind = (k-1)*4*n3+(j-1)*(4*nx+4)+4*i-3+n4
                     ind2 = (k-1)*2*n1+(j-1)*(nx+2)+2*i-1+n2
                     i100 = d100+ind  
                     i18 = ind2+(/d9,d9+n1/)
                      u2(I100) = matmul(L3,u1(i18))
                   else
                     ind = (k-1)*4*n6+(j-1)*(4*nx+4)+4*i-3
                     ind2 = (k-1)*2*n5+(j-1)*(nx+2)+2*i-1
                      i100 = d100+ind
                      i18 = ind2+(/d9,d9+n5/)
                      u2(I100) = matmul(L3,u1(i18))
                  endif
!                   if(k.eq.1) then
!                      if(opt.eq.1)then
!                        i27 = ind2+(/d9,d9+n1,d9+2*n1/)
!                      else
!                        i27 = ind2+(/d9,d9+n5,d9+2*n5/)
!                      endif
!                        u2(I100) = matmul(L1,u1(i27))
!                   elseif(k.eq.nz/4)then
!                      if(opt.eq.1)then
!                        i27 = ind2+(/d9-n1,d9,d9+n1/)
!                      else
!                        i27 = ind2+(/d9-n5,d9,d9+n5/)
!                      endif
!                       u2(I100) = matmul(L2,u1(i27))
!                    else
!                      if(opt.eq.1)then
!                       i27 = ind2+(/d9,d9+n1,d9+2*n1/)
!                       g27 = ind2+(/d9-n1,d9,d9+n1/)
!                      else
!                       i27 = ind2+(/d9,d9+n5,d9+2*n5/)
!                       g27 = ind2+(/d9-n5,d9,d9+n5/)
!                      endif
!                       u2(i100) =(matmul(L1,u1(i27))+matmul(L2,u1(g27)))/2
!                    endif
               enddo
            enddo
        enddo

         contains
              function cal_co100(tp)
                  integer  tp,i,j,k,ind,ijk
                   real(8) coor(100,3), cal_co100(100,27)

                   do k = 1,4
                       do j= 1,5
                           do i = 1,5
                               ijk = (k-1)*25+(j-1)*5+i
                               coor(ijk,1) = (i-3)*0.5
                               coor(ijk,2) = (j-3)*0.5
                               if(tp.eq.1)then
                                  coor(ijk,3)  = (k-1)*0.5-1.25
                               else
                                  coor(ijk,3)  = (k-1)*0.5-0.25
                               endif
                           enddo
                        enddo
                   enddo

                 do ijk = 1,100
                    do k  = 1,3
                       do j = 1,3
                          do i = 1,3
                             ind = (k-1)*9+(j-1)*3+i
                             cal_co100(ijk,ind) =l2(i,coor(ijk,1))*l2(j,coor(ijk,2))*&
                                          l2(k,coor(ijk,3))
                           enddo
                        enddo
                    enddo
                 enddo

               end function

               function cal_co18()
                   integer  i,j,k,ind,ijk
                   real(8) coor(100,3), cal_co18(100,18)

                   do k = 1,4
                       do j= 1,5
                           do i = 1,5
                               ijk = (k-1)*25+(j-1)*5+i
                               coor(ijk,1) = (i-3)*0.5
                               coor(ijk,2) = (j-3)*0.5
                               coor(ijk,3)  = (k-1)*0.5-0.75

                           enddo
                        enddo
                   enddo

                 do ijk = 1,100
                    do k  = 1,2
                       do j = 1,3
                          do i = 1,3
                             ind = (k-1)*9+(j-1)*3+i
                             cal_co18(ijk,ind)=l2(i,coor(ijk,1))*l2(j,coor(ijk,2))*&
                                          l1(k,coor(ijk,3))
                           enddo
                        enddo
                      enddo
                 enddo
              end function
     end subroutine     

     subroutine cubic_interp(nx,ny,nz,nl,u1,u2,cod_sch)
       implicit none
       integer nx,ny,nz,nl,cod_sch       
       complex(8)::u1(*),u2(nl)
        integer i,j,k,m,ind,ind2,nd,n1,n2,n3,n4,m1,m2,m3,i25(25),i100(100),d9(9),&
            d81(81),d648(648,3),d25(25),d100(100,3),i648(648),xnl1,ynl1,&
            xnl2,ynl2,ind3,ind4
       real(8) L1(648,100)

        ! print *,nx,ny
         l1 = cal_co648()
        ! print *,nx,ny
         nd =(nx+1)*(ny+1)*(nz+1)
         xnl1 = nx/2*(ny/2+1)*(nz/2+1); xnl2 = nx*(ny+1)*(nz+1) 
         ynl1 = ny/2*(nx/2+1)*(nz/2+1); ynl2 = ny*(nx+1)*(nz+1)
       if(cod_sch.eq.1)then 
       !! first coding scheme
           m1 = (nx/2+1)*ny/2+(ny/2+1)*nx/2+(nx/2+1)*(ny/2+1)         
           n2 = (nx/2+1)*ny/2+(ny/2+1)*nx/2
           n3 = (nx+1)*ny+(ny+1)*nx+(nx+1)*(ny+1)
           n4 = (nx+1)*ny+(ny+1)*nx
        ! print *,m1,n2,n3,n4 
        do i = 1,5
           d25((i-1)*5+1:i*5) = (i-1)*m1+(/0,nx+1,2*nx+2,3*nx+3,4*nx+4/)
        enddo
        do i = 1,4
           d100((i-1)*25+1:i*25,1) = d25+i-1
        enddo
        do i = 1,5
           d25((i-1)*5+1:i*5) = (i-1)*m1+(/0,1,2,3,4/)+nx/2
        enddo
        do i = 1,4
           d100((i-1)*25+1:i*25,2) = d25+(i-1)*(nx+1)
        enddo
        do i = 1,5
           d25((i-1)*5+1:i*5) = (i-1)*(nx/2+1)+(/0,1,2,3,4/)
        enddo
        do i = 1,4
           d100((i-1)*25+1:i*25,3) = d25+(i-1)*m1+n2
        enddo

        do i = 1,9
           d9(i) = (i-1)*(2*nx+1)
        enddo
        do i = 1,9
           d81((i-1)*9+1:i*9) = (i-1)*n3 +d9
        enddo
        do i = 1,8
           d648((i-1)*81+1:i*81,1) = d81+i-1
        enddo 
        
        do i = 1,9
           d9(i) = i-1+nx
        enddo
        do i =1,9
           d81((i-1)*9+1:i*9) = (i-1)*n3 +d9
        enddo
        do i =1,8
           d648((i-1)*81+1:i*81,2) = d81+(i-1)*(2*nx+1)
        enddo
        d9 = (/0:8:1/)
        do i =1,9
           d81((i-1)*9+1:i*9) = (i-1)*(nx+1) +d9
        enddo
        do i =1,8
           d648((i-1)*81+1:i*81,3) = d81+(i-1)*n3+n4
        enddo
      else
       !! for second coding scheme
       m1 = nx/2*(ny/2+1); m2 = ny/2*(nx/2+1); m3 =(nx/2+1)*(ny/2+1)
       n1 = nx*(ny+1); n2 = ny*(nx+1); n3 = (nx+1)*(ny+1)
       do i = 1,5
           d25((i-1)*5+1:i*5) = (i-1)*m1+(/0,nx/2,nx,3*nx/2,2*nx/)
        enddo
        do i = 1,4
           d100((i-1)*25+1:i*25,1) = d25+i-1
        enddo
        do i = 1,5
           d25((i-1)*5+1:i*5) = (i-1)*m2+(/0,1,2,3,4/)
        enddo
        do i = 1,4
           d100((i-1)*25+1:i*25,2) = d25+(i-1)*(nx/2+1)
        enddo
        do i = 1,5
           d25((i-1)*5+1:i*5) = (i-1)*(nx/2+1)+(/0,1,2,3,4/)
        enddo
        do i = 1,4
           d100((i-1)*25+1:i*25,3) = d25+(i-1)*m3
        enddo

        do i = 1,9
           d9(i) = (i-1)*nx
        enddo
        do i = 1,9
           d81((i-1)*9+1:i*9) = (i-1)*n1 +d9
        enddo
        do i = 1,8
           d648((i-1)*81+1:i*81,1) = d81+i-1
        enddo

        do i = 1,9
           d9(i) = i-1
        enddo
        do i =1,9
           d81((i-1)*9+1:i*9) = (i-1)*n2 +d9
        enddo
        do i =1,8
           d648((i-1)*81+1:i*81,2) = d81+(i-1)*(nx+1)
        enddo
      
       d9 = (/0:8:1/)
        do i =1,9
           d81((i-1)*9+1:i*9) = (i-1)*(nx+1) +d9
        enddo
        do i =1,8
           d648((i-1)*81+1:i*81,3) = d81+(i-1)*n3
        enddo

      endif
      ! l1 = cal_co648()
        
        do k = 1,nz/8
          do j = 1,ny/8
            do i = 1,nx/8
              !! ex
              if(cod_sch.eq.1)then
               ind = (k-1)*4*m1+(j-1)*(4*nx+4)+4*i-3 
               ind2 = (k-1)*8*n3+(j-1)*(16*nx+8)+8*i-7
              else
               ind = (k-1)*4*m1+(j-1)*nx*2+4*i-3
               ind2 = (k-1)*8*n1+(j-1)*(8*nx)+8*i-7  
              endif
               i100 = d100(:,1)+ind
               i648 = d648(:,1)+ind2
               u2(i648) = matmul(l1,u1(i100))    
              !! ey
              if(cod_sch.eq.2)then
               ind = (k-1)*4*m2+(j-1)*(2*nx+4)+4*i-3+xnl1
               ind2 = (k-1)*8*n2+(j-1)*(8*nx+8)+8*i-7+xnl2
              endif
               i100 = d100(:,2)+ind
               i648 = d648(:,2)+ind2
               u2(i648) = matmul(l1,u1(i100))
              !! ez
              if(cod_sch.eq.1)then
               ind3 =  (k-1)*4*m1+(j-1)*(2*nx+4)+4*i-3
               ind4 =  (k-1)*8*n3+(j-1)*(8*nx+8)+8*i-7
              else
               ind3 = (k-1)*4*m3+(j-1)*(2*nx+4)+4*i-3+xnl1+ynl1
               ind4 = (k-1)*8*n3+(j-1)*(8*nx+8)+8*i-7+xnl2+ynl2
              endif
               i100 = d100(:,3)+ind3
               i648 = d648(:,3)+ind4
               u2(i648) = matmul(l1,u1(i100))
             !  do m = 1,648
             !    if(i648(m).eq.155704) print *,i,j,k
             !  enddo   
            enddo
          enddo
        enddo

   contains

       function cal_co648()
        integer i,j,k,ind,ijk
        real(8) cx1(5),cx2(4),coor(648,3), cal_co648(648,100)

         cal_co648 = 0
         cx1 = (/-1.d0,-.5d0,0.d0,.5d0,1.d0/)
         cx2 = (/-.75d0,-.25d0,.25d0,.75d0/) !!inner interp
        ! cx2 = (/-1.d0,-0.33333d0,0.33333d0,1.d0/)    !! outer interp
            do k = 1,8
               do j= 1,9
                   do i = 1,9
                       ijk = (k-1)*81+(j-1)*9+i
                       coor(ijk,1) = (i-5)*0.25
                       coor(ijk,2) = (j-5)*0.25
                       coor(ijk,3)  = (k-1)*0.25-0.875
                    !   coor(ijk,3) = -7.d0/6+(k-1)/3.d0 
                   enddo
                enddo
           enddo

         do ijk = 1,648
            do k  = 1,4
               do j = 1,5
                  do i = 1,5
                     ind = (k-1)*25+(j-1)*5+i
                     cal_co648(ijk,ind)=lbf(4,cx1,i,coor(ijk,1))*lbf(4,cx1,j,coor(ijk,2))*&
                                  lbf(3,cx2,k,coor(ijk,3))
                   enddo
                enddo
            enddo
        enddo
 
       end function
      end subroutine

          function  lbf(l_ord,cx,id,x)
            integer::id,i,l_ord
            double precision:: x,cx(l_ord+1),num,den, lbf

          !  do i = 1,l_ord+1
          !     cx(i) = -1.d0+(i-1)*2.d0/l_ord
          !  enddo

            num = 1.d0; den = 1.d0
            do i = 1,l_ord+1
              if(i.ne.id)then
               num = num*(x-cx(i))
               den = den*(cx(id)-cx(i))
              endif
            enddo
           
            lbf=num/den
        end function


          function l2(id,x)
            integer::id
            double precision:: x, l2

            if (id==1) then
              l2=x*(x-1)/2
            elseif (id==2) then
              l2=-(x+1)*(x-1)
            else
              l2=x*(x+1)/2
            endif
          end function
   
        function l1(id,x)
            integer::id
            double precision:: x, l1

             if (id==1) then
                 l1=0.5-x
             else
                 l1=x+0.5

             endif
         end function

   subroutine lin_interp_n(Nx,Ny,Nz,np,u1,u2)
    !! subroutine for linear interpolation of nodes
    !! reference: Lu Kangmei et al(2011)
        implicit none
        integer,intent(in) :: nx,ny,nz,np
        complex(8),allocatable:: u1(:)
        complex(8):: u2(np)
        integer i,j,k,jN,Kn,N1,N2,N3,N4
        integer D8(8),I8(8),D27(27),I27(27),tmp1(12),tmp2(9)
        real*8 L0(27,8)

        L0 = cal_co8()

        N1 = (NX/2+1)*(NY/2+1)
        N2 = (NX+1)*(NY+1)
        
        D8 = (/1,2,2+NX/2,3+NX/2,1+N1,2+N1,2+N1+NX/2,3+N1+NX/2/)
            
        tmp2 = (/1,2,3,nx+2,nx+3,nx+4,2*nx+3,2*nx+4,2*nx+5/)
        D27 = (/TMP2, TMP2+N2, TMP2+2*N2/)
        !print *,u1(281),d12(9)

        do k = 1,nz/2
           do j = 1,ny/2
              do i = 1,nx/2
                jN = 2*N2*(K-1)+(2*NX+2)*(J-1)+2*(I-1)
                kN = N1*(K-1)+(NX/2+1)*(J-1)+I-1
                I27 = D27+JN
                I8 = D8+KN
                U2(I27) = matmul(L0,U1(I8))
              !  if(kn==0) print *,u1(i12)
              enddo
            enddo
        enddo
    end subroutine

      function  cal_co8()  
        double precision::cal_co8(27,8)
        integer i, j, k, ijk(8,7),ind

        cal_co8 = 0
        ijk = reshape((/1,3,7,9,19,21,25,27,&
                        2,2,4,6,10,12,16,18,&
                        4,6,8,8,20,20,22,24,&
                        10,12,16,18,22,24,26,26,&
                        5,5,5,5,23,23,23,23,&
                        11,11,17,17,11,11,17,17,&
                        13,15,13,15,13,15,13,15/),(/8,7/))
                                    

        do i =1,8
           cal_co8(ijk(i,1),i) = 1.d0
           cal_co8(ijk(i,2:4),i) = 0.5d0
           cal_co8(ijk(i,5:7),i) = 0.25d0
           cal_co8(14,i) = 0.125d0
        enddo

      end function
 
     subroutine lin_interp_e(Nx,Ny,nz,nl,u1,u2,cod_sch)
        !! subroutine for linear interpolation in vector field
        !! reference: Lu Kangmei et al(2011)	
        implicit none
        integer,intent(in) :: nx,ny,nz,nl,cod_sch
       ! complex(8),allocatable:: u1(:)
        complex(8):: u1(*),u2(nl)
        integer i,j,k,jN,Kn,N1,N2,N3,N4,n5,n6,n7,jn1,kn1,jn2,kn2,jn3,kn3
        integer D12(12),I12(12),D54(54),I54(54),tmp1(12),tmp2(9),ids(9),&
          ids1(18),ids2(18),ids3(18),ids4(4),ids5(4)
        real*8 L0(54,12)

        L0 = cal_co12() 
        !write(109,*) L0    
       ids4 = (/1,4,9,12/); ids5 = (/2,3,10,11/)
     if(cod_sch.eq.1)then 
        N1 = (NX/2+1)*NY/2+(NY/2+1)*NX/2
        N2 = N1+(NX/2+1)*(NY/2+1)
        N3 = (NX+1)*NY+(NY+1)*NX
        N4 = N3+(NX+1)*(NY+1)
        D12 = (/1,1+NX/2,2+NX/2,2+NX,1+N1,2+N1,2+N1+NX/2,3+N1+NX/2,&
                1+N2, 1+n2+nx/2,2+n2+nx/2,2+nx+n2/)
        tmp1 = (/1,2,nx+1,nx+2,nx+3,2*nx+2,2*nx+3,3*nx+2,3*nx+3,3*nx+4,4*nx+3,4*nx+4/)  
        tmp2 = (/1,2,3,nx+2,nx+3,nx+4,2*nx+3,2*nx+4,2*nx+5/)
        D54 = (/TMP1, TMP2+N3, TMP1+N4, TMP2+N3+N4, TMP1+2*N4/)  
      else
       n1 = (ny/2+1)*nx/2*(nz/2+1)
       n2 = (ny/2+1)*nx/2*(nz/2+1)+(nx/2+1)*ny/2*(nz/2+1) 
       n3 = (ny+1)*nx*(nz+1)
       n4 = (ny+1)*nx*(nz+1)+(nx+1)*ny*(nz+1)
       n5 = nx*(ny+1); n6 = ny*(nx+1); n7=(nx+1)*(ny+1)
       d12 =(/1,n1+1,n1+2,1+nx/2, n2+1,n2+2,n2+nx+1,n2+nx+2,  &
           1+nx/2*(ny/2+1),n1+ny/2*(nx/2+1)+1,n1+ny/2*(nx/2+1)+2,1+nx/2*(ny/2+1)+nx/2/) 
      ! tmp1 = ()
        tmp2=(/1,nx+1,2*nx+1,1+n5,1+n5+nx,1+n5+2*nx,1+2*n5,1+2*n5+nx,1+2*n5+2*nx/)
        ids=(/1,6,11,22,27,32,43,48,53/)
        d54(ids) = tmp2; d54(ids+1) =tmp2+1
        ids1= (/ids,ids+1/)
        tmp2 =n3+(/1,2,3,n6+1,n6+2,n6+3,2*n6+1,2*n6+2,2*n6+3/)
        ids =(/3,4,5,24,25,26,45,46,47/)
        ids2 = (/ids,ids+5/)   
        d54(ids) = tmp2; d54(ids+5) =tmp2+nx+1
        tmp2 =n4+(/1,2,3,nx+2,nx+3,nx+4, 2*nx+3,2*nx+4,2*nx+5/)
        ids =(/13:21/)
        ids3 =(/ids,ids+21/)
        d54(ids) = tmp2; d54(ids+21) =tmp2+n7

      endif
!print *,u1(281),d12(9)

        do k = 1,nz/2
        do j = 1,ny/2
        do i = 1,nx/2
        if(cod_sch.eq.1)then
         jN = 2*N4*(K-1)+(4*NX+2)*(J-1)+2*(I-1)
         kN = N2*(K-1)+(NX+1)*(J-1)+I-1
         i54 = d54+ jn; i12 =d12+kn
       else
         jn1 = 2*(k-1)*n5+nx*(j-1)*2 +2*(i-1)        
         jn2 = 2*(k-1)*n6+(nx+1)*(j-1)*2+2*(i-1) 
         jn3 = 2*(k-1)*n7+(nx+1)*(j-1)*2+2*(i-1) 
         i54(ids1) = d54(ids1)+jn1
         i54(ids2) = d54(ids2)+jn2
         i54(ids3) = d54(ids3)+jn3
         kn1 =(k-1)*nx/2*(ny/2+1)+(j-1)*nx/2+i-1
         i12(ids4) = D12(ids4)+kn1
         kn2 =(k-1)*ny/2*(nx/2+1)+(j-1)*(nx/2+1)+i-1
         i12(ids5) = D12(ids5)+kn2
         kn3 =(k-1)*(nx/2+1)*(ny/2+1)+(j-1)*(nx/2+1)+i-1
         i12(5:8) = D12(5:8)+kn3
       endif
        U2(I54) = matmul(L0,U1(I12))
!  if(kn==0) print *,u1(i12)
        enddo
        enddo
        enddo
  end subroutine


  function cal_co12()

        double precision::cal_co12(54,12)
        integer i, j, k, ijk(12,8),ind

        cal_co12 = 0
        ijk = transpose(reshape((/1,2,6,7,22,23,27,28, &
                       3,8,4,9,24,29,25,30, &
                       5,10,4,9,26,31,25,30, &
                       11,12,6,7,32,33,27,28, &
                       13,34,14,35,16,37,17,38,&
                       15,36,14,35,18,39,17,38,&
                       19,40,16,37,20,41,17,38,&
                       21,42,18,39,20,41,17,38,&
                       43,44,48,49,22,23,27,28,&
                       45,50,46,51,24,29,25,30,&
                       47,52,46,51,26,31,25,30,&
                       53,54,48,49,32,33,27,28/),(/8,12/)))

          do i = 1,12
             cal_co12(ijk(i,1:2),i) = 1.d0  !!outer surface edges(24)
             cal_co12(ijk(i,3:6),i) = 0.5d0  !!interior surface edges(24)
             cal_co12(ijk(i,7:8),i) = 0.25d0  !! pure interior edges(6)
          enddo
        return

 end function

    subroutine quad_interp_n(U3,nk,Nx,Ny,Nz,ND,abu)
    integer::Nx,Ny,Nz,ND,nk,abu(6)
    complex(8)::U3(ND)
 
    integer i,j,k,jN,N1,N2,Nx4,Ny4,Nz4,idc(3)
    integer D27(27),D125(125),temp5(5),temp25(25),I27(27),I125(125)
    real*8 L0(125,27)

    N1=(Nx+1);      N2=(Nx+1)*(Ny+1)
    Nx4=Nx/4;       Ny4=Ny/4;       Nz4=Nz/4
    D27=(/1, 3, 5, 2*N1+1, 2*N1+3, 2*N1+5, 4*N1+1, 4*N1+3, 4*N1+5, &
        &(/1, 3, 5, 2*N1+1, 2*N1+3, 2*N1+5, 4*N1+1, 4*N1+3, 4*N1+5/)+2*N2, &
        &(/1, 3, 5, 2*N1+1, 2*N1+3, 2*N1+5, 4*N1+1, 4*N1+3, 4*N1+5/)+4*N2 /)
    temp5=(/1:5/)
    temp25=(/temp5,  temp5+N1,  temp5+2*N1,  temp5+3*N1,  temp5+4*N1/)
    D125 =(/temp25, temp25+N2, temp25+2*N2, temp25+3*N2, temp25+4*N2/)
   ! print *,u3(d27)

    do k=1,Nz4
      do j=1,Ny4
        do i=1,Nx4
         if(4*i.lt.abu(1))then
             idc(1) =1
           elseif(4*i.ge.abu(1).and.4*i.le.abu(2))then
             idc(1) = 2
           else
             idc(1) = 3
           endif
           if(4*j.lt.abu(3))then
             idc(2) = 1              
           elseif(4*j.ge.abu(3).and.4*j.le.abu(4))then
             idc(2) = 2                  
           else
             idc(2) = 3
           endif
           if(4*k.lt.abu(5))then
             idc(3) = 1     
           elseif(4*k.ge.abu(5).and.4*k.le.abu(6))then
             idc(3) = 2          
           else
             idc(3) = 3
          endif
          L0 = co27(idc,nk)
          jN=4*(k-1)*(Nx+1)*(Ny+1)+4*(j-1)*(Nx+1)+4*(i-1)
          I27=D27+jN
          I125=D125+jN
          U3(I125)=matmul(L0,U3(I27))
         ! stop
         ! if(i.eq.1.and.j.eq.1.and.k.eq.1) print *,u3(i27),l0
        enddo
      enddo
    enddo

   contains 
    function co27(idc,nk)
      !! qudratic interpolation engaged with expansion coef 
      integer idc(3)        
      double precision co27(125,27)
      integer i, j, k, ijk,nk,ind
      double precision coor(125,3),cox(5),coz(5),lamh,lamv,bas1,bas2

      lamh = lambda_h**(nk-1)
      lamv = lambda_v**(nk-1)
      bas1 = 2.d0/(1+lamh+lamh**2+lamh**3)
      bas2 = 2.d0/(1+lamv+lamv**2+lamv**3)
      cox(1) = -1;coz(1) = -1
      cox(2) = -1+bas1; coz(2) = -1+bas2;
      cox(3) = -1+bas1*(1+lamh); coz(3) = -1+bas2*(1+lamv);
      cox(4) = -1+bas1*(1+lamh+lamh**2); coz(4) = -1+bas2*(1+lamv+lamv**2);
      cox(5) = 1; coz(5)  =1

      do k=1,5
        do j=1,5
            do i=1,5
                ijk=(k-1)*25+(j-1)*5+i
                if(idc(1).eq.1)then
                  coor(ijk,1) = -cox(6-i)
                elseif(idc(1).eq.2)then
                  coor(ijk,1) = (i-3)*0.5
                else
                  coor(ijk,1) = cox(i)
                endif
               if(idc(2).eq.1)then
                  coor(ijk,2) = -cox(6-j)
                elseif(idc(2).eq.2)then
                  coor(ijk,2) = (j-3)*0.5
                else
                  coor(ijk,2) = cox(j)
                endif
               if(idc(3).eq.1)then
                  coor(ijk,3) = -coz(6-k)
                elseif(idc(3).eq.2)then
                  coor(ijk,3) = (k-3)*0.5
                else
                  coor(ijk,3) = coz(k)
                endif
 
            enddo
        enddo
      enddo
    !  print *,coor(1,:)
      do ijk=1,125
        do k=1,3
            do j=1,3
                do i=1,3
                    ind = (k-1)*9+(j-1)*3+i
                    co27(ijk,ind) =ll2(i,1,idc(1),nk,coor(ijk,1))*&
                     ll2(j,1,idc(2),nk,coor(ijk,2))*ll2(k,2,idc(3),nk,coor(ijk,3))
                    ! if(ijk.eq.1) print *,ll2(i,1,idc(1),nk,coor(ijk,1))
                enddo
            enddo
        enddo
      enddo

    end function
 
    function ll2(id,ig,im,nk,x)
        
    integer::id,ig,im,nk
    double precision:: x, midx,ll2
 
    select case(im)
      case(1)
        if(ig.eq.1)then 
          midx = -1+2*lambda_h**nk/(1+lambda_h**nk)
        else
          midx = -1+2*lambda_v**nk/(1+lambda_v**nk)
        endif
      case(2)
        midx = 0
      case(3)
        if(ig.eq.1)then
          midx = -1+2/(1+lambda_h**nk)
        else
          midx = -1+2/(1+lambda_v**nk) 
        endif
    end select

    if (id==1) then
        ll2=(x-midx)*(x-1)/(-1-midx)/(-1-1)
    elseif (id==2) then
        ll2=(x+1)*(x-1)/(midx+1)/(midx-1)
    else
        ll2=(x-midx)*(x+1)/(1-midx)/(1+1)
    endif
 
 
    end
    end subroutine
 
end module
