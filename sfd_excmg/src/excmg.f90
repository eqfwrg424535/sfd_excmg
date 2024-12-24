module excmg
 !use global_parameter
 implicit none 
  !real(8), parameter::pi=3.1592653589793238462643d0
   real(8), parameter::order = 2.0d0

  contains

    function extint_v(nx,ny,nz,nl,u0,u1,a1,b1,c1,nair)
     !! prolongation for VFEM  
      integer :: nx,ny,nz,nl,nair
      real(8),allocatable::a1(:),b1(:),c1(:),a2(:),b2(:),c2(:)
      complex(8) :: u0(*),u1(*),extint_v(nl)
      complex(8),allocatable:: utmp(:),w0(:),w1(:)
      integer i,j,k,nl_new,n1,nl2,nl3,xnl,ynl,xnl2,xnl3,ynl2,ynl3
      real(8) co1,co2,co3

      !! check the length of a1,b1,and c1, if they cannot be devisible by 4, stop
      if(mod(size(a1),4).ne.0.or.mod(size(b1),4).ne.0.or.mod(size(c1),4).ne.0)then
         print *,'error in interpolation, please check the input mesh'
        stop
      endif 
      
      allocate(a2(nx/2),b2(ny/2),c2(nz/2))
      do i = 1,nx/2
        a2(i) = a1(i*2-1)+a1(i*2)
      enddo
      do i = 1,ny/2
        b2(i) = b1(i*2-1)+b1(i*2)
      enddo
      do i = 1,nz/2
        c2(i) = c1(i*2-1)+c1(i*2)
      enddo
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
   
     ! call lin_interp_e(Nx/2,Ny/2,Nz/2,n1,u0,utmp,a1,b1,c1,nair) 
      call quad_interp_e(nx/2,ny/2,nz/2,n1,u0,utmp,a2,b2,c2,nair/2)
      deallocate(a2,b2,c2)
      
       do i = 1,n1
          utmp(i) = co1*u1(i)-co2*utmp(i)
         ! if(abs(utmp(i)).le.1.d-50.and.utmp(i).ne.0) print *,i,w1(i),utmp(i) 
       enddo
   
       call cubic_interp(nx,ny,nz,nl,utmp,extint_v,a1,b1,c1,nair)
    
    ! print *,'prolong completed'
     !   extint_v(1:xnl) = extint_v(1:xnl)/nx
     !   extint_v(xnl+1:xnl+ynl) = extint_v(xnl+1:xnl+ynl)/ny
     !   extint_v(xnl+ynl+1:nl) = extint_v(xnl+ynl+1:nl)/nz
    ! extint_v = cmplx(real(extint_v)/2,imag(extint_v)/4)
    ! extint_v =  extint_v/2
     do i = 1,nl
         if(isnan(abs(extint_v(i)))) print *,'NaN found in prolongation',i
     enddo
      deallocate(utmp) 
    end function

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


    subroutine quad_interp_e(nx,ny,nz,nl,u1,u2,a1,b1,c1,nair)
        implicit none
        integer,intent(in) :: nx,ny,nz,nl,nair
        complex(8):: u1(*),u2(nl)
        integer i,j,k,ind,ind2,nd,N1,N2,N3,N4,n5,n6,ind3,ind4, d9(9,3),d27(27),&
                d5(5),d25(25),d100(100,3),i27(27),g27(27),i100(100),xnl1,ynl1,&
                xnl2,ynl2,i18(18)
        ! real(8) L1(100,27),L2(100,27),l3(100,18)
        real(8),allocatable :: a1(:),b1(:),c1(:),aax(:),bby(:),ccz(:),l1(:,:,:),tmp(:)
        
        allocate(aax(nx+1),bby(ny+1),ccz(nz+1),tmp(nair))
       tmp = c1(1:nair)
       aax(1) = -sum(a1)/2; bby(1) =-sum(b1)/2; ccz(1) = -sum(tmp)
       do i = 1,nx
          aax(i+1) = aax(i)+a1(i)    
       enddo
       do i = 1,ny
          bby(i+1) = bby(i)+b1(i)    
       enddo
       do i = 1,nz
          ccz(i+1) = ccz(i)+c1(i)    
       enddo
        ! L1 = cal_co100(1)
        ! L2 = cal_co100(2)
        !l3 = cal_co18()
        nd  = (nx+1)*(ny+1)*(nz+1)
        xnl1 = nx/2*(ny/2+1)*(nz/2+1); ynl1 = (nx/2+1)*ny/2*(nz/2+1);
        xnl2 = nx*(ny+1)*(nz+1); ynl2 = (nx+1)*ny*(nz+1)
        !allocate(utmp1(nd),utmp2(nd),utmp3(nd))
        !utmp1 = 0; utmp2= 0; utmp3 =0
     
         n1 = (ny/2+1)*nx/2; n2 = nx*(ny+1)
         n3 = (nx/2+1)*ny/2; n4 = ny*(nx+1)
         n5 = (nx/2+1)*(ny/2+1); n6 = (nx+1)*(ny+1)
     
         !! x-edges
         d9(:,1) = (/0,nx/2,nx, n1,n1+nx/2, n1+nx, 2*n1,2*n1+nx/2,2*n1+nx/)
        do i= 1,5
            d25((i-1)*5+1:5*i) = (i-1)*n2+(/0,nx,2*nx,3*nx,4*nx/)
        enddo
        do i= 1,4
           d100((i-1)*25+1:25*i,1) = d25+i-1
        enddo
        
          !! y-edges
          d9(:,2) = xnl1+(/0,1,2,n3,n3+1,n3+2, 2*n3,2*n3+1,2*n3+2/)
        do i = 1,5
            d25((i-1)*5+1:5*i) = xnl2+(i-1)*n4+(/0:4/)
        enddo
        do i =1,4
           d100((i-1)*25+1:25*i,2) = d25+(i-1)*(nx+1)
        enddo
        
         !! z- edges
         d9(:,3) = xnl1+ynl1+(/0,1,2, nx/2+1,nx/2+2,nx/2+3, nx+2,nx+3,nx+4/)
        do i = 1,5
           d25((i-1)*5+1:5*i) = (i-1)*(nx+1)+(/0:4/)
        enddo
         do i = 1,4
           d100((i-1)*25+1:25*i,3) = d25+(i-1)*n6+xnl2+ynl2
        enddo
        
        
        !! interpolate
        do k = 1,Nz/4
            do j = 1,Ny/4
                do i = 1,Nx/4
                    l1 = cal_coef(aax,bby,ccz,2,i,j,k)
                    ! ex
                    ind = (k-1)*4*n2+(j-1)*(4*nx)+4*i-3
                    ind2 = (k-1)*2*n1+(j-1)*nx+2*i-1
                    i100 = d100(:,1)+ind
                    i18 = ind2+(/d9(:,1),d9(:,1)+1/)
                    u2(I100) = matmul(L1(:,:,1),u1(i18))
                    ! ey
                      ind = (k-1)*4*n4+(j-1)*(4*nx+4)+4*i-3
                      ind2 = (k-1)*2*n3+(j-1)*(nx+2)+2*i-1
                     i100 = d100(:,2)+ind 
                      i18 = ind2+(/d9(:,2),d9(:,2)+nx/2+1/)
                      u2(i100) = matmul(L1(:,:,2),u1(i18))
                    ! ez
                      ind = (k-1)*4*n6+(j-1)*(4*nx+4)+4*i-3
                     ind2 = (k-1)*2*n5+(j-1)*(nx+2)+2*i-1
                      i100 = d100(:,3)+ind
                      i18 = ind2+(/d9(:,3),d9(:,3)+n5/)
                      u2(I100) = matmul(L1(:,:,3),u1(i18))
                      deallocate(l1)
                enddo
             enddo
        enddo
    
        deallocate(aax,bby,ccz,tmp)

         !contains
         !     function cal_co100(tp)
         !         integer  tp,i,j,k,ind,ijk
         !          real(8) coor(100,3), cal_co100(100,27)
         !
         !          do k = 1,4
         !              do j= 1,5
         !                  do i = 1,5
         !                      ijk = (k-1)*25+(j-1)*5+i
         !                      coor(ijk,1) = (i-3)*0.5
         !                      coor(ijk,2) = (j-3)*0.5
         !                      if(tp.eq.1)then
         !                         coor(ijk,3)  = (k-1)*0.5-1.25
         !                      else
         !                         coor(ijk,3)  = (k-1)*0.5-0.25
         !                      endif
         !                  enddo
         !               enddo
         !          enddo
         !
         !        do ijk = 1,100
         !           do k  = 1,3
         !              do j = 1,3
         !                 do i = 1,3
         !                    ind = (k-1)*9+(j-1)*3+i
         !                    cal_co100(ijk,ind) =l2(i,coor(ijk,1))*l2(j,coor(ijk,2))*&
         !                                 l2(k,coor(ijk,3))
         !                  enddo
         !               enddo
         !           enddo
         !        enddo
         !
         !      end function
         !
         !      function cal_co18()
         !          integer  i,j,k,ind,ijk
         !          real(8) coor(100,3), cal_co18(100,18)
         !
         !          do k = 1,4
         !              do j= 1,5
         !                  do i = 1,5
         !                      ijk = (k-1)*25+(j-1)*5+i
         !                      coor(ijk,1) = (i-3)*0.5
         !                      coor(ijk,2) = (j-3)*0.5
         !                      coor(ijk,3)  = (k-1)*0.5-0.75
         !
         !                  enddo
         !               enddo
         !          enddo
         !
         !        do ijk = 1,100
         !           do k  = 1,2
         !              do j = 1,3
         !                 do i = 1,3
         !                    ind = (k-1)*9+(j-1)*3+i
         !                    cal_co18(ijk,ind)=l2(i,coor(ijk,1))*l2(j,coor(ijk,2))*&
         !                                 l1(k,coor(ijk,3))
         !                  enddo
         !               enddo
         !             enddo
         !        enddo
         !     end function
     end subroutine     

     subroutine cubic_interp(nx,ny,nz,nl,u1,u2,a1,b1,c1,nair)
       implicit none
       integer nx,ny,nz,nl,nair
       complex(8)::u1(*),u2(nl)
        integer i,j,k,m,ind,ind2,nd,n1,n2,n3,n4,m1,m2,m3,i25(25),i100(100),d9(9),&
            d81(81),d648(648,3),d25(25),d100(100,3),i648(648),xnl1,ynl1,&
            xnl2,ynl2,ind3,ind4
       real(8) a1(nx),b1(ny),c1(nz)
       real(8),allocatable :: aax(:),bby(:),ccz(:),l1(:,:,:),tmp(:)
       
       allocate(aax(nx+1),bby(ny+1),ccz(nz+1),tmp(nair))
       tmp = c1(1:nair)
       aax(1) = -sum(a1)/2; bby(1) =-sum(b1)/2; ccz(1) = -sum(tmp)
       do i = 1,nx
          aax(i+1) = aax(i)+a1(i)    
       enddo
       do i = 1,ny
          bby(i+1) = bby(i)+b1(i)    
       enddo
       do i = 1,nz
          ccz(i+1) = ccz(i)+c1(i)    
       enddo
        ! print *,nx,ny
        ! l1 = cal_co648() !!l1 should related with the mesh
       
        ! print *,nx,ny
         nd =(nx+1)*(ny+1)*(nz+1)
         xnl1 = nx/2*(ny/2+1)*(nz/2+1); xnl2 = nx*(ny+1)*(nz+1) 
         ynl1 = ny/2*(nx/2+1)*(nz/2+1); ynl2 = ny*(nx+1)*(nz+1)
     
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
      ! l1 = cal_co648()
        
        do k = 1,nz/8
          do j = 1,ny/8
            do i = 1,nx/8
              l1 = cal_coef(aax,bby,ccz,4,i,j,k) ! 4 indicates the interpolation cell size
              !! ex
               ind = (k-1)*4*m1+(j-1)*nx*2+4*i-3
               ind2 = (k-1)*8*n1+(j-1)*(8*nx)+8*i-7  
               i100 = d100(:,1)+ind
               i648 = d648(:,1)+ind2
               u2(i648) = matmul(l1(:,:,1),u1(i100))    
              !! ey
               ind = (k-1)*4*m2+(j-1)*(2*nx+4)+4*i-3+xnl1
               ind2 = (k-1)*8*n2+(j-1)*(8*nx+8)+8*i-7+xnl2
               i100 = d100(:,2)+ind
               i648 = d648(:,2)+ind2
               u2(i648) = matmul(l1(:,:,2),u1(i100))
              !! ez
               ind = (k-1)*4*m3+(j-1)*(2*nx+4)+4*i-3+xnl1+ynl1
               ind2 = (k-1)*8*n3+(j-1)*(8*nx+8)+8*i-7+xnl2+ynl2
               i100 = d100(:,3)+ind
               i648 = d648(:,3)+ind2
               u2(i648) = matmul(l1(:,:,3),u1(i100))
             !  do m = 1,648
             !    if(i648(m).eq.155704) print *,i,j,k
               deallocate(l1)
             !  enddo   
            enddo
          enddo
        enddo

   !contains
   !
   !    function cal_co648()
   !     integer i,j,k,ind,ijk
   !     real(8) cx1(5),cx2(4),coor(648,3), cal_co648(648,100)
   !
   !      cal_co648 = 0
   !      cx1 = (/-1.d0,-.5d0,0.d0,.5d0,1.d0/)
   !      cx2 = (/-.75d0,-.25d0,.25d0,.75d0/) !!inner interp
   !     ! cx2 = (/-1.d0,-0.33333d0,0.33333d0,1.d0/)    !! outer interp
   !         do k = 1,8
   !            do j= 1,9
   !                do i = 1,9
   !                    ijk = (k-1)*81+(j-1)*9+i
   !                    coor(ijk,1) = (i-5)*0.25
   !                    coor(ijk,2) = (j-5)*0.25
   !                    coor(ijk,3)  = (k-1)*0.25-0.875
   !                 !   coor(ijk,3) = -7.d0/6+(k-1)/3.d0 
   !                enddo
   !             enddo
   !        enddo
   !
   !      do ijk = 1,648
   !         do k  = 1,4
   !            do j = 1,5
   !               do i = 1,5
   !                  ind = (k-1)*25+(j-1)*5+i
   !                  cal_co648(ijk,ind)=lbf(4,cx1,i,coor(ijk,1))*lbf(4,cx1,j,coor(ijk,2))*&
   !                               lbf(3,cx2,k,coor(ijk,3))
   !                enddo
   !             enddo
   !         enddo
   !     enddo
   !
   !    end function
        deallocate(aax,bby,ccz,tmp)
   end subroutine
   
    function cal_coef(aax,bby,ccz,merge,n1,n2,n3) result(coef)  
        ! merge:2--linear; 3--quadratic; 4--cubic
        real(8) aax(:),bby(:),ccz(:),cx1(merge+1),cx2(merge+1),cx3(merge)
        real(8),allocatable :: coef(:,:,:),coor(:,:)
        integer merge,i,j,k,ind,ijk,n1,n2,n3,nx,ny,nz
    
        nx = size(aax); ny = size(bby); nz = size(ccz) 
        allocate(coor((merge*2+1)**2*merge*2,3),coef((merge*2+1)**2*merge*2,(merge+1)**2*merge,3))
       ! allocate(cx1(merge+1),cx2(merge+1),cx3(merge)) !cx1--nodes; cx2--edge midpoints
        !! for ex
        do k = 1,merge*2
           do j= 1,merge*2+1
              do i = 1,merge*2+1
                 ijk = (k-1)*(merge*2+1)**2+(j-1)*(merge*2+1)+i 
                 coor(ijk,1) = (aax((n1-1)*merge*2+k)+aax((n1-1)*merge*2+k+1))/2
                 coor(ijk,2) = bby((n2-1)*merge*2+i)
                 coor(ijk,3) = ccz((n3-1)*merge*2+j)
              enddo
           enddo
        enddo
        
        cx1 = bby((n2-1)*merge*2+1:n2*merge*2+1:2)
        cx2 = ccz((n3-1)*merge*2+1:n3*merge*2+1:2)
        cx3 = (aax((n1-1)*merge*2+1:n1*merge*2:2)+aax((n1-1)*merge*2+2:n1*merge*2+1:2))/2
         do ijk = 1,(merge*2+1)**2*merge*2
            do k  = 1,merge
               do j = 1,merge+1
                  do i = 1,merge+1
                     ind = (k-1)*(merge+1)**2+(j-1)*(merge+1)+i             
                     coef(ijk,ind,1)=lbf(merge,cx1,i,coor(ijk,2))*lbf(merge,cx2,j,coor(ijk,3))*&
                                  lbf(merge-1,cx3,k,coor(ijk,1))
                   enddo
                enddo
            enddo
        enddo
       !! for ey
        do k = 1,merge*2
           do j= 1,merge*2+1
              do i = 1,merge*2+1
                 ijk = (k-1)*(merge*2+1)**2+(j-1)*(merge*2+1)+i 
                 coor(ijk,1) =  aax((n1-1)*merge*2+i)
                 coor(ijk,2) = (bby((n2-1)*merge*2+k)+bby((n2-1)*merge*2+k+1))/2
                coor(ijk,3) = ccz((n3-1)*merge*2+j)
              enddo
           enddo
        enddo
        
        cx1 = aax((n1-1)*merge*2+1:n1*merge*2+1:2)
        cx2 = ccz((n3-1)*merge*2+1:n3*merge*2+1:2)
        cx3 =(bby((n2-1)*merge*2+1:n2*merge*2:2)+bby((n2-1)*merge*2+2:n2*merge*2+1:2))/2
        
         do ijk = 1,(merge*2+1)**2*merge*2
            do k  = 1,merge
               do j = 1,merge+1
                  do i = 1,merge+1
                     ind = (k-1)*(merge+1)**2+(j-1)*(merge+1)+i
                     coef(ijk,ind,2)=lbf(merge,cx1,i,coor(ijk,1))*lbf(merge,cx2,j,coor(ijk,3))*&
                                  lbf(merge-1,cx3,k,coor(ijk,2))
                   enddo
                enddo
            enddo
         enddo
        !! for ez
        do k = 1,merge*2
           do j= 1,merge*2+1
              do i = 1,merge*2+1
                 ijk = (k-1)*(merge*2+1)**2+(j-1)*(merge*2+1)+i 
                 coor(ijk,1) = aax((n1-1)*merge*2+i)
                 coor(ijk,2) = bby((n2-1)*merge*2+j)
                coor(ijk,3) = (ccz((n3-1)*merge*2+k)+ccz((n3-1)*merge*2+k+1))/2
              enddo
           enddo
        enddo
        cx1 = aax((n1-1)*merge*2+1:n1*merge*2+1:2)
        cx2 = bby((n2-1)*merge*2+1:n2*merge*2+1:2)
        cx3 =(ccz((n3-1)*merge*2+1:n3*merge*2:2)+ccz((n3-1)*merge*2+2:n3*merge*2+1:2))/2
        
         do ijk = 1,(merge*2+1)**2*merge*2
            do k  = 1,merge
               do j = 1,merge+1
                  do i = 1,merge+1
                     ind = (k-1)*(merge+1)**2+(j-1)*(merge+1)+i
                     coef(ijk,ind,3)=lbf(merge,cx1,i,coor(ijk,1))*lbf(merge,cx2,j,coor(ijk,2))*&
                                  lbf(merge-1,cx3,k,coor(ijk,3))
                   enddo
                enddo
            enddo
         enddo
    
         deallocate(coor)
    end function

          function  lbf(l_ord,cx,id,x)
          ! Lagrange basis function (l_ord is the order)
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
 
!     subroutine lin_interp_e(Nx,Ny,nz,nl,u1,u2,a1,b1,c1,nair)
!        !! subroutine for linear interpolation in vector field
!        !! reference: Lu Kangmei et al(2011)	
!        implicit none
!        integer,intent(in) :: nx,ny,nz,nl,nair
!        real*8,allocatable :: aax(:),bby(:),ccz(:),L0(:,:,:)
!       ! complex(8),allocatable:: u1(:)
!        complex(8):: u1(*),u2(nl)
!        integer i,j,k,jN,Kn,N1,N2,N3,N4,n5,n6,n7,jn1,kn1,jn2,kn2,jn3,kn3
!        integer D12(12),I12(12),D54(54),I54(54),tmp1(12),tmp2(9),ids(9),&
!          ids1(18),ids2(18),ids3(18),ids4(4),ids5(4)
!        
!        !L0 = cal_co12() 
!        !write(109,*) L0    
!       ids4 = (/1,4,9,12/); ids5 = (/2,3,10,11/)
!  
!       n1 = (ny/2+1)*nx/2*(nz/2+1)
!       n2 = (ny/2+1)*nx/2*(nz/2+1)+(nx/2+1)*ny/2*(nz/2+1) 
!       n3 = (ny+1)*nx*(nz+1)
!       n4 = (ny+1)*nx*(nz+1)+(nx+1)*ny*(nz+1)
!       n5 = nx*(ny+1); n6 = ny*(nx+1); n7=(nx+1)*(ny+1)
!       d12 =(/1,n1+1,n1+2,1+nx/2, n2+1,n2+2,n2+nx+1,n2+nx+2,  &
!           1+nx/2*(ny/2+1),n1+ny/2*(nx/2+1)+1,n1+ny/2*(nx/2+1)+2,1+nx/2*(ny/2+1)+nx/2/) 
!      ! tmp1 = ()
!        tmp2=(/1,nx+1,2*nx+1,1+n5,1+n5+nx,1+n5+2*nx,1+2*n5,1+2*n5+nx,1+2*n5+2*nx/)
!        ids=(/1,6,11,22,27,32,43,48,53/)
!        d54(ids) = tmp2; d54(ids+1) =tmp2+1
!        ids1= (/ids,ids+1/)
!        tmp2 =n3+(/1,2,3,n6+1,n6+2,n6+3,2*n6+1,2*n6+2,2*n6+3/)
!        ids =(/3,4,5,24,25,26,45,46,47/)
!        ids2 = (/ids,ids+5/)   
!        d54(ids) = tmp2; d54(ids+5) =tmp2+nx+1
!        tmp2 =n4+(/1,2,3,nx+2,nx+3,nx+4, 2*nx+3,2*nx+4,2*nx+5/)
!        ids =(/13:21/)
!        ids3 =(/ids,ids+21/)
!        d54(ids) = tmp2; d54(ids+21) =tmp2+n7
!
!
!        do k = 1,nz/2
!          do j = 1,ny/2
!            do i = 1,nx/2
!              l0= cal_coef(aax,bby,ccz,2,i,j,k)
!              jn1 = 2*(k-1)*n5+nx*(j-1)*2 +2*(i-1)        
!              jn2 = 2*(k-1)*n6+(nx+1)*(j-1)*2+2*(i-1) 
!              jn3 = 2*(k-1)*n7+(nx+1)*(j-1)*2+2*(i-1) 
!              i54(ids1) = d54(ids1)+jn1
!              i54(ids2) = d54(ids2)+jn2
!         i54(ids3) = d54(ids3)+jn3
!         kn1 =(k-1)*nx/2*(ny/2+1)+(j-1)*nx/2+i-1
!         i12(ids4) = D12(ids4)+kn1
!         kn2 =(k-1)*ny/2*(nx/2+1)+(j-1)*(nx/2+1)+i-1
!         i12(ids5) = D12(ids5)+kn2
!         kn3 =(k-1)*(nx/2+1)*(ny/2+1)+(j-1)*(nx/2+1)+i-1
!         i12(5:8) = D12(5:8)+kn3
!    
!         U2(I54) = matmul(L0,U1(I12))
!         deallocate(l0)
!!  if(kn==0) print *,u1(i12)
!        enddo
!        enddo
!        enddo
!        deallocate(aax,bby,ccz)
!  end subroutine


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

   
 
end module
