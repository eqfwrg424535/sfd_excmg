module  divergence_corr
   use global_parameter
    use model_selection
   implicit none
   !character,parameter::pre_type = 'S'
  
    contains
     subroutine div_correc(nl,a,b,c,x,ep,sigma,el,emap,dc,d_mat,d_ia,d_ja,t_mat,t_ia,&
            t_ja,tol_cor,icount)
           implicit none
           real*8,allocatable:: el(:),sigma(:,:),a(:),b(:),c(:)
           integer(8),allocatable:: d_ia(:),t_ia(:)!,nnz,nnz2
           integer,allocatable::t_ja(:),d_ja(:),emap(:,:)!,edgeid(:),nid(:),dian(:)
           integer::nl,icount
           complex(8) x(:),tmp
           complex(8),allocatable :: d_mat(:),t_mat(:),ep(:)
           complex(8),allocatable :: t(:),t1(:),t2(:),xold(:),phi(:),ele_slice(:,:)
         !  real(8),allocatable::drmat(:), tr(:),ti(:),phir(:),phii(:)
           logical,allocatable::dc(:)
           real*8 tol_cor,w,nm1(12,8),nm2(12,8),nm3(12,8),te(12,8),hx,hy,hz,ele_r,ele_i
           real iter
           integer i,j,k,h,np,ie,nx,ny,nz,n1,n2,n3,itmax,ind,ind2,ind3,ids(12),idn(8),&
               iter1
           character(len=5) str

!         NM1 = reshape((/-4,-2,-2,-1,-2,-2,-1,-1,-2,-1,-2,-1,&
!                         4,2,2,1,-2,-2,-1,-1,-2,-1,-2,-1,&
!                         2,4,1,2,2,2,1,1,-1,-2,-1,-2, &
!                        -2,-4,-1,-2,2,2,1,1,-1,-2,-1,-2,&
!                        -2,-1,-4,-2,-1,-1,-2,-2,2,1,2,1,&
!                         2,1,4,2,-1,-1,-2,-2,2,1,2,1,&
!                         1,2,2,4,1,1,2,2,1,2,1,2,&
!                        -1,-2,-2,-4,1,1,2,2,1,2,1,2/),(/12,8/))*1.d0
!
!        NM2 = reshape((/-2,-1,-2,-1,-4,-2,-2,-1,-2,-2,-1,-1,&
!                         2,1,2,1,-2,-1,-4,-2,-1,-1,-2,-2,&
!                         2,1,2,1,2,1,4,2,-1,-1,-2,-2,&
!                        -2,-1,-2,-1,4,2,2,1,-2,-2,-1,-1,&
!                        -1,-2,-1,-2,-2,-4,-1,-2,2,2,1,1,&
!                         1,2,1,2,-1,-2,-2,-4,1,1,2,2,&
!                         1,2,1,2,1,2,2,4,1,1,2,2,&
!                        -1,-2,-1,-2,2,4,1,2,2,2,1,1/),(/12,8/))*1.d0
!
!        NM3 = reshape((/-2,-2,-1,-1,-2,-1,-2,-1,-4,-2,-2,-1,&
!                         2,2,1,1,-1,-2,-1,-2,-2,-4,-1,-2,&
!                         1,1,2,2,1,2,1,2,-1,-2,-2,-4,&
!                        -1,-1,-2,-2,2,1,2,1,-2,-1,-4,-2,&
!                        -2,-2,-1,-1,-2,-1,-2,-1,4,2,2,1,&
!                         2,2,1,1,-1,-2,-1,-2,2,4,1,2,&
!                         1,1,2,2,1,2,1,2,1,2,2,4,&
!                        -1,-1,-2,-2,2,1,2,1,2,1,4,2/),(/12,8/))*1.d0

           
           itmax = 2000
           np = size(d_ia)-1
           nx = size(a);ny=size(b);nz = size(c)
           n1 = nx*(ny+1)*(nz+1)
           n2 = (nx+1)*ny*(nz+1)
           n3 =(nx+1)*(ny+1)
           allocate(t(np),t1(np),t2(np),phi(np))
      
          do i =1,size(t_ja)
           !  xold(i) = x(i)/el(i)**2
              if(t_ja(i).gt.nl)then
                print *,i,t_ja(i)
                stop
              endif
          enddo
        !  print *,np,nl
           call zamux_new(nl,np, x, t2, t_mat, t_Ja, t_Ia)
           do i =1,np
              if(dc(i)) t2(i) = 0
!              if(abs(t2(i)).ge.1.d100)then
!                 print *,'t2 error',i
!                  stop
!              endif
           enddo
      
          t1  =0
          t = t2
          do i = 1,np
              if(dc(i)) t(i) = 0
          enddo

   
         phi = 0
         call pcg(d_mat,d_ia,d_ja,np,t,phi,tol_cor,itmax,iter1) 
         print *,'iters of pcg',iter1
        ! print *,np,maxval(emap(:,1)),maxval(emap(:,2))   
     
           do i = 1,size(x)        
             if(.not.(dc(emap(i,1)).and.dc(emap(i,2))))then
                x(i) = x(i)+(phi(emap(i,2))-phi(emap(i,1)))/el(i)
              endif
           !  x(i) = x(i)*el(i)**2
           enddo

           !x = xold
           if(output_div)then
           !! calculate the divergence of corrected electric field
            call zamux(np, x, t2, t_mat, t_Ja, t_Ia)
            print *,'after correction',vnorm(t2,np)
           !   stop
            write(str,'(i5)') icount
             open(54,file='divj_cc_'//trim(adjustl(str))//'.dat') 
            do i = 1,73
              do j = 1,73
               ind = 24*73**2+(i-1)*73+j 
               write(54,120) real(t2(ind)), imag(t2(ind))
              enddo           
            enddo
          !  do i = 1,73
          !     ind = 24*73**2+(i-1)*73+1
          !     write(54,120) imag(t2(ind:ind+72))
          !  enddo  
          close(54)
            allocate(ele_slice(72,72))
            !! calculate the electric field at the cell centers of the xoy-slice    
            do i = 1,72 
              do j = 1,72
                 ind = 24*72*73+(i-1)*72+j
                 ind2 = 72*73**2+24*72*73+(i-1)*73+j
                 ele_r = (real(x(ind))+real(x(ind+72)))/2
                 ele_i = (imag(x(ind))+imag(x(ind+72)))/2
                 ele_slice(i,j) = cmplx(ele_r,ele_i)
              enddo
            enddo
         !  print *,ele_slice(1,1)
         !  stop
            open(55,file = 'ele_ccdc_'//trim(adjustl(str))//'.dat')
            do i = 1,72
              write(55,121) real(ele_slice(i,:))
            enddo
            do i = 1,72
              write(55,121) imag(ele_slice(i,:))
            enddo
            close(55)  
            deallocate(ele_slice)
           endif
121 format(72(f9.2,2x))
120 format(2(e12.5,2x))
           deallocate(phi,t,t1,t2)
          ! deallocate(phir,phii,tr,ti,drmat)
    end subroutine
            


   subroutine  diagcsr(nnz,n, K_mat, Ia, Ja, M_j)
    !! subroutine for calculating the diagonal array of a CSR matrix
    implicit none
    integer(8), intent(in) :: nnz,Ia(n+1)
    integer, intent(in) ::  n, Ja(nnz)
    complex(8), intent(in) :: K_mat(nnz)
    complex(8),intent(out) :: M_j(n)
    integer(8) i,j,k,leni

   ! leni = size(ia)
   ! if(leni==n+1)then
      do i = 1,n
        do j = Ia(i),Ia(i+1)-1
            if (Ja(j)==i)then
                M_j(i) = K_mat(j)
            endif
        enddo
      enddo
    end subroutine

 subroutine ssorfac(nnz,n, A, Ia, Ja,omega, ssor, jlu, ju, ierr)
    ! subroutine for SSOR factorization of complex matrix in CSR format
    ! the output matrix is the same size as A, and omega is the relaxation
    ! factor(0~2)
    implicit none
    integer(8):: nnz,ia(*)
    real(8):: omega,omm
    complex(8) :: a(*),iw(n), ssor(*)
    integer(8) j,jlu(*),ju(*),ju0,js,jcol
    integer ii, ierr,i,n,ja(*)

    !omm =sqrt(omega*(2-omega))
    ju0 = n+2
   ! print *,ju0
    jlu(1) = ju0

!c main loop

       do 500 ii = 1, n
!c generating row number ii of L and U.

           do 100 j=ia(ii),ia(ii+1)-1

!c     copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.

              jcol = ja(j)
              if (jcol .eq. ii) then
                  ssor(ii) = a(j)
                  ju(ii)  = ju0
              else
                  ssor(ju0) = a(j)
                  jlu(ju0) = ja(j)
                  ju0 = ju0 +1
              endif
100        continue
          jlu(ii+1) = ju0

           if (ssor(ii) .eq. 0.0d0) goto 600
           ssor(ii) = 1.0d0/ssor(ii)
          !  print *,ii,ssor(ii)
500    continue
       ierr = 0
       return

!c     zero pivot :

600    ierr = ii

       return
   end subroutine


 subroutine zilu0(n, idd,a, ja, ia, alu, jlu, ju, iw, ierr)
    implicit none
    complex(8) a(*), alu(*), tl
    integer ja(*),ii,n,ierr,idd(*)
    integer(8) ia(*),ju(*),jlu(*),iw(*),i,ju0,js,jf,jm,jcol,j,jj,jrow,jw

    !c------------------ right preconditioner ------------------------------*
    !c                    ***   ilu(0) preconditioner.   ***                *
    !c----------------------------------------------------------------------*
    !c Note that this has been coded in such a way that it can be used
    !c with pgmres. Normally, since the data structure of the L+U matrix is
    !c the same as that the A matrix, savings can be made. In fact with
    !c some definitions (not correct for general sparse matrices) all we
    !c need in addition to a, ja, ia is an additional diagonal.
    !c ILU0 is not recommended for serious problems. It is only provided
    !c here for comparison purposes.

    ju0 = n+2
    jlu(1) = ju0
    !c initialize work vector to zero's

    do 31 i=1, n
        iw(i) = 0
31  continue
    !c main loop
    do 500 ii = 1, n
        js = ju0

        !c generating row number ii of L and U.

        do 100 j=ia(ii),ia(ii+1)-1
        !c     copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.

            jcol = ja(j)
           !! for omp
           IF(abs(jcol-ii-idd(1)).le.idd(2).and.jcol.le.idd(3).and.jcol.ge.idd(1))then
            if (jcol .eq. ii) then
           ! IF(abs(jcol-ii).le.idd(4)-idd(3))then ! for serial
           !   if(jcol.eq.ii)then
                alu(ii) = a(j)
              !  print *,alu(ii)
                iw(jcol) = ii
                ju(ii)  = ju0
            !    if(a(j).eq.0) print *,ii
              else
                alu(ju0) = a(j)
               ! jlu(ju0) = ja(j)
                jlu(ju0) = ja(j)
                iw(jcol) = ju0
                ju0 = ju0+1
            !    if(ju0.ge.33*idd(3)) print *,ii,jcol
             endif
           endif
100     continue
        jlu(ii+1) = ju0
        jf = ju0-1
        jm = ju(ii)-1
        
        !c     exit if diagonal element is reached.
        do 150 j=js, jm
            jrow = jlu(j)
            tl = alu(j)*alu(jrow)
            alu(j) = tl
           ! if(ii.eq.10) print *,alu(j)            
            !c     perform  linear combination

            do 140 jj = ju(jrow), jlu(jrow+1)-1
                jw = iw(jlu(jj))
                if (jw .ne. 0) alu(jw) = alu(jw) - tl*alu(jj)
140         continue
150     continue

        !c     invert  and store diagonal element.

        if (alu(ii) .eq. 0.0d0) then
         !   print *,ii,alu(ii)          
            goto 600
        endif
        alu(ii) = 1.0d0/alu(ii)
       ! if(ii.eq.2162) print *,alu(ii)
        !c     reset pointer iw to zero

        iw(ii) = 0
        do 201 i = js, jf
201     iw(jlu(i)) = 0
500 continue
    ierr = 0
    !ju(idd(3)) = ju0
    return

    !c     zero pivot :

600 ierr = ii
    return

  end subroutine

   subroutine slusol(n,y,x,alu,jlu,ju)
      complex(8) :: x(n), y(n), alu(*),sum
      integer(8) :: k, jlu(*), ju(*)
      integer i,n

    do 40 i = 1, n
        x(i) = y(i)
       do 41 k = jlu(i),ju(i)-1
          x(i) = x(i)- alu(k)* x(jlu(k))
41     continue
40  continue

    do 90 i = n, 1, -1
          do 91 k=ju(i),jlu(i+1)-1
             x(i) = x(i)- alu(k)*x(jlu(k))
91        continue
          x(i) = alu(i)*x(i)
90     continue
         return

   end subroutine


   subroutine pcg(k_mat,ia,ja,nl,b,u,tol,maxit,iter)
     !! Jacobi-CG
     implicit none
     integer(8),intent(in)::ia(nl+1)
     integer, intent(in) :: Ja(*)
     integer, intent(in) :: nl, maxit
     complex(8), intent(in) :: K_mat(*)
     complex(8):: r(nl),rold(nl),p(nl),v(nl),b(nl),u(nl)
     real(8), intent(in) :: tol
     integer,intent(out) ::iter
     integer i,j,k, n_less, ian(nl+1)
     integer(8) nnz
     complex(8) alpha,beta,m_j(nl)
     real(8):: error,babs,res

     nnz = ia(nl+1)-1
    ! print *,nnz
    ! n_less = nnz; ian= ia
     call  diagcsr(nnz,nl, K_mat, Ia, Ja, M_j)
   ! diagcsr(nnz,n, K_mat, Ia, Ja, M_j) 
    !call ddiag(nnz,nl, K_mat, Ia, Ja, M_j)
     !m_j = 1.d0
   !  do i = 1,nl
   !    if(m_j(i).eq.0)then
   !       print *,'in pcg',i
   !      stop
   !    endif
   !  enddo

     call zamux(nl,u,r,k_mat,ja,ia)
     r = b-r
     babs = vnorm(b,nl)
     error = vnorm(r,nl)/babs
     k = 0
     p = r/m_j

     do i = 1,maxit
       rold = r
      ! if(mod(i,100).eq.0)
      ! print *,i,error
       if(error<tol.and.error.gt.0)then
          iter = k
          res = error
         ! print *,iter
          exit
       endif
       call zamux(nl,p,v,k_mat,ja,ia)
       alpha = dot_product(r/m_J,r)/dot_product(p,v)
       u = u+alpha*p
       r = r-alpha*v
       beta = dot_product(r/m_j,r)/dot_product(rold/m_J,rold)
       p = r/m_j + beta*p
       error = vnorm(r,nl)/babs
       k = k+1
     !  print *,K,ERROR
     !  stop
     end do
     iter = k

    end subroutine


    subroutine damux(transp,n,x,y,a,ja,ia)
      real*8 a(*),  x(n), y(n), tmp
      integer i,n,ja(*),transp,k,ia(*)

     ! print *,size(x),size(y)
     ! print *,n
      if(transp.eq.0)then
        do i = 1,n
           tmp = 0
           do k = ia(i),ia(i+1)-1
              tmp =tmp+a(k)*x(ja(k))
           enddo
           y(i) =tmp
        enddo
      else
        ! y = 0
        do i = 1,n
           do k=ia(i), ia(i+1)-1
            y(ja(k)) = y(ja(k)) + x(i)*a(k)            
           enddo
         !  print *,i,ia(i),ja(k)
        enddo
      endif
       
    end subroutine
 
   subroutine zamux (n, x, y, a,ja,ia)
      complex(8)  x(*), y(*), a(*), t
      integer n, ja(*)      
      integer i
      integer(8) k, ia(*) 

      do 100 i = 1,n

!     compute the inner product of row i with vector x
         t = 0.0d0
         do 99 k=ia(i), ia(i+1)-1
            t = t + a(k)*x(ja(k))
 99      continue

!     store result in y(i)
         y(i) = t
 100  continue
      return
    end subroutine

   subroutine zamux_new(n1,n2, x, y, a, Ja, Ia)
       complex(8)  x(n1), y(n2), a(*), t
      integer n1,n2, ja(*)
      integer i
      integer(8) k, ia(n2+1)

      do 100 i = 1,n2

!     compute the inner product of row i with vector x
         t = 0.0d0
         do 99 k=ia(i), ia(i+1)-1
            t = t + a(k)*x(ja(k))
 99      continue

!     store result in y(i)
         y(i) = t
 100  continue
      return

   end subroutine

   subroutine zatmux (n, x, y, a, ja, ia)
      complex(8) x(*), y(*), a(*)
      integer n, ja(*)
      integer i
      integer(8) k, ia(*)
!c-----------------------------------------------------------------------
!c     zero out output vector

      do 1 i=1,n
         y(i) = 0.0
 1    continue

!c loop over the rows

      do 100 i = 1,n
         do 99 k=ia(i), ia(i+1)-1
            y(ja(k)) = y(ja(k)) + x(i)*conjg(a(k))
 99      continue
 100  continue

      return
   end subroutine

   subroutine fddc_assemb(a,b,c,nair,sigma,sigma_anm,pmat,pia,pja,tmat,tia,tja,el,&
          emap,dc,ep)
    implicit none
    real(8) a(:),b(:),c(:)
    complex(8),allocatable:: v(:),pmat(:),tmat(:),ep(:),pp(:),v2(:),tmat2(:)
    real*8,allocatable::sigma(:,:),el(:),sigma_anm(:,:) !length of each edge
    integer(8),allocatable::tia(:),pia(:)
    integer,allocatable ::vr(:),vc(:),tja(:),pja(:),me(:,:),me2(:,:),emap(:,:)
    !real*8
    !::k11(8,8),k22(8,8),k33(8,8),t1(8,4),t2(8,4),t3(8,4),ke(8,8),te(8,12)
    real*8 hx,hy,hz,sigmas(6),ve
    integer i,j,k,nx,ny,nz,nl,np,ne,xnl,ynl,n1,nj,nk,lmap(12,2),ie,ind,nair,idx(6),eid(12),&
          nid1(12),nid2(12)
    integer(8) nnz,nnz2,ii
    logical,allocatable :: dc(:)   

    idx = (/3,4,2,5,1,6/)
    nx =size(a); ny =size(b); nz = size(c)
    np = (nx+1)*(ny+1)*(nz+1)
    xnl = nx*(ny+1)*(nz+1)
    ynl = (nx+1)*ny*(nz+1)
    nl = xnl+ynl+(nx+1)*(ny+1)*nz
    n1 = nx*ny
    allocate(el(nl),emap(nl,2))
    allocate(v(np*7),vr(np*7),vc(np*7),dc(np))
    v = 0; vr = 0; vc = 0
    dc = .false.
    nnz = 0
    do k = 1,nz+1
      do j =1,ny+1
        do i =1,nx+1
           ind = (k-1)*(nx+1)*(ny+1)+(j-1)*(nx+1)+I
           ie = (k-1)*nx*ny+(j-1)*nx+i
           if(k.ne.1.and.k.ne.nz+1.and.j.ne.1.and.j.ne.ny+1.and.i.ne.1.and.i.ne.nx+1)then
             !! calculate the average conductivity first
           !  sigmas =(sigma(ie,1)*a(i)*b(j)*c(k)+sigma(ie-1,1)*a(i-1)*b(j)*c(k)+sigma(ie-nx,1)*a(i)*b(j-1)*c(k)+&
           !  sigma(ie-nx-1,1)*a(i-1)*b(j-1)*c(k)+sigma(ie-n1,1)*a(i)*b(j)*c(k-1)+sigma(ie-n1-1,1)*a(i-1)*b(j)*c(k-1)+&
           !  sigma(ie-n1-nx,1)*a(i)*b(j-1)*c(k-1)+sigma(ie-n1-nx-1,1)*a(i-1)*b(j-1)*c(k-1))/((a(i-1)+a(i))*&
            ! (b(j-1)+b(j))*(c(k-1)+c(k)))
            !order: bottom-front-left-right-back-top
           ! if(k.le.nair+1)then           
           !   sigmas = (/(a(i-1)+a(i))*(b(j-1)+b(j))/c(k-1), (a(i-1)+a(i))*(c(k-1)+c(k))/b(j-1),&
           !            (b(j-1)+b(j))*(c(k-1)+c(k))/a(i-1), (b(j-1)+b(j))*(c(k-1)+c(k))/a(i), &
           !            (a(i-1)+a(i))*(c(k-1)+c(k))/b(j),(a(i-1)+a(i))*(b(j-1)+b(j))/c(k)/)/4
           ! else
             ve = ((a(i-1)+a(i))*(b(j-1)+b(j))*(c(k-1)+c(k)))
             sigmas = (/(sigma(ie-n1-nx-1,3)*a(i-1)*b(j-1)+sigma(ie-n1-nx,3)*a(i)*b(j-1)+&
                     sigma(ie-n1-1,3)*a(i-1)*b(j)+sigma(ie-n1,3)*a(i)*b(j))/c(k-1)*2,&
                     (sigma(ie-n1-nx-1,2)*a(i-1)*c(k-1)+sigma(ie-n1-nx,2)*a(i)*c(k-1)+&
                      sigma(ie-nx-1,2)*a(i-1)*c(k)+sigma(ie-nx,2)*a(i)*c(k))/b(j-1)*2,&
                     (sigma(ie-n1-nx-1,1)*b(j-1)*c(k-1)+sigma(ie-n1-1,1)*b(j)*c(k-1)+&
                      sigma(ie-nx-1,1)*b(j-1)*c(k)+sigma(ie-1,1)*b(j)*c(k))/a(i-1)*2,&
                     (sigma(ie-n1-nx,1)*b(j-1)*c(k-1)+sigma(ie-n1,1)*b(j)*c(k-1)+&
                      sigma(ie-nx,1)*b(j-1)*c(k)+sigma(ie,1)*b(j)*c(k))/a(i)*2,&
                     (sigma(ie-n1-1,2)*a(i-1)*c(k-1)+sigma(ie-n1,2)*a(i)*c(k-1)+&
                      sigma(ie-1,2)*a(i-1)*c(k)+sigma(ie,2)*a(i)*c(k))/b(j)*2,&
                     (sigma(ie-nx-1,3)*a(i-1)*b(j-1)+sigma(ie-nx,3)*a(i)*b(j-1)+&
                     sigma(ie-1,3)*a(i-1)*b(j)+sigma(ie,3)*a(i)*b(j))/c(k)*2  /) 
           ! endif
           !  v((ind-1)*7+1:ind*7) = (/-2/(c(k-1)+c(k))/c(k-1),-2/(b(j-1)+b(j))/b(j-1), -2/(a(i-1)+a(i))/a(i-1),&
           !                       2/(c(k-1)*c(k))+2/(b(j-1)*b(j))+2/(a(i-1)*a(i)),&
           !                       -2/(a(i-1)+a(i))/a(i),-2/(b(j-1)+b(j))/b(j),-2/(c(k-1)+c(k))/c(k)/)*sigmas*ve     
             v((ind-1)*7+1:(ind-1)*7+3) = -sigmas(1:3);v((ind-1)*7+5:(ind-1)*7+7) = -sigmas(4:6);
             v((ind-1)*7+4) = sum(sigmas)
             vr((ind-1)*7+1:ind*7) = ind
             vc((ind-1)*7+1:ind*7) =(/ind-(nx+1)*(ny+1),ind-nx-1,ind-1,ind,ind+1,ind+nx+1,ind+(nx+1)*(ny+1)/)
            nnz = nnz + 7
         else
            nnz = nnz+1
            v((ind-1)*7+1) = 1.d0
            vr((ind-1)*7+1) = ind
            vc((ind-1)*7+1) = ind
            dc(ind) = .true.
           endif
        ENDDO
        enddo
    enddo

    allocate(pmat(nnz),pia(np+1),pja(nnz))
    k =0
    pia(1) = 1
  !  open(34,file='pmat.dat')
    do i = 1,np
          do j = 1,7
               ind = 7*int((i-1),8)+j
              ! print *,ind
               if(vc(ind).ne.0.and.vr(ind).ne.0)then
                    k = k +1
                    pmat(k) = v(ind)
                    pja(k) = vc(ind)
   !                 write(34,659) i,vc(ind),real(v(ind)),imag(v(ind)) 
                endif
          enddo
          pia(i+1) = int(k,8)+1
     enddo
    deallocate(v,vr,vc)
  !  close(34)
  !   stop 
    !! discretization of the divergence operator
    allocate(v(np*6),vr(np*6),vc(np*6))
    nnz = 0
    v = 0; vr = 0; vc = 0
    do k = 1,nz+1
      do j =1,ny+1
        do i =1,nx+1         
            ind = (k-1)*(nx+1)*(ny+1)+(j-1)*(nx+1)+I
            ie = (k-1)*nx*ny+(j-1)*nx+i
           if(k.ne.1.and.k.ne.nz+1.and.j.ne.1.and.j.ne.ny+1.and.i.ne.1.and.i.ne.nx+1)then     
             !! calculate the average conductivity first
           !   sigmas =(sigma(ie,1)*a(i)*b(j)*c(k)+sigma(ie-1,1)*a(i-1)*b(j)*c(k)+sigma(ie-nx,1)*a(i)*b(j-1)*c(k)+&
           !   sigma(ie-nx-1,1)*a(i-1)*b(j-1)*c(k)+sigma(ie-n1,1)*a(i)*b(j)*c(k-1)+sigma(ie-n1-1,1)*a(i-1)*b(j)*c(k-1)+&
           !   sigma(ie-n1-nx,1)*a(i)*b(j-1)*c(k-1)+sigma(ie-n1-nx-1,1)*a(i-1)*b(j-1)*c(k-1))/((a(i-1)+a(i))*&
           !   (b(j-1)+b(j))*(c(k-1)+c(k)))
            ve = (a(i-1)+a(i))*(b(j-1)+b(j))*(c(k-1)+c(k))
          !  v((ind-1)*6+1:ind*6) =(/-2/(a(i-1)+a(i)),2/(a(i-1)+a(i)),-2/(b(j-1)+b(j)),2/(b(j-1)+b(j)),&
           ! -2/(c(k-1)+c(k)),2/(c(k-1)+c(k))/)*sigmas*ve  
           ! if(k.le.nair+1)then
           !   sigmas = (/-(a(i-1)+a(i))*(b(j-1)+b(j)),  -(a(i-1)+a(i))*(c(k-1)+c(k)),&
           !            -(b(j-1)+b(j))*(c(k-1)+c(k)),  (b(j-1)+b(j))*(c(k-1)+c(k)), &
           !            (a(i-1)+a(i))*(c(k-1)+c(k)), (a(i-1)+a(i))*(b(j-1)+b(j))/)/4
           ! else   
        !!the divergence of secondary field
             sigmas =(/-(sigma(ie-n1-nx-1,3)*a(i-1)*b(j-1)+sigma(ie-n1-nx,3)*a(i)*b(j-1)+&
                     sigma(ie-n1-1,3)*a(i-1)*b(j)+sigma(ie-n1,3)*a(i)*b(j))*2,&
                     -(sigma(ie-n1-nx-1,2)*a(i-1)*c(k-1)+sigma(ie-n1-nx,2)*a(i)*c(k-1)+&
                      sigma(ie-nx-1,2)*a(i-1)*c(k)+sigma(ie-nx,2)*a(i)*c(k))*2,&
                     -(sigma(ie-n1-nx-1,1)*b(j-1)*c(k-1)+sigma(ie-n1-1,1)*b(j)*c(k-1)+&
                      sigma(ie-nx-1,1)*b(j-1)*c(k)+sigma(ie-1,1)*b(j)*c(k))*2,&
                     (sigma(ie-n1-nx,1)*b(j-1)*c(k-1)+sigma(ie-n1,1)*b(j)*c(k-1)+&
                      sigma(ie-nx,1)*b(j-1)*c(k)+sigma(ie,1)*b(j)*c(k))*2,&
                     (sigma(ie-n1-1,2)*a(i-1)*c(k-1)+sigma(ie-n1,2)*a(i)*c(k-1)+&
                      sigma(ie-1,2)*a(i-1)*c(k)+sigma(ie,2)*a(i)*c(k))*2,&
                     (sigma(ie-nx-1,3)*a(i-1)*b(j-1)+sigma(ie-nx,3)*a(i)*b(j-1)+&
                     sigma(ie-1,3)*a(i-1)*b(j)+sigma(ie,3)*a(i)*b(j))*2  /)
           ! endif     
            v((ind-1)*6+1:ind*6) = sigmas(idx)
           ! if(ind.eq.37690) print *,ie,sigmas!sigma(ie-n1-nx-1,3),sigma(ie-nx-1,3)
          !!the divergence of primary field
!            sigmas =(/-(sigma_anm(ie-n1-nx-1,3)*a(i-1)*b(j-1)+sigma_anm(ie-n1-nx,3)*a(i)*b(j-1)+&
!                     sigma_anm(ie-n1-1,3)*a(i-1)*b(j)+sigma_anm(ie-n1,3)*a(i)*b(j))*2,&
!                     -(sigma_anm(ie-n1-nx-1,2)*a(i-1)*c(k-1)+sigma_anm(ie-n1-nx,2)*a(i)*c(k-1)+&
!                      sigma_anm(ie-nx-1,2)*a(i-1)*c(k)+sigma_anm(ie-nx,2)*a(i)*c(k))*2,&
!                     -(sigma_anm(ie-n1-nx-1,1)*b(j-1)*c(k-1)+sigma_anm(ie-n1-1,1)*b(j)*c(k-1)+&
!                      sigma_anm(ie-nx-1,1)*b(j-1)*c(k)+sigma_anm(ie-1,1)*b(j)*c(k))*2,&
!                     (sigma_anm(ie-n1-nx,1)*b(j-1)*c(k-1)+sigma_anm(ie-n1,1)*b(j)*c(k-1)+&
!                      sigma_anm(ie-nx,1)*b(j-1)*c(k)+sigma_anm(ie,1)*b(j)*c(k))*2,&
!                     (sigma_anm(ie-n1-1,2)*a(i-1)*c(k-1)+sigma_anm(ie-n1,2)*a(i)*c(k-1)+&
!                      sigma_anm(ie-1,2)*a(i-1)*c(k)+sigma_anm(ie,2)*a(i)*c(k))*2,&
!                     (sigma_anm(ie-nx-1,3)*a(i-1)*b(j-1)+sigma_anm(ie-nx,3)*a(i)*b(j-1)+&
!                     sigma_anm(ie-1,3)*a(i-1)*b(j)+sigma_anm(ie,3)*a(i)*b(j))*2  /)
!            v2((ind-1)*6+1:ind*6) = sigmas(idx)
 
            vr((ind-1)*6+1:ind*6) = ind
            vc((ind-1)*6+1:ind*6) =(/(k-1)*nx*(ny+1)+(j-1)*nx+i-1,(k-1)*nx*(ny+1)+(j-1)*nx+i,&
                                 xnl+(k-1)*(nx+1)*ny+(j-2)*(nx+1)+i,xnl+(k-1)*(nx+1)*ny+(j-1)*(nx+1)+i,&
                                xnl+ynl+(k-2)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i,xnl+ynl+(k-1)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i/)
          !  el(vc((ind-1)*6+1:ind*6)) = (/a(i-1),a(i),b(j-1),b(j),c(k-1),c(k)/)
          !  emap(vc((ind-1)*6+1:ind*6),1)=(/ind-1,ind,ind-nx-1,ind,ind-(nx+1)*(ny+1),ind/)
          !  emap(vc((ind-1)*6+1:ind*6),2)=(/ind,ind+1,ind,ind+nx+1,ind,ind+(nx+1)*(ny+1)/)
           nnz = nnz+6
          else
            nnz = nnz+1
            v((ind-1)*6+1) = 0
            vr((ind-1)*6+1) = ind
            vc((ind-1)*6+1) = ind
         endif
        enddo
      enddo
    enddo

    ! assemble emap and el
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
          ! index of the edges
             eid = (/(k-1)*nx*(ny+1)+(j-1)*nx+i,(k-1)*nx*(ny+1)+j*nx+i,k*nx*(ny+1)+(j-1)*nx+i,k*nx*(ny+1)+j*nx+i, &
             xnl+(k-1)*(nx+1)*ny+(j-1)*(nx+1)+i,xnl+k*(nx+1)*ny+(j-1)*(nx+1)+i,xnl+(k-1)*(nx+1)*ny+(j-1)*(nx+1)+i+1,&
             xnl+k*(nx+1)*ny+(j-1)*(nx+1)+i+1,xnl+ynl+(k-1)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i,xnl+ynl+(k-1)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i+1,&
             xnl+ynl+(k-1)*(nx+1)*(ny+1)+j*(nx+1)+i,xnl+ynl+(k-1)*(nx+1)*(ny+1)+j*(nx+1)+i+1/) 
          ! index of the edge starting nodes
             nid1 =(/(k-1)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i,(k-1)*(nx+1)*(ny+1)+j*(nx+1)+i,k*(nx+1)*(ny+1)+(j-1)*(nx+1)+i,&
              k*(nx+1)*(ny+1)+j*(nx+1)+i,(k-1)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i,k*(nx+1)*(ny+1)+(j-1)*(nx+1)+i,&
              (k-1)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i+1, k*(nx+1)*(ny+1)+(j-1)*(nx+1)+i+1,(k-1)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i,&
              (k-1)*(nx+1)*(ny+1)+(j-1)*(nx+1)+i+1,(k-1)*(nx+1)*(ny+1)+j*(nx+1)+i,(k-1)*(nx+1)*(ny+1)+j*(nx+1)+i+1/)
          ! index of the edge ending nodes
             nid2(1:4) = nid1(1:4)+1
             nid2(5:8) = nid1(5:8)+nx+1
             nid2(9:12) = nid1(9:12)+(nx+1)*(ny+1)
             emap(eid,1) = nid1
             emap(eid,2) = nid2
             el(eid(1:4)) = a(i)
             el(eid(5:8)) = b(j) 
             el(eid(9:12)) = c(k)
          enddo
      enddo
    enddo
   
    allocate(tmat(nnz),tia(np+1),tja(nnz))
    k =0
    tia(1) = 1
   ! open(75,file='sfd_divmat.dat')
    do i = 1,np
          do j = 1,6
               ind = 6*int((i-1),8)+j
              ! print *,ind
               if(vc(ind).ne.0.and.vr(ind).ne.0)then
                    k = k +1
                    tmat(k) = v(ind)
     !               tmat2(k) = v2(ind)
     !               write(75,659) i,vc(ind),real(v(ind)),imag(v(ind))
                    tja(k) = vc(ind)
                endif
          enddo
          tia(i+1) = int(k,8)+1
     enddo
    ! print *,tmat(1)
    ! allocate(pp(np))
    ! deallocate(tmat2,v2)
    ! close(75)
    ! stop
    deallocate(v,vr,vc)
659 format(2(i9,2x),2(e9.2,2x))

  end subroutine
        
      subroutine cal_divj(kmat,ia,ja,n1, a,b,c,x,nzid)
      !! calculating div(J) for A-Phi method
      !! although div(A) = 0 is enforced, div(J) may not equal to 0 (across interfaces)
      !! can be called at any stage of the solving process
      !! the FE matrix is reused
      !! nzid decides the depth of the output data
        integer nx,ny,nz,nzid
        real*8,allocatable:: a(:),b(:),c(:)
        complex(8),allocatable :: kmat(:),t(:),divj(:,:)
        complex(8) x(:)
        integer(8) :: ia(:)
        integer :: ja(:)
        integer i,j,k,np,nnz,n1,n2,n3,ind

        nx =size(a); ny = size(b); nz =size(c)
        n2 = (nx-1)*(ny-1)
        n3 = size(ia)-1
        allocate(t(n1),divj(nx+1,ny+1))
        divj = 0
        k = ia(n3-n1+1); nnz=ia(n3+1)-1
        call zamux(n1,x(n3-n1+1:n3),t,kmat(k:nnz),ja(k:nnz),ia(n3-n1+1:n3+1)-k)

        open(23,file='aphi_3d1_divj.dat') 
        do i = 1,nx+1
           do j= 1,ny+1
              if(i.ne.1.and.i.ne.nx+1.and.j.ne.1.and.j.ne.ny+1)then
                 ind = nzid*n2+(j-2)*(nx-1)+i-1
                 divj(i,j) = t(ind)
                 write(23,1001) real(divj(i,j)), imag(divj(i,j))
              endif
           enddo
        enddo
1001 format(2(e12.5,2x))       
        close(23)
      stop
       deallocate(t,divj)

      end subroutine
      
end module
