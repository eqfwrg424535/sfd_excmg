! annotate the first two lines if you are building in Windows
! also use mkl_blas95 instead 
include 'mkl_pardiso.f90'
include 'mkl_spblas.f90'
module solvers
    use ifport
    use global_parameter
    use divergence_corr
    use mkl_pardiso
    use mkl_spblas
    
    implicit none
    
    real*8,parameter :: ssor_omg = 1.d0 ! relaxation factor for the SSOR factor, (tuning range: 0~2)
    contains
    ! subroutine for inspecting the memory usage, not recommended, one can use some 
    ! expertise tools like GridView for this purpose
      subroutine system_mem_usage(valueRSS)
      implicit none
      integer, intent(out) :: valueRSS

      character(len=200):: filename=' '
      character(len=80) :: line
      character(len=8)  :: pid_char=' '
      integer :: pid
      logical :: ifxst

      valueRSS=-1    ! return negative number if not found

       !--- get process ID

       pid=getpid()
       write(pid_char,'(I8)') pid
       filename='/proc/'//trim(adjustl(pid_char))//'/status'

       !--- read system file

       inquire (file=filename,exist=ifxst)
       if (.not.ifxst) then
         write (*,*) 'system file does not exist'
         return
       endif

       open(unit=100, file=filename, action='read')
       do
         read (100,'(a)',end=120) line
         if (line(1:6).eq.'VmRSS:') then
          read (line(7:),*) valueRSS
          exit
         endif
       enddo
120    continue
       close(100)

      return
   end subroutine system_mem_usage

     subroutine pcocr(nnz,k_mat,ia,ja,n,idd, x,b,eps,itmax,iter)
  !! Sogabe and Zhang, 2007, CAM
  implicit none
  complex(8),allocatable::k_mat(:),m_s(:),m_i(:),m_j(:)
  complex(8),allocatable::x(:),b(:),r(:),r1(:),z(:),s(:),p(:),a(:),w(:)
  complex(8) ::rho,rho_old,alpha,beta,aw

  real(8) eps,tolb,normr,err_rec(itmax)
  integer(8),allocatable :: jlu(:),ju(:),iw(:)
  integer(8)::nnz,ia(n+1)
  integer ja(nnz),n, itmax, i, k,ierr, idd(7)
  real iter(3)

     allocate(jlu(nnz),ju(n))
     allocate(r(n),p(n),r1(n),s(n))
     allocate(z(n),a(n),w(n))
    !  allocate(p(n),s(n),d(n))

    if (pre_type=='I') then
        allocate(M_i(nnz),iw(n*3))
        call zilu0(n, idd,k_mat, ja, ia, M_i, jlu, ju, iw, ierr)
        deallocate(iw)
     elseif(pre_type=='J') then
        allocate(M_j(n))
        call diagcsr(nnz,n, K_mat, Ia, Ja, M_j)
     elseif(pre_type=='S') then
        allocate(M_s(nnz))
        call ssorfac(nnz,n, K_mat, Ia, Ja, ssor_omg, M_s, jlu, ju, ierr)
     endif

       write(*,*) 'preconditioner LOADED',ierr
      call zamux (n, x, r, K_mat, Ja, Ia)
      r = b - r
      normr = vnorm(r,n)
      tolb =eps*vnorm(b,n)
      if(normr.le.tolb) return
       p = 0; a=0;w=0;beta = 0

      do i = 1,itmax
          if(normr.le.tolb)then
             iter(1)=i; iter(2) =vnorm(r,n)/vnorm(b,n)
             exit
          endif
          if(mod(i,100).eq.0) print *,i,vnorm(r,n)/vnorm(b,n)

          r1 = r
          p = r1+beta*p
        !  a = z+beta*a
          call zamux(n,p/m_j,s,K_mat, Ja, Ia)
          call zamux(n,r1/m_j,z,K_mat, Ja, Ia)
          rho = dot_product(r1,z)
          rho_old = rho
          if(pre_type == 'I') then
            call lusol(n,a,w,M_I,jlu,ju,ssor_omg)
          elseif(pre_type=='J') then
            w =a/m_j
          elseif(pre_type == 'S') then
             call lusol(n,a,w,M_s,jlu,ju,ssor_omg)
          endif
          aw = dot_product(s,s)
          alpha = rho/aw
          x = x+alpha*p
          r = r-alpha*s
        !  s = s-alpha*w
          call  zamux (n, r/m_j, z, K_mat, Ja, Ia)
          rho = dot_product(r,z)
          beta = rho/rho_old
          normr = vnorm(r,n)
      enddo


       if(allocated(M_i)) deallocate(M_i)
       if(allocated(M_S)) deallocate(M_S)
       if(allocated(M_J)) deallocate(M_J)
       deallocate(jlu,ju)
   deallocate(r,p,r1,w,z,s,a)

  end subroutine

subroutine pgpbicg(modelname,ff,polar,nnz, n, ax,by,cz,sig,sig_a,ep,idd, k_mat,&
 ia, ja, x, b, tol,maxit,iters)
!! precontioned gpbicg algorithm from Itoh(2022. JJIAM)
!!    gp stands for generalized product-type
!! several preconditioning techniques were proposed, the right
!! preconditioning is used in this subroutine, divergence correction
!! is also supported
  IMPLICIT NONE

  integer(8),intent(in) :: nnz,ia(n+1)
  integer, intent(in) :: Ja(nnz)
  integer(8),allocatable :: jlu(:),ju(:),iw(:),d_ia(:),t_ia(:)
  integer,allocatable::ia_sp(:),ja_sp(:),itmp(:),jtmp(:),ian(:),d_ja(:),t_ja(:),emap(:,:),eid(:)
  !integer :: Ibilut(n+1), Jbilut(nnz), maxfil
 integer, intent(in) :: n, maxit
  complex(8), intent(in) :: K_mat(nnz), b(n)
  complex(8),allocatable :: M_j(:),M_I(:),M_S(:),ep(:),pp(:),d_mat(:),t_mat(:)
  real*8,allocatable ::ax(:),by(:),cz(:),sig(:,:),sig_a(:,:),el(:)
  real(8), intent(in) :: ff,tol
  complex(8), intent(inout) :: x(n)
  real,intent(out) :: iters(3)

  integer i, J,m,k,half, ierr, info, vrss, stag, morestep, maxstagstep,&
              maxmstep, flag,nnz2,idd(7)
  integer(8) nmi
  integer,dimension(128) :: ipar
  real(8), dimension(128) :: dpar
  real(8), dimension(maxit) :: err_rec
  complex(8), allocatable :: r(:),rp(:),y(:),up(:),p(:),w(:), t(:),tp(:),&
      ap(:), at(:),u(:), z(:),r0(:),tmp(:),tpold(:),rold(:),jw(:)
  complex(8) rho, rhoold, alpha, beta, eta,  omega, one
  parameter (one = (1.d0,0.d0))
  real(8) babs, error, tolb, tau,tau0, theta, cmp1
  logical even,fag(n)
  logical,allocatable::dc(:) 
  character(len=30) modelname,str1,str2
  character(len=2) polar
  
  if(pre_type/='T') allocate(ju(n),jlu(nnz*20))
  allocate(r(n),rp(n),y(n),p(n),up(n),w(n),t(n), tp(n))
  allocate(ap(n),at(n),u(n),z(n),tmp(n),r0(n),rold(n))

  err_rec = 0
  if(div_cor)then
    call fddc_assemb(ax,by,cz,32,sig,sig_a,d_mat,d_ia,d_ja,t_mat,t_ia,t_ja,el,emap,dc,ep)
 endif

 if (pre_type=='I') then
    allocate(M_i(nnz*20),iw(n*2),jw(n+1))
    call zilu0(n,idd, k_mat, ja, ia, M_i, jlu, ju, iw, ierr)
 !   call zilut(n,idd,k_mat,ja,ia,n,1.d-4,m_i,jlu,ju,nnz*20,jw,iw,ierr)
    nmi = jlu(n+1)-1
    print *,'nnz in p',nmi
     deallocate(iw,jw)
elseif(pre_type=='J') then
 allocate(M_j(n))
    call diagcsr(nnz,n, K_mat, Ia, Ja, M_j)
elseif(pre_type=='S') then
    allocate(M_s(nnz))
    call ssorfac(nnz,n, K_mat, Ia, Ja, ssor_omg, M_s, jlu, ju, ierr)
endif

write(*,*) 'preconditioner LOADED',ierr
call system_mem_usage(vrss)
write(*,*) 'memory inspection in solver:', vrss

call zamux (n, x, r, K_mat, Ja, Ia)
r = b - r
r0 = r
!r = b - matmul(A, x)
! write(149,*) r
 !! calculate M^-1*u_m
   if(pre_type == 'I') then
        call lusol(n,r,rp,M_I(1:nmi),jlu(1:nmi),ju,ssor_omg)
    elseif(pre_type=='J') then
        rp = r/M_j
    elseif(pre_type == 'S') then
        call lusol(n,r,rp,M_s,jlu,ju,ssor_omg)
  endif

 t = 0; w = 0; beta = 0
 p = 0; u = 0; tp = 0; z= 0
 print *,"initial residual",vnorm(r0,n)
 babs = vnorm(b,n)
 error = vnorm(r0,n)/vnorm(b,n)
  do i = 1,maxit
    if(mod(i,corr_intv).eq.0)then
       ! print *,i,error
       if(div_cor)then
          call div_correc(n,ax,by,cz,x,pp,sig_a,el,emap,dc,d_mat,d_ia,d_ja,t_mat,t_ia,&
          t_ja,corr_tol,i/corr_intv)
          !! update r,r0,rp accordingly
         call zamux (n, x, r,K_mat, Ja, Ia)
         r = b - r
         r0 = r
       !  print *,'after dc:',vnorm(r,n)/babs

       endif
     endif

     tpold = tp
     rold = r
     p = rp+beta*(p-u)
     call zamux(n, p, ap, K_mat, Ja, Ia)
     alpha = dot_product(r0,r)/dot_product(r0,ap)
     y = t-r-alpha*w+alpha*ap
     t = r- alpha*ap
    if(pre_type == 'I') then
        call lusol(n,ap,tmp,M_I(1:nmi),jlu(1:nmi),ju,ssor_omg)
    elseif(pre_type=='J') then
        tmp = ap/M_j
    elseif(pre_type == 'S') then
        call lusol(n,ap,tmp,M_s,jlu,ju,ssor_omg)
    endif
     tp= rp - alpha*tmp
     call zamux(n, tp, at, K_mat, Ja, Ia)
     if(i==1)then
       omega = dot_product(at,t)/dot_product(at,at)
       eta = 0
     else
       omega = (dot_product(y,y)*dot_product(at,t)-dot_product(y,t)*dot_product(at,y))/&
             (dot_product(at,at)*dot_product(y,y)-dot_product(y,at)*dot_product(at,y))
       eta = (dot_product(at,at)*dot_product(y,t)-dot_product(y,at)*dot_product(at,t))/&
             (dot_product(at,at)*dot_product(y,y)-dot_product(y,at)*dot_product(at,y))
     endif
     u = omega*tmp+eta*(tpold-rp+beta*u)
     z = omega*rp+eta*z-alpha*U
     x = x + alpha*p + z
     !! convergence check
    ! call zamux (n, x, r,K_mat, Ja, Ia)
    ! r = b - r
    ! err_rec(i) = vnorm(r,n)/babs
     r = t - eta*y - omega*at
     err_rec(i) = vnorm(r,n)/babs
    ! print *,i,err_rec(i)
     if(vnorm(r,n)/babs.le.tol)then
         print *,'solver convergence achieved',I
         iters(1) = i
         exit
     endif

     if(pre_type == 'I') then
        call lusol(n,r,rp,M_I(1:nmi),jlu(1:nmi),ju,ssor_omg)
    elseif(pre_type=='J') then
        rp = r/M_j
    elseif(pre_type == 'S') then
        call lusol(n,r,rp,M_s,jlu,ju,ssor_omg)
    endif
    beta = alpha/omega*dot_product(r0,r)/dot_product(r0,rold)
    w = at+beta*ap
  enddo

 if(i.eq.maxit) iters(1) = maxit
 if(output_ccv)then
      write(str1,'(e9.2)') ff
     if(fwd_acc.ne.1)then
        write(str2,*) 'gpbicg' 
     else
        write(str2,*) 'excmg'   
    endif
    open(15,file=trim(adjustl(modelname))//'_'//trim(adjustl(str1))//'_'//polar//'_'//trim(adjustl(str2))//'.ccv')
    write(15,*)  iters(1)
      do k = 1,maxit
         if(err_rec(k).ne.0)then
           write(15,'(e9.2)') err_rec(k)
         endif
      enddo
     close(15)
endif

if(allocated(M_i)) deallocate(M_i)
if(allocated(M_S)) deallocate(M_S)
if(allocated(M_J)) deallocate(M_J)

deallocate(r,r0,rp,t,tp,u,z)
deallocate(w,y,p,ap,at,tmp,tpold)
if(pre_type/='T') deallocate(jlu,ju)

 end subroutine

  subroutine ptfqmr(modelname,ff,polar,nnz,ax,by,cz,sig,sig_a,ep,idd,K_mat, Ia, Ja, b, x, eps, n, itmax, iters)
    !! preconditioned transpose-free QMR algorithm (Freund, 1993, SISC)
    implicit none
    integer(8), intent(in) :: nnz
    !character :: pre_type
    integer, intent(in) :: Ja(nnz)
    integer(8),allocatable :: ia(:),jlu(:),ju(:),iw(:),d_ia(:),t_ia(:)
    integer,allocatable ::ia_sp(:),ja_sp(:),itmp(:),jtmp(:),ian(:),d_ja(:),t_ja(:),emap(:,:),eid(:)
    !integer :: Ibilut(n+1), Jbilut(nnz), maxfil
    integer, intent(in) :: n, itmax
    complex(8), intent(in) :: K_mat(nnz), b(n)
    complex(8),allocatable :: M_j(:),M_I(:),M_S(:),m_sp(:),mtmp(:),ep(:),pp(:),d_mat(:),t_mat(:)
    real*8,allocatable ::ax(:),by(:),cz(:),sig(:,:),sig_a(:,:),el(:)
    real(8), intent(in) :: ff,eps
    complex(8), intent(inout) :: x(n)
    real,intent(out) :: iters(3)

    integer i, J,m,k,half, ierr, info, vrss, stag, morestep, maxstagstep,&
                maxmstep, flag,nnz2,idd(7),dccount
    integer(8) nmi
    integer,dimension(128) :: ipar
    real(8), dimension(128) :: dpar
    real(8), dimension(itmax) :: err_rec
    complex(8), allocatable :: r(:),y(:),u(:),w(:), d(:),az(:),&
        ay(:), ad(:),v(:),r0(:),tmp(:),jw(:)
    complex(8) rho, rhoold, alpha, beta, eta,  sigma, one
    parameter (one = (1.d0,0.d0))
    real(8) babs, error, tolb, tau,tau0, theta, cmp1
    logical even,fag(n)
    logical,allocatable::dc(:)
    character(len=30) modelname,str1,str2
    character(len=2) polar
    
    if(pre_type/='T') allocate(ju(n),jlu(nnz*20))
    allocate(r(n),y(n),u(n),ay(n),d(n),w(n))
    allocate(az(n),v(n),tmp(n),ad(n),r0(n))
    idd(1) = 1; idd(3) = n
    idd(2) =idd(4)

    err_rec = 0
     if(div_cor)then
       call fddc_assemb(ax,by,cz,32,sig,sig_a,d_mat,d_ia,d_ja,t_mat,t_ia,t_ja,el,emap,dc,ep)
    endif

    !print *,n,ia(n+1)
    ! calculate the precondition matrix M
    if (pre_type=='I') then
        allocate(M_i(nnz*20),iw(n*2),jw(n+1))
        call zilu0(n,idd, k_mat, ja, ia, M_i, jlu, ju, iw, ierr)
     !   call zilut(n,idd,k_mat,ja,ia,n,1.d-4,m_i,jlu,ju,nnz*20,jw,iw,ierr)
        nmi = jlu(n+1)-1
        print *,'nnz in p',nmi
         deallocate(iw,jw)
    elseif(pre_type=='J') then
        allocate(M_j(n))
        call diagcsr(nnz,n, K_mat, Ia, Ja, M_j)
    elseif(pre_type=='S') then
        allocate(M_s(nnz))
        call ssorfac(nnz,n, K_mat, Ia, Ja, ssor_omg, M_s, jlu, ju, ierr)
    endif

    write(*,*) 'preconditioner LOADED',ierr
    call system_mem_usage(vrss)
    write(*,*) 'memory inspection in solver:', vrss
    dccount = 0
110  call zamux (n, x, r, K_mat, Ja, Ia)
    r = b - r
    !r = b - matmul(A, x)
    ! write(149,*) r
     !! calculate M^-1*u_m
  if(pre_type == 'I') then
     call lusol(n,r,y,M_I(1:nmi),jlu(1:nmi),ju,ssor_omg)
  elseif(pre_type=='J') then
     y = r/M_j
  elseif(pre_type == 'S') then
     call lusol(n,tmp,v,M_s,jlu,ju,ssor_omg)
  endif

    w = v
    theta = 0; eta = 0
    d = 0 ; ad = d
    rho = dot_product(r0,u)
    tau =vnorm(u,n)
    tau0 = tau
    if(pre_type == 'I') then
        call lusol(n,tmp,v,M_I(1:nmi),jlu(1:nmi),ju,ssor_omg)
    elseif(pre_type=='J') then
         v = tmp/M_j
    elseif(pre_type == 'S') then
        call lusol(n,tmp,v,M_s,jlu,ju,ssor_omg)
    endif

    babs = vnorm(b,n)
    tolb = eps*babs
  !  print *,tolb,eps
     error = vnorm(r,n)/babs
    !print *,error
    do i = 1,itmax
        rhoold = rho
        sigma = dot_product(r0,v)
        alpha = rho/sigma
        ay = y-alpha*v
        if(mod(i,corr_intv).eq.0)then
             print *,i,error
            dccount = dccount+1
            iters(1) =dccount*corr_intv
            if(dccount.ge.itmax/corr_intv) goto 124 !in case of endless loop
            if(div_cor)then
               call div_correc(n,ax,by,cz,x,pp,sig_a,el,emap,dc,d_mat,d_ia,d_ja,t_mat,t_ia,&
            t_ja,corr_tol,i/corr_intv)
               goto 110 ! restart
            endif
        endif
        ! each iteration contains two steps
        do j = 1,2
           m = 2*(i-1)+j
           if(j.eq.1)then
              u = u-alpha*w
              d = y+(theta**2*eta/alpha)*d
             ! ad =w+sigma*ad
           else
             call zamux (n, ay, tmp, K_mat, Ja, Ia)
             if(pre_type == 'I') then
                call lusol(n,tmp,az,M_I(1:nmi),jlu(1:nmi),ju,ssor_omg)
             elseif(pre_type=='J') then
                az = tmp/M_j
             elseif(pre_type == 'S') then
                call lusol(n,tmp,az,M_s,jlu,ju,ssor_omg)
             endif
              u = u-alpha*az
              d = ay+(theta**2*eta/alpha)*d

           endif
           theta = vnorm(u,n)/tau
           cmp1 = 1/sqrt(1+theta**2)
           tau =tau*theta*cmp1
           eta = cmp1**2*alpha

           x = x + eta*d
         !  r = r - eta*ad
          if(j.eq.2)then
            call zamux (n, x, r, K_mat, Ja, Ia)
            r = b - r
           error = vnorm(r,n)/babs
           err_rec(i) =error
          endif
         ! if(tau*sqrt((m+1)*1.d0)<tolb)then
         if(error<eps)then
             iters(1) = dccount*corr_intv+i
             call zamux (n, x, r, K_mat, Ja, Ia)
             r = b-r
             error = vnorm(r,n)/vnorm(b,n)
             goto 124
          endif
      enddo

        rho = dot_product(r0,u)
        beta = rho/rhoold
        y = u+beta*ay
        call zamux (n, y, tmp, K_mat, Ja, Ia)
       if(pre_type == 'I') then
            call lusol(n,tmp,w,M_I(1:nmi),jlu(1:nmi),ju,ssor_omg)
        elseif(pre_type=='J') then
            w = tmp/M_j
        elseif(pre_type == 'S') then
            call lusol(n,tmp,w,M_s,jlu,ju,ssor_omg)
        endif
        v= w+beta*(az+beta*v)

    enddo

124 if(output_ccv)then
    write(str1,'(e9.2)') ff
     if(fwd_acc.ne.1)then
        write(str2,*) 'tfqmr' 
     else
        write(str2,*) 'excmg'   
    endif
  open(15,file=trim(adjustl(modelname))//'_'//trim(adjustl(str1))//'_'//polar//'_'//trim(adjustl(str2))//'.ccv')
     write(15,*) iters(1)
      do k = 1,itmax
         if(err_rec(k).ne.0)then
           write(15,'(e9.2)') err_rec(k)
         endif
      enddo
     close(15)
    endif

    if(allocated(M_i)) deallocate(M_i)
    if(allocated(M_S)) deallocate(M_S)
    if(allocated(M_J)) deallocate(M_J)
    if(allocated(M_sp)) deallocate(M_sp)
    if(allocated(ja_sp)) deallocate(ja_sp)
    if(allocated(ia_sp)) deallocate(ia_sp)
    if(allocated(mtmp)) deallocate(mtmp)
    if(allocated(jtmp)) deallocate(jtmp)
    if(allocated(itmp)) deallocate(itmp)

    deallocate(r,r0,y,u,d,ad)
    deallocate(v,ay,w,az,tmp)
    if(pre_type/='T') deallocate(jlu,ju)

    end subroutine

    
      subroutine Bicg_stab_p(modelname,ff,polar,ax,by,cz,sig,sig_a,ep,nnz, K_mat, Ia, Ja, &
       b, x, eps, n,idd, itmax, iters,l_mat)
    !---------------------------------------------------------------------------------------!
    !   this code is based on the algo given in following paper:
    !   H. A. VAN DER VORST,J.SCI.STAT.COMPUT, 1992
    !---------------------------------------------------------------------------------------!
    !  preconditioner has been added for fast convergence in MT3D

    implicit none
    integer(8), intent(in) :: nnz,Ia(n+1)
    !character :: pre_type
    integer, intent(in) :: Ja(nnz)
    integer(8),allocatable :: jlu(:),ju(:),iw(:),d_ia(:),t_ia(:),ivp(:),ilp(:),&
                      jlh(:),jh(:),li(:)
    integer,allocatable :: ian(:),ja_sp(:),ia_sp(:),jtmp(:),itmp(:),d_ja(:),&
        t_ja(:),emap(:,:),jmp(:),imp(:),me(:,:),me2(:,:),jlp(:),jvp(:),eid(:),lj(:)
    !integer :: Ibilut(n+1), Jbilut(nnz), maxfil
    integer, intent(in) :: n, itmax
    !complex(8), intent(in) :: K_mat(nnz)
    complex(8), allocatable,optional::l_mat(:)
    complex(8),allocatable ::k_mat(:),M_j(:),M_I(:),M_S(:),m_sp(:),ep(:),d_mat(:),t_mat(:),&
        mlp(:),mvp(:),mh(:),pp(:),w(:)
    real(8), intent(in) :: eps,ff
    real(8),allocatable ::el(:),ax(:),by(:),cz(:),sig(:,:),sig_a(:,:),mp(:)
   ! complex(8), intent(inout) :: x(n)
    real,intent(out) :: iters(3)

    integer i, k,half, ierr, info, vrss,nnsp,tb,te,cnt,nnz2,cr,tid,nth,&
            n1,n2,n3,idd(7),np,nl,ng(5)
    integer(8) j,asd,nmi,nmh,k1,k2
    integer,dimension(128) :: ipar
    real(8), dimension(128) :: dpar
    real(8), dimension(itmax) :: err_rec
    complex(8),dimension(n) :: b,r,rhat,p,y,s,z,vk,t,xmin,x
   ! complex(8), allocatable :: r(:),rhat(:),p(:),y(:),v(:), s(:),z(:),&
   !     t(:), trvec(:),xmin(:),w(:)
    complex(8) rho, rhoold, alpha, beta, omega,alpha_in,beta_in,cone,aii
    parameter(cone=(1.d0,0.d0))
    real(8) babs, error, tol, t_amux, t_lusol
    logical fag(n)
    logical,allocatable::dc(:)
    character(2) polar
    character(30) str1,modelname,str2
    character(1) gid

    np = idd(4); !nl = n-np

   allocate(ju(n),jlu(nnz))
   ! allocate(xmin())
   ! allocate(r(n),rhat(n),p(n),y(n),xmin(n))
   ! allocate(v(n),s(n),z(n),t(n))
    if(div_cor)then
       call fddc_assemb(ax,by,cz,32,sig,sig_a,d_mat,d_ia,d_ja,t_mat,t_ia,t_ja,el,emap,dc,ep)
    endif

      err_rec = 0
    !print *,nnz,n*3
   tol = eps
    call system_clock(count_rate = cnt)
    call system_clock(tb)

    ! calculate the precondition matrix M
    if (pre_type=='I') then
        allocate(M_i(nnz),iw(n*3))
        call zilu0(n,idd,k_mat,ja,ia,m_i,jlu,ju,iw,ierr)
        nmi = jlu(n+1)-2
        print *,ierr,nmi
     !   nmi = nnz
       deallocate(iw)
   elseif(pre_type=='J') then
        allocate(M_j(n))
        call diagcsr(nnz,n, K_mat, Ia, Ja, M_j)
    elseif(pre_type=='S') then
        allocate(M_s(nnz))
        call ssorfac(nnz,n, K_mat, Ia, Ja, ssor_omg, M_s, jlu, ju, ierr)
    endif

    call system_clock(te)
    write(*,*) 'preconditioner LOADED',ierr,(te-tb)/real(cnt)
    call system_mem_usage(vrss)
    write(*,*) 'memory inspection in solver:', vrss
 !   stop
     call zamux(n, x, r,K_mat, Ja, Ia)

     r = b - r

    !r = b - matmul(A, x)
    rhat = r
    rho = 1.d0
    alpha = 1.d0
    omega = 1.d0
    vk = 0
    p = 0
    babs = vnorm(b,n)
   ! write(*,*) babs
    error = vnorm(r,n)/babs
    write(*,*) 'initial relative residual', vnorm(r,n)/babs
    if(isnan(error))  return
    !print *,itmax
    do i=1, itmax
        iters(1) = i
        err_rec(i)=error
      !  print *,i,error
       if (error < eps) then
            if(i.eq.itmax) x = xmin 
            exit
        endif
        !! recording the best approximate solution
        if(mod(i,corr_intv)==0) then
           ! print *,i,error
           ! call cal_divj(k_mat,ia,ja,idd(4), ax,by,cz,x,16)
            if(div_cor) then
           !   continue
              call div_correc(n,ax,by,cz,x,ep,sig_a,el,emap,dc,d_mat,d_ia,d_ja,t_mat,t_ia,&
            t_ja,corr_tol,i/corr_intv)
              call zamux (n, x, r,K_mat, Ja, Ia)
              r = b - r
              rhat = r
           endif
           ! if(i/100.eq.1) xmin = x
            !err_rec(i/100) = error
           ! if(i/100.gt.1.and.err_rec(i/100).le.err_rec(i/100-1))then
            xmin = x
           ! endif
        endif
        rhoold = rho
        rho = dot_product(rhat, r)
        beta = (rho / rhoold) * (alpha / omega)
        p = r + beta * (p - omega * vk)
        ! call system_clock(tb)
        !! here add preconditioner
        if(pre_type == 'I') then
            call lusol(n,p,y,M_I(1:nmi),jlu(1:nmi),ju,ssor_omg)
        elseif(pre_type=='J') then
            y = p/M_j
        elseif(pre_type == 'S') then
            call lusol(n,p,y,M_s,jlu,ju,ssor_omg)
        endif
    
         call zamux(n, y, vk,K_mat, Ja, Ia)
        alpha = rho / dot_product(rhat, vk)
        s = r - alpha * vk

        if(pre_type == 'I') then
            call lusol(n,s,z,M_I(1:nmi),jlu(1:nmi),ju,ssor_omg)
        elseif(pre_type=='J') then
            z = s/M_j
        elseif(pre_type == 'S') then
            call lusol(n,s,z,M_s,jlu,ju,ssor_omg)
        endif
    
       call zamux(n,z,t,k_mat,ja,ia)
       ! print *,vnorm(s,n),vnorm(z,n)
        omega = dot_product(t, s) / dot_product(t, t)
        x = x + alpha * y  + omega * z
        r = s - omega * t
        !call hdpre_asd(ng,ff,ax,by,cz,sigma,me,me2,r,r,dc,mp,&
        !     imp,jmp,mlp,ilp,jlp,mvp,ivp,jvp)
       error = vnorm(r,n) /babs
      ! print *,vnorm(r,n),babs 
      ! stop 
      if(isnan(error)) exit
    enddo
   ! write(*,*) 'In bicgstab',i-1,'err_final', error

   if(output_ccv)then
     write(str1,'(e9.2)') ff
      if(fwd_acc.ne.1)then
        write(str2,*) 'bicgstab' 
     else
        write(str2,*) 'excmg'   
    endif
    ! write(gid,'(i1)') size(ax)
     open(14,file=trim(adjustl(modelname))//'_'//trim(adjustl(str1))//'_'//polar//'_'//trim(adjustl(str2))//'.ccv')
    ! open(14,file='sfd_dctest_none.dat')
     write(14,*) iters(1)
     do j = 1,itmax
        if(err_rec(j).ne.0)then
         write(14,'(e9.2)') err_rec(j)
        endif
      enddo
    close(14)
    endif
    !if(pre_type/='T') info = MKL_SPARSE_DESTROY(A)
    if(allocated(M_i)) deallocate(M_i)
    if(allocated(M_S)) deallocate(M_S)
    if(allocated(M_J)) deallocate(M_J)
    if(allocated(M_sp)) deallocate(M_sp)
    if(allocated(ja_sp)) deallocate(ja_sp)
    if(allocated(ia_sp)) deallocate(ia_sp)
   if(allocated(mh)) deallocate(mh,jlh,jh)
    if(allocated(eid)) deallocate(eid)
  !   print *,allocated(dc),allocated(d_mat)
   if(div_cor) deallocate(el,emap,dc,d_mat,d_ia,d_ja,t_mat,t_ia,t_ja)
   if(allocated(Mp)) deallocate(mp,imp,jmp,mvp,ivp,jvp,mlp,ilp,jlp,me,me2,dc)
    if(allocated(ju)) deallocate(jlu,ju)

  end subroutine

     subroutine lusol(n, y, x, alu, jlu, ju,w)
      complex(8) :: x(n), y(n), alu(*),sum
      integer(8) :: k, jlu(*), ju(*)
      integer i,n
      real(8) w
     ! character,optional::pty

    !c This routine solves the system (LU) x = y,
    !c given an LU decomposition of a matrix stored in (alu, jlu, ju)
    !c modified sparse row format
       do 40 i = 1, n
           x(i) = y(i)
         !  sum = 0
          do 41 k=jlu(i),ju(i)-1
            !print *,k,jlu(k)
            x(i) = x(i)- alu(k)* x(jlu(k))
41        continue
         ! x(i) = x(i)-sum
         if(pre_type=='S')  x(i) = x(i)*alu(i)*w
       !  if(present(pty)) x(i)=x(i)/alu(i)/w
40     continue

       if(pre_type=='S')then
           do i =1,n
             x(i) = -x(i)/alu(i)*(2-w)/w
           enddo
       endif

       do 90 i = n, 1, -1
          do 91 k=ju(i),jlu(i+1)-1
             x(i) = x(i)- alu(k)*x(jlu(k))
91        continue
        ! x(i) = x(i)-sum
         if(pre_type=='S')then
           x(i) = alu(i)*x(i)*w
         else
           x(i) = alu(i)*x(i)
         endif
90     continue
         return

     end subroutine

     
       subroutine pardiso_spd(a,ia,ja,b,nrhs,x,n,nnz)
    !use mkl_pardiso
    implicit none
    integer :: ia(n+1), ja(nnz)
    complex(8) :: a(nnz),b(n),x(n)
    integer::n,nrhs,nnz,vrss

    ! nrhs = 1 in this situation
    !.. internal solver memory pointer
    type(mkl_pardiso_handle), allocatable  :: pt(:)
    !.. all other variables
    integer maxfct, mnum, mtype, phase, error, msglvl
    integer error1
    integer, allocatable :: iparm( : )

    integer i, idum(1)
    complex(8) ddum(1)
    !.. fill all arrays containing matrix data.
    maxfct = 1
    mnum = 1

    !.. set up pardiso control parameter
    !..
    allocate( iparm ( 64 ) )

    iparm = 0

    iparm(1) = 1 ! no solver default
    iparm(2) = 3 ! fill-in reordering from metis
    iparm(4) = 0 ! no iterative-direct algorithm
    iparm(5) = 0 ! no user fill-in reducing permutation
    iparm(6) = 0 ! solution on the first n compoments of x
    iparm(8) = 9 ! numbers of iterative refinement steps
    iparm(10) = 13 ! perturbe the pivot elements with 1e-13
    iparm(11) = 1 ! use nonsymmetric permutation and scaling mps
 iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). try iparm(13) = 1 in case of inappropriate accuracy
    iparm(14) = 0 ! output: number of perturbed pivots
    iparm(18) = -1 ! output: number of nonzeros in the factor lu
    iparm(19) = -1 ! output: mflops for lu factorization
    iparm(20) = 0 ! output: numbers of cg iterations
    !iparm(60) = 2 !! OOC MODE

    error  = 0 ! initialize error flag
    msglvl = 0 ! 1 for print statistical information
    mtype  =  13 ! complex and nonsymmetric
    !mtype  = -2 ! real and symmetric indefinite
    !mtype  = 3 ! complex and structually symmetric

    !.. initiliaze the internal solver memory pointer. this is only
    ! necessary for the first call of the pardiso solver.

    allocate ( pt ( 64 ) )
    do i = 1, 64
        pt( i )%dummy =  0
    end do

    !.. reordering and symbolic factorization, this step also allocates
    ! all memory that is necessary for the factorization

    phase = 11 ! only reordering and symbolic factorization

    call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
        idum, nrhs, iparm, msglvl, ddum, ddum, error)

    !write(*,*) 'reordering completed ... '
    if (error /= 0) then
        write(*,*) 'the following error was detected: ', error
        goto 1000
    end if
    !write(*,*) 'number of nonzeros in factors = ',iparm(18)
    !write(*,*) 'number of factorization mflops = ',iparm(19)

    !.. factorization.
    phase = 22 ! only factorization
    call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
        idum, nrhs, iparm, msglvl, ddum, ddum, error)
    write(*,*) 'factorization completed ... ',max(iparm(15),iparm(16)+iparm(63))
    call system_mem_usage(vrss)
    write(*,*) 'memory inspection in pardiso:', vrss

    if (error /= 0) then
        write(*,*) 'the following error was detected: ', error
        goto 1000
    endif

    !.. back substitution and iterative refinement
    iparm(8) = 2 ! max numbers of iterative refinement steps
    phase = 33 ! only factorization
 !   do i=1, 2
     call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
            idum, 1, iparm, msglvl, b, x, error)
        write(*,*) 'solve completed ... '
        if (error /= 0) then
            write(*,*) 'the following error was detected: ', error
            goto 1000
        endif
 !   enddo

1000 continue
    !.. termination and release of memory
    phase = -1 ! release internal memory
    call pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
        idum, nrhs, iparm, msglvl, ddum, ddum, error1)

    if ( allocated( iparm ) )   deallocate( iparm )

    if (error1 /= 0) then
        write(*,*) 'the following error on release stage was detected: ', error1
    endif

    end subroutine

    
     subroutine pbicg(modelname,ff,polar,idd,ax,by,cz,sig,sig_a,e_dc,k_mat,ia,ja,b,&
     x,eps,n,itmax,iter)
    !! preconditioned bi-conjugate gradient method
    !! Jaccobi and ILU preconditioner are on option
     
    implicit none
    integer(8), intent(in) :: Ia(n+1)
    !character :: pre_type
  !  integer, intent(in) :: Ja(nnz)
    integer, intent(in) :: n, itmax
  !  complex(8), intent(in) :: K_mat(nnz), b(n)
  !  complex(8):: M_j(n),yt(n),y(n),z(n),zt(n),p(n),pt(n),q(n),qt(n),r(n),rt(n)
    complex(8),allocatable::k_mat(:),b(:),m_i(:),m_j(:),yt(:),y(:),z(:),zt(:),p(:),&
       pt(:),q(:),qt(:),r(:),rt(:),x(:),e_dc(:),pp(:),d_mat(:),t_mat(:)
    real(8), intent(in) :: ff,eps
     real(8),allocatable:: ax(:),by(:),cz(:),sig(:,:),sig_a(:,:),el(:)
  !  complex(8),intent(inout) :: x(n)
    complex(8):: rho,rho1,alpha,beta,ptq
    integer(8) nnz,nmi
    integer,allocatable:: ja(:),d_ja(:),t_ja(:),emap(:,:)
    integer(8),allocatable:: iw(:),ju(:),jlu(:),d_ia(:),t_ia(:)
    real,intent(out):: iter
    integer :: i,j,info,vrss,idd(7),ierr
    real(8), dimension(itmax) :: err_rec
    real(8) :: babs, error
    character(len=30) str,modelname
    character(len=2) polar
    !real(8),external :: vnorm
    logical,allocatable:: dc(:)

   nnz =ia(n+1)-1
   allocate(M_j(nnz),yt(n),y(n),z(n),zt(n),p(n))
   allocate(pt(n),q(n),qt(n),r(n),rt(n))
   ! err_rec = 0
    print *,itmax,nnz,n
    if (div_cor)then 
       call fddc_assemb(ax,by,cz,32,sig,sig_a,d_mat,d_ia,d_ja,t_mat,t_ia,t_ja,el,emap,dc,e_dc)
    endif

    !! loading the Jaccobi preconditioner
    if(pre_type=='J')then
       call diagcsr(nnz,n,k_mat,ia,ja,m_j)
    elseif(pre_type=='I')then  
       allocate(M_i(nnz),iw(n*3),ju(n),jlu(nnz))
      !  idd = (/1,n,n,0,0,0,0/)
        call zilu0(n,idd,k_mat,ja,ia,m_i,jlu,ju,iw,ierr)
        nmi = jlu(n+1)-1
    else
      allocate(M_i(nnz),ju(n),jlu(nnz))
       m_i = 0
       call ssorfac(nnz,n, K_mat, Ia, Ja, ssor_omg, M_i, jlu, ju, ierr) 
       nmi = jlu(n+1)-1
      ! print *,nmi
    endif
   ! call system_mem_usage(vrss)
   ! write(*,*) 'memory inspection in solver:', vrss

    !rho = 1
    call zamux(n, x, r, K_mat, Ja, Ia)
    r = b - r
    !r = b - matmul(A, x)
    rt = r
    rho = 1.d0
    alpha = 1.d0
    beta = 1.d0

    p = 0
    babs = vnorm(b,n)
   ! write(*,*) babs
    error = vnorm(r,n)
  ! write(*,*) 'initial error:',error/babs   

    do i = 1,itmax
       iter = i
       err_rec(i) = error/babs
      if(error<eps*babs)then
           print *,iter,error/babs
           exit
      endif
       if(mod(i,corr_intv)==0)then
         print *,i,error/babs
          if(div_cor) then
           !   continue
              call div_correc(n,ax,by,cz,x,pp,sig_a,el,emap,dc,d_mat,d_ia,d_ja,t_mat,t_ia,&
            t_ja,corr_tol,i/corr_intv)
              call zamux (n, x, r,K_mat, Ja, Ia)
              r = b - r
              rt = r
           endif
       endif

   !       print *,i,error/sqrt(babs)
      !    err_rec(i/100) = error
   !    endif

      !! the first several steps will change if another precond is used
       if(pre_type=='J')then
         y = r/m_j(1:n)
       else
        call lusol(n,r,y,M_i(1:nmi),jlu(1:nmi),ju,ssor_omg)
       endif
       z = y
       yt = rt
       if(pre_type=='J')then
        zt = yt/conjg(m_j(1:n))
       else
        call lusol(n,r,y,M_i(1:nmi),jlu(1:nmi),ju,ssor_omg)
     !   m_j = conjg(m_i(1:nmi))
     !   call slusol(n,yt,zt,M_j,jlu(1:nmi),ju)
       endif
       rho1= rho
       rho = dot_product(rt,z)
       beta = rho/rho1
       if(i==1) then
         p = z; pt =zt;
       else
         p = z+beta*p
         pt = zt + conjg(beta)*pt
       endif

       call zamux(n, p,q,k_mat,ja,ia)
       call zatmux(n,pt,qt,k_mat,ja,ia)
       ptq = dot_product(pt,q)
       alpha = rho/ptq
       x = x + alpha*p
       r = r - alpha*q
       rt = rt- conjg(alpha)*qt
       error = vnorm(r,n)

    enddo

    ! open(42,file = 'mcsem_sfd_ccv_')
     ! write(str,'(e9.2)') ff
  !   write(gid,'(i1)') size(ax)/16
    if(output_ccv)then
   open(14,file=trim(adjustl(modelname))//'_'//trim(adjustl(str))//'_'//polar//'.ccv')
   write(14,*) iter   
   do j = 1,itmax
        if(err_rec(j).ne.0)then
         write(14,'(e9.2)') err_rec(j)
        endif
      enddo
    close(14)
    endif
     if(pre_type=='I') deallocate(m_i,jlu,ju,iw)
     if(pre_type=='S') deallocate(m_i,jlu,ju)
     deallocate(m_j,z,zt,yt,y,p,pt,q,qt,r,rt)    
    end subroutine
    
end module
