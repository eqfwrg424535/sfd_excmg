! annotate the first line if you are building in Windows
include 'lapack.f90'
module model_selection
    !! module for calculating background field and generating resistivity models
    !! 1--COMMEMI3D-1(1':anis,han,2018); 2--COMMEMI3D-2; 3--DTM3.0(2014); 
    !! 4--DTM1.0(2007); 5--MT Junggar basin(2019); 6--the layered anis(Bai et al,2022)
     !currently, model 1' and 5 only suit for NEFEM
    use global_parameter
    use lapack95
    implicit none
   
    !external zcgesv
    contains
    
!! function for some simple topographies
    !function topo_func(x,y)
    !   real*8 x,y,topo_func,ll1,ll2,l3,sl,h,aas,bbs
    !
    !  select case(bath_type)
    !    case(1) !! pyramid
    !       ll1= bath(1); ll2 = bath(2);l3=bath(3); h = bath(4)
    !       aas = (ll2-ll1)/2; bbs =(l3-ll1)/2; sl =l3/ll2
    !      if(abs(x).le.ll1/2.and.abs(y).le.ll1/2)then
    !         topo_func = h
    !      elseif((abs(x).gt.ll1/2.or.abs(y).gt.ll1/2).and.&
    !              abs(x).le.ll2/2.and.abs(y).le.l3/2)then
    !           if(y.ge.-sl*x.and.y.le.sl*x)then
    !               topo_func = -h/aas*(x-ll1/2)+h
    !            elseif(sl*x.ge.-y.and.sl*x.le.y)then
    !               topo_func = -h/bbs*(y-ll1/2)+h
    !            elseif(y.ge.sl*x.and.y.le.-sl*x)then
    !               topo_func = h/aas*(x+ll1/2)+h
    !            else
    !               topo_func = h/bbs*(y+ll1/2)+h
    !            endif
    !      else
    !         topo_func = 0
    !      endif        
    !    case(2) !! smooth hill
    !      ll1= bath(1); ll2 = bath(1); h = bath(3)
    !      if(abs(x).le.ll1/2.and.abs(y).le.ll1/2)then
    !         topo_func = h*cos(pi*x/ll1)*cos(pi*y/ll1)
    !      else
    !         topo_func = 0
    !      endif
    !    case(3) !! volcano
    !       ll1= bath(1); ll2 = bath(2); h = bath(3)
    !       aas = (ll2-ll1)/2
    !      if(abs(x).le.ll1/2.and.abs(y).le.ll1/2)then
    !         topo_func = h
    !         if(sqrt(x**2+y**2).le.160.d0)then
    !           topo_func = h-sqrt(160.d0**2-x**2-y**2)
    !         endif
    !      elseif((abs(x).gt.ll1/2.or.abs(y).gt.ll1/2).and.&
    !              abs(x).le.ll2/2.and.abs(y).le.ll2/2)then
    !           if(y.ge.-x.and.y.le.x)then
    !               topo_func = -h/aas*(x-ll1/2)+h
    !            elseif(x.ge.-y.and.x.le.y)then
    !               topo_func = -h/aas*(y-ll1/2)+h
    !            elseif(y.ge.x.and.y.le.-x)then
    !               topo_func = h/aas*(x+ll1/2)+h
    !            else
    !               topo_func = h/aas*(y+ll1/2)+h
    !            endif
    !      else
    !           topo_func = 0
    !      endif
    !    case(4) !! half sphere
    !      topo_func = sqrt(50.d0**2-x**2-y**2)
    !  end select
    !end function

    subroutine read_cond(mname,a,b,c,n_air,conv)
      implicit none
      integer nx,ny,nz,n_air,i
      real*8,allocatable ::a(:),b(:),c(:),conv(:,:)      
      character(*) mname
      character(len=200) firstline

         open(35,file=mname//'.ws')
          read(35,*) firstline !not used
         read(35,*) nx,ny,nz,n_air
       !  n_air = 60
         allocate(a(nx),b(ny),c(nz+n_air))
         read(35,*) a
         read(35,*) b
         read(35,*) c(n_air+1:n_air+nz)
       !  print *,sum(a)
         allocate(conv(6,nx*ny*nz))
     
          do i = 1,nx*ny*nz
             read(35,*) conv(:,i)
          enddo
         close(35)
       !! add air layers  
        do i = 1,n_air
           c(i) = c(n_air+1)*lambda_v**(n_air-i)
        enddo

    end subroutine

   ! subroutine bgcond(conv,nx,ny,nz,patm)
   !!! subroutine for extracting background model from input
   !!! or defining the model by user
   !  real*8,allocatable::conv(:,:),patm(:)
   !  integer nx,ny,nz,ne,i,j,k,ie
   !
   !  allocate(patm(6))
   !  if(bgmsd)then
   !    patm =(/rho_air,1.d3,1.d0,2.d0,2.d0,5.d2/) !!from top to bottom
   !  else
   !    patm(1) = rho_air
   !    j  =1 
   !    do k = 1,nz
   !       ie = (k-1)*nx*ny+1        
   !       do i = 1,j
   !        if(conv(1,ie).ne.patm(i))then
   !          j = j+1
   !          patm(j) = conv(1,ie)
   !        endif
   !       enddo
   !    enddo
   !  endif
   !
   ! end subroutine

    subroutine load_sigma(PaTM,a,b,c,N_air,abr,conv,sigma,topo)
    implicit none
    integer,intent(in)::N_air,abr(:)
    real(8),allocatable,optional::topo(:)
    real(8) :: PaTM(:),conv(:,:),condc(3,3),asx,asy,asz
    INTEGER nx,ny,nz,i,j,h,M,L,N,K,NE,nab,ide,neair,sumsf
    real(8),allocatable :: a(:),b(:),c(:),aax(:),bby(:),ccz(:)
    !DATA abr /12,22,12,22,16,25/ !location of the anomalous body
    !real(8,intent(in),optional ::sig_in(:)
    real(8),allocatable,intent(out) :: sigma(:,:)

    nx= size(a); ny =size(b); nz=size(c)
    NE=Nx*Ny*Nz; neair = nx*ny*n_air
    nab = size(conv,1)
    allocate(sigma(ne,6))
    allocate(aax(nx),bby(ny),ccz(nz))
    do i = 1,nx
      aax(i) = -sum(a)/2+sum(a(1:i))-a(i)/2
    enddo
    do i = 1,ny
      bby(i) = -sum(b)/2+sum(b(1:i))-b(i)/2
    enddo
    do i = 1,nz
      ccz(i) = -sum(c(1:n_air))+sum(c(1:i))-c(i)/2
    enddo

    print *, nab,n_air,ne
    sigma = 0
         do k = 1,n_air
            do j = 1,ny
                do i = 1,nx
                    ide = (k-1)*nx*ny+(j-1)*nx+i
                    sigma(ide,1:3) = 1.d0/patm(1)
                enddo
            enddo
         enddo

         do k = n_air+1,nz
              do j = 1,ny
                 do i = 1,nx
                    ide = (k-1)*nx*ny+(j-1)*nx+i
                    sigma(ide,1:3) = 1.d0/patm(2)
                enddo
              enddo
          enddo
        ! correction for the Salt model
        ! if(model.eq.10) sigma(1064961:1081344,1:3) = 1.d0/500

        ! if(model.ne.10.and.model.ne.11)then
        !  do n = 1,nab
        !       condc = conten_cal(conv(n,:))
        !      ! print *,condc  
        !      ! if(n.eq.1.or.n.eq.3) print *,'rhozz',condc(3,3) 
        !        do k = abr((n-1)*6+5),abr(n*6)
        !          do j = abr((n-1)*6+3),abr((n-1)*6+4)
        !             do i = abr((n-1)*6+1),abr((n-1)*6+2)
        !               ide = (k-1)*nx*ny+(j-1)*nx+i
        !                sigma(ide,1) = condc(1,1) !sigma_xx
        !                sigma(ide,2) = condc(2,2) !sigma_yy
        !                sigma(ide,3) = condc(3,3) !sigma_zz
        !                sigma(ide,4) = condc(1,2) !sigma_xy
        !                sigma(ide,5) = condc(1,3) !sigma_xz
        !                sigma(ide,6) = condc(2,3) !sigma_yz
        !              ! if(ide.eq.3072) print *,sigma(ide,:)
        !              if(n.ne.1.and.disk)then
        !                asx = aax(i); asy = bby(j); asz = ccz(k) 
        !              if(sqrt((asx-850)**2+(asy)**2).gt.1.d3)then
        !                 sigma(ide,1:3) = 1.d0/patm(3)
        !               endif
        !              endif 
        !            enddo
        !          enddo
        !        enddo
        !  enddo
        !endif
        !   print *,sigma(90,1)
        !  stop
        sumsf = 0
          !! recover some elements for the bathymetry model
          if(present(topo))then
              n =1 !n=nab
              do k = abr((n-1)*6+5),abr(n*6)
                  do j = abr((n-1)*6+3),abr((n-1)*6+4)
                     do i = abr((n-1)*6+1),abr((n-1)*6+2)
                        ide = (k-1)*nx*ny+(j-1)*nx+i
                        asx = aax(i); asy = bby(j); asz = ccz(k)
                       ! if(asz.le.-topo_func(asx,asy))then
                      if(asz.ge.-topo((j-abr(3))*(abr(4)-abr(3)+1)+i+1-abr(1)))then
                        sigma(ide,1:3) = 1.d0/patm(1)
                        !  print *,ide,asx,asy,topo_func(asx,asy)
                         sumsf = sumsf+1
                        !  if( sigma(ide,1)) print *,asx,asy
                       !   print *,asz,sigma(ide,1)
                       ! else
                    !       print *,asx,asy,asz
                      endif
                     enddo
                  enddo
              enddo
          endif 
      !   print *,abr(1:6),sumsf
      !   stop 
        deallocate(aax,bby,ccz)
    end subroutine
    
  ! function for calculating the conductivity tensor(also contained in the python script)  
  !  function conten_cal(conv)
  !   implicit none
  !   real(8) conv(6),conten_cal(3,3)
  !   real(8) d1(3,3),d2(3,3),d3(3,3),dtmp(3,3),sigma_c(3,3)
  !
  !   sigma_c = 0.d0
  !  ! d1 = 0; d2 = 0; d3 = 0;  
  !   sigma_c(1,1) = 1.d0/conv(1)
  !   sigma_c(2,2) = 1.d0/conv(2)
  !   sigma_c(3,3) = 1.d0/conv(3)
  !   d1 = reshape((/cos(conv(4)),sin(conv(4)), 0.d0, &
  !          -sin(conv(4)),cos(conv(4)), 0.d0, &
  !          0.d0, 0.d0, 1.d0/),(/3,3/))
  !  if(euler_mode.eq.1)then
  !    d2 = reshape((/1.d0, 0.d0, 0.d0, & 
  !          0.d0, cos(conv(5)), sin(conv(5)), &
  !          0.d0, -sin(conv(5)),cos(conv(5))/),(/3,3/))
  !  else
  !     d2 = reshape((/cos(conv(5)),0.d0,sin(conv(5)), &
  !                  0.d0,1.d0,0.d0, &            
  !                 -sin(conv(5)),0.d0,cos(conv(5))/),(/3,3/))
  !   endif
  !    d3 = reshape((/cos(conv(6)),sin(conv(6)), 0.d0, &
  !          -sin(conv(6)),cos(conv(6)), 0.d0, &
  !          0.d0, 0.d0, 1.d0/),(/3,3/))       
  !   dtmp = matmul(d1,matmul(d2,d3))
  !!   print *,sigma_c
  !!   print *,matmul(sigma_c,dtmp)
  !   conten_cal = matmul(dtmp,matmul(sigma_c,transpose(dtmp)))
  !
  !  end function


   subroutine mt1d_anl(ff,conv,c,ex,hx,imp)
     !!calculating background field through analtical formulas
     implicit none
     real(8) ::ff,c(:),conv(:,:)     
     real(8),allocatable::aon(:,:),qp(:),qm(:)
    ! complex(8),allocatable::ex(:)
     complex(8),allocatable::ex(:),hx(:),sprod(:,:,:),xip(:),xim(:)
     real(8) w,ada,add,da12,a1,a2,axx,ayy,axy,qdq,dq
     complex(8) kp,km,sp,sm,cp,cm,S(4,4),s_m(4,4),g(2,2),&
           e0,m(4,4),detg,cbot(4),mc(4),f(4),tmp,imp(4),atmp(8)
     integer i,j,k,nz

     m = 0
     nz=size(c)
     allocate(ex(2*nz+2),hx(2*nz+2))
     allocate(sprod(4,4,nz+1),xip(nz+1),xim(nz+1))
     allocate(qp(nz+1),qm(nz+1))
     w = 2*pi*ff
    ! print *,ff
     sprod  = 0; 
     do i =1,4
        sprod(i,i,nz+1)=(1.d0,0.d0)
     enddo

     do i =nz,1,-1
        if(i.ne.nz+1)then
         axx = conv(i,1)-conv(i,5)**2/conv(i,3)
         ayy = conv(i,2)-conv(i,6)**2/conv(i,3)
         axy = conv(i,4)-conv(i,5)*conv(i,6)/conv(i,3)
        else
         axx = conv(nz,1)-conv(nz,5)**2/conv(nz,3)
         ayy = conv(nz,2)-conv(nz,6)**2/conv(nz,3)
         axy = conv(nz,4)-conv(nz,5)*conv(nz,6)/conv(nz,3)
        endif
      !  print *,i,axx,ayy,axy
         ada=axx+ayy; add = axx-ayy
         if(add.ge.0)then
            da12 = sqrt(add**2+4*axy**2)
         else
            da12 = -sqrt(add**2+4*axy**2)
         endif
        a1 = (ada+da12)/2; a2 = (ada-da12)/2
       
        kp = sqrt(-cmplx(0,1,8)*w*mu)*sqrt(a1)
        km = sqrt(-cmplx(0,1,8)*w*mu)*sqrt(a2)
        xip(i) = -kp/(cmplx(0,1,8)*w*mu)
        xim(i) = -km/(cmplx(0,1,8)*w*mu)
        
        if(axy.eq.0)then
           qp(i) =0; qm(i) =0  
        else
          qp(i) = 2*axy/(add+da12)
          qm(i) = (add-da12)/axy/2
        endif
        qdq = qp(i)*qm(i)
        dq = 1-qdq
         
        if(i.eq.nz)then
           m(1,2) = (1.d0,0.d0)
           m(2,2) =cmplx(qp(i),0.d0)
           m(3,2) =-xip(i)*qp(i)
           m(4,2) =xip(i)
           m(1,4) = cmplx(qm(i),0.d0)
           m(2,4) =(1.d0,0.d0)
           m(3,4) = -xim(i)
           m(4,4)=xim(i)*qm(i)
       endif
 
       tmp = kp*c(i)
       sp = (exp(tmp)-exp(-tmp))/2.d0 
       cp = (exp(tmp)+exp(-tmp))/2.d0
       tmp = km*c(i)
       sm = (exp(tmp)-exp(-tmp))/2.d0
       cm = (exp(tmp)+exp(-tmp))/2.d0
      ! print *,exp(tmp)
       
       s(1,:) = (/cp-qdq*cm,-qm(i)*(cp-cm),qm(i)*(sp/xip(i)-sm/xim(i)),&
                sp/xip(i)-qdq*sm/xim(i)/)
       s(2,:) =(/qp(i)*(cp-cm),-qdq*cp+cm,qdq*sp/xip(i)-sm/xim(i),&
                qp(i)*(sp/xip(i)-sm/xim(i))/)
       s(3,:)=(/-qp(i)*(xip(i)*sp-xim(i)*sm),qdq*xip(i)*sp-xim(i)*sm,&
                 -qdq*cp+cm,-qp(i)*(cp-cm)/)
       s(4,:) =(/xip(i)*sp-qdq*xim(i)*sm,-qm(i)*(xip(i)*sp-xim(I)*sm),&
                qm(i)*(cp-cm),cp-qdq*cm/)
      
      s = s/dq
      sprod(:,:,i) = matmul(s,sprod(:,:,i+1))
    enddo

    s_m = matmul(sprod(:,:,1),m)
    !!xy mode     
    e0 = (1.d0,0.d0)
    g = reshape((/s_m(1,2),s_m(2,2),s_m(1,4),s_m(2,4)/),(/2,2/))

    detg =g(1,1)*g(2,2)-g(1,2)*g(2,1)
    cbot = 0
    cbot(2) = (g(2,2)*real(e0)-g(1,2)*imag(e0))/detg
    cbot(4) = (-g(2,1)*real(e0)+g(1,1)*imag(e0))/detg
    mc = matmul(m,cbot)
   ! print *,mc
    do i = 1,nz+1
        f = matmul(sprod(:,:,i),mc)
      ! print *,f
       ex(i) = f(1)
     !  ex(i+nz+1) = f(2)
       hx(i) = f(4)
       if(i.eq.nz+1) atmp(1:4) = f
    enddo
   ! write(22,*) ex
   ! stop
   !!yx mode
    e0 = (0.d0,1.d0)
    g = reshape((/s_m(1,2),s_m(2,2),s_m(1,4),s_m(2,4)/),(/2,2/))

    detg =g(1,1)*g(2,2)-g(1,2)*g(2,1)
    cbot = 0
    cbot(2) = (g(2,2)*real(e0)-g(1,2)*imag(e0))/detg
    cbot(4) = (-g(2,1)*real(e0)+g(1,1)*imag(e0))/detg
    mc = matmul(m,cbot)

    do i = 1,nz+1
        f = matmul(sprod(:,:,i),mc)
        ex(i+nz+1) = f(2)
        hx(i+nz+1) = f(3)
        if(i.eq.nz+1) atmp(5:8) = f
    enddo
   !! impeadance on the bottom
    tmp = atmp(3)*atmp(8)-atmp(4)*atmp(7)  
    imp(1) =(atmp(1)*atmp(8)-atmp(5)*atmp(4))/tmp 
    imp(2) =(atmp(5)*atmp(3)-atmp(1)*atmp(7))/tmp 
    imp(3) =(atmp(2)*atmp(8)-atmp(6)*atmp(4))/tmp
    imp(4) =(atmp(6)*atmp(3)-atmp(2)*atmp(7))/tmp
    deallocate(sprod,xip,xim,qp,qm)
  ! write(82,*) ex
  ! stop
   end subroutine

subroutine mt1dfem_n(f,c_z,nz,nair,patm,cont,e_dc)
!! this subroutine is desinged for 1d arbitary anisotropic problem
!! --support ANY layers of anisotropic material
!! --3rd order element is used to discretize the original equation
!! --apparent resistivities and phases are calculated for comparison, stored
!! in asp as (rho_xy,rho_yx,phi_xy,phi_yx)
!! --this routine employs mkl direct solver zcgesv due to the small size
!! of the matrix, iterative solvers like BiCG are also applicable
implicit none
integer,intent(in) :: nz,nair
integer(kind=4) :: mt1dne, np
integer(kind=4) :: i,h,nj,me(4,nz),j,nk,n,dj,nab
real(kind=8) :: patm(:),cont(:,:)
real(kind=8):: rho(3*nz+1),c_z(*),f,r1,r2,r3,r4,w,k1(4,4),k2(4,4),asp(4)
complex(kind=8) :: v(6*nz+2,6*nz+2),b1(6*nz+2),b2(6*nz+2),m,a,e1(6*nz+2),&
                 v1(6*nz+2,6*nz+2),v2(6*nz+2,6*nz+2),e2(6*nz+2)
complex(kind=8),intent(out) :: e_dc(2*nz+2)
logical :: aoi(nz)
integer(kind=4) :: iter,info
!real(kind=8) :: tol
!real(kind=8),external :: cnorm
!complex(kind=8) :: ax(3*(nz+1)),p(3*(nz+1)),r(3*(nz+1))
complex(kind=8) :: ex,ey,hx,hy,zxy,zyx
real(kind=8) :: pc,phase
complex(8),allocatable:: work(:)
complex,allocatable::swork(:)
real(8),allocatable::rwork(:),aont(:,:)
integer,allocatable ::ipiv(:)

aoi = .false.;!array for judging a layer is anisotropic or not(true=anis)
!nab = size(abr)/6
allocate(aont(nz,4))
do i = 1,nz
    if(.not.(all(cont(i,4:6)==0.d0).and.cont(i,1)==cont(i,2).and.cont(i,1)==cont(i,3))) aoi(i) = .true.
    aont(i,1) = cont(i,1)-cont(i,5)**2/cont(i,3)
    aont(i,2) = cont(i,4)-cont(i,5)*cont(i,6)/cont(i,3)
    aont(i,3) = cont(i,4)-cont(i,6)*cont(i,5)/cont(i,3)
    aont(i,4) = cont(i,2)-cont(i,6)**2/cont(i,3)
enddo
!print *,aont
mt1dne = nz; np = 3*nz+1
allocate(work(2*np),swork(2*np*(2*np+1)),rwork(2*np))
allocate(ipiv(2*np))

!print *,'in mt1d',abr

do i=1,np
    rho(i)=patm(2)
enddo
do i=1,nair
    rho(i)=patm(1)
enddo

v = 0; b1 = 0; b2 = 0;

do i=1,mt1dne
    me(1,i)=1+(i-1)*3
    me(2,i)=2+(i-1)*3
    me(3,i)=3+(i-1)*3
    me(4,i)=4+(i-1)*3
enddo

data((k1(i,j),j=1,4),i=1,4) /148.d0, -189.d0, 54.d0,-13.d0, &
                            -189.d0,432.d0,-297.d0,54.d0,  &
                             54.d0,  -297.d0,432.d0,-189.d0,&
                             -13.d0,54.d0,-189.d0,148.d0/
data((k2(i,j),j=1,4),i=1,4)/512.d0,396.d0,-144.d0,76.d0,&
                            396.d0,2592.d0,-324.d0,-144.d0,&
                            -144.d0,-324.d0,2592.d0,396.d0,&
                             76.d0, -144.d0,396.d0,512.d0/

    w=2.d0*pi*f
    m=cmplx(0.d0,1.d0)*w*mu
    do h=1,mt1dne
        do j=1,4
            nj=me(j,h)
            do n=1,4
                nk=me(n,h)
              if(.not.aoi(h))then
                  v(nj,nk)=v(nj,nk)+k1(j,n)/40.d0/c_z(h)-k2(j,n)*m*c_z(h)/6720.d0/rho(h)
                  v(nj+np,nk+np) =  v(nj+np,nk+np)+k1(j,n)/40.d0/c_z(h)-&
                                  k2(j,n)*m*c_z(h)/6720.d0/rho(h)
              else
                  v(nj,nk)=v(nj,nk)+k1(j,n)/40.d0/c_z(h)-aont(h,1)*k2(j,n)*m*c_z(h)/6720.d0
                  v(nj,nk+np) = v(nj,nk+np)-aont(h,2)*k2(j,n)*m*c_z(h)/6720.d0
                  v(nj+np,nk) = v(nj+np,nk)-aont(h,3)*k2(j,n)*m*c_z(h)/6720.d0
                  v(nj+np,nk+np) = v(nj+np,nk+np)+k1(j,n)/40.d0/c_z(h)- &
                                 aont(h,4)*k2(j,n)*m*c_z(h)/6720.d0
             endif
           enddo
        enddo
    enddo

    !! force Neumann boundary on the bottom
    a=sqrt(-m/rho(np))
    v(np,np)=v(np,np)+a
    v(2*np,2*np) = v(2*np,2*np)+a
  !! xy-mode
    b1(1)=(1.d0,0.d0)
    do i=2,np*2
        b1(i)=b1(i)-v(i,1)
    enddo
    v1 = v
    do i=1,np*2
        v1(1,i)=(0.d0,0.d0)
        v1(i,1)=(0.d0,0.d0)
    enddo
    v1(1,1)=(1.d0,0.d0)
    e1 = 0
   ! print *,'before solve',tol,itmax
   !  call bcgmall(v1,b1,e1,np*2,NP*2,TOL,ITMAX)
   ! following direct solver is recommended, since bicg might not 
   ! converge in maximum steps
     call zcgesv(2*np,1,v1,2*np,ipiv,b1,2*np,e1,2*np,work,swork,&
       rwork,iter,info)

 !! yx-mode
    b2(np+1)=(1.d0,0.d0)
    do i=1,np*2
      if(i.ne.np+1)then
        b2(i)=b2(i)-v(i,np+1)
      endif
    enddo
    v2 = v
    do i=1,np*2
        v2(np+1,i)=(0.d0,0.d0)
        v2(i,np+1)=(0.d0,0.d0)
    enddo
    v2(np+1,np+1)=(1.d0,0.d0)
    e2 = 0
   ! call bcgmall(v2,b2,e2,np*2,NP*2,TOL,ITMAX)
      
    call zcgesv(2*np,1,v2,2*np,ipiv,b2,2*np,e2,2*np,work,swork,&
       rwork,iter,info)
    !write(123,*) e1
    e_dc(1:nz+1) = e1(1:np:3)
    e_dc(nz+2:2*nz+2) = e2(np+1:2*np:3)
    !ex =  e1(3*nair+1); ey = e2(3*nair+1+np)
    !
    !hy = (-11*e1(3*nair+1)+18*e1(3*nair+2)-9*e1(3*nair+3)+&
    !    2*e1(3*nair+4)) /2.d0/c_z(nair+1)/m
    !hx = (-11*e2(3*nair+1+np)+18*e2(3*nair+2+np)-9*e2(3*nair+3+np)+&
    !    2*e2(3*nair+4+np)) /2.d0/c_z(nair+1)/m
    !zxy = -ex/hy
    !zyx = ey/hx
    ! r1 =  -atan(imag(zxy)/real(zxy))*180.d0/pi
    ! r2 =  -atan(imag(zyx)/real(zyx))*180.d0/pi
    ! asp = (/abs(zxy)**2/w/mu,abs(zyx)**2/w/mu,r1,r2/)
    !  print *,abr
    deallocate(work,swork,rwork,ipiv)
   deallocate(aont)
end subroutine

 SUBROUTINE BCGMALL(A,B,X,N,NP,TOL,ITMAX)
   ! use solvers, only: vnorm
    IMPLICIT NONE
    INTEGER(KIND=4) :: ITER,I,ITMAX,N,NP
    REAL(KIND=8) :: VNRM,RSS,TOL
    COMPLEX(KIND=8) :: A(NP,N),B(N),X(N),AX(N),P(N),R(N)
    COMPLEX(KIND=8) :: ALPHA,GAMMA,RNORM
   ! REAL(KIND=8),EXTERNAL :: VNORM
  !  COMPLEX(KIND=8),EXTERNAL :: CNORM
    iter = 0
    VNRM=VNORM(B,N)
    DO I=1,N
        X(I)=0.
    ENDDO
    CALL PRODUCTALL(A,X,AX,N,NP)
    DO I=1,N
        R(I)=B(I)-AX(I)
        P(I)=R(I)
    ENDDO
    RSS=1.D0
    DO WHILE(ITER<=ITMAX.AND.RSS>TOL)
        RNORM=CNORM(R,R,N)
        CALL PRODUCTALL(A,P,AX,N,NP)
        ALPHA=RNORM/CNORM(AX,P,N)
        DO I=1,N
            X(I)=X(I)+ALPHA*P(I)
            R(I)=R(I)-ALPHA*AX(I)
        ENDDO
        GAMMA=CNORM(R,R,N)/RNORM
        RNORM=GAMMA*RNORM
        DO I=1,N
            P(I)=R(I)+GAMMA*P(I)
        ENDDO
        ITER=ITER+1
        RSS=VNORM(R,N)/VNRM
    ENDDO
    WRITE(*,*) 'MT1D_ITER=',ITER,'MT1D_RSS=',RSS
!    IF(RSS<TOL)THEN
!        WRITE(*,*) 'CONVERGENCE ACHIEVED.'
!    ELSE
!        IF(ITER>ITMAX)THEN
!            WRITE(*,*) 'MT1D: ITMAX EXCEEDED, NO CONVERGENCE.'
!        ENDIF
!    ENDIF
    RETURN
end

 SUBROUTINE PRODUCTALL(A,X,AX,IDN,NP)
    IMPLICIT NONE
    INTEGER(KIND=4) :: I,J,IDN,NP
    COMPLEX(KIND=8) :: A(NP,IDN),X(IDN),AX(IDN)
    DO I=1,IDN
        AX(I)=(0.,0.)
    ENDDO
    DO I=1,IDN
        DO J=1,IDN
            AX(I)=AX(I)+A(I,J)*X(J)
        ENDDO
    ENDDO
    RETURN
    END

    FUNCTION CNORM(A,B,IDN)
    IMPLICIT NONE
    COMPLEX(KIND=8) :: CNORM
    INTEGER(KIND=4) :: I,IDN
    COMPLEX(KIND=8) :: A(IDN),B(IDN),SUM
    SUM=(0.,0.)
    DO I=1,IDN
        SUM=SUM+A(I)*B(I)
    ENDDO
    CNORM=SUM
    RETURN
    END
  
    function vnorm(z,idn)
    implicit none
    real(kind=8) :: vnorm
    integer(kind=4) :: i,idn
    real(kind=8) :: sum
    complex(kind=8) :: z(idn)
    sum=0.
    do i=1,idn
        sum=sum+abs(z(i))**2
    enddo
    vnorm=sqrt(sum)
    return
    end function

 
    subroutine model_set(model, a, b, c,abr,abu, PaTM,xdet,conv)
      !! model settings 
      !! contains: COMMEMI3D-1, COMMEMI3D-2, DTM1, DTM3, and several CSEM models
       implicit none	   
       real(8),allocatable,intent(out):: a(:), b(:), c(:),PaTM(:),xdet(:),&
              conv(:,:)
       integer,allocatable,intent(out):: abr(:),abu(:)
       integer :: model,i,j,Nx,Ny,Nz,N_air
       real(8):: aa(7),ab(21)
       character (LEN=275) ::firstline	
	 !  N_air = 1 !! same for all models
	   allocate(abu(6))
       select case (model)
         case (1) !!COMMEMI3D-1
         allocate(a(144),b(144),c(160),abr(6),PaTM(3))
          PaTM = (/1.d10,100.d0,0.5d0/)  
        
          a(23:122) = 100.d0
          do i = 1,22
           a(i+122) = 100.d0*lambda_h**i
         enddo
          a(22:1:-1) = a(123:144) 
         ! b = 2*a;! c = a
          b = a
    
         c(21:140) = 50.d0
          do i = 1,20
            c(140+i) = 50.d0*lambda_v**i
          enddo
          c(20:1:-1) = c(141:160) 
       !   abr = (/4,6,4,6,4,6/)
          allocate(conv(1,6))
          conv= 0; ! conv(1,1:3) = 1.d1; !conv(2,1:3) = 1.d3
       !    abr = (/6,7,8,9,8,9, 10,11,8,9,8,9/)
          conv(1,1:3) = patm(3)
        !  print *,sum(a),sum(c)
          abr = (/68,77,63,82,46,85/)
          xdet =(/-3000:3000:100/)*1.d0 ! the x-coord of the receivers
				 
       case (2) !!COMMEMI3D-2
         allocate(a(12),b(12),c(13),abr(24),PaTM(5))
         !allocate(a(16),b(16),c(16),abr(24),PaTM(5))
         allocate(xdet(41),conv(4,6))
          PaTM = (/1.d9,1.d-1,1.d2,1.d-1,1.d0/) 

          c = (/58,32,18,5, 2,3,5, 10,10, 20,32,58,105/)*1.d3
           a = (/15,13,12,10, 10,10,10,10, 10,12,13,15/)*1.d3
           b = a
      
          abr = (/1,size(a),1,size(b),5,7, 1,size(a),1,size(b),8,9, &
                 5,6,5,8,5,7, 7,8,5,8,5,7/)
     
          conv = 0
          conv(1,1:3) = patm(2)
          conv(2,1:3) = patm(3); conv(3,1:3) = patm(5)
          conv(4,:) = conv(2,:)
      
          xdet = (/-4000:4000:200/)*1.d0
       case (3)	   
        !Dublin model 3.0(2016)
          allocate(a(24),b(27),c(24),abr(12),PaTM(4))
          allocate(xdet(2),conv(2,6))
         PaTM = (/1.d12,1.d3,3.d0,0.3d0/) 
     
         a = (/207,125,80,50,30,4, 4,4, 4,4,4,4,4,4,4,4,4,4, 4,14,30,75,120,213/)*1.d3	
         b = (/220,102,40,20,10,4, 4,8,8,8,8,8,8,8,8,8,8,8,8,4,4, 4,30,50,85,131,196/)*1.d3
         c = (/175.0,75.0,30.0,12.0,4.9,1.8,0.83,0.37,0.1,0.1,&
              0.1,0.45,1.3,2.15,3.0,8.0,23.0,51.0,113.0,195.0,300.0,600.0,1000.0,2000.0 /)*1.d3

         abr = (/7,8,20,21,11,14, 9,18,7,21,11,14/)
         conv = 0
         conv(1,1:3) = 3.d0
         conv(2,1:3) =0.3d0
         xdet  =0
       case (4) 
       !! dtm 1.0
        !! There are 3 anomolous bodies in this model, whose conductivity are 
        !! extremely contrast
       !  allocate(a(24),b(24),c(24),abr(18),PaTM(5))
         allocate(a(16),b(16),c(16),abr(18),PaTM(5))
         allocate(xdet(17),conv(3,6))
         PaTM = (/1.d9,1.d2,1.d1,1.d0,1.d4/)
         do i = 1,4
            b(12+i) = 5.d3*1.8**i
         enddo
         b(1:4) = b(16:13:-1) 
	     b(5:12) = 5.d3; !b(12:13) = 2.5d3
        ! a(20:24) = b(19:23)
        ! a(1:5) = a(24:20:-1)
        ! a(6:19) = 5.d3
         A = B; b(4) = 5.d3; b(13) = 5.d3  
         do i = 1,6
            c(10+i) = 5.d3*2.0**i
         enddo
         c(1:4) = c(14:11:-1)
         c(5:10) = 5.d3; c(12) = 1.d4
       !  print *,a
       !  print *,c
       !  stop
       !  abr = (/9,16,12,13,11,13, 10,12,12,17,14,14, 13,15,8,13,14,19/)
         abr = (/5,12,8,8,6,8, 9,11,8,12,9,9, 6,8,4,8,9,12/)
        ! xdet = (/-37.43,-30.0,-26.0,-23.43,-21.43,-18.86,-14.86,&
        !          -10.0,-6.0,-3.71,-1.43,0.0,1.43,3.71,6.0,10.0,&
        !          14.86,18.86,21.43,23.43,26.0,30.0,37.43/)*1.d3
         xdet =(/-37.5, -32.5, -27.5, -22.5, -17.5, -12.5, -7.5, -2.5, 0.0,&
                2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5/)*1.d3
          conv(1,1:3) = patm(3)
         conv(2,1:3) = patm(4)
         conv(3,1:3) = patm(5)
      !   abu =(/8,17,8,17,10,20/)
      case(5)
      !! complex anisotropic model from Liu et al(2019,JGR Solid Earth),
      !! related with the real geoelectric structure in the Junggar Basin
        allocate(a(32),b(26),c(24),patm(3),abr(42))
        patm = (/1.d9,5.d2,1.d2/)
        allocate(xdet(41),conv(7,6))         
      !  a(8:57) = 2.d3  
        a(5:28) = 4.d3
        do i =1,4
           a(28+i) = 4.d3*1.7d0**i
        enddo
        a(1:4) = a(32:29:-1)
        a(27) =2.d3; a(6) =2.d3; a(29) =4.d3; a(4)=4.d3
        b(1:4) = a(1:4); b(23:26) = a(29:32)
        b(5:22) = 4.d3 
      
        c = (/20000,8000,3000,800,200,&
              200,200,1000,2600,3000,3000,&
              3000,3000,4000,  4000,4000,5000,  5000,8000, &
              10000,20000,30000,50000,80000/)*1.d0
        !! order:A1,A2,A3,C1,C2,C2',C3
       abr = (/17,23,14,18,6,14, 22,26,9,13,6,14, 17,29,9,14,15,19,&
              27,29,9,18,6,11, 15,16,9,18,6,14, &
              5,11,9,18,6,14,  4,16,9,18,15,17/) 
       conv =0
       conv(1,:) = (/200.d0,2000.d0,200.d0,pi/2,0.d0,0.d0/)
       conv(2,:) = (/10.d0,500.d0,10.d0,pi/2,0.d0,0.d0/)
       conv(3,:) = conv(2,:)
       conv(3,4) = pi/9
       conv(4:6,1:3) = 10.d0
       conv(7,1:3) = 30.d0
       xdet = (/-20000:20000:1000/)*1.d0 
       
      case(6) 
         ! anisotropic layered model by Bai et al(2022,CAGEO)
         ! there are four layers in total, the 2nd and 3rd are anisotropic  
         allocate(a(55),b(55),c(55),abr(12),PaTM(3))     
         PaTM = (/1.d9,500.d0,1000.d0/) 
!        !  aa(1)=7.d3
          do i = 1,19
             a(i+36) = 5.d1*1.5d0**i
          enddo
          a(1:19) = a(55:37:-1)
          a(20:36) =5.d1
      
          b = a
           print *,sum(a)
        !   stop
           c(1:16) = a(4:19)
           c(11:16) = (/15,10,10,10,10,10/)*1.d0
           c(17:42) = (/10,35,125,240,400,450,350,250,100,&
                       100,140,240,500,360,240,100,&  
                       100,150,225,325,200,200,325,225,150,100/)*1.d0
        !  c(5:14) = 1.d3
        !  c(17:30) = 1.d2; c(22:25) = 4.2d2
          do i = 1,13
            c(42+i) =1.d2*1.8d0**i
          enddo
          print *,sum(c(1:16)),sum(c)
          abr=(/1,size(a),1,size(b),26,32,1,size(a),1,size(b),33,42/)
          allocate(conv(2,6),xdet(61)) ! anisotropic parameters
      
        ! Bai et al,2022,COGSC
          conv(1,:) = (/2.d2,1.d3,2.d2,pi/18,pi/6,pi/9/)
          conv(2,:) = (/1.d3,2.d2,1.d3,pi/18,pi/9,pi/6/)
      end select
    end subroutine
    
end module
