! COPYRIGHT (c) 2012-2018 Yvan Notay - ULB
!
! This file is part of AGMG software package
! Release 3.3.5-aca built on "May  9 2019" by Yvan Notay
!
! ALL USAGE OF AGMG IS SUBJECT TO LICENSE. PLEASE REFER TO THE FILE "LICENSE".
! IF YOU OBTAINED A ZCOPY OF THIS SOFTWARE WITHOUT THIS FILE,
! PLEASE CONTACT info@agmg.eu
!
! In particular, if you have a free academic license:
!
! (1) You must be a member of an educational, academic or research institution.
!     The license agreement automatically terminates once you no longer fulfill
!     this requirement.
!
! (2) You are obliged to cite AGMG in any publication or report as:
!     "Yvan Notay, AGMG software and documentation;
!      see http://agmg.eu".
!
! (3) You may not make available to others the software in any form, either
!     as source or as a precompiled object.
!
! (4) You may not use AGMG for the benefit of any third party or for any
!     commercial purposes. Note that this excludes the use within the
!     framework of a contract with an industrial partner.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! See the accompanying userguide for more details on how to use the software,
! and the README file for installation instructions.
!
! See the Web pages <http://agmg.eu> for
! release information and possible upgrade.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DICLAIMER:
!    AGMG is provided on an "AS IS" basis, without any explicit or implied
!    WARRANTY; see the file "LICENSE" for more details.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   If you use AGMG for research, please observe that your work benefits
!   our past research efforts that allowed the development of AGMG.
!   Hence, even if you do not see it directly, the results obtained thanks
!   to the use of AGMG depend on the results in publications [1-3] below,
!   where the main algorithms used in AGMG are presented and justified.
!   It is then a normal duty to cite these publications (besides citing
!   AGMG itself) in any scientific work depending on the usage of AGMG,
!   as you would do with any former research result you are using.
!
! [1] Y. Notay, An aggregation-based algebraic multigrid method,
!    Electronic Transactions on Numerical Analysis, vol. 37, pp. 123-146, 2010
!
! [2] A. Napov and Y. Notay, An algebraic multigrid method with guaranteed
!    convergence rate, SIAM J. Sci. Comput., vol. 34, pp. A1079-A1109, 2012.
!
! [3] Y. Notay, Aggregation-based algebraic multigrid for convection-diffusion
!    equations, SIAM J. Sci. Comput., vol. 34, pp. A2288-A2316, 2012.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file provides zagmg and dependencies; zagmg is a
! sequential implementation for complex matrices in double precision
! of the method presented in [1], where the used algorithms are described
! in detail. From realease 3.x, the coarsening has been modified according
! to the results in [2,3], see the release notes in the README file
! for more details.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! SUBROUTINE zagmg (MAIN DRIVER): see bottom of this file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! PARAMETERS DEFINITON -  may be tuned by expert users
  MODULE zagmg_mem
    SAVE
!
! INTEGER
!
!  maxlev  is the maximal number of levels
!          (should be large enough - much larger than needed is armless).
!  complex_len is the length of 1 COMPLEX(kind(0.0d0)) in byte
!        (used only to display information on memory usage).
!  nsmooth  is the number of pre- and post-smoothing steps
!       (see smoothtype for details).
!
!  smoothtype indicates which smoother is used:
!
!    if smoothtype==1, the smoother is based on Gauss-Seidel sweeps;
!    if smoothtype==0, the smoother is based on SOR sweeps with automatic
!       estimation of the relaxation parameter (often reduces to Gauss-Seidel);
!    if smoothtype==-1, the smoother is based on SOR sweeps with relaxation
!       parameter specified in constant omega defined below;
!
!     Scheme used in all three cases: nsmooth sweeps for both pre- and
!     post-smoothing, with:
!        pre-smoothing: Forward sweep, then Backward sweep, then Forward, etc
!       post-smoothing: Backward sweep, then Forward sweep, then Backward, etc;
!
!    if smoothtype==2, the smoother is ILU(0), with
!      (nsmooth/2) pre- and (nsmooth/2+mod(nsmooth,1)) post-smoothing steps.
!
!  nstep   is the maximum number of coarsening steps;
!          nstep<0 means that coarsening is stopped only according to
!          the matrix size, see parameter maxcoarsesize.
!  nlvcyc  is the number of coarse levels (from bottom) on which V-cycle
!          formulation is enforced (Rmk: K-cycle always allowed at first
!          coarse level).
!  npass   is the maximal number of pairwise aggregation passes for each
!          coarsening step, according to the algorithms in [2,3].
!  maxcoarsesize: when the size of the coarse grid matrix is less than or
!                 equal to maxcoarsesize*N^(1/3),  it is factorized exactly
!                 and coarsening is stopped;
!         maxcoarsesizeslow: in case of slow coarsening,
!                 exact factorization can be performed when the size of
!                 the coarse grid matrix is less than or equal to
!                 maxcoarsesizeslow*N^(1/3).
!         (N is the number of rows of the input matrix)
    INTEGER, PARAMETER :: maxlev=40, complex_len=16
    INTEGER, PARAMETER :: nsmooth=1, smoothtype=0, nstep=-1, nlvcyc=0
    INTEGER, PARAMETER :: npass=2,maxcoarsesize=40,maxcoarsesizeslow=400
!
! REAL
!
!  omega is the relaxation parameter for SOR sweeps (smoothtype=-1)
!  resi is the threshold for the relative residual error in inner FCG & GCR
!       iterations, see Algorithm 3.2 in [1]
!  trspos is a threshold: if a row has a positive offdiagonal entry larger
!         than trspos times the diagonal entry, the corresponding node is
!         transferred unaggregated to the coarse grid
!  kaptg_ is the threshold used to accept or not a tentative aggregate
!         when applying the coarsening algorithms from [2,3];
!         kaptg_blocdia is used for control based on bloc diagonal smoother [2];
!         kaptg_dampJac is used for control based on Jacobi smoother [3].
!  checkdd is the threshold to keep outside aggregation nodes where
!         the matrix is strongly diagonally dominant (based on mean of row
!         and column);
!         In fact, AGMG uses the maximum of |checkdd| and of the default value
!            according to kaptg_ as indicated in [2,3]
!            (hence |checkdd| < 1 ensures that one uses this default value)
!         checkdd <0 : consider |checkdd|, but base the test on minus the
!               sum of offdiagonal elements, without taking the absolute value.
!  targetcoarsefac is the target coarsening factor (parameter tau in the main
!         coarsening algorithm in [2,3]): further pairwise aggregation passes
!         are omitted once the number of nonzero entries has been reduced by a
!         factor of at least targetcoarsefac.
!  fracnegrcsum: if, at some level, more than fracnegrcsum*nl nodes,
!         where nl is the total number of nodes at that level, have
!         negative mean row and column sum, then the aggregation algorithm
!         of [2,3] is modified, exchanging all diagonal entries for the mean
!         row and column sum (that is, the algorithm is applied to a
!         modified matrix with mean row and column sum enforced to be zero).
! K_ItGramS is the K factor in the rule to decide if a second iteration of
!         Gram Schmidt orthogonalization is needed (used within GCR)
    REAL(kind(0.0d0)), PARAMETER :: omega=0.8, resi=0.2, trspos=0.45
    REAL(kind(0.0d0)), PARAMETER :: kaptg_blocdia=8.0, kaptg_dampJac=10.0
    REAL(kind(0.0d0)), PARAMETER :: checkdd=0.5
    REAL(kind(0.0d0)), PARAMETER :: targetcoarsefac=2.0**npass
    REAL(kind(0.0d0)), PARAMETER :: fracnegrcsum=0.25
    REAL(kind(0.0d0)), PARAMETER :: K_ItGramS=10.0
!!!!!!!!!!!!!!!!!!!! END of PARAMETERS DEFINITION -----------------
!!!!!!!!!!!!!!!!!!! Internal variables declaration
!
    TYPE InnerData
       COMPLEX(kind(0.0d0)), DIMENSION(:), POINTER :: a
       INTEGER, DIMENSION(:), POINTER :: ja
       INTEGER, DIMENSION(:), POINTER :: ia
       INTEGER, DIMENSION(:), POINTER :: il
       INTEGER, DIMENSION(:), POINTER :: iu
       COMPLEX(kind(0.0d0)), DIMENSION(:), POINTER :: p
       INTEGER, DIMENSION(:), POINTER :: idiag
       INTEGER, DIMENSION(:), POINTER :: ind
       INTEGER, DIMENSION(:), POINTER :: iext
       INTEGER, DIMENSION(:), POINTER :: ilstout
       INTEGER, DIMENSION(:), POINTER :: lstout
       INTEGER, DIMENSION(:), POINTER :: ilstin
       INTEGER, DIMENSION(:), POINTER :: lstin
       INTEGER, DIMENSION(:), POINTER :: iblockl
    END TYPE InnerData
!
    TYPE(InnerData) :: dt(maxlev)
    COMPLEX(kind(0.0d0)), ALLOCATABLE :: scald(:)
    INTEGER :: nn(maxlev),kstat(2,maxlev)=0,innermax(maxlev)
    INTEGER :: nlev,nwrkcum,iout,nrst,nbblock
    INTEGER :: maxcoarset,maxcoarseslowt,smoothtp
    REAL(kind(0.0d0)) :: wcplex(4),fracnz(maxlev),omeg=1.0d0,flop=0.0d0
    LOGICAL :: spd,wfo,wff,woo,allzero,trans,transint,zerors,gcrin
    LOGICAL :: coasing,coarsit,spdcoarse
    REAL(kind(0.0d0)), PARAMETER :: cplxmax=3.0, xsi=0.6
    REAL(kind(0.0d0)), PARAMETER :: repsmach=SQRT(EPSILON(1.0d0)),epsmach=EPSILON(1.0d0)
    REAL(kind(0.0d0)), PARAMETER :: UNP=1.0d0+10*epsmach
    INTEGER :: nlc(2),nlcp(2),nlc1(2),icum,npassr
    INTEGER :: imult
    REAL(kind(0.0d0)) :: ngl(2),nglp(2),nlctot(2),ngl1(2),ngltot(2),RELRESL1
    INTEGER, DIMENSION(:), POINTER :: iextL1,ilstoutL1,lstoutL1,ilstinL1
    INTEGER, DIMENSION(:), POINTER :: lstinL1,ipermL1
    INTEGER, PARAMETER :: IRANK=-9999
    REAL(kind(0.0d0)) :: cputm=0.0d0,eltm=0.0d0,eltmp=0.0d0
    LOGICAL :: preprocessed=.FALSE., firstcall=.TRUE., solve=.FALSE.
    REAL(kind(0.0d0)), PARAMETER ::    &
         checkddJ=MAX(ABS(checkdd),kaptg_dampJac/(kaptg_dampJac-2))
    REAL(kind(0.0d0)), PARAMETER ::    &
         checkddB=MAX(ABS(checkdd),(kaptg_blocdia+1)/(kaptg_blocdia-1))
  END MODULE zagmg_mem
!!!!!!!!!!!!!!!!!!! END of Internal variables declaration
!!!!!!!!!!!!!!!!!!! TIMING
  SUBROUTINE zagmg_mestime(id,cputm,eltm)
    IMPLICIT NONE
    INTEGER, SAVE :: cpt_init(10)=-1
    REAL,    SAVE :: t1(10)
    INTEGER    :: cpt_fin,cpt_max,freq,cpt
    REAL       :: t2
    REAL(kind(0.0d0)) :: cputm,eltm
    INTEGER :: id
    IF (id.GT.0) THEN
       !Next line may be uncommented if FORTRAN 95 function
       !CPU_TIME is implemented
       !   CALL CPU_TIME(t2)
       CALL SYSTEM_CLOCK(cpt_fin,freq,cpt_max)
       cpt = cpt_fin - cpt_init(id)
       IF (cpt_fin.LT.cpt_init(id)) cpt = cpt + cpt_max
       eltm = dble(cpt) / freq
       cputm = dble(t2 - t1(id))
    ELSE
       CALL SYSTEM_CLOCK(cpt_init(-id),freq,cpt_max)
       !Next line may be uncommented if FORTRAN 95 function
       !CPU_TIME is implemented
       !   CALL CPU_TIME(t1(-id))
    END IF
    RETURN
  END SUBROUTINE zagmg_mestime
!!!!!!!!!!!!!!!!!!! END of TIMING
  MODULE zagmg_PARDISSO
    TYPE PARDISO_HANDLE
       INTEGER(KIND=8) DUMMY
    END TYPE PARDISO_HANDLE
  END MODULE zagmg_PARDISSO
  MODULE zagmg_PARDISO
    USE zagmg_PARDISSO
    INTERFACE
       SUBROUTINE PARDISO( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA  &
            , PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
         USE zagmg_PARDISSO
         TYPE(PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
         INTEGER,          INTENT(IN)    :: MAXFCT
         INTEGER,          INTENT(IN)    :: MNUM
         INTEGER,          INTENT(IN)    :: MTYPE
         INTEGER,          INTENT(IN)    :: PHASE
         INTEGER,          INTENT(IN)    :: N
         INTEGER,          INTENT(IN)    :: IA(*)
         INTEGER,          INTENT(IN)    :: JA(*)
         INTEGER,          INTENT(IN)    :: PERM(*)
         INTEGER,          INTENT(IN)    :: NRHS
         INTEGER,          INTENT(INOUT) :: IPARM(*)
         INTEGER,          INTENT(IN)    :: MSGLVL
         INTEGER,          INTENT(OUT)   :: ERROR
         COMPLEX(kind(0.0d0)), INTENT(IN)    :: A(*)
         COMPLEX(kind(0.0d0)), INTENT(INOUT) :: B(*)
         COMPLEX(kind(0.0d0)), INTENT(OUT)   :: X(*)
       END SUBROUTINE PARDISO
    END INTERFACE
  END MODULE zagmg_PARDISO
MODULE zagmg_ALLROUTINES
PRIVATE
PUBLIC :: zagmg_setupL1, zagmg_relmem,  zagmg_smoothsetup
PUBLIC :: zagmg_applyprec, zagmg_GCR, zagmg_FlexCG
PUBLIC :: zagmg_partroword
PUBLIC :: zagmg_diagscal
CONTAINS
    SUBROUTINE PARDISO( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA    &
         , PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
      USE zagmg_PARDISSO
      TYPE(PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
      INTEGER,          INTENT(IN)    :: MAXFCT
      INTEGER,          INTENT(IN)    :: MNUM
      INTEGER,          INTENT(IN)    :: MTYPE
      INTEGER,          INTENT(IN)    :: PHASE
      INTEGER,          INTENT(IN)    :: N
      INTEGER,          INTENT(IN)    :: IA(*)
      INTEGER,          INTENT(IN)    :: JA(*)
      INTEGER,          INTENT(IN)    :: PERM(*)
      INTEGER,          INTENT(IN)    :: NRHS
      INTEGER,          INTENT(INOUT) :: IPARM(*)
      INTEGER,          INTENT(IN)    :: MSGLVL
      INTEGER,          INTENT(OUT)   :: ERROR
      COMPLEX(kind(0.0d0)), INTENT(IN)    :: A(*)
      COMPLEX(kind(0.0d0)), INTENT(INOUT) :: B(*)
      COMPLEX(kind(0.0d0)), INTENT(OUT)   :: X(*)
      ERROR=1
      X(1:N)=cmplx(0,0,kind(0.0d0))
      STOP
    END SUBROUTINE PARDISO
  SUBROUTINE zagmg_relmem(ijb)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l, ier, ijb
    COMPLEX(kind(0.0d0)) :: fdum(1)
    DEALLOCATE(scald)
    DO l=1,nlev-1
       IF(nn(l).GT.0) THEN
            IF (nn(l+1).GT.0)     &
                      DEALLOCATE(dt(l)%ind)
            DEALLOCATE(dt(l)%a,dt(l)%ja,dt(l)%il,dt(l)%iu)
            IF (smoothtp.NE.1) DEALLOCATE(dt(l)%p)
       END IF
    END DO
    IF ( (nn(nlev).GT.0.AND.nlev.GT.1) .OR. (nlev.EQ.1.AND.ijb.LT.0) )   &
              DEALLOCATE(dt(nlev)%a,dt(nlev)%ja,dt(nlev)%ia)
    IF (.NOT.allzero)                              &
            CALL zagmg_DIRseq(nn(nlev),fdum,-2)
    RETURN
  END SUBROUTINE zagmg_relmem
  SUBROUTINE zagmg_applyprec( N,f,X,ijb,a,ja,ia)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: N
    COMPLEX(kind(0.0d0)), OPTIONAL :: f(N), X(N)
    COMPLEX(kind(0.0d0)), OPTIONAL :: a(*)
    INTEGER, OPTIONAL :: ja(*), ia(N+1), ijb
    COMPLEX(kind(0.0d0)), ALLOCATABLE, SAVE :: S(:)
    REAL(kind(0.0d0)) :: tau=1.35
    INTEGER :: iter, nvec, i
       iter=ijb-2
       nvec=min(2,iter)+1
       ALLOCATE (S(nwrkcum+nvec*N))
    S(1:N)=f(1:n)*scald(1:N)
    IF (iter.GT.1) THEN
       CALL zagmg_prec_matv(1,S,X,S(2*N+1),S(3*N+1),.TRUE.,a,ja,ia)
       X(1:N)=tau*X(1:N)
       DO i=2,iter-1
          CALL zagmg_prec_matv(1,S,S(N+1),S(2*N+1),S(3*N+1),.TRUE.,a,ja,ia)
          X(1:N)=X(1:N)+tau*S(N+1:2*N)
          S(1:N)=S(1:N)-tau*S(2*N+1:3*N)
       END DO
       CALL zagmg_prec_matv(1,S,S(N+1),S(2*N+1),S(3*N+1),.FALSE.,a,ja,ia)
       X(1:N)=X(1:N)+tau*S(N+1:2*N)
    ELSE
       CALL zagmg_prec_matv(1,S,X,S(N+1),S(2*N+1),.FALSE.,a,ja,ia)
    END IF
    X(1:N)=X(1:n)*scald(1:N)
    kstat(2,1)=kstat(2,1)+iter
    DEALLOCATE(S)
    RETURN
  END SUBROUTINE zagmg_applyprec
  SUBROUTINE zagmg_partroword(n, a, ja, ia, idiag, iext, nzext)
    IMPLICIT NONE
    INTEGER :: n, ia(n+1), idiag(n)
    INTEGER, OPTIONAL :: iext(*), nzext
    INTEGER :: ja(*)
    INTEGER, ALLOCATABLE, TARGET, DIMENSION(:)  :: iw, jcol
    COMPLEX(kind(0.0d0)) :: a(*), p
    COMPLEX(kind(0.0d0)), ALLOCATABLE, TARGET :: w(:)
    INTEGER :: i, j, k, ipos, nzu, nze, ila, mal, mil, ifi
    ALLOCATE (w(n),iw(n),jcol(n))
    jcol=0
    ipos=ia(1)
    DO i=1,n
       nzu=0
       nze= 0
       p=cmplx(0,0,kind(0.0d0))
       ifi=ia(i)
       ia(i)=ipos
       DO k=ifi,ia(i+1)-1
          j=ja(k)
          IF (j.GT.i) THEN
             IF (jcol(j).EQ.0) THEN
               nzu=nzu+1
               w(nzu)=a(k)
               iw(nzu)=j
               jcol(j)=nzu
             ELSE
               w(jcol(j))=w(jcol(j))+a(k)
             END IF
          ELSE IF (j.LT.i) THEN
            IF (jcol(j).EQ.0) THEN
               a(ipos)=a(k)
               ja(ipos)=j
               jcol(j)=ipos
               ipos=ipos+1
            ELSE
               a(jcol(j))=a(jcol(j))+a(k)
            END IF
          ELSE
             p=p+a(k)
          END IF
       END DO
       DO k=ia(i),ipos-1
          jcol(ja(k))=0
       END DO
       idiag(i)=ipos
       a(ipos)=p
       ja(ipos)=i
       ipos=ipos+1
       DO k=1,nzu
          a(ipos)=w(k)
          ja(ipos)=iw(k)
          jcol(iw(k))=0
          ipos=ipos+1
       END DO
    END DO
    ia(n+1)=ipos
    DEALLOCATE (w,iw,jcol)
    RETURN
  END SUBROUTINE zagmg_partroword
  SUBROUTINE zagmg_diagscal(n,a,ja,ia,idiag,as)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n, ja(*), ia(n+1), idiag(n)
    COMPLEX(kind(0.0d0)) :: a(*), val
    COMPLEX(kind(0.0d0)) :: as(*)
    INTEGER :: i,k,k1,idum(1),ier
    DO i=1,n
       val=a(idiag(i))
       scald(i)=sqrt(abs(val)/val)
    END DO
    DO i=1,n
       DO k=ia(i),ia(i+1)-1
          as(k)=a(k)*scald(i)*scald(ja(k))
       END DO
    END DO
    RETURN
  END SUBROUTINE zagmg_diagscal
  SUBROUTINE zagmg_csrdlu(n,a,ja,ia,idiag,ao,jao,iup,iext,iextop)
    IMPLICIT NONE
    INTEGER :: n, ja(*), jao(n+1:*)
    INTEGER, TARGET :: ia(n+1), idiag(n+1)
    INTEGER, OPTIONAL, TARGET :: iup(n+1), iext(n+1),iextop(n+1)
    COMPLEX(kind(0.0d0)) :: a(*),ao(*)
    INTEGER :: i,j,ili,iui,iei,ipos,ili0,iui0,iei0,nl,nu,k,next
    INTEGER, POINTER :: il(:), iu(:), iexto(:)
    nl=0
    nu=0
    il => idiag(1:n+1)
    IF (PRESENT(iup)) THEN
       iu => iup(1:n+1)
       DO i=1,n
          ili=idiag(i)-ia(i)
          iui=ia(i+1)-idiag(i)-1
          iup(i)=iui
          idiag(i)=ili
          nl=nl+ili
          nu=nu+iui
       END DO
    ELSE
       iu => ia(1:n+1)
       DO i=1,n
          ili=idiag(i)-ia(i)
          iui=ia(i+1)-idiag(i)-1
          ia(i)=iui
          idiag(i)=ili
          nl=nl+ili
          nu=nu+iui
       END DO
    END IF
    ipos=0
    ili=n+1
    iui=ili+nl
    iei=iui+nu
       DO i=1,n
          ili0=ili
          iui0=iui
          DO j=1,il(i)
             ipos=ipos+1
             ao(ili)=a(ipos)
             jao(ili)=ja(ipos)
             ili=ili+1
          END DO
          ipos=ipos+1
          ao(i)=a(ipos)
          DO j=1,iu(i)
             ipos=ipos+1
             ao(iui)=a(ipos)
             jao(iui)=ja(ipos)
             iui=iui+1
          END DO
          iu(i)=iui0
          il(i)=ili0
       END DO
       iu(n+1)=iui
       il(n+1)=ili
    RETURN
  END SUBROUTINE zagmg_csrdlu
  SUBROUTINE zagmg_csrdluT(n,a,ja,ia,idiag,ao,jao,il,iu)
    IMPLICIT NONE
    INTEGER :: n, ja(*), ia(n+1), idiag(n), jao(n+1:*), il(n+1), iu(n+1)
    COMPLEX(kind(0.0d0)) :: a(*),ao(*)
    INTEGER :: i,j,ili,iui,iei,ipos,ili0,iui0,iei0,nl,nu,k,next
    nl=0
    nu=0
       il(2:n+1)=0
       iu(2:n+1)=0
       DO i=1,n
          iui=idiag(i)-ia(i)
          ili=ia(i+1)-idiag(i)-1
          nl=nl+ili
          nu=nu+iui
          DO k=ia(i),idiag(i)-1
             iu(ja(k)+1)=iu(ja(k)+1)+1
          END DO
          DO k=idiag(i)+1,ia(i+1)-1
             il(ja(k)+1)=il(ja(k)+1)+1
          END DO
       END DO
    ipos=0
    ili=n+1
    iui=ili+nl
    iei=iui+nu
       il(1)=ili
       iu(1)=iui
       DO i=1,n
          il(i+1) = il(i) + il(i+1)
          iu(i+1) = iu(i) + iu(i+1)
       END DO
       DO i=1,n
          DO k=ia(i),idiag(i)-1
             j = ja(k)
             next = iu(j)
             ao(next) = a(k)
             jao(next) = i
             iu(j) = next+1
          END DO
          ao(i)=a(idiag(i))
          DO k=idiag(i)+1,ia(i+1)-1
             j = ja(k)
             next = il(j)
             ao(next) = a(k)
             jao(next) = i
             il(j) = next+1
          END DO
       END DO
       DO i=n,1,-1
          il(i+1) = il(i)
          iu(i+1) = iu(i)
       END DO
       il(1)=ili
       iu(1)=iui
    RETURN
  END SUBROUTINE zagmg_csrdluT
  SUBROUTINE zagmg_csrmv(n, a, ja, ifja, ia, x, ifx, y, iad, flop, lstin)
    IMPLICIT NONE
    INTEGER :: n, ifja, ja(ifja:*), ia(n+1), iad, ifx, i, k, j
    INTEGER, OPTIONAL :: lstin(0:*)
    COMPLEX(kind(0.0d0)) :: a(*), x(ifx:*), y(n), t
    REAL(kind(0.0d0)) :: flop
    IF (n.LE.0) RETURN
       IF (iad.LT.-1) THEN
          DO i=1,n
             t=y(i)
             DO k=ia(i),ia(i+1)-1
                t = t - a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       ELSE IF (iad.GT.1) THEN
          DO i=1,n
             t=y(i)
             DO k=ia(i),ia(i+1)-1
                t = t + a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       ELSE IF (iad.EQ.-1) THEN
          DO i=1,n
             t=cmplx(0,0,kind(0.0d0))
             DO k=ia(i),ia(i+1)-1
                t = t - a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       ELSE
          DO i=1,n
             t=cmplx(0,0,kind(0.0d0))
             DO k=ia(i),ia(i+1)-1
                t = t + a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       END IF
    flop=flop+dble(2*(ia(n+1)-ia(1)))
    RETURN
  END SUBROUTINE zagmg_csrmv
  SUBROUTINE zagmg_mcsrmv(n, a, ja, ifja, ia, x, ifx, y, iad, flop)
    IMPLICIT NONE
    INTEGER :: n, ifja, ja(ifja:*), ia(n+1), iad, ifx, i, k, j
    COMPLEX(kind(0.0d0)) :: a(*), x(ifx:*), y(n), t
    REAL(kind(0.0d0)) :: flop
    IF (n.LE.0) RETURN
       IF (iad.LT.-1) THEN
          DO i=1,n
             t=y(i)-a(i)*x(i)
             DO k=ia(i),ia(i+1)-1
                t = t - a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
          flop=flop+dble(2*(ia(n+1)-ia(1)+n))
       ELSE IF (iad.GT.1) THEN
          DO i=1,n
             t=y(i)+a(i)*x(i)
             DO k=ia(i),ia(i+1)-1
                t = t + a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
          flop=flop+dble(2*(ia(n+1)-ia(1)+n))
       ELSE IF (iad.EQ.-1) THEN
          DO i=1,n
             t=-a(i)*x(i)
             DO k=ia(i),ia(i+1)-1
                t = t - a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
          flop=flop+dble(2*(ia(n+1)-ia(1))+n)
       ELSE
          DO i=1,n
             t=a(i)*x(i)
             DO k=ia(i),ia(i+1)-1
                t = t + a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
          flop=flop+dble(2*(ia(n+1)-ia(1))+n)
       END IF
    RETURN
  END SUBROUTINE zagmg_mcsrmv
  SUBROUTINE zagmg_csrlsolve(n, a, ja, ifja, ia, p, x, y, iunit, flop)
    IMPLICIT NONE
    INTEGER :: n, ifja, ja(ifja:*), ia(n+1), iunit, i, k
    COMPLEX(kind(0.0d0)) :: a(*), p(n), x(n), y(n), t
    REAL(kind(0.0d0)) :: flop
    IF (n.LE.0) RETURN
    IF (iunit.LT.0) THEN
       x(1) = y(1)
       DO i=2,n
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))*p(ja(k))
          END DO
          x(i) = t
       END DO
       flop=flop+dble(n+3*(ia(n+1)-ia(1)))
    ELSE IF (iunit.GT.0) THEN
       x(1) = p(1)*y(1)
       DO i=2,n
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = p(i)*t
       END DO
       flop=flop+dble(n+2*(ia(n+1)-ia(1)))
    ELSE
       x(1) = y(1)
       DO i=2,n
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = t
       END DO
       flop=flop+dble(2*(ia(n+1)-ia(1)))
    END IF
    RETURN
  END SUBROUTINE zagmg_csrlsolve
  SUBROUTINE zagmg_csrusolve(n, a, ja, ifja, ia, p, x, y, iunit, flop)
    IMPLICIT NONE
    INTEGER :: n, ifja, ja(ifja:*), ia(n+1), iunit, i, k
    COMPLEX(kind(0.0d0)) :: a(*), p(n), x(n), y(n), t
    REAL(kind(0.0d0)) :: flop
    IF (n.LE.0) RETURN
    IF (iunit.LT.0) THEN
       x(n) = y(n)
       DO i=n-1,1,-1
          t = cmplx(0,0,kind(0.0d0))
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = y(i) + p(i)*t
       END DO
       flop=flop+dble(n+2*(ia(n+1)-ia(1)))
    ELSE IF (iunit.GT.0) THEN
       x(n) = p(n)*y(n)
       DO i=n-1,1,-1
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = p(i)*t
       END DO
       flop=flop+dble(n+2*(ia(n+1)-ia(1)))
    ELSE
       x(n) = y(n)
       DO i=n-1,1,-1
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = t
       END DO
       flop=flop+dble(2*(ia(n+1)-ia(1)))
    END IF
    RETURN
  END SUBROUTINE zagmg_csrusolve
  SUBROUTINE zagmg_FlexCG(N,f,X,ITER,RESID,a,ja,ia,ijb)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER       :: N, ITER, ijb
    COMPLEX(kind(0.0d0)) :: f(N), X(N)
    INTEGER       :: ja(*), ia(N+1)
    COMPLEX(kind(0.0d0)) :: a(*)
    INTEGER       :: MAXIT, ierr, kk, k, idum(1), i
    REAL(kind(0.0d0)) :: TOL, BNORM, RESID, dum0, RESID0, FREL, RHO
    COMPLEX(kind(0.0d0)) :: ALPHA, BET0, BET1, BET, dum(4), td
    COMPLEX(kind(0.0d0)), ALLOCATABLE :: SD(:),S(:),fsc(:)
    LOGICAL :: init,intmv,intmv2
    COMPLEX(kind(0.0d0)), external :: ZDOTCM
    REAL(kind(0.0d0)), external :: DZNRM2
    INTEGER , parameter :: IONE=1
    intmv=.FALSE.
    IF ( (MOD(ijb,10).EQ.0 .OR. (ijb.GT.201)) .AND. nlev.GT.1) intmv=.TRUE.
    intmv2=intmv .OR. ijb.GT.201
    IF (intmv2) THEN
       ALLOCATE (S(2*N+1:nwrkcum+5*N),SD(N))
    ELSE
       ALLOCATE (S(2*N+1:nwrkcum+5*N),SD(N),fsc(N))
    END IF
    flop=0.0d0
    kstat=0
    init=.FALSE.
    IF (MOD(ijb,100).GE.10) init=.TRUE.
    IF (wfo) THEN
       WRITE(iout,940) IRANK
    END IF
    IF (wff) THEN
       IF (  trans) THEN
          WRITE(iout,941)
       END IF
       IF (nlev.GT.1) THEN
          IF (smoothtp.EQ.-1) THEN
             WRITE(iout,943) nsmooth
             WRITE(iout,944) omeg
          ELSE IF (smoothtp.EQ.1) THEN
             WRITE(iout,945) nsmooth
             WRITE(iout,947)
          ELSE
             WRITE(iout,946) nsmooth+mod(nsmooth,2),nsmooth-mod(nsmooth,2)
             WRITE(iout,947)
          END IF
       END IF
    END IF
    TOL = RESID
    MAXIT = ITER
    RESID = 1.0d0
    ITER = 0
    dum(3) = DZNRM2(N, f, IONE)**2
    flop=flop+dble(2*N)
    IF (intmv2) THEN
       f(1:N)=f(1:N)*scald(1:N)
       flop=flop+dble(N)
    END IF
    IF (init) THEN
       IF (intmv2) THEN
          x(1:N)=x(1:N)/scald(1:N)
          flop=flop+dble(2*N)
       END IF
       IF (intmv) THEN
          IF (smoothtp.EQ.1) THEN
             f(1:n)=f(1:n)-x(1:n)/dt(1)%a(1:n)
             CALL zagmg_csrmv(n,dt(1)%a,dt(1)%ja,n+1,dt(1)%il,x,1,f,-2,flop)
          ELSE
             CALL zagmg_mcsrmv(n,dt(1)%a,dt(1)%ja,n+1,dt(1)%il,x,1,f,-2,flop)
             f(1:n)=f(1:n)-x(1:n)/(dt(1)%p(1:n))
          END IF
          CALL zagmg_csrmv(n,dt(1)%a,dt(1)%ja,n+1,dt(1)%iu,x,1,f,-2,flop)
       ELSE
          IF (intmv2) THEN
             CALL zagmg_matv(N,x,SD,dt(1)%a,dt(1)%ja,dt(1)%ia,transint )
          ELSE
             CALL zagmg_matv(N,x,SD,a,ja,ia,trans )
          END IF
          f(1:n)=f(1:n)-SD(1:n)
       END IF
       dum(2) = DZNRM2(N, f, IONE)**2
       flop=flop+dble(2*N)
       BNORM=SQRT(dble(dum(3)))
       RESID=SQRT(dble(dum(2)))
       RESID0=RESID
       IF (BNORM.EQ.0.0d0) THEN
          IF(wff) THEN
             WRITE(iout,998)
             WRITE(iout,999)
          END IF
          X(1:N)=cmplx(0,0,kind(0.0d0))
          DEALLOCATE(S,SD)
          IF(ALLOCATED(fsc)) DEALLOCATE(fsc)
          RETURN
       END IF
       RESID=RESID/BNORM
       IF(wff.AND. (MAXIT.LE.0 .OR. RESID.LE.TOL)) THEN
          WRITE(iout,900) 0, resid*bnorm, resid
       END IF
    END IF
    DO WHILE ( ITER.LT.MAXIT .AND. RESID.GT.TOL )
       ITER = ITER + 1
       IF (.NOT.intmv2) THEN
          fsc(1:N)=f(1:N)*scald(1:N)
          CALL zagmg_prec_matv(1,fsc,S(1+3*N),S(1+4*N),S(1+5*N),intmv,a,ja,ia)
          S(1+3*N:4*N)=S(1+3*N:4*N)*scald(1:N)
          flop=flop+dble(2*N)
       ELSE IF (ijb.GT.201) THEN
          CALL zagmg_prec_matv(1,f,S(1+3*N),S(1+4*N),S(1+5*N),intmv)
       ELSE
          CALL zagmg_prec_matv(1,f,S(1+3*N),S(1+4*N),S(1+5*N),intmv,a,ja,ia)
       END IF
       kk=3
       IF ( ITER.GT.1 ) THEN
          kk=2
          dum(1) = - ZDOTCM(N, SD(1:n), IONE, S(1+3*N:4*N), IONE)
          BET0=dum(1)
          BET1=BET0/RHO
          CALL ZAXPY(N, BET1, S(1+2*N:3*N), IONE, S(1+3*N:4*N), IONE)
          flop=flop+dble(4*N)
          IF (intmv) THEN
             CALL ZAXPY(N, BET1, SD(1:n), IONE, S(1+4*N:5*N), IONE)
             flop=flop+dble(2*N)
          END IF
       END IF
       CALL ZCOPY(N, S(1+3*N:4*N), IONE, S(1+2*N:3*N), IONE)
       IF (intmv)  THEN
          CALL ZCOPY(N, S(1+4*N:5*N), IONE, SD(1:n), IONE)
       ELSE
          IF (intmv2) THEN
             CALL zagmg_matv(N,S(1+2*N),SD,dt(1)%a,dt(1)%ja,dt(1)%ia,transint )
          ELSE
             CALL zagmg_matv(N,S(1+2*N),SD,a,ja,ia,trans )
          END IF
       END IF
       dum(1) =  ZDOTCM(N, S(1+2*N:3*N), IONE, SD(1:n), IONE)
       dum(2) =  ZDOTCM(N,S(1+2*N:3*N),IONE,f,IONE)
       flop=flop+dble(4*N)
       IF (ITER.EQ.1) THEN
          IF (.NOT.init) THEN
             BNORM=SQRT(dble(dum(3)))
             RESID0=BNORM
             IF (BNORM.EQ.0.0d0) THEN
                IF(wff) THEN
                   WRITE(iout,998)
                   WRITE(iout,999)
                END IF
                X(1:N)=cmplx(0,0,kind(0.0d0))
                ITER=0
                DEALLOCATE(S,SD)
                IF(ALLOCATED(fsc)) DEALLOCATE(fsc)
                RETURN
             END IF
             IF(wff) THEN
                WRITE(iout,900) 0, resid*bnorm, resid
             END IF
             IF ( RESID.LE.TOL ) THEN
                ITER=0
                EXIT
             END IF
          END IF
       END IF
       RHO=dble(dum(1))
       ALPHA=dum(2)/RHO
       IF (RHO.EQ.0.0d0) EXIT
       IF (ITER.EQ.1 .AND. (.NOT.init)) THEN
          CALL ZCOPY(N,S(1+2*N:3*N),IONE,X,IONE)
          CALL ZSCAL(N,ALPHA,X,IONE)
          flop=flop+dble(N)
       ELSE
          CALL ZAXPY(N, ALPHA, S(1+2*N:3*N), IONE, X, IONE)
          flop=flop+dble(2*N)
       END IF
       CALL ZAXPY(N, -ALPHA, SD(1:n), IONE, f, IONE)
       dum0 = DZNRM2(N,f,IONE)**2
       flop=flop+dble(4*N)
       RESID=dum0
       RESID=SQRT(RESID)
       RESID=RESID/BNORM
       IF (wff) THEN
          WRITE(iout,900) iter, resid*bnorm, resid
       END IF
    END DO
    IF (intmv2) THEN
       x(1:N)=x(1:N)*scald(1:N)
       flop=flop+dble(N)
    END IF
    IF( resid.GT.tol ) THEN
       IF (woo) THEN
          WRITE(iout,'()')
          WRITE(iout,950) iter
          WRITE(iout,951)
          WRITE(iout,'()')
       END IF
       iter=-iter
    ELSE
       IF (wff) THEN
          WRITE(iout,952) iter
          WRITE(iout,'()')
       END IF
    END IF
    RESID=RESID*BNORM/RESID0
    kstat(2,1)=ABS(iter)
    DEALLOCATE(S,SD)
    IF(ALLOCATED(fsc)) DEALLOCATE(fsc)
    RETURN
900 FORMAT('****  Iter=',i5,'        Resid=',e9.2,                 &
         '        Relat. res.=',e9.2)
940 FORMAT(i3,                                                     &
         '*SOLUTION: flexible conjugate gradient iterations (FCG(1))')
941 FORMAT('****     Rmk: solve system with the transpose of the input matrix')
943 FORMAT('****     AMG preconditioner with',i2,             &
         ' SOR pre- and post-smoothing sweeps')
944 FORMAT('****         at each level  (relaxation factor: omega=',f5.2,')')
945 FORMAT('****     AMG preconditioner with',i2,             &
         ' Gauss-Seidel pre- and post-smoothing sweeps')
946 FORMAT('****     AMG preconditioner with',i2,             &
         ' ILU(0) pre- plus',i2,' ILU(0) post-smoothing step(s)')
947 FORMAT('****         at each level')
950 FORMAT('**** !!!   WARNING:',I5,' ITERATIONS WERE')
951 FORMAT('**** INSUFFICIENT TO ACHIEVE THE DESIRED ACCURACY')
952 FORMAT('****  - Convergence reached in',I5,' iterations -')
998 FORMAT('**** The norm of the right hand side is zero:')
999 FORMAT('****     set x equal to the zero vector and exit')
  END SUBROUTINE zagmg_FlexCG
  SUBROUTINE zagmg_GCR(N,f,X,ITER,RESID,a,ja,ia,ijb)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER       :: N, ITER, ijb
    COMPLEX(kind(0.0d0)) ::  f(N), X(N)
    INTEGER       :: ja(*), ia(N+1)
    COMPLEX(kind(0.0d0)) :: a(*)
    INTEGER   :: MAXIT,i,irst,ierr,j,k,kk,idum(1),ii
    COMPLEX(kind(0.0d0)) ::  ALPHA, BET0, tmp2
    REAL(kind(0.0d0)) ::  RESID, RESID2, BNORM2, RESID0
    REAL(kind(0.0d0)) ::  TOL, TOL2BNORM2, TOLT, dum0, t1, t2, rho, tmp1
    COMPLEX(kind(0.0d0)) :: dum(3),xd,t
    COMPLEX(kind(0.0d0)), ALLOCATABLE :: Sc(:),Siv(:),SiR(:)
    COMPLEX(kind(0.0d0)), ALLOCATABLE :: Su(:),R(:),fsc(:)
    LOGICAL :: init,intmv,intmv2,scndp
    COMPLEX(kind(0.0d0)), external :: ZDOTCM
    REAL(kind(0.0d0)), external :: DZNRM2
    INTEGER , parameter :: IONE=1
    INTEGER  :: itm1, m, info, itm
    CHARACTER*1 :: CHA
    intmv=.FALSE.
    IF ( (MOD(ijb,10).EQ.0 .OR. (ijb.GT.201)) .AND. nlev.GT.1) intmv=.TRUE.
    intmv2=intmv .OR. ijb.GT.201
    flop=0.0d0
    kstat=0
    init=.FALSE.
    IF (MOD(ijb,100).GE.10) init=.TRUE.
    IF (wfo) THEN
       WRITE(iout,938) IRANK,nrst
    END IF
    IF (wff) THEN
       IF (  trans) THEN
          WRITE(iout,941)
       END IF
       IF (nlev.GT.1) THEN
          IF (smoothtp.EQ.-1) THEN
             WRITE(iout,943) nsmooth
             WRITE(iout,944) omeg
          ELSE IF (smoothtp.EQ.1) THEN
             WRITE(iout,945) nsmooth
             WRITE(iout,947)
          ELSE
             WRITE(iout,946) nsmooth/2,(nsmooth+1)/2
             WRITE(iout,947)
          END IF
       END IF
    END IF
    TOL = RESID
    TOL2BNORM2 = TOL
    RESID2= 1.0d0
    MAXIT = ITER
    ITER = 0
    m=MIN(nrst,MAXIT)
    m=MAX(m,2)
    ALLOCATE (Su(N*m),Sc(N*m),SiR(((m+1)*m)/2),Siv(2*m+1),R(nwrkcum))
    dum(3) = DZNRM2(N, f, IONE)**2
    flop=flop+dble(2*N)
    IF (.NOT.intmv2) ALLOCATE(fsc(N))
    IF (intmv2) THEN
       f(1:N)=f(1:N)*scald(1:N)
       flop=flop+dble(N)
    END IF
    IF (init) THEN
       IF (intmv2) THEN
          x(1:N)=x(1:N)/scald(1:N)
          flop=flop+dble(2*N)
       END IF
       IF (intmv) THEN
          IF (smoothtp.EQ.1) THEN
             f(1:n)=f(1:n)-x(1:n)/dt(1)%a(1:n)
             CALL zagmg_csrmv(n,dt(1)%a,dt(1)%ja,n+1,dt(1)%il,x,1,f,-2,flop)
          ELSE
             CALL zagmg_mcsrmv(n,dt(1)%a,dt(1)%ja,n+1,dt(1)%il,x,1,f,-2,flop)
             f(1:n)=f(1:n)-x(1:n)/(dt(1)%p(1:n))
          END IF
          CALL zagmg_csrmv(n,dt(1)%a,dt(1)%ja,n+1,dt(1)%iu,x,1,f,-2,flop)
       ELSE
          IF (intmv2) THEN
             CALL zagmg_matv(N,x,Sc,dt(1)%a,dt(1)%ja,dt(1)%ia,transint )
          ELSE
             CALL zagmg_matv(N,x,Sc,a,ja,ia,trans )
          END IF
          f(1:n)=f(1:n)-Sc(1:N)
       END IF
       dum(2) = DZNRM2(N, f, IONE)**2
       flop=flop+dble(2*N)
       BNORM2=dble(dum(3))
       RESID2=dble(dum(2))
       IF (BNORM2.EQ.0.0d0) THEN
          IF(wff) THEN
             WRITE(iout,998)
             WRITE(iout,999)
          END IF
          X(1:N)=cmplx(0,0,kind(0.0d0))
          DEALLOCATE (Su,Sc,R,Siv,SiR)
          IF(ALLOCATED(fsc)) DEALLOCATE(fsc)
          RETURN
       END IF
       TOL2BNORM2 = TOL*TOL*BNORM2
       IF (wff.AND. (MAXIT.LE.0 .OR. RESID2.LE.TOL2BNORM2)) THEN
          WRITE(iout,900) 0, 0,' ',SQRT(resid2),SQRT(resid2/bnorm2)
       END IF
       RESID0=SQRT(RESID2)
    END IF
    itm  = -1
    irst = 0
    CHA=' '
    DO WHILE ( ITER.LT.MAXIT .AND. RESID2.GT.TOL2BNORM2 )
       itm  = itm  + 1
       ITER = ITER + 1
       IF (itm.EQ.m) THEN
          CALL ZTPTRS('U','N','U',m,IONE,SiR,Siv,m,info)
          IF (irst.EQ.0 .AND. (.NOT.init)) THEN
             CALL ZGEMV('N',N,m,cmplx(1,0,kind(0.0d0)),Su,       &
                  N,Siv,IONE,cmplx(0,0,kind(0.0d0)),X,IONE)
             flop=flop+dble(2*m*N+m*(m+1))
          ELSE
             CALL ZGEMV('N',N,m,cmplx(1,0,kind(0.0d0)),Su,        &
                  N,Siv,IONE,cmplx(1,0,kind(0.0d0)),X,IONE)
             flop=flop+dble((2*m+1)*N+m*(m+1))
          END IF
          itm=0
          irst=irst+1
       END IF
       IF (.NOT.intmv2) THEN
          fsc(1:N)=f(1:N)*scald(1:N)
          CALL zagmg_prec_matv(1,fsc,Su(1+itm*N),Sc(1+itm*N),R,intmv,a,ja,ia)
          Su(1+itm*N:(1+itm)*N)=Su(1+itm*N:(1+itm)*N)*scald(1:N)
          flop=flop+dble(2*N)
       ELSE IF (ijb.GT.201) THEN
          CALL zagmg_prec_matv(1,f,Su(1+itm*N),Sc(1+itm*N),R,intmv)
       ELSE
          CALL zagmg_prec_matv(1,f,Su(1+itm*N),Sc(1+itm*N),R,intmv,a,ja,ia)
       END IF
       IF (.NOT.intmv) THEN
          IF (intmv2) THEN
             CALL zagmg_matv(N, Su(1+itm*N), Sc(1+itm*N)             &
                  , dt(1)%a, dt(1)%ja, dt(1)%ia, trans )
          ELSE
             CALL zagmg_matv(N, Su(1+itm*N), Sc(1+itm*N)             &
                  , a, ja, ia, trans )
          END IF
       END IF
       CHA=' '
       scndp=.FALSE.
       flop=flop+dble(2*N)
10     CONTINUE
       IF (itm.GT.0)       &
          CALL ZGEMV('C',N,itm,cmplx(1,0,kind(0.0d0)),Sc,N,Sc(1+itm*N:N+itm*N), &
               IONE,cmplx(0,0,kind(0.0d0)),Siv(m+2:m+itm+1),IONE)
       Siv(m+1)=DZNRM2(N, Sc(1+itm*N:N+itm*N), IONE)**2
       Siv(m+itm+2)=ZDOTCM(N, Sc(1+itm*N:(itm+1)*N), IONE, f, IONE)
       flop=flop+dble((4+4*itm)*N+4*itm)
       kk=0
       IF (ITER.EQ.1) THEN
             Siv(m)=dum(3)
             kk=1
       END IF
       IF (ITER.EQ.1) THEN
        IF (.NOT.init) THEN
          BNORM2=dble(Siv(m))
          RESID2=BNORM2
          IF (BNORM2.EQ.0.0d0) THEN
             IF(wff) THEN
                WRITE(iout,998)
                WRITE(iout,999)
             END IF
             X(1:N)=cmplx(0,0,kind(0.0d0))
             ITER=0
             DEALLOCATE (Su,Sc,R,Siv,SiR)
             IF(ALLOCATED(fsc)) DEALLOCATE(fsc)
             RETURN
          END IF
          TOL2BNORM2 = TOL*TOL*BNORM2
          RESID0=SQRT(RESID2)
          IF ( RESID2.LE.TOL2BNORM2 ) THEN
             ITER=0
             EXIT
          END IF
        END IF
        IF (wff) THEN
           WRITE(iout,900) 0, 0,' ',SQRT(resid2),SQRT(resid2/bnorm2)
        END IF
        TOLT = repsmach*RESID2
       END IF
       t1=dble(Siv(m+1))
       t2=0.0d0
       DO i=0,itm-1
          t2=t2+conjg(Siv(m+2+i))*Siv(m+2+i)/SiR(1+i+(i*(i+1))/2)
          Siv(m+2+i)=-Siv(m+2+i)/SiR(1+i+(i*(i+1))/2)
       END DO
       rho=t1-t2
       alpha=Siv(m+itm+2)
       IF (itm.GT.0) THEN
          CALL ZGEMV('N',N,itm,cmplx(1,0,kind(0.0d0)),Sc,N,Siv(m+2:m+itm+1),IONE,    &
                 cmplx(1,0,kind(0.0d0)),Sc(1+itm*N:N+itm*N),IONE)
          IF (scndp) THEN
             SiR(1+(itm*(itm+1))/2:itm+(itm*(itm+1))/2)=                   &
               SiR(1+(itm*(itm+1))/2:itm+(itm*(itm+1))/2)-Siv(m+2:m+itm+1)
          ELSE
             SiR(1+(itm*(itm+1))/2:itm+(itm*(itm+1))/2)=-Siv(m+2:m+itm+1)
             IF ( rho*(K_ItGramS**2).LE.t1 ) THEN
                scndp=.TRUE.
                GOTO 10
             END IF
          END IF
       END IF
       IF (rho.LE.0.0d0) THEN
          itm=itm-1
          EXIT
       END IF
       SiR(1+itm+(itm*(itm+1))/2)=rho
       bet0=alpha/rho
       Siv(1+itm)=bet0
       CALL ZAXPY(N, -bet0, Sc(1+itm*N:(itm+1)*N), IONE, f, IONE)
       RESID2 = RESID2 - (ABS(alpha)**2)/rho
       IF (RESID2.LE.TOLT) THEN
          RESID2 = DZNRM2(N,f,IONE)**2
          flop=flop+dble(2*N)
          TOLT = repsmach*RESID2
       END IF
       IF (wff)THEN
          WRITE(iout,900) iter,irst,CHA,SQRT(resid2),SQRT(resid2/bnorm2)
       END IF
    END DO
    IF (itm.GE.0) THEN
       itm1=itm+1
       CALL ZTPTRS('U','N','U',itm1, IONE,SiR,Siv,m,info)
       IF (irst.EQ.0 .AND. (.NOT.init)) THEN
          CALL ZGEMV('N',N,itm1,cmplx(1,0,kind(0.0d0)),Su,        &
               N,Siv,IONE,cmplx(0,0,kind(0.0d0)),X,IONE)
          flop=flop+dble(2*(itm+1)*N+(itm+1)*(itm+2))
       ELSE
          CALL ZGEMV('N',N,itm1,cmplx(1,0,kind(0.0d0)),Su,        &
               N,Siv,IONE,cmplx(1,0,kind(0.0d0)),X,IONE)
          flop=flop+dble((2*(itm+1)+1)*N+(itm+1)*(itm+2))
       END IF
    END IF
    IF (intmv2) THEN
       x(1:N)=x(1:N)*scald(1:N)
       flop=flop+dble(N)
    END IF
    RESID=SQRT(RESID2/BNORM2)
    IF( resid.GT.tol ) THEN
       IF (woo) THEN
          WRITE(iout,'()')
          WRITE(iout,950) iter
          WRITE(iout,951)
          WRITE(iout,'()')
       END IF
       iter=-iter
    ELSE
       IF (wff) THEN
          WRITE(iout,952) iter
          WRITE(iout,'()')
       END IF
    END IF
    RESID=RESID*SQRT(BNORM2)/RESID0
    DEALLOCATE (Su,Sc,R,Siv,SiR)
    IF(ALLOCATED(fsc)) DEALLOCATE(fsc)
    kstat(2,1)=ABS(iter)
    RETURN
900 FORMAT('****  Iter=',i5,' (',i2,' rest.)',A1,'       Resid=',e9.2,    &
         '        Relat. res.=',e9.2)
938 FORMAT(i3,'*SOLUTION: GCR iterations (GCR(',i2,'))')
941 FORMAT('****     Rmk: solve system with the transpose of the input matrix')
943 FORMAT('****     AMG preconditioner with',i2,             &
         ' SOR pre- and post-smoothing sweeps')
944 FORMAT('****         at each level  (relaxation factor: omega=',f5.2,')')
945 FORMAT('****     AMG preconditioner with',i2,             &
         ' Gauss-Seidel pre- and post-smoothing sweeps')
946 FORMAT('****     AMG preconditioner with',i2,             &
         ' ILU(0) pre- plus',i2,' ILU(0) post-smoothing step(s)')
947 FORMAT(  '****         at each level')
950 FORMAT('**** !!!   WARNING:',I5,' ITERATIONS WERE')
951 FORMAT('**** INSUFFICIENT TO ACHIEVE THE DESIRED ACCURACY')
952 FORMAT('****  - Convergence reached in',I5,' iterations -')
998 FORMAT('**** The norm of the right hand side is zero:')
999 FORMAT('****     set x equal to the zero vector and exit')
  END SUBROUTINE zagmg_GCR
  SUBROUTINE zagmg_matv(n, x, y, a, ja, ia, transpose,         &
       iext, lstin, l)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n, ja(*), ia(n+1), i, kk, k1, k2, ier, idum(1)
    INTEGER, OPTIONAL :: iext(n), lstin(0:*), l
    COMPLEX(kind(0.0d0)) :: y(n), a(*), t, xx
    COMPLEX(kind(0.0d0)) :: x(n)
    LOGICAL :: transpose
    IF (.NOT.transpose) THEN
       DO i = 1, n
          k1 = ia(i)
          xx = x(ja(k1))
          t = a(k1) * xx
          k2 = ia(i+1)-1
          DO kk = k1+1, k2
             xx = x(ja(kk))
             t = t + a(kk)*xx
          END DO
          y(i) = t
       END DO
    ELSE
       y(1:n)=cmplx(0,0,kind(0.0d0))
       DO i = 1, n
          xx = x(i)
          DO kk = ia(i), ia(i+1)-1
             y(ja(kk)) = y(ja(kk)) + a(kk)*xx
          END DO
       END DO
    END IF
    flop = flop + dble(2 *(ia(n+1)-ia(1))-n)
    RETURN
  END SUBROUTINE zagmg_matv
  RECURSIVE SUBROUTINE zagmg_prec_matv(l,B,X,AX,R,matv,a,ja,ia)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l
    COMPLEX(kind(0.0d0)), OPTIONAL ::  B(*), X(*), AX(*), R(*), A(*)
    LOGICAL, OPTIONAL :: matv
    INTEGER, OPTIONAL :: ja(*), ia(*)
    COMPLEX(kind(0.0d0)) ::  dum(1)
    LOGICAL :: update
    INTEGER :: is,n,nnext,XN,BN,AXN,RN,idum(1),i,j,XNe,BNe,nnexte,nnext1
    INTEGER, PARAMETER :: smoothilu=max(smoothtype,1)-1
    INTEGER, PARAMETER :: npostsmooth=nsmooth+smoothilu*mod(nsmooth,2)
    INTEGER, PARAMETER ::  npresmooth=nsmooth-smoothilu*mod(nsmooth,2)
    n=nn(l)
    nnext=nn(l+1)
    IF (l.EQ.nlev) THEN
       X(1:N)=B(1:N)
       IF (PRESENT(a)) THEN
          CALL zagmg_DIRseq(n,X,2,a,ja,ia)
       ELSE
          CALL zagmg_DIRseq(n,X,2)
       END IF
       RETURN
    END IF
    IF (npresmooth.EQ.0) THEN
       X(1:n)=0.0d0
    ELSE IF (npresmooth.EQ.1) THEN
       CALL zagmg_fwGS(l,B,X,AX,.FALSE.,R)
    ELSE
       update=.FALSE.
       DO is=2,npresmooth,2
          CALL zagmg_fwbwsolve1(l,B,X,AX,update,R)
          update=.TRUE.
       END DO
       IF (MOD(npresmooth,2).EQ.1) THEN
          CALL zagmg_fwGS(l,B,X,AX,.TRUE.,R)
       END IF
    END IF
    IF (nnext.GT.0) THEN
       XN=1
       BN=XN+nnext
       IF (l+1.EQ.nlev) BN=XN
       IF (npresmooth.EQ.0) THEN
         CALL zagmg_restag(N,nnext,B,R(BN),dt(l)%ind,flop)
       ELSE
         CALL zagmg_restag(N,nnext,AX,R(BN),dt(l)%ind,flop)
       END IF
       IF (l+1.EQ.nlev) THEN
          CALL zagmg_DIRseq(nnext,R(BN),2)
       ELSE IF ( innermax(l+1).LE.1 ) THEN
          AXN=BN+nnext
          RN=AXN+nnext
          CALL zagmg_prec_matv(l+1,R(BN),R(XN),R(AXN),R(RN),.FALSE.)
          kstat(1,l+1)=MAX(kstat(1,l+1),1)
          kstat(2,l+1)=kstat(2,l+1)+1
       ELSE
          CALL zagmg_inner_iter(nnext,R(XN),R(BN),l+1)
       END IF
       CALL zagmg_prolag(N,nnext,X,R(XN),dt(l)%ind,flop)
    END IF
    IF (npostsmooth.EQ.0) THEN
    ELSE IF (npostsmooth.EQ.1) THEN
       CALL zagmg_bwGS(l,B,X,AX,matv)
    ELSE
       IF (MOD(npostsmooth,2).EQ.1) THEN
          CALL zagmg_bwGS(l,B,X,AX,.FALSE.)
       END IF
       update=.FALSE.
       DO is=2,npostsmooth,2
          IF (is.GE.npostsmooth-1) update=matv
          CALL zagmg_fwbwsolve2(l,B,X,AX,update,R)
       END DO
    END IF
    RETURN
  END SUBROUTINE zagmg_prec_matv
  SUBROUTINE zagmg_fwGS(l,B,X,AX,update,R)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    COMPLEX(kind(0.0d0)) ::  B(*), X(*), AX(*), R(*)
    LOGICAL :: update
    n=nn(l)
    IF (update) THEN
       IF (smoothtp.EQ.1) THEN
         CALL zagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a   &
                              ,R,AX,1,flop)
       ELSE
         CALL zagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%p   &
                              ,R,AX,1,flop)
       END IF
       X(1:n)=X(1:n)+R(1:n)
       flop=flop+dble(n)
       IF (smoothtp.EQ.1) THEN
         CALL  zagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,R,1,AX,-1,flop)
       ELSE
         CALL  zagmg_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,R,1,AX,-1,flop)
       END IF
    ELSE
       IF (smoothtp.EQ.1) THEN
         CALL zagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a,  &
                               X,B,1,flop)
       ELSE
         CALL zagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%p,  &
                               X,B,1,flop)
       END IF
       IF (smoothtp.EQ.1) THEN
         CALL  zagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,X,1,AX,-1,flop)
       ELSE
         CALL  zagmg_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,X,1,AX,-1,flop)
       END IF
    END IF
    RETURN
  END SUBROUTINE zagmg_fwGS
  SUBROUTINE zagmg_bwGS(l,B,X,AX,matv)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    COMPLEX(kind(0.0d0)) ::  B(*), X(*), AX(*)
    LOGICAL :: matv
    n=nn(l)
    AX(1:n)=B(1:n)
    IF (smoothtp.EQ.1) THEN
      CALL  zagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,-2,flop)
    ELSE
      CALL  zagmg_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,-2,flop)
    END IF
    IF (smoothtp.EQ.1) THEN
      CALL zagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%a,X,AX,1,flop)
    ELSE
      CALL zagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%p,X,AX,1,flop)
    END IF
    IF (.NOT.matv) RETURN
    IF (smoothtp.EQ.1) THEN
      CALL  zagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,2,flop)
    ELSE
      CALL  zagmg_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,2,flop)
    END IF
    RETURN
  END SUBROUTINE zagmg_bwGS
  SUBROUTINE zagmg_fwbwsolve1(l,B,X,AX,update,R)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    COMPLEX(kind(0.0d0)) ::  B(*), X(*), AX(*), R(*)
    LOGICAL :: update
    n=nn(l)
    IF (.NOT.update) THEN
       IF (smoothtp.EQ.1) THEN
         CALL zagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a    &
                              ,R,B,-1,flop)
       ELSE
         CALL zagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%p    &
                              ,R,B,-1,flop)
       END IF
       IF (smoothtp.EQ.1) THEN
         CALL zagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%a    &
                              ,X,R,1,flop)
       ELSE
         CALL zagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%p    &
                              ,X,R,1,flop)
       END IF
       AX(1:n)=B(1:n)-R(1:n)
       flop=flop+dble(n)
       IF (smoothtp.EQ.1) THEN
         CALL  zagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,-2,flop)
       ELSE
         CALL  zagmg_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,-2,flop)
       END IF
    ELSE
       IF (smoothtp.EQ.1) THEN
         CALL zagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a    &
                              ,R,AX,-1,flop)
       ELSE
         CALL zagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%p    &
                              ,R,AX,-1,flop)
       END IF
       IF (smoothtp.EQ.1) THEN
         CALL zagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%a    &
                              ,R(n+1),R,1,flop)
       ELSE
         CALL zagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%p    &
                              ,R(n+1),R,1,flop)
       END IF
       X(1:n)=X(1:n)+R(n+1:2*n)
       AX(1:n)=AX(1:n)-R(1:n)
       flop=flop+dble(2*n)
       IF (smoothtp.EQ.1) THEN
         CALL  zagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,R(n+1),1,AX  &
                           ,-2,flop)
       ELSE
         CALL  zagmg_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,R(n+1),1,AX  &
                           ,-2,flop)
       END IF
    END IF
    RETURN
  END SUBROUTINE zagmg_fwbwsolve1
  SUBROUTINE zagmg_fwbwsolve2(l,B,X,AX,matv,R)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    COMPLEX(kind(0.0d0)) ::  B(*), X(*), AX(*), R(*)
    LOGICAL :: matv
    n=nn(l)
    IF (smoothtp.EQ.1) THEN
      CALL  zagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,X,1,AX,0,flop)
    ELSE
      CALL  zagmg_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,X,1,AX,0,flop)
    END IF
    R(1:n)=B(1:n)-AX(1:n)
    IF (smoothtp.EQ.1) THEN
      CALL zagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a    &
                           ,R(n+1),R,1,flop)
    ELSE
      CALL zagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%p    &
                           ,R(n+1),R,1,flop)
    END IF
    IF (matv) THEN
      IF (smoothtp.EQ.1) THEN
         AX(1:n)=AX(1:n)+R(n+1:2*n)/dt(l)%a(1:n)
      ELSE
         AX(1:n)=AX(1:n)+R(n+1:2*n)/dt(l)%p(1:n)
      END IF
      flop=flop+dble(2*n)
    END IF
    R(n+1:2*n) = R(n+1:2*n) - X(1:n)
    IF (smoothtp.EQ.1) THEN
      CALL zagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%a    &
                           ,R,R(n+1),-1,flop)
    ELSE
      CALL zagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%p    &
                           ,R,R(n+1),-1,flop)
    END IF
    X(1:n)=X(1:n) + R(1:n)
    flop=flop+dble(3*n)
    IF (.NOT.matv) RETURN
    CALL  zagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,2,flop)
    IF (smoothtp.NE.1) AX(1:n)=AX(1:n)+dt(l)%a(1:n)*R(1:n)
    RETURN
  END SUBROUTINE zagmg_fwbwsolve2
  RECURSIVE SUBROUTINE zagmg_inner_iter(n,X,R,l)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER   :: N, ITER, l, ierr, k, i
    REAL(kind(0.0d0)) ::  RESID, BNORM, det
    COMPLEX(kind(0.0d0)) :: X(N), R(MAX(N,1),*)
    COMPLEX(kind(0.0d0)) :: alpha1,alpha2,bet0,bet1,bet2,rho1,rho2,gamm0,gamm1,zeta
    COMPLEX(kind(0.0d0)) :: dum(5),y1,y2
    COMPLEX(kind(0.0d0)), external :: ZDOTCM
    REAL(kind(0.0d0)), external :: DZNRM2
    INTEGER , parameter :: IONE=1
    dum(3)=DZNRM2(N, R(1:N,1), IONE)**2
    IF (dble(dum(3)).EQ.0.0d0) THEN
       X(1:N)=cmplx(0,0,kind(0.0d0))
       flop=flop+dble(2*N)
       RETURN
    END IF
    ITER = 1
    CALL zagmg_prec_matv(l,R,X,R(1,2),R(1,3),.TRUE.)
    dum(1) = ZDOTCM(N,X,IONE,R(1:N,2),IONE)
    dum(2) = ZDOTCM(N,X,IONE,R(1:N,1),IONE)
    IF (ABS(dum(1)).LE.repsmach*ABS(dum(2))) THEN
       flop=flop+dble(2*N)
       GOTO 100
    END IF
    BNORM = dble(dum(3))
    rho1=dum(1)
    alpha1=dum(2)
    bet0=alpha1/rho1
    CALL ZAXPY(N, -bet0, R(1:N,2), IONE, R(1:N,1), IONE)
    RESID = DZNRM2(N,R(1:N,1),IONE)**2
    IF (RESID.LE.resi*resi*BNORM) THEN
       CALL ZSCAL( N, bet0, X, IONE )
       flop=flop+dble(11*N)
       kstat(1,l)=MAX(kstat(1,l),iter)
       kstat(2,l)=kstat(2,l)+iter
       RETURN
    END IF
    ITER = 2
    CALL zagmg_prec_matv(l,R,R(1,3),R(1,4),R(1,5),.TRUE.)
    dum(1) = ZDOTCM(N,R(1:N,3),IONE,R(1:N,2),IONE)
    dum(2) = ZDOTCM(N,R(1:N,3),IONE,R(1:N,1),IONE)
    dum(3) = ZDOTCM(N,R(1:N,3),IONE,R(1:N,4),IONE)
    IF (.NOT.spd) dum(4) = ZDOTCM(N,X,IONE,R(1:N,4),IONE)
    gamm0 =dum(1)
    alpha2=dum(2)
    rho2  =dum(3)
    IF (spd) THEN
       gamm1 = conjg(gamm0)
    ELSE
       gamm1 = dum(4)
    END IF
    bet1=rho2-gamm0*gamm1/rho1
    bet2=(alpha1-alpha2*gamm1/bet1)/rho1
    zeta=alpha2/bet1
    IF (ABS(bet1).LE.repsmach*ABS(alpha2) .OR.   &
        ABS(bet2)*repsmach.GE.1.0d0)            THEN
       flop=flop+dble(6*N)
       IF (.NOT.spd) flop=flop+dble(2*N)
      GOTO 200
    END IF
    CALL ZSCAL(N, bet2, X, IONE)
    CALL ZAXPY(N, zeta, R(1:N,3), IONE, X, IONE)
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    flop=flop+dble(19*N)
    IF (.NOT.spd) flop=flop+dble(2*N)
    RETURN
100 CONTINUE
    dum(1)=DZNRM2(N, R(1:N,2), IONE)**2
    dum(2)=ZDOTCM(N, R(1:N,2), IONE, R(1:N,1), IONE )
    BNORM = dble(dum(3))
110 CONTINUE
    rho1=dum(1)
    alpha1=dum(2)
    bet0=alpha1/rho1
    RESID=BNORM-dble(conjg(alpha1)*bet0)
    IF (RESID.LE.resi*resi*BNORM) THEN
       CALL ZSCAL( N, bet0, X, IONE )
       flop=flop+dble(7*N)
       kstat(1,l)=MAX(kstat(1,l),iter)
       kstat(2,l)=kstat(2,l)+iter
       RETURN
    END IF
    CALL ZAXPY(N, -bet0, R(1:N,2), IONE, R(1:N,1), IONE)
    ITER = 2
    CALL zagmg_prec_matv(l,R,R(1,3),R(1,4),R(1,5),.TRUE.)
    dum(1) = ZDOTCM(N,R(1:N,4),IONE,R(1:N,2),IONE)
    dum(2) = ZDOTCM(N,R(1:N,4),IONE,R(1:N,1),IONE)
    dum(3) = DZNRM2(N,R(1:N,4),IONE)**2
    gamm0 =dum(1)
    alpha2=dum(2)
    rho2  =dum(3)
    gamm1 = conjg(gamm0)
    bet1=rho2-gamm0*gamm1/rho1
    bet2=(alpha1-alpha2*gamm1/bet1)/rho1
    zeta=alpha2/bet1
    CALL ZSCAL(N, bet2, X, IONE)
    CALL ZAXPY(N, zeta, R(1:N,3), IONE, X, IONE)
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    flop=flop+dble(19*N)
    RETURN
200 CONTINUE
    dum(1) = DZNRM2(N,R(1:N,2),IONE)**2
    dum(2) = ZDOTCM(N,R(1:N,4),IONE,R(1:N,2),IONE)
    dum(3) = DZNRM2(N,R(1:N,4),IONE)**2
    dum(4) = ZDOTCM(N,R(1:N,2),IONE,R(1:N,1),IONE)
    dum(5) = ZDOTCM(N,R(1:N,4),IONE,R(1:N,1),IONE)
     dum(4) = dum(4)+bet0*dum(1)
     dum(5) = dum(5)+bet0*dum(2)
     det = dble(dum(1)*dum(3)-dum(2)*conjg(dum(2)))
     y1  = (dum(3)*dum(4)-conjg(dum(2))*dum(5))/det
     y2  = (-dum(2)*dum(4)+dum(1)*dum(5))/det
     CALL ZSCAL( N, y1, X, IONE )
     CALL ZAXPY( N, y2, R(1:N,3), IONE, X, IONE )
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    flop=flop+dble(23*N)
    RETURN
  END SUBROUTINE zagmg_inner_iter
  SUBROUTINE zagmg_prolag(n, nc, V, B, ind, flop)
    IMPLICIT NONE
    INTEGER :: n, nc, ind(n), k, i
    COMPLEX(kind(0.0d0)) :: V(n), B(nc)
    REAL(kind(0.0d0)) :: flop
    IF (nc.LE.0) RETURN
    DO i = 1, n
       k = ind(i)
       IF (k.GT.0) THEN
          V(i) = V(i)+B(k)
       END IF
    END DO
    flop = flop + dble(n)
    RETURN
  END SUBROUTINE zagmg_prolag
  SUBROUTINE zagmg_restag(n, nc, V, B, ind, flop)
    IMPLICIT NONE
    INTEGER :: n, nc, ind(n), k, i
    COMPLEX(kind(0.0d0)) :: V(n), B(nc)
    REAL(kind(0.0d0)) :: flop
    IF (nc.LE.0) RETURN
    B(1:nc)=cmplx(0,0,kind(0.0d0))
    DO i = 1, n
       k = ind(i)
       IF (k.GT.0) B(k) = B(k) + V(i)
    END DO
    flop = flop + dble(n)
    RETURN
  END SUBROUTINE zagmg_restag
  SUBROUTINE zagmg_setupL1(n,a,ja,ia)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n
    INTEGER :: ja(*),ia(n+1)
    COMPLEX(kind(0.0d0)) :: a(*)
    INTEGER :: nc,ierr,i,j,k,nz,nza,icmm
    LOGICAL :: slcoarse
    INTEGER, POINTER, DIMENSION(:) :: jap
    COMPLEX(kind(0.0d0)), POINTER, DIMENSION(:) :: ap
    REAL(kind(0.0d0)) :: eta,dum(3)
    CHARACTER(len=13) :: prtint
    COMPLEX (kind(0.0d0)) :: fff(1)
    nn(1)=n
    nlc(1)=n
    nlc(2)=ia(n+1)-ia(1)
    nza=nlc(2)
    ngl=nlc
    maxcoarset=FLOOR(maxcoarsesize*(n**(1.0d0/3)))
    maxcoarseslowt=FLOOR(maxcoarsesizeslow*(n**(1.0d0/3)))
    IF ( 0.EQ.nstep  .OR. 1.EQ.maxlev .OR. n.LE.maxcoarset ) nlev=1
       coasing=.FALSE.
       wcplex=1.0d0
       wcplex(4)=0.0d0
       nlctot=nlc
       ngltot=ngl
       ngl1=ngl
       nlc1=nlc
       fracnz(1)=1.0d0
       allzero=.FALSE.
       nglp=0.0d0
       IF (wfo) THEN
          IF (n.GT.9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') dble(n)
          ELSE
             WRITE(prtint(1:12),'(i12)') n
          END IF
          WRITE(iout,'()')
          WRITE(iout,918) prtint(1:12)
          IF (nlc(2).GT.9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') dble(nlc(2))
          ELSE
             WRITE(prtint(1:12),'(i12)') nlc(2)
          END IF
          WRITE(iout,919) prtint(1:12),dble(nlc(2))/dble(n)
          WRITE(iout,'()')
       END IF
    nlcp=nlc
    nglp=ngl
    IF (1.NE.nlev) THEN
       CALL zagmg_aggregation(1,n,a,ja,ia,nc)
       CALL zagmg_setup(2,nc)
    ELSE
       DEALLOCATE(dt(1)%idiag)
       CALL zagmg_DIRseq(n,fff,1,a,ja,ia)
       nwrkcum=1
       IF (ngl(1).GT.maxcoarset) wcplex(4)=log10(dble(ngl(1))/maxcoarset)
       RETURN
    END IF
       nz=ia(n+1)-ia(1)
       IF (transint) THEN
          ALLOCATE(ap(nz),jap(nz-n),dt(1)%il(n+1),dt(1)%iu(n+1))
          CALL zagmg_csrdluT(n,a,ja,ia,dt(1)%idiag              &
               ,ap,jap,dt(1)%il,dt(1)%iu)
          DEALLOCATE(dt(1)%idiag)
       ELSE
          ALLOCATE(ap(nz),jap(nz-n),dt(1)%iu(n+1))
          dt(1)%il => dt(1)%idiag
          CALL zagmg_csrdlu(n,a,ja,ia,dt(1)%idiag              &
               ,ap,jap,dt(1)%iu )
          NULLIFY(dt(1)%idiag)
       END IF
       dt(1)%a  => ap
       dt(1)%ja => jap
       NULLIFY(ap,jap)
       innermax(nlev)=0
       innermax(1)=1
       eta=xsi/((1-xsi)*(cplxmax-1))
       icmm=1
       DO i=2,nlev-1
          innermax(i)=min(2,floor(xsi**(i-1)/(eta*fracnz(i)*icmm)))
          IF (nlev-i.LE.nlvcyc .AND. i.GT.2) innermax(i)=1
          icmm=icmm*innermax(i)
          wcplex(2)=wcplex(2)+icmm*fracnz(i)
          wcplex(3)=wcplex(3)+(2**(i-1))*fracnz(i)
       END DO
       wcplex(2)=wcplex(2)+(2**(nlev-1))*fracnz(nlev)
       wcplex(3)=wcplex(3)+(2**(nlev-1))*fracnz(nlev)
       IF (nsmooth.GT.1 .OR. smoothtype.GT.1)  THEN
          nwrkcum=2*MAX(nn(nlev-1),1)
       ELSE
          nwrkcum=MAX(nn(nlev),1)
       END IF
       DO i=nlev-2,1,-1
          nwrkcum=nwrkcum+3*MAX(nn(i+1),1)
          IF (innermax(i+1).GT.1) nwrkcum=nwrkcum+2*MAX(nn(i+1),1)
          IF (nsmooth.GT.1 .OR. smoothtype.GT.1)  nwrkcum=max(2*nn(i),nwrkcum)
       END DO
       IF (wfo) THEN
          WRITE(iout,'()')
          WRITE(iout,954) nlctot(1)/dble(n)
          WRITE(iout,955) nlctot(2)/dble(nza)
       END IF
       IF (wff) THEN
          WRITE(iout,956) wcplex(3)
          WRITE(iout,957) wcplex(2)
          WRITE(iout,'()')
       END IF
    RETURN
918 FORMAT('****','        Number of unknowns:', A12)
919 FORMAT('****','                 Nonzeros :', A12,                &
         ' (per row:',f7.2,')')
920 FORMAT('****','        Number of variables:',A12,              &
         '          (reduction ratio:',f5.2,')')
921 FORMAT('****','                   Nonzeros:',A12,              &
         ' (per row:',f4.1,    &
         '; red. ratio:',f5.2,')')
954 FORMAT('****','                  Grid complexity:',f9.2)
955 FORMAT('****','              Operator complexity:',f9.2)
956 FORMAT('****','  Theoretical Weighted complexity:',f9.2, &
                ' (K-cycle at each level)' )
957 FORMAT('****','    Effective Weighted complexity:',f9.2, &
                ' (V-cycle enforced where needed)' )
958 FORMAT('****','              Weighted complexity:',f9.2)
  END SUBROUTINE zagmg_setupL1
  RECURSIVE SUBROUTINE zagmg_setup(l,n)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n
    INTEGER :: nc,ierr,i,j,k,nz
    LOGICAL :: slcoarse
    INTEGER, POINTER, DIMENSION(:) :: jap
    COMPLEX(kind(0.0d0)), POINTER, DIMENSION(:) :: ap
    LOGICAL, SAVE :: slowcoarse
    REAL(kind(0.0d0)) :: fw,eta,dum(3),aa
    CHARACTER(len=13) :: prtint
    COMPLEX (kind(0.0d0)) :: fff(1)
    nn(l)=n
    nlc(1)=n
    nlc(2)=0
    IF (n.GT.0) THEN
          nlc(2)=dt(l)%ia(n+1)-dt(l)%ia(1)
    END IF
    ngl=nlc
    IF (l.EQ.2) slowcoarse=.FALSE.
    IF (l.EQ.2) wcplex(1)=ngl1(2)/ngl(2)
    slcoarse = 2*nlcp(1).LT.3*nlc(1) .AND. 2*nlcp(2).LT.3*nlc(2)
    IF( l.EQ.nstep+1  .OR. l.EQ.maxlev                          &
       .OR. ( nlc(1).LE.maxcoarset )                            &
       .OR. ( nlcp(1).LT.2*nlc(1) .AND. nlcp(2).LT.2*nlc(2)       &
            .AND. nlc(1).LE.maxcoarseslowt )                    &
       .OR. ( slowcoarse .AND. slcoarse )                       &
       .OR. nlcp(1).EQ.nlc(1) )                       THEN
       nlev=l
    END IF
    slowcoarse=slcoarse
    fracnz(l)=ngl(2)/ngl1(2)
    nlctot=nlctot+nlc
    ngltot=ngltot+ngl
    IF (wfo) THEN
       WRITE(iout,'()')
       WRITE(iout,914) IRANK,l
       IF (n.GT.0) THEN
          IF (n.GT.9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') dble(n)
          ELSE
             WRITE(prtint(1:12),'(i12)') n
          END IF
          WRITE(iout,920) prtint(1:12),dble(nlcp(1))/dble(n)
          IF (nlc(2).GT.9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') dble(nlc(2))
          ELSE
             WRITE(prtint(1:12),'(i12)') nlc(2)
          END IF
          WRITE(iout,921) prtint(1:12),dble(nlc(2))/dble(n),        &
               dble(nlcp(2))/dble(nlc(2))
       ELSE
          WRITE(iout,924) n
       END IF
    END IF
    nlcp=nlc
    nglp=ngl
    IF (n.EQ.0) THEN
       nc=0
       allzero=.TRUE.
       RETURN
    END IF
    IF (l /= nlev) THEN
       CALL zagmg_aggregation(l,n,dt(l)%a,dt(l)%ja,dt(l)%ia,nc)
       CALL zagmg_setup(l+1,nc)
    ELSE
       IF (ngl(1).GT.maxcoarset) wcplex(4)=log10(dble(ngl(1))/maxcoarset)
          DEALLOCATE(dt(l)%idiag)
          CALL zagmg_DIRseq(n,fff,1,dt(l)%a,dt(l)%ja,dt(l)%ia)
          IF (coasing .AND. wfo) THEN
             WRITE(iout,912) IRANK
          END IF
          RETURN
    END IF
    IF (n.EQ.0) RETURN
    nz=dt(l)%ia(n+1)-dt(l)%ia(1)
    IF (transint) THEN
       ALLOCATE(ap(nz),jap(nz-n),dt(l)%il(n+1),dt(l)%iu(n+1))
       CALL zagmg_csrdluT(n,dt(l)%a,dt(l)%ja,dt(l)%ia,dt(l)%idiag &
            ,ap,jap,dt(l)%il,dt(l)%iu)
       DEALLOCATE(dt(l)%a,dt(l)%ja)
       DEALLOCATE(dt(l)%idiag,dt(l)%ia)
    ELSE
       ALLOCATE(ap(nz),jap(nz-n))
       dt(l)%iu => dt(l)%ia
       dt(l)%il => dt(l)%idiag
       CALL zagmg_csrdlu(n,dt(l)%a,dt(l)%ja,dt(l)%ia,dt(l)%idiag &
            ,ap,jap )
       DEALLOCATE(dt(l)%a,dt(l)%ja)
       NULLIFY(dt(l)%idiag,dt(l)%ia)
    END IF
    dt(l)%a  => ap
    dt(l)%ja => jap
    NULLIFY(ap,jap)
    RETURN
912 FORMAT(i3,'*','        Warning: coarsest grid matrix treated as singular')
914 FORMAT(i3,'*','                      Level:',I12)
918 FORMAT('****','        Number of unknowns:', A12)
919 FORMAT('****','                 Nonzeros :', A12,                &
         ' (per row:',f7.2,')')
920 FORMAT('****','        Number of variables:',A12,              &
         '          (reduction ratio:',f5.2,')')
921 FORMAT('****','                   Nonzeros:',A12,              &
         ' (per row:',f4.1,    &
         '; red. ratio:',f5.2,')')
924 FORMAT('****','        Number of variables:',i12)
  END SUBROUTINE zagmg_setup
  SUBROUTINE zagmg_smoothsetup
    USE zagmg_mem
    IMPLICIT NONE
    REAL(kind(0.0d0)) :: unmominv
    INTEGER :: l,i,llev
    smoothtp=smoothtype
    IF (smoothtype.EQ.0) THEN
       IF (spd) THEN
          smoothtp=1
       ELSE IF (omeg.LE.0.97) THEN
          smoothtp=-1
       ELSE
          smoothtp=1
       END IF
    END IF
    llev=nlev-1
    IF (smoothtp.EQ.-1) THEN
       IF (smoothtype.EQ.-1) omeg=omega
       unmominv=1.0d0-1.0d0/omeg
       DO l=1,llev
          IF (nn(l).GT.0) THEN
            ALLOCATE(dt(l)%p(nn(l)))
            DO i=1,nn(l)
               dt(l)%p(i)=omeg/dt(l)%a(i)
               dt(l)%a(i)=unmominv*dt(l)%a(i)
            END DO
          END IF
       END DO
    ELSE IF (smoothtp.EQ.1) THEN
       DO l=1,llev
          IF (nn(l).GT.0) THEN
            DO i=1,nn(l)
               dt(l)%a(i)=1.0d0/dt(l)%a(i)
            END DO
          END IF
       END DO
    ELSE IF (smoothtp.EQ.2) THEN
       DO l=1,llev
          IF (nn(l).GT.0) THEN
            ALLOCATE(dt(l)%p(nn(l)))
            CALL zagmg_ilu0(nn(l),dt(l)%a,dt(l)%ja,dt(l)%iu,dt(l)%il,dt(l)%p)
          END IF
       END DO
    END IF
    RETURN
  END SUBROUTINE zagmg_smoothsetup
  SUBROUTINE zagmg_ilu0(n,a,ja,iu,il,p)
    IMPLICIT NONE
    INTEGER :: n, ja(n+1:*), il(n+1), iu(n+1)
    COMPLEX(kind(0.0d0)) :: a(*), p(n), t
    INTEGER :: i,j,k,l
    DO i=1,n
      t=a(i)
      DO j=il(i),il(i+1)-1
         k=ja(j)
         DO l=iu(k),iu(k+1)-1
            IF (ja(l).EQ.i) THEN
               t=t-a(j)*a(l)*p(k)
               EXIT
            END IF
         END DO
      END DO
      a(i)=a(i)-t
      p(i)=1.0d0/t
    END DO
    RETURN
  END SUBROUTINE zagmg_ilu0
  SUBROUTINE zagmg_omegaest(n,a,ja,ia,idiag,si,w,iext)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n)
    INTEGER, OPTIONAL :: iext(*)
    COMPLEX(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n),w(2*n),t,tm
    INTEGER :: i,j,j1,j2,jd,ll,l2,kl,kk
    INTEGER :: ifirst,ilast,kb
    w(1:2*n)=0.0d0
    ifirst=1
    ilast=n
     DO i=ilast,ifirst,-1
       j1 =ia(i)
       jd=idiag(i)
       j2=ia(i+1)-1
       DO kk=j1,jd-1
        j=ja(kk)
          l2=ia(j+1)-1
          kl=0
          DO ll=idiag(j)+1,l2
            IF (ja(ll).EQ.i) THEN
              kl=ll
              EXIT
            END IF
          END DO
          IF (kl.EQ.0) THEN
             w(j)=w(j)+ABS(a(kk))
          ELSE
             w(j)=w(j)+ABS(a(kk)-conjg(a(kl)))-ABS(a(kl))
          END IF
       END DO
       t=w(i)
       DO kk=jd+1,j2
          t=t+ABS(a(kk))
       END DO
       w(i)=t/(max(dble(a(idiag(i))),si(i)))
     END DO
     tm=0.0d0
    ifirst=1
    ilast=n
     DO i=ifirst,ilast
       j1 =ia(i)
       jd=idiag(i)
       j2=ia(i+1)-1
       DO kk=jd+1,j2
        j=ja(kk)
          kl=0
          DO ll=ia(j),idiag(j)-1
            IF (ja(ll).EQ.i) THEN
               kl=ll
               EXIT
            END IF
          END DO
          IF (kl.EQ.0) THEN
             w(n+j)=w(n+j)+ABS(a(kk))*w(i)
          ELSE
             w(n+j)=w(n+j)+(ABS(a(kk)-conjg(a(kl)))-ABS(a(kl)))*w(i)
          END IF
       END DO
       t=w(n+i)
       DO kk=j1,jd-1
          t=t+ABS(a(kk))*w(ja(kk))
       END DO
       t=t/(max(dble(a(idiag(i))),si(i)))
       tm=max(t,tm)
     END DO
    IF (tm.LE.1.0d0) THEN
      omeg=1.0d0
    ELSE
      omeg=2/(1.0d0+SQRT(tm))
    END IF
    RETURN
  END SUBROUTINE zagmg_omegaest
  SUBROUTINE zagmg_aggregation(l,n,a,ja,ia,nc)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n,nc
    INTEGER :: ja(*),ia(n+1)
    COMPLEX (kind(0.0d0)) :: a(*), dum
    INTEGER :: ier,i,j,k,maxdg,np,kpass,nzc,m1,ndd,nzp,isize,nddp,npass1,nz,i0
    LOGICAL :: skipass
    INTEGER, POINTER, DIMENSION(:) :: jan,ian,idiagn,iextn,ind2,lcg,lcgn
    COMPLEX(kind(0.0d0)), POINTER, DIMENSION(:) :: an
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: sinn
    INTEGER, POINTER, DIMENSION(:) :: jap,iap,idiagp,iextp,lcg1
    COMPLEX(kind(0.0d0)), POINTER, DIMENSION(:) :: ap
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: sip
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ldd,iw,iperm,riperm
    REAL(kind(0.0d0)), ALLOCATABLE, DIMENSION(:) :: si1,w
    COMPLEX(kind(0.0d0)), ALLOCATABLE, DIMENSION(:) :: wc
    IF (l.EQ.1) THEN
       IF (wfo) THEN
          WRITE (iout,901) IRANK
       END IF
       IF (wff) THEN
          IF (spd) THEN
             WRITE (iout,906)
          ELSE IF (  transint) THEN
             WRITE (iout,908)
          END IF
          IF (.not.spd) then
             WRITE (iout,902) 'Jacobi',kaptg_dampJac,checkddJ
          ELSE
             WRITE (iout,902) 'BlockD',kaptg_blocdia,checkddB
          END IF
          IF (checkdd.LT.0) THEN
             WRITE (iout,904)
          END IF
          WRITE(iout,903) npass,targetcoarsefac
          WRITE (iout,905) trspos
       END IF
    END IF
    IF (l.EQ.1) THEN
       ALLOCATE(si1(n),ind2(n),iperm(n),riperm(n))
       CALL zagmg_setCMK(n,ja,ia,dt(l)%idiag,riperm,iperm)
    ELSE
       ALLOCATE(si1(n),ind2(n),iperm(n))
       iperm(1:n)=1
    END IF
    CALL zagmg_prepareagg(n,a,ja,ia,dt(l)%idiag,ind2,iperm,si1,ndd,l )
    IF (ndd.EQ.n) THEN
       nc=0
       nzc=0
       DEALLOCATE(si1,iperm,ind2)
       IF (l.EQ.1) THEN
          DEALLOCATE(riperm)
          IF (smoothtype.EQ.0) omeg=1.0d0
       END IF
       GOTO 999
    END IF
    IF (smoothtype.EQ.0 .AND. l.EQ.1 .AND. (.NOT.spd)) THEN
        ALLOCATE(w(2*n))
        CALL zagmg_omegaest(n,a,ja,ia,dt(l)%idiag,si1,w )
        DEALLOCATE(w)
    END IF
    ALLOCATE(ldd(ndd),lcg(2*(n-ndd)))
    IF (dble(n).GT.targetcoarsefac*(n-ndd)) THEN
       skipass=.TRUE.
       npass1=npass+1
    ELSE
       skipass=.FALSE.
       npass1=npass
    END IF
    IF (l.GT.1) THEN
       IF (spd) THEN
          CALL zagmg_findpairs_SI(l,n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,iperm )
       ELSE
          CALL zagmg_findpairs_GI(l,n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,iperm )
       END IF
       DEALLOCATE(iperm)
    ELSE
       IF (spd) THEN
          CALL zagmg_findpairs_SI1(l,n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,riperm,iperm )
       ELSE
          CALL zagmg_findpairs_GI1(l,n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,riperm,iperm )
       END IF
       DEALLOCATE(iperm,riperm)
    END IF
10  CONTINUE
    nz=ia(n+1)-ia(1)
    IF (npass1.GT.1) THEN
       isize=nc
    ELSE
       isize=1
    END IF
    ALLOCATE(an(nz-2*(n-nc)+ndd),jan(nz-2*(n-nc)+ndd)             &
         ,ian(nc+1),idiagn(nc+1),sinn(isize),wc(nc),iw(2*nc))
    CALL zagmg_setcg(n,a,ja,ia,dt(l)%idiag,si1,ind2,lcg,nc,an,jan,ian      &
         ,idiagn,sinn,npass1.GT.1,maxdg,iw,wc )
    DEALLOCATE(wc,iw)
    nzc=ian(nc+1)-1
    IF (dble(nz).GT.targetcoarsefac*nzc .OR. npass1.LE.1) THEN
       DEALLOCATE(si1,sinn,lcg)
       IF(ALLOCATED(ldd)) DEALLOCATE(ldd)
       dt(l)%ind  => ind2
       dt(l+1)%a    => an
       dt(l+1)%ja   => jan
       dt(l+1)%ia   => ian
       dt(l+1)%idiag=> idiagn
       NULLIFY(ind2,an,jan,ian,idiagn)
       GOTO 999
    END IF
    DEALLOCATE(ind2)
    lcg1 => lcg
    NULLIFY(lcg)
    m1=1
    DO kpass=2,npass1
       m1=2*m1
       np  = nc
       nzp = nzc
       ap     => an
       jap    => jan
       iap    => ian
       idiagp => idiagn
       sip    => sinn
       NULLIFY(an,jan,ian,idiagn,sinn)
       ALLOCATE(lcg(2*np),ind2(np),w(maxdg),iw(maxdg))
       ind2(1:np)=-1
       IF (spd) THEN
          CALL zagmg_findpairs_SF(l,np,ap,jap,iap,idiagp,sip  &
               ,ind2,lcg,nc                                  &
               ,m1,lcg1,a,ja,ia,dt(l)%idiag,si1,w,iw )
       ELSE
          CALL zagmg_findpairs_GF(l,np,ap,jap,iap,idiagp,sip     &
               ,ind2,lcg,nc                                     &
               ,m1,lcg1,a,ja,ia,dt(l)%idiag,si1,w,iw )
       END IF
       DEALLOCATE(w,iw)
       IF (kpass.NE.npass1) THEN
          isize=nc
       ELSE
          isize=1
          DEALLOCATE(si1)
       END IF
       ALLOCATE(an(nzp-2*(np-nc)),jan(nzp-2*(np-nc))                 &
            ,ian(nc+1),idiagn(nc+1),sinn(isize),wc(nc),iw(2*nc))
       CALL zagmg_setcg(np,ap,jap,iap,idiagp,sip,ind2,lcg,nc,an     &
            ,jan,ian,idiagn,sinn,kpass.NE.npass1,maxdg,iw,wc )
       DEALLOCATE(ap,jap,iap,idiagp,sip,ind2,wc,iw)
       ALLOCATE(lcgn(2*m1*nc))
       CALL zagmg_lcgmix(nc,m1,lcg1,lcg,lcgn)
       DEALLOCATE(lcg,lcg1)
       lcg1 => lcgn
       NULLIFY(lcgn)
       nzc=ian(nc+1)-1
       IF ( kpass.NE.npass1 .AND. dble(nz).GT.targetcoarsefac*nzc ) THEN
          DEALLOCATE(si1)
          EXIT
       END IF
    END DO
    DEALLOCATE(sinn)
    ALLOCATE(dt(l)%ind(n))
    CALL zagmg_setind(nc,ndd,ldd,lcg1,2*m1,dt(l)%ind)
    DEALLOCATE(lcg1,ldd)
       dt(l+1)%a    => an
       dt(l+1)%ja   => jan
       dt(l+1)%ia   => ian
       dt(l+1)%idiag=> idiagn
       NULLIFY(an,jan,ian,idiagn)
999 CONTINUE
    RETURN
901 FORMAT(i3,'*SETUP: Coarsening by multiple pairwise aggregation')
902 FORMAT('****       Quality threshold (',A6,'):',f6.2, &
         ' ;  Strong diag. dom. trs:',f5.2)
903 FORMAT('****         Maximal number of passes:',i3,     &
         '  ; Target coarsening factor:',f5.2)
904 FORMAT('****           Diag. dom. checked w.r.t. sum of offdiag', &
          ' (no absolute vaues)')
905 FORMAT('****',22x,'Threshold for rows with large pos. offdiag.:',f5.2)
906 FORMAT('****  Rmk: Setup performed assuming the matrix symmetric')
908 FORMAT('****  Rmk: Setup performed for the transpose of the input matrix')
  END SUBROUTINE zagmg_aggregation
  SUBROUTINE zagmg_setCMK(n,ja,ia,idiag,riperm,iperm,iext)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),riperm(*),iperm(n)
    INTEGER, OPTIONAL :: iext(*)
    LOGICAL :: exc
    INTEGER :: i,j,jj,jk,jd,k,kk,j1,j2,i1,i2,ijs,ijs1,ijs2,dg,mindg,kdim
    INTEGER :: ifirst,ilast,kb,i0,i2rec,i2rec0
    ifirst=1
    ilast=n
    i2=ifirst
    mindg=n+1
    DO i = ifirst,ilast
       dg=ia(i+1)-ia(i)
       IF (dg.GT.1) THEN
          iperm(i)=-dg
          IF (dg.LT.mindg) THEN
             mindg=dg
             jj=i
          END IF
       ELSE
          riperm(i2)=i
          iperm(i)=i2
          i2=i2+1
       END IF
    ENDDO
    ijs=ifirst-1
    i1=i2
15  CONTINUE
    IF (i2.LE.ilast) THEN
      riperm(i2)=jj
      iperm(jj)=i2
    END IF
    DO WHILE (i1.LE.i2 .AND. i2.LT.ilast)
       i=riperm(i1)
       ijs1=i2+1
       j1 =ia(i)
       jd=idiag(i)
       j2 = ia (i+1)-1
       DO kk = j1,jd-1
          j=ja(kk)
          IF (iperm(j).LT.0) THEN
             i2=i2+1
             riperm(i2)=j
          END IF
       ENDDO
       DO kk = jd+1,j2
          j=ja(kk)
          IF (iperm(j).LT.0) THEN
             i2=i2+1
             riperm(i2)=j
          END IF
       ENDDO
       ijs2=i2
       exc=.TRUE..AND. ijs2.GT.ijs1
       DO WHILE(exc)
          exc=.FALSE.
          DO kk=ijs1+1,ijs2
             IF( iperm(riperm(kk)).GT.iperm(riperm(kk-1)) )THEN
                j=riperm(kk)
                riperm(kk)=riperm(kk-1)
                riperm(kk-1)=j
                exc=.TRUE.
             END IF
          END DO
       END DO
       DO kk=ijs1,ijs2
          iperm(riperm(kk))=kk
       END DO
       i1=i1+1
    END DO
    IF (i2.LT.ilast) THEN
       jj=0
       DO WHILE (jj.EQ.0)
          ijs=ijs+1
          IF (ijs.GT.ilast) THEN
             mindg=mindg+1
             ijs=ifirst
          END IF
          ijs1=ijs
          IF (iperm(ijs1).LT.0 .AND. ia(ijs1+1)-ia(ijs1).EQ.mindg) &
               jj=ijs1
       END DO
       i2=i2+1
       GOTO 15
    END IF
    RETURN
  END SUBROUTINE zagmg_setCMK
  SUBROUTINE zagmg_prepareagg(n,a,ja,ia,idiag,ind2,iperm,si,ndd,l,iext)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind2(n),iperm(n)
    INTEGER, OPTIONAL :: iext(*)
    INTEGER :: ndd, l
    COMPLEX(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) , TARGET :: si(n)
    REAL(kind(0.0d0)) :: checkddl,oda,odm,ods,vald
    INTEGER :: i,j,jj,jk,jd,k,kk,j1,j2,i1,i2,ijs,ijs1,ijs2,dg,kdim,nnegrcs
    INTEGER :: ifirst,ilast,i0,kb
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: odmax,odabs,osi
    IF (.NOT.spd) THEN
       checkddl=checkddJ
    ELSE
       checkddl=checkddB
    END IF
    IF (spd .AND. n.EQ.0) THEN
         ndd=0
         RETURN
    END IF
    IF (.NOT.spd) THEN
       IF (checkdd.GT.0) THEN
          ALLOCATE(odmax(n),odabs(n))
          odabs(1:n)=0.0d0
       ELSE
          ALLOCATE(odmax(n))
       END IF
       osi => si(1:n)
       si(1:n)=0.0d0
       odmax(1:n)=0.0d0
       ifirst=1
       ilast=n
        DO i=ilast,ifirst,-1
          j =ia(i)
          jd=idiag(i)
          jj=ia(i+1)-1
          DO k = j,jd-1
            jk=ja(k)
             osi(jk)=osi(jk)+dble(a(k))
             odmax(jk)=max(odmax(jk),dble(a(k)))
             IF (checkdd.GT.0) odabs(jk)=odabs(jk)+abs(a(k))
          ENDDO
          DO k = jd+1,jj
            jk=ja(k)
             osi(jk)=osi(jk)+dble(a(k))
             odmax(jk)=max(odmax(jk),dble(a(k)))
             IF (checkdd.GT.0) odabs(jk)=odabs(jk)+abs(a(k))
          ENDDO
        ENDDO
    END IF
    ndd=0
    nnegrcs=0
    ifirst=1
    ilast=n
     DO i=ifirst,ilast
       j1 =ia(i)
       jd=idiag(i)
       j2 = ia (i+1)-1
       vald = dble(a(jd))
       odm=0.0d0
       oda=0.0d0
       ods=0.0d0
       DO kk = j1,jd-1
          ods=ods+dble(a(kk))
          odm=max(odm,dble(a(kk)))
          IF (checkdd.GT.0) oda=oda+abs(a(kk))
       ENDDO
       DO kk = jd+1,j2
          ods=ods+dble(a(kk))
          odm=max(odm,dble(a(kk)))
          IF (checkdd.GT.0) oda=oda+abs(a(kk))
       ENDDO
       IF (.NOT.spd) THEN
          ods=(osi(i)+ods)/2
          odm=max(odm,odmax(i))
          IF (checkdd.GT.0) oda=(oda+odabs(i))/2
       END IF
       IF ((vald+ods).LT.-repsmach*ABS(vald)) nnegrcs=nnegrcs+1
       si(i)=-ods
       IF ( (checkdd.GT.0 .AND. vald.GT.checkddl*oda)     &
            .OR. (checkdd.LT.0 .AND. vald.GT.checkddl*abs(ods)) ) THEN
          ind2(i)=0
          ndd=ndd+1
       ELSE
          ind2(i)=-1
          IF (odm.GT.trspos*vald) iperm(i)=0
       ENDIF
     END DO
    zerors=.FALSE.
    IF (nnegrcs.GT.fracnegrcsum*n) THEN
       zerors=.TRUE.
       ndd=0
       ind2(1:n)=-1
    END IF
100 CONTINUE
    IF (.NOT.spd) THEN
       IF (checkdd.GT.0) THEN
          DEALLOCATE(odmax,odabs)
       ELSE
          DEALLOCATE(odmax)
       END IF
    END IF
    RETURN
  END SUBROUTINE zagmg_prepareagg
  SUBROUTINE zagmg_findpairs_GF(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,m1,lcg1,a1,ja1,ia1,idiag1,si1,rtent,jtent )
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    COMPLEX(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: m1,ja1(*),ia1(*),idiag1(*),jtent(*),lcg1(m1,n)
    COMPLEX(kind(0.0d0)) :: a1(*)
    REAL(kind(0.0d0)) :: si1(*),rtent(*)
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
    ifirst=.TRUE.
  1 CONTINUE
    kaptg=imult*kaptg_dampJac
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
      DO WHILE (nmark.LT.nt)
       isel=ijs
       ijs=ijs+1
       IF (ind(isel).GE.0) CYCLE
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       IF (lcg1(2,isel).EQ.0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       ntentleft=0
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i.EQ.idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j).GE.0) CYCLE
          IF(lcg1(2,j).EQ.0) CYCLE
          kk=0
          IF (i.LT.idiag(isel)) THEN
             j2=ia(j+1)-1
             DO jk=idiag(j)+1,j2
                IF (ja(jk).EQ.isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ELSE
             DO jk=ia(j),idiag(j)-1
                IF (ja(jk).EQ.isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ENDIF
          vals=-dble(a(i))/2
          IF(kk.NE.0) vals=vals-dble(a(kk))/2
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
             eta1=2*si(isel)
             eta2=2*si(j)
          ELSE
             rsi=-si(isel)+dble(a(idiag(isel)))
             rsj=-si(j)+dble(a(idiag(j)))
             eta1=2*dble(a(idiag(isel)))
             eta2=2*dble(a(idiag(j)))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          IF (sig1.GT.0.0d0) THEN
             del1=rsi
          ELSE
             del1=rsi+2*sig1
          END IF
          IF (sig2.GT.0.0d0) THEN
             del2=rsj
          ELSE
             del2=rsj+2*sig2
          END IF
          IF (vals.GT.0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1).LT.epsr .AND. ABS(del2).LT.epsr) THEN
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1).LT.epsr) THEN
                IF (del2.LT.-epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2).LT.epsr) THEN
                IF (del1.LT.-epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12.LT.-epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp.LT.0.0d0) CYCLE
                valp=((eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1.LE.0.0d0 .OR. del2.LE.0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp.LT.0.0d0) CYCLE
             vals=(eta1*eta2)/(eta1+eta2)
             valp=vals/valp
          END IF
          IF (valp.GT.kaptg) CYCLE
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          ntentleft=ntentleft+1
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          rtent(ntentleft)=tent
          jtent(ntentleft)=j
          CYCLE
9         CONTINUE
          rtent(ntentleft)=val
          jtent(ntentleft)=ipair
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair.EQ.0) GOTO 25
20     CONTINUE
       CALL zagmg_checktentagg_GF
       IF (.NOT.acc) THEN
          ipair = 0
          IF (ntentleft.GT.0) THEN
             i=1
             j=1
             DO WHILE (i.LE.ntentleft)
                IF (jtent(j).GT.0) THEN
                   tent=rtent(j)
                   IF (ipair.EQ.0) GOTO 22
                   IF (16*(tent-val).LT.-1) GOTO 22
                   IF (16*(tent-val).LT.1 .AND. j.LT.ipair) GOTO 22
                   GOTO 23
22                 CONTINUE
                   val=tent
                   ipair=jtent(j)
                   ijtent=j
23                 CONTINUE
                   i=i+1
                END IF
                j=j+1
             END DO
             ntentleft=ntentleft-1
             jtent(ijtent)=0
             GOTO 20
          END IF
       END IF
25     CONTINUE
       IF (ipair.EQ.0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    RETURN
  CONTAINS
  SUBROUTINE zagmg_checktentagg_GF
    INTEGER , PARAMETER :: mm=max(2**(npass+1),8)
    REAL(kind(0.0d0)) :: W(mm,mm), sig(mm), AGe(mm), v(mm)
    REAL(kind(0.0d0)) :: alpha, alp, tmp, beta, beta0, f1, f2
    INTEGER :: j,jj,k,l,m, setdim, l2, k2
    INTEGER  :: info, setdim1
    INTEGER :: set(mm), l1, wdthT
    REAL(kind(0.0d0)) :: T
    LOGICAL :: exc
    IF (m1.eq.2) THEN
       IF (lcg1(2,isel).LT.0) THEN
          IF (lcg1(2,ipair).LT.0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             setdim=2
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             set(3)=lcg1(2,ipair)
             setdim=3
          END IF
          l1=1
       ELSE
          IF (lcg1(2,ipair).LT.0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             setdim=3
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             set(4)=lcg1(2,ipair)
             setdim=4
          END IF
          l1=2
       END IF
    ELSE
       l1=m1
       IF (lcg1(m1,isel).LT.0) l1=-lcg1(m1,isel)
       set(1:l1)=lcg1(1:l1,isel)
       l2=m1
       IF (lcg1(m1,ipair).LT.0) l2=-lcg1(m1,ipair)
       set(l1+1:l1+l2)=lcg1(1:l2,ipair)
       setdim=l1+l2
    END IF
    exc=.TRUE.
    DO WHILE(exc)
       exc=.FALSE.
       DO l=2,SetDim
          IF( set(l).LT.set(l-1) )THEN
             jj=set(l)
             set(l)=set(l-1)
             set(l-1)=jj
             exc=.TRUE.
          END IF
       END DO
    END DO
    DO j=1,SetDim
       jj=Set(j)
       sig(j)=si1(jj)
       IF (zerors) THEN
          W(j,j)=sig(j)
          AGe(j)=0.0d0
       ELSE
          W(j,j)=dble(a1(idiag1(jj)))
          AGe(j)=W(j,j)-sig(j)
       END IF
       l2=j+1
       DO l=l2,SetDim
          W(j,l)=0.0d0
          W(l,j)=0.0d0
       END DO
       k2=ia1(jj+1)-1
       DO k=idiag1(jj)+1,k2
          DO l=l2,SetDim
             m=Set(l)
             IF(ja1(k).EQ.m)THEN
                alpha=dble(a1(k))/2
                W(j,l)=alpha
                W(l,j)=alpha
                EXIT
             END IF
          END DO
       END DO
       DO k=ia1(jj),idiag1(jj)-1
          DO l=1,j-1
             m=Set(l)
             IF(ja1(k).EQ.m)THEN
                alpha=dble(a1(k))/2
                W(j,l)=W(j,l)+alpha
                W(l,j)=W(j,l)
                EXIT
             END IF
          END DO
       END DO
    END DO
    DO j=1,SetDim
       DO k=1,SetDim
          IF (j.ne.k) THEN
             sig(j)=sig(j)+W(j,k)
          END IF
       ENDDO
       IF (sig(j).LT.0.0d0)  AGe(j)=AGe(j)+2*sig(j)
       v(j)=W(j,j)
       W(j,j)=umdbndmum1*W(j,j)-abs(sig(j))
       IF (j.EQ.1) THEN
          beta0=v(j)
          alp=abs(AGe(j))
       ELSE
          beta0=beta0+v(j)
          alp=max(alp,abs(AGe(j)))
       END IF
    END DO
    beta=dbndmum1/beta0
    DO j=1,SetDim
       DO k=1,SetDim
          W(j,k)=W(j,k)+beta*v(j)*v(k)
       END DO
    END DO
    IF (alp.LT.repsmach*beta0) THEN
       SetDim1=SetDim-1
    ELSE
       SetDim1=SetDim
    END IF
    acc=.FALSE.
    SELECT CASE (SetDim1)
    CASE (1)
       GOTO 11
    CASE (2)
       GOTO 12
    CASE (3)
       GOTO 13
    CASE (4)
       GOTO 14
    CASE (5)
       GOTO 15
    CASE (6)
       GOTO 16
    CASE (7)
       GOTO 17
    CASE (8)
       GOTO 18
    CASE DEFAULT
       CALL DPOTRF('U',SetDim1,W,mm,info)
       IF (info.NE.0) RETURN
       GOTO 10
    END SELECT
18  CONTINUE
    IF (W(8,8).LE.0.0d0) RETURN
    W(7,7) = W(7,7) - (W(7,8)/W(8,8)) * W(7,8)
    T = W(6,8)/W(8,8)
    W(6,7) = W(6,7) - T * W(7,8)
    W(6,6) = W(6,6) - T * W(6,8)
    T = W(5,8)/W(8,8)
    W(5,7) = W(5,7) - T * W(7,8)
    W(5,6) = W(5,6) - T * W(6,8)
    W(5,5) = W(5,5) - T * W(5,8)
    T = W(4,8)/W(8,8)
    W(4,7) = W(4,7) - T * W(7,8)
    W(4,6) = W(4,6) - T * W(6,8)
    W(4,5) = W(4,5) - T * W(5,8)
    W(4,4) = W(4,4) - T * W(4,8)
    T = W(3,8)/W(8,8)
    W(3,7) = W(3,7) - T * W(7,8)
    W(3,6) = W(3,6) - T * W(6,8)
    W(3,5) = W(3,5) - T * W(5,8)
    W(3,4) = W(3,4) - T * W(4,8)
    W(3,3) = W(3,3) - T * W(3,8)
    T = W(2,8)/W(8,8)
    W(2,7) = W(2,7) - T * W(7,8)
    W(2,6) = W(2,6) - T * W(6,8)
    W(2,5) = W(2,5) - T * W(5,8)
    W(2,4) = W(2,4) - T * W(4,8)
    W(2,3) = W(2,3) - T * W(3,8)
    W(2,2) = W(2,2) - T * W(2,8)
    T = W(1,8)/W(8,8)
    W(1,7) = W(1,7) - T * W(7,8)
    W(1,6) = W(1,6) - T * W(6,8)
    W(1,5) = W(1,5) - T * W(5,8)
    W(1,4) = W(1,4) - T * W(4,8)
    W(1,3) = W(1,3) - T * W(3,8)
    W(1,2) = W(1,2) - T * W(2,8)
    W(1,1) = W(1,1) - T * W(1,8)
17  CONTINUE
    IF (W(7,7).LE.0.0d0) RETURN
    W(6,6) = W(6,6) - (W(6,7)/W(7,7)) * W(6,7)
    T = W(5,7)/W(7,7)
    W(5,6) = W(5,6) - T * W(6,7)
    W(5,5) = W(5,5) - T * W(5,7)
    T = W(4,7)/W(7,7)
    W(4,6) = W(4,6) - T * W(6,7)
    W(4,5) = W(4,5) - T * W(5,7)
    W(4,4) = W(4,4) - T * W(4,7)
    T = W(3,7)/W(7,7)
    W(3,6) = W(3,6) - T * W(6,7)
    W(3,5) = W(3,5) - T * W(5,7)
    W(3,4) = W(3,4) - T * W(4,7)
    W(3,3) = W(3,3) - T * W(3,7)
    T = W(2,7)/W(7,7)
    W(2,6) = W(2,6) - T * W(6,7)
    W(2,5) = W(2,5) - T * W(5,7)
    W(2,4) = W(2,4) - T * W(4,7)
    W(2,3) = W(2,3) - T * W(3,7)
    W(2,2) = W(2,2) - T * W(2,7)
    T = W(1,7)/W(7,7)
    W(1,6) = W(1,6) - T * W(6,7)
    W(1,5) = W(1,5) - T * W(5,7)
    W(1,4) = W(1,4) - T * W(4,7)
    W(1,3) = W(1,3) - T * W(3,7)
    W(1,2) = W(1,2) - T * W(2,7)
    W(1,1) = W(1,1) - T * W(1,7)
16  CONTINUE
    IF (W(6,6).LE.0.0d0) RETURN
    W(5,5) = W(5,5) - (W(5,6)/W(6,6)) * W(5,6)
    T = W(4,6)/W(6,6)
    W(4,5) = W(4,5) - T * W(5,6)
    W(4,4) = W(4,4) - T * W(4,6)
    T = W(3,6)/W(6,6)
    W(3,5) = W(3,5) - T * W(5,6)
    W(3,4) = W(3,4) - T * W(4,6)
    W(3,3) = W(3,3) - T * W(3,6)
    T = W(2,6)/W(6,6)
    W(2,5) = W(2,5) - T * W(5,6)
    W(2,4) = W(2,4) - T * W(4,6)
    W(2,3) = W(2,3) - T * W(3,6)
    W(2,2) = W(2,2) - T * W(2,6)
    T = W(1,6)/W(6,6)
    W(1,5) = W(1,5) - T * W(5,6)
    W(1,4) = W(1,4) - T * W(4,6)
    W(1,3) = W(1,3) - T * W(3,6)
    W(1,2) = W(1,2) - T * W(2,6)
    W(1,1) = W(1,1) - T * W(1,6)
15  CONTINUE
    IF (W(5,5).LE.0.0d0) RETURN
    W(4,4) = W(4,4) - (W(4,5)/W(5,5)) * W(4,5)
    T = W(3,5)/W(5,5)
    W(3,4) = W(3,4) - T * W(4,5)
    W(3,3) = W(3,3) - T * W(3,5)
    T = W(2,5)/W(5,5)
    W(2,4) = W(2,4) - T * W(4,5)
    W(2,3) = W(2,3) - T * W(3,5)
    W(2,2) = W(2,2) - T * W(2,5)
    T = W(1,5)/W(5,5)
    W(1,4) = W(1,4) - T * W(4,5)
    W(1,3) = W(1,3) - T * W(3,5)
    W(1,2) = W(1,2) - T * W(2,5)
    W(1,1) = W(1,1) - T * W(1,5)
14  CONTINUE
    IF (W(4,4).LE.0.0d0) RETURN
    W(3,3) = W(3,3) - (W(3,4)/W(4,4)) * W(3,4)
    T = W(2,4)/W(4,4)
    W(2,3) = W(2,3) - T * W(3,4)
    W(2,2) = W(2,2) - T * W(2,4)
    T = W(1,4)/W(4,4)
    W(1,3) = W(1,3) - T * W(3,4)
    W(1,2) = W(1,2) - T * W(2,4)
    W(1,1) = W(1,1) - T * W(1,4)
13  CONTINUE
    IF (W(3,3).LE.0.0d0) RETURN
    W(2,2) = W(2,2) - (W(2,3)/W(3,3)) * W(2,3)
    T = W(1,3)/W(3,3)
    W(1,2) = W(1,2) - T * W(2,3)
    W(1,1) = W(1,1) - T * W(1,3)
12  CONTINUE
    IF (W(2,2).LE.0.0d0) RETURN
    W(1,1) = W(1,1) - (W(1,2)/W(2,2)) * W(1,2)
11  CONTINUE
    IF (W(1,1).LE.0.0d0) RETURN
10  CONTINUE
    acc=.TRUE.
    RETURN
  END SUBROUTINE zagmg_checktentagg_GF
  END SUBROUTINE zagmg_findpairs_GF
  SUBROUTINE zagmg_findpairs_SF(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,m1,lcg1,a1,ja1,ia1,idiag1,si1,rtent,jtent )
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    COMPLEX(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: m1,ja1(*),ia1(*),idiag1(*),jtent(*),lcg1(m1,n)
    COMPLEX(kind(0.0d0)) :: a1(*)
    REAL(kind(0.0d0)) :: si1(*),rtent(*)
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
    ifirst=.TRUE.
  1 CONTINUE
    kaptg=imult*kaptg_blocdia
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
      DO WHILE (nmark.LT.nt)
       isel=ijs
       ijs=ijs+1
       IF (ind(isel).GE.0) CYCLE
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       IF (lcg1(2,isel).EQ.0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       ntentleft=0
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i.EQ.idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j).GE.0) CYCLE
          IF(lcg1(2,j).EQ.0) CYCLE
          vals=-dble(a(i))
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
          ELSE
             rsi=-si(isel)+dble(a(idiag(isel)))
             rsj=-si(j)+dble(a(idiag(j)))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          epsr=repsmach*ABS(vals)
          IF (sig1.GT.0.0d0) THEN
             del1=rsi
             eta1=rsi+2*sig1
          ELSE
             del1=rsi+2*sig1
             eta1=rsi
          END IF
          IF (eta1.LT.-epsr) CYCLE
          IF (sig2.GT.0.0d0) THEN
             del2=rsj
             eta2=rsj+2*sig2
          ELSE
             del2=rsj+2*sig2
             eta2=rsj
          END IF
          IF (eta2.LT.-epsr) CYCLE
          IF (vals.GT.0.0d0) THEN
             IF (ABS(del1).LT.epsr .AND. ABS(del2).LT.epsr) THEN
               valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1).LT.epsr) THEN
                IF (del2.LT.-epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2).LT.epsr) THEN
                IF (del1.LT.-epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12.LT.-epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp.LT.0.0d0) CYCLE
                valp=(vals+(eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1.LE.0.0d0 .OR. del2.LE.0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp.LT.0.0d0) CYCLE
             vals=vals+(eta1*eta2)/(eta1+eta2)
             IF (vals.LT.0.0d0) CYCLE
             valp=vals/valp
          END IF
          IF (valp.GT.kaptg) CYCLE
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          ntentleft=ntentleft+1
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          rtent(ntentleft)=tent
          jtent(ntentleft)=j
          CYCLE
9         CONTINUE
          rtent(ntentleft)=val
          jtent(ntentleft)=ipair
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair.EQ.0) GOTO 25
20     CONTINUE
       CALL zagmg_checktentagg_SF
       IF (.NOT.acc) THEN
          ipair = 0
          IF (ntentleft.GT.0) THEN
             i=1
             j=1
             DO WHILE (i.LE.ntentleft)
                IF (jtent(j).GT.0) THEN
                   tent=rtent(j)
                   IF (ipair.EQ.0) GOTO 22
                   IF (16*(tent-val).LT.-1) GOTO 22
                   IF (16*(tent-val).LT.1 .AND. j.LT.ipair) GOTO 22
                   GOTO 23
22                 CONTINUE
                   val=tent
                   ipair=jtent(j)
                   ijtent=j
23                 CONTINUE
                   i=i+1
                END IF
                j=j+1
             END DO
             ntentleft=ntentleft-1
             jtent(ijtent)=0
             GOTO 20
          END IF
       END IF
25     CONTINUE
       IF (ipair.EQ.0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    RETURN
  CONTAINS
  SUBROUTINE zagmg_checktentagg_SF
    INTEGER , PARAMETER :: mm=max(2**(npass+1),8)
    REAL(kind(0.0d0)) :: W(mm,mm), sig(mm), AGe(mm), v(mm)
    REAL(kind(0.0d0)) :: alpha, alp, tmp, beta, beta0, f1, f2
    INTEGER :: j,jj,k,l,m, setdim, l2, k2
    INTEGER  :: info, setdim1
    INTEGER :: set(mm), l1, wdthT
    REAL(kind(0.0d0)) :: T
    LOGICAL :: exc
    IF (m1.eq.2) THEN
       IF (lcg1(2,isel).LT.0) THEN
          IF (lcg1(2,ipair).LT.0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             setdim=2
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             set(3)=lcg1(2,ipair)
             setdim=3
          END IF
          l1=1
       ELSE
          IF (lcg1(2,ipair).LT.0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             setdim=3
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             set(4)=lcg1(2,ipair)
             setdim=4
          END IF
          l1=2
       END IF
    ELSE
       l1=m1
       IF (lcg1(m1,isel).LT.0) l1=-lcg1(m1,isel)
       set(1:l1)=lcg1(1:l1,isel)
       l2=m1
       IF (lcg1(m1,ipair).LT.0) l2=-lcg1(m1,ipair)
       set(l1+1:l1+l2)=lcg1(1:l2,ipair)
       setdim=l1+l2
    END IF
    exc=.TRUE.
    DO WHILE(exc)
       exc=.FALSE.
       DO l=2,SetDim
          IF( set(l).LT.set(l-1) )THEN
             jj=set(l)
             set(l)=set(l-1)
             set(l-1)=jj
             exc=.TRUE.
          END IF
       END DO
    END DO
    DO j=1,SetDim
       jj=Set(j)
       sig(j)=si1(jj)
       IF (zerors) THEN
          W(j,j)=sig(j)
          AGe(j)=0.0d0
       ELSE
          W(j,j)=dble(a1(idiag1(jj)))
          AGe(j)=W(j,j)-sig(j)
       END IF
       l2=j+1
       DO l=l2,SetDim
          W(j,l)=0.0d0
          W(l,j)=0.0d0
       END DO
       k2=ia1(jj+1)-1
       DO k=idiag1(jj)+1,k2
          DO l=l2,SetDim
             m=Set(l)
             IF(ja1(k).EQ.m)THEN
                alpha=dble(a1(k))
                W(j,l)=alpha
                W(l,j)=alpha
                EXIT
             END IF
          END DO
       END DO
    END DO
    DO j=1,SetDim
       DO k=1,SetDim
          IF (j.ne.k) THEN
             sig(j)=sig(j)+W(j,k)
          END IF
       ENDDO
       IF (sig(j).LT.0.0d0)  AGe(j)=AGe(j)+2*sig(j)
       W(j,j)=W(j,j)-abs(sig(j))
       tmp=2*abs(sig(j))
       W(j,j)=W(j,j)-bndmum1m1*tmp
       v(j)=tmp+AGe(j)
       IF (j.EQ.1) THEN
          beta0=v(j)
          alp=abs(AGe(j))
       ELSE
          beta0=beta0+v(j)
          alp=max(alp,abs(AGe(j)))
       END IF
    END DO
    beta=bndmum1m1/beta0
    DO j=1,SetDim
       DO k=1,SetDim
          W(j,k)=W(j,k)+beta*v(j)*v(k)
       END DO
    END DO
    IF (alp.LT.repsmach*beta0) THEN
       SetDim1=SetDim-1
    ELSE
       SetDim1=SetDim
    END IF
    acc=.FALSE.
    SELECT CASE (SetDim1)
    CASE (1)
       GOTO 11
    CASE (2)
       GOTO 12
    CASE (3)
       GOTO 13
    CASE (4)
       GOTO 14
    CASE (5)
       GOTO 15
    CASE (6)
       GOTO 16
    CASE (7)
       GOTO 17
    CASE (8)
       GOTO 18
    CASE DEFAULT
       CALL DPOTRF('U',SetDim1,W,mm,info)
       IF (info.NE.0) RETURN
       GOTO 10
    END SELECT
18  CONTINUE
    IF (W(8,8).LE.0.0d0) RETURN
    W(7,7) = W(7,7) - (W(7,8)/W(8,8)) * W(7,8)
    T = W(6,8)/W(8,8)
    W(6,7) = W(6,7) - T * W(7,8)
    W(6,6) = W(6,6) - T * W(6,8)
    T = W(5,8)/W(8,8)
    W(5,7) = W(5,7) - T * W(7,8)
    W(5,6) = W(5,6) - T * W(6,8)
    W(5,5) = W(5,5) - T * W(5,8)
    T = W(4,8)/W(8,8)
    W(4,7) = W(4,7) - T * W(7,8)
    W(4,6) = W(4,6) - T * W(6,8)
    W(4,5) = W(4,5) - T * W(5,8)
    W(4,4) = W(4,4) - T * W(4,8)
    T = W(3,8)/W(8,8)
    W(3,7) = W(3,7) - T * W(7,8)
    W(3,6) = W(3,6) - T * W(6,8)
    W(3,5) = W(3,5) - T * W(5,8)
    W(3,4) = W(3,4) - T * W(4,8)
    W(3,3) = W(3,3) - T * W(3,8)
    T = W(2,8)/W(8,8)
    W(2,7) = W(2,7) - T * W(7,8)
    W(2,6) = W(2,6) - T * W(6,8)
    W(2,5) = W(2,5) - T * W(5,8)
    W(2,4) = W(2,4) - T * W(4,8)
    W(2,3) = W(2,3) - T * W(3,8)
    W(2,2) = W(2,2) - T * W(2,8)
    T = W(1,8)/W(8,8)
    W(1,7) = W(1,7) - T * W(7,8)
    W(1,6) = W(1,6) - T * W(6,8)
    W(1,5) = W(1,5) - T * W(5,8)
    W(1,4) = W(1,4) - T * W(4,8)
    W(1,3) = W(1,3) - T * W(3,8)
    W(1,2) = W(1,2) - T * W(2,8)
    W(1,1) = W(1,1) - T * W(1,8)
17  CONTINUE
    IF (W(7,7).LE.0.0d0) RETURN
    W(6,6) = W(6,6) - (W(6,7)/W(7,7)) * W(6,7)
    T = W(5,7)/W(7,7)
    W(5,6) = W(5,6) - T * W(6,7)
    W(5,5) = W(5,5) - T * W(5,7)
    T = W(4,7)/W(7,7)
    W(4,6) = W(4,6) - T * W(6,7)
    W(4,5) = W(4,5) - T * W(5,7)
    W(4,4) = W(4,4) - T * W(4,7)
    T = W(3,7)/W(7,7)
    W(3,6) = W(3,6) - T * W(6,7)
    W(3,5) = W(3,5) - T * W(5,7)
    W(3,4) = W(3,4) - T * W(4,7)
    W(3,3) = W(3,3) - T * W(3,7)
    T = W(2,7)/W(7,7)
    W(2,6) = W(2,6) - T * W(6,7)
    W(2,5) = W(2,5) - T * W(5,7)
    W(2,4) = W(2,4) - T * W(4,7)
    W(2,3) = W(2,3) - T * W(3,7)
    W(2,2) = W(2,2) - T * W(2,7)
    T = W(1,7)/W(7,7)
    W(1,6) = W(1,6) - T * W(6,7)
    W(1,5) = W(1,5) - T * W(5,7)
    W(1,4) = W(1,4) - T * W(4,7)
    W(1,3) = W(1,3) - T * W(3,7)
    W(1,2) = W(1,2) - T * W(2,7)
    W(1,1) = W(1,1) - T * W(1,7)
16  CONTINUE
    IF (W(6,6).LE.0.0d0) RETURN
    W(5,5) = W(5,5) - (W(5,6)/W(6,6)) * W(5,6)
    T = W(4,6)/W(6,6)
    W(4,5) = W(4,5) - T * W(5,6)
    W(4,4) = W(4,4) - T * W(4,6)
    T = W(3,6)/W(6,6)
    W(3,5) = W(3,5) - T * W(5,6)
    W(3,4) = W(3,4) - T * W(4,6)
    W(3,3) = W(3,3) - T * W(3,6)
    T = W(2,6)/W(6,6)
    W(2,5) = W(2,5) - T * W(5,6)
    W(2,4) = W(2,4) - T * W(4,6)
    W(2,3) = W(2,3) - T * W(3,6)
    W(2,2) = W(2,2) - T * W(2,6)
    T = W(1,6)/W(6,6)
    W(1,5) = W(1,5) - T * W(5,6)
    W(1,4) = W(1,4) - T * W(4,6)
    W(1,3) = W(1,3) - T * W(3,6)
    W(1,2) = W(1,2) - T * W(2,6)
    W(1,1) = W(1,1) - T * W(1,6)
15  CONTINUE
    IF (W(5,5).LE.0.0d0) RETURN
    W(4,4) = W(4,4) - (W(4,5)/W(5,5)) * W(4,5)
    T = W(3,5)/W(5,5)
    W(3,4) = W(3,4) - T * W(4,5)
    W(3,3) = W(3,3) - T * W(3,5)
    T = W(2,5)/W(5,5)
    W(2,4) = W(2,4) - T * W(4,5)
    W(2,3) = W(2,3) - T * W(3,5)
    W(2,2) = W(2,2) - T * W(2,5)
    T = W(1,5)/W(5,5)
    W(1,4) = W(1,4) - T * W(4,5)
    W(1,3) = W(1,3) - T * W(3,5)
    W(1,2) = W(1,2) - T * W(2,5)
    W(1,1) = W(1,1) - T * W(1,5)
14  CONTINUE
    IF (W(4,4).LE.0.0d0) RETURN
    W(3,3) = W(3,3) - (W(3,4)/W(4,4)) * W(3,4)
    T = W(2,4)/W(4,4)
    W(2,3) = W(2,3) - T * W(3,4)
    W(2,2) = W(2,2) - T * W(2,4)
    T = W(1,4)/W(4,4)
    W(1,3) = W(1,3) - T * W(3,4)
    W(1,2) = W(1,2) - T * W(2,4)
    W(1,1) = W(1,1) - T * W(1,4)
13  CONTINUE
    IF (W(3,3).LE.0.0d0) RETURN
    W(2,2) = W(2,2) - (W(2,3)/W(3,3)) * W(2,3)
    T = W(1,3)/W(3,3)
    W(1,2) = W(1,2) - T * W(2,3)
    W(1,1) = W(1,1) - T * W(1,3)
12  CONTINUE
    IF (W(2,2).LE.0.0d0) RETURN
    W(1,1) = W(1,1) - (W(1,2)/W(2,2)) * W(1,2)
11  CONTINUE
    IF (W(1,1).LE.0.0d0) RETURN
10  CONTINUE
    acc=.TRUE.
    RETURN
  END SUBROUTINE zagmg_checktentagg_SF
  END SUBROUTINE zagmg_findpairs_SF
  SUBROUTINE zagmg_findpairs_GI(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,ipc )
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    COMPLEX(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: ipc(n)
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
    ifirst=.TRUE.
  1 CONTINUE
    kaptg=imult*kaptg_dampJac
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
      DO WHILE (nmark.LT.nt)
       isel=ijs
       ijs=ijs+1
       IF (ind(isel).EQ.0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       IF (ind(isel).GE.0) CYCLE
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       IF (ipc(isel).EQ.0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i.EQ.idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j).GE.0) CYCLE
          IF(ipc(j).EQ.0) CYCLE
          kk=0
          IF (i.LT.idiag(isel)) THEN
             j2=ia(j+1)-1
             DO jk=idiag(j)+1,j2
                IF (ja(jk).EQ.isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ELSE
             DO jk=ia(j),idiag(j)-1
                IF (ja(jk).EQ.isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ENDIF
          vals=-dble(a(i))/2
          IF(kk.NE.0) vals=vals-dble(a(kk))/2
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
             eta1=2*si(isel)
             eta2=2*si(j)
          ELSE
             rsi=-si(isel)+dble(a(idiag(isel)))
             rsj=-si(j)+dble(a(idiag(j)))
             eta1=2*dble(a(idiag(isel)))
             eta2=2*dble(a(idiag(j)))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          IF (sig1.GT.0.0d0) THEN
             del1=rsi
          ELSE
             del1=rsi+2*sig1
          END IF
          IF (sig2.GT.0.0d0) THEN
             del2=rsj
          ELSE
             del2=rsj+2*sig2
          END IF
          IF (vals.GT.0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1).LT.epsr .AND. ABS(del2).LT.epsr) THEN
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1).LT.epsr) THEN
                IF (del2.LT.-epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2).LT.epsr) THEN
                IF (del1.LT.-epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12.LT.-epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp.LT.0.0d0) CYCLE
                valp=((eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1.LE.0.0d0 .OR. del2.LE.0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp.LT.0.0d0) CYCLE
             vals=(eta1*eta2)/(eta1+eta2)
             valp=vals/valp
          END IF
          IF (valp.GT.kaptg) CYCLE
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair.EQ.0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    IF (ifirst .AND. 3*nc.GT.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       ifirst=.FALSE.
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_dampJac
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',es9.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE zagmg_findpairs_GI
  SUBROUTINE zagmg_findpairs_SI(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,ipc )
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    COMPLEX(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: ipc(n)
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
    ifirst=.TRUE.
  1 CONTINUE
    kaptg=imult*kaptg_blocdia
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
      DO WHILE (nmark.LT.nt)
       isel=ijs
       ijs=ijs+1
       IF (ind(isel).EQ.0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       IF (ind(isel).GE.0) CYCLE
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       IF (ipc(isel).EQ.0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i.EQ.idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j).GE.0) CYCLE
          IF(ipc(j).EQ.0) CYCLE
          vals=-dble(a(i))
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
          ELSE
             rsi=-si(isel)+dble(a(idiag(isel)))
             rsj=-si(j)+dble(a(idiag(j)))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          epsr=repsmach*ABS(vals)
          IF (sig1.GT.0.0d0) THEN
             del1=rsi
             eta1=rsi+2*sig1
          ELSE
             del1=rsi+2*sig1
             eta1=rsi
          END IF
          IF (eta1.LT.-epsr) CYCLE
          IF (sig2.GT.0.0d0) THEN
             del2=rsj
             eta2=rsj+2*sig2
          ELSE
             del2=rsj+2*sig2
             eta2=rsj
          END IF
          IF (eta2.LT.-epsr) CYCLE
          IF (vals.GT.0.0d0) THEN
             IF (ABS(del1).LT.epsr .AND. ABS(del2).LT.epsr) THEN
               valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1).LT.epsr) THEN
                IF (del2.LT.-epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2).LT.epsr) THEN
                IF (del1.LT.-epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12.LT.-epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp.LT.0.0d0) CYCLE
                valp=(vals+(eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1.LE.0.0d0 .OR. del2.LE.0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp.LT.0.0d0) CYCLE
             vals=vals+(eta1*eta2)/(eta1+eta2)
             IF (vals.LT.0.0d0) CYCLE
             valp=vals/valp
          END IF
          IF (valp.GT.kaptg) CYCLE
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair.EQ.0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    IF (ifirst .AND. 3*nc.GT.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       ifirst=.FALSE.
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_blocdia
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',es9.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE zagmg_findpairs_SI
  SUBROUTINE zagmg_findpairs_GI1(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,riperm,iperm )
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    COMPLEX(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: iperm(n),riperm(n)
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
    ifirst=.TRUE.
  1 CONTINUE
    kaptg=imult*kaptg_dampJac
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
      DO WHILE (nmark.LT.nt)
       isel=ijs
       isel=riperm(ijs)
       ijs=ijs+1
       IF (ind(isel).EQ.0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       IF (ind(isel).GE.0) CYCLE
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       IF (iperm(isel).EQ.0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i.EQ.idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j).GE.0) CYCLE
          IF(iperm(j).EQ.0) CYCLE
          kk=0
          IF (i.LT.idiag(isel)) THEN
             j2=ia(j+1)-1
             DO jk=idiag(j)+1,j2
                IF (ja(jk).EQ.isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ELSE
             DO jk=ia(j),idiag(j)-1
                IF (ja(jk).EQ.isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ENDIF
          vals=-dble(a(i))/2
          IF(kk.NE.0) vals=vals-dble(a(kk))/2
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
             eta1=2*si(isel)
             eta2=2*si(j)
          ELSE
             rsi=-si(isel)+dble(a(idiag(isel)))
             rsj=-si(j)+dble(a(idiag(j)))
             eta1=2*dble(a(idiag(isel)))
             eta2=2*dble(a(idiag(j)))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          IF (sig1.GT.0.0d0) THEN
             del1=rsi
          ELSE
             del1=rsi+2*sig1
          END IF
          IF (sig2.GT.0.0d0) THEN
             del2=rsj
          ELSE
             del2=rsj+2*sig2
          END IF
          IF (vals.GT.0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1).LT.epsr .AND. ABS(del2).LT.epsr) THEN
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1).LT.epsr) THEN
                IF (del2.LT.-epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2).LT.epsr) THEN
                IF (del1.LT.-epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12.LT.-epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp.LT.0.0d0) CYCLE
                valp=((eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1.LE.0.0d0 .OR. del2.LE.0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp.LT.0.0d0) CYCLE
             vals=(eta1*eta2)/(eta1+eta2)
             valp=vals/valp
          END IF
          IF (valp.GT.kaptg) CYCLE
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. iperm(j).LT.iperm(ipair))  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair.EQ.0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    IF (ifirst .AND. 3*nc.GT.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       ifirst=.FALSE.
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_dampJac
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',es9.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE zagmg_findpairs_GI1
  SUBROUTINE zagmg_findpairs_SI1(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,riperm,iperm )
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    COMPLEX(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: iperm(n),riperm(n)
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
    ifirst=.TRUE.
  1 CONTINUE
    kaptg=imult*kaptg_blocdia
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
      DO WHILE (nmark.LT.nt)
       isel=ijs
       isel=riperm(ijs)
       ijs=ijs+1
       IF (ind(isel).EQ.0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       IF (ind(isel).GE.0) CYCLE
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       IF (iperm(isel).EQ.0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i.EQ.idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j).GE.0) CYCLE
          IF(iperm(j).EQ.0) CYCLE
          vals=-dble(a(i))
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
          ELSE
             rsi=-si(isel)+dble(a(idiag(isel)))
             rsj=-si(j)+dble(a(idiag(j)))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          epsr=repsmach*ABS(vals)
          IF (sig1.GT.0.0d0) THEN
             del1=rsi
             eta1=rsi+2*sig1
          ELSE
             del1=rsi+2*sig1
             eta1=rsi
          END IF
          IF (eta1.LT.-epsr) CYCLE
          IF (sig2.GT.0.0d0) THEN
             del2=rsj
             eta2=rsj+2*sig2
          ELSE
             del2=rsj+2*sig2
             eta2=rsj
          END IF
          IF (eta2.LT.-epsr) CYCLE
          IF (vals.GT.0.0d0) THEN
             IF (ABS(del1).LT.epsr .AND. ABS(del2).LT.epsr) THEN
               valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1).LT.epsr) THEN
                IF (del2.LT.-epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2).LT.epsr) THEN
                IF (del1.LT.-epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12.LT.-epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp.LT.0.0d0) CYCLE
                valp=(vals+(eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1.LE.0.0d0 .OR. del2.LE.0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp.LT.0.0d0) CYCLE
             vals=vals+(eta1*eta2)/(eta1+eta2)
             IF (vals.LT.0.0d0) CYCLE
             valp=vals/valp
          END IF
          IF (valp.GT.kaptg) CYCLE
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. iperm(j).LT.iperm(ipair))  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair.EQ.0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    IF (ifirst .AND. 3*nc.GT.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       ifirst=.FALSE.
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_blocdia
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',es9.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE zagmg_findpairs_SI1
  SUBROUTINE zagmg_lcgmix(nc,m,lcg1,lcg,lcgn)
    INTEGER :: nc,m,lcg1(m,*),lcg(2,*),lcgn(2*m,*),i,l,l1,l2
    IF (m.eq.2) THEN
       DO i=1,nc
          IF(lcg(2,i).EQ.0) THEN
             lcgn(1,i)=lcg1(1,lcg(1,i))
             lcgn(2,i)=0
             lcgn(4,i)=-1
          ELSE IF(lcg(2,i).LT.0) THEN
             IF (lcg1(2,lcg(1,i)).LT.0) THEN
                lcgn(1,i)=lcg1(1,lcg(1,i))
                lcgn(2,i)=-1
                lcgn(4,i)=-1
             ELSE
                lcgn(1,i)=lcg1(1,lcg(1,i))
                lcgn(2,i)=lcg1(2,lcg(1,i))
                lcgn(4,i)=-2
             END IF
          ELSE
             IF (lcg1(2,lcg(1,i)).LT.0) THEN
                IF (lcg1(2,lcg(2,i)).LT.0) THEN
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(1,lcg(2,i))
                   lcgn(4,i)=-2
                ELSE
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(1,lcg(2,i))
                   lcgn(3,i)=lcg1(2,lcg(2,i))
                   lcgn(4,i)=-3
                END IF
             ELSE
                IF (lcg1(2,lcg(2,i)).LT.0) THEN
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(2,lcg(1,i))
                   lcgn(3,i)=lcg1(1,lcg(2,i))
                   lcgn(4,i)=-3
                ELSE
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(2,lcg(1,i))
                   lcgn(3,i)=lcg1(1,lcg(2,i))
                   lcgn(4,i)=lcg1(2,lcg(2,i))
                END IF
             END IF
          END IF
       END DO
    ELSE
       DO i=1,nc
          IF(lcg(2,i).EQ.0) THEN
             lcgn(1,i)=lcg1(1,lcg(1,i))
             lcgn(2,i)=0
             lcgn(2*m,i)=-1
          ELSE
             lcgn(2,i)=-1
             l1=m
             IF (lcg1(m,lcg(1,i)).LT.0) l1=-lcg1(m,lcg(1,i))
             lcgn(1:l1,i)=lcg1(1:l1,lcg(1,i))
             IF(lcg(2,i).LT.0) THEN
                l=l1
             ELSE
                l2=m
                IF (lcg1(m,lcg(2,i)).LT.0) l2=-lcg1(m,lcg(2,i))
                lcgn(l1+1:l1+l2,i)=lcg1(1:l2,lcg(2,i))
                l=l1+l2
             END IF
             IF(l.LT.2*m) lcgn(2*m,i)=-l
          END IF
       END DO
    END IF
    RETURN
  END SUBROUTINE zagmg_lcgmix
  SUBROUTINE zagmg_setind(nc,ndd,ldd,lcg,m,ind)
    INTEGER :: nc,m,lcg(m,*),nll,ldd(ndd),ind(*),i,k,l
    DO i=1,ndd
       ind(ldd(i))=0
    END DO
    DO i=1,nc
       l=m
       IF (lcg(m,i).LT.0) l=-lcg(m,i)
       DO k=1,l
          ind(lcg(k,i))=i
       END DO
    END DO
    RETURN
  END SUBROUTINE zagmg_setind
  SUBROUTINE zagmg_setcg(n,a,ja,ia,idiag,si,ind,lcg           &
       ,nc,a2,ja2,ia2,idiag2,si2,ysi,maxdg,iw,w,iext,iext2)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind(n),nc,lcg(2,nc)
    INTEGER :: ja2(*),ia2(nc+1),idiag2(nc),maxdg
    INTEGER, TARGET :: iw(2*nc)
    INTEGER, OPTIONAL :: iext(*),iext2(*)
    COMPLEX(kind(0.0d0)) :: a(*),a2(*),w(nc),vald
    REAL(kind(0.0d0)) :: si(n),si2(*),sii
    LOGICAL :: ysi
    INTEGER :: nz,nzu,i,jj,jc,jcol,ki,kb,kf,jpos
    INTEGER, POINTER, DIMENSION(:) :: iw1, iw2
    iw1 => iw(1:nc)
    iw2 => iw(nc+1:2*nc)
    nz = 0
    iw1(1:nc)=0
    maxdg=0
    ia2(1)=1
    DO i = 1,nc
       sii=0.0d0
       vald=cmplx(0,0,kind(0.0d0))
       nzu=0
       DO ki= 1,2
          jj = lcg(ki,i)
          IF (ki.EQ.1 .OR. jj.GT.0) THEN
             IF (ysi) sii=sii+dble(si(jj))
             kf = ia(jj+1) - 1
             DO kb = ia(jj),kf
                jc = ja(kb)
                jcol = ind(jc)
                IF (jcol.GT.0) THEN
                   IF (jcol.LT.i) THEN
                      jpos = iw1(jcol)
                      IF (jpos.EQ.0) THEN
                         nz = nz+1
                         ja2(nz) = jcol
                         iw1(jcol) = nz
                         a2(nz) = a(kb)
                      ELSE
                         a2(jpos) = a2(jpos) + a(kb)
                      ENDIF
                   ELSE IF (jcol.GT.i) THEN
                      jpos = iw1(jcol)
                      IF (jpos.EQ.0) THEN
                         nzu = nzu+1
                         iw2(nzu) = jcol
                         iw1(jcol) = nzu
                         w(nzu) = a(kb)
                      ELSE
                         w(jpos) = w(jpos) + a(kb)
                      ENDIF
                   ELSE
                      vald=vald+a(kb)
                      IF (ysi .AND. jc.NE.jj) sii=sii+dble(a(kb))
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       nz=nz+1
       a2(nz)=vald
       idiag2(i)=nz
       ja2(nz)=i
       a2(nz+1:nz+nzu)=w(1:nzu)
       ja2(nz+1:nz+nzu)=iw2(1:nzu)
       nz=nz+nzu
       maxdg=max(maxdg,nz-ia2(i))
       DO kb = ia2(i), nz
          iw1(ja2(kb))=0
       ENDDO
       IF (ysi) si2(i)=sii
       ia2(i+1)=nz+1
    ENDDO
    RETURN
  END SUBROUTINE zagmg_setcg
  SUBROUTINE zagmg_DIRseq(n,f,ijob,a,ja,ia)
    USE zagmg_mem
    USE zagmg_PARDISO
    IMPLICIT NONE
    INTEGER :: n,ijob
    COMPLEX(kind(0.0d0)), TARGET :: f(n)
    COMPLEX(kind(0.0d0)), OPTIONAL, TARGET :: a(*)
    INTEGER, OPTIONAL :: ia(n+1)
    INTEGER, OPTIONAL, TARGET :: ja(*)
    REAL(kind(0.0d0)) , SAVE :: iflop
    COMPLEX(kind(0.0d0)), SAVE :: aaa
    TYPE(PARDISO_HANDLE), ALLOCATABLE, SAVE :: PT(:)
    INTEGER, SAVE :: IPARM(64)
    COMPLEX(kind(0.0d0)), ALLOCATABLE, SAVE :: X(:)
    INTEGER, SAVE :: MAXFCT, MNUM, MTYPE, PHASE, NRHS, MSGLVL, ERROR
!
!
      TYPE ZMUMPS_ROOT_STRUC
        SEQUENCE
        INTEGER MBLOCK, NBLOCK, NPROW, NPCOL
        INTEGER MYROW, MYCOL
        INTEGER ROOT_SIZE, TOT_ROOT_SIZE
        INTEGER :: CNTXT_BLACS, truc
        INTEGER, DIMENSION(:), POINTER :: RG2L_ROW
        INTEGER, DIMENSION(:), POINTER :: RG2L_COL
        INTEGER , DIMENSION(:), POINTER :: IPIV
        INTEGER, DIMENSION( 9 ) :: DESCRIPTOR, DESCB
        LOGICAL yes, gridinit_done
        INTEGER LPIV, brol
!       Used to access Schur easily from root structure
        COMPLEX*16, DIMENSION(:), POINTER :: SCHUR_POINTER
        INTEGER SCHUR_MLOC, SCHUR_NLOC, SCHUR_LLD, machin
!
!      Data for nullspace/QR
!
        COMPLEX*16, DIMENSION(:), POINTER :: QR_TAU
        DOUBLE PRECISION     QR_RCOND
!
!      Givens rotations
!
        INTEGER MAXG, GIND
        COMPLEX*16, DIMENSION(:),POINTER::GROW, GCOS, GSIN
!
!      RRRLU data
!
        INTEGER ELG_MAX,NULL_MAX
        INTEGER ELIND,EUIND,NLUPDATE,NUUPDATE
        INTEGER,DIMENSION(:),POINTER::PERM_ROW,PERM_COL
        INTEGER,DIMENSION(:),POINTER::ELROW, EUROW, PTREL, PTREU
        COMPLEX*16, DIMENSION(:), POINTER :: ELELG, EUELG, DL
!
      END TYPE ZMUMPS_ROOT_STRUC
      TYPE ZMUMPS_STRUC
        SEQUENCE
!
! This structure contains all parameters
! for the interface to the user, plus internal
! information
!
! *****************
! INPUT PARAMETERS
! *****************
!    -----------------
!    MPI Communicator
!    -----------------
        INTEGER COMM
!    ------------------
!    Problem definition
!    ------------------
!    Solver (SYM=0 unsymmetric,SYM=1 symmetric Positive Definite,
!        SYM=2 general symmetric)
!    Type of parallelism (PAR=1 host working, PAR=0 host not working)
        INTEGER SYM, PAR
        INTEGER JOB
!    --------------------
!    Order of Input matrix
!    --------------------
        INTEGER N
!
!    ----------------------------------------
!    Assembled input matrix : User interface
!    ----------------------------------------
        INTEGER NZ
        COMPLEX*16, DIMENSION(:), POINTER :: A
        INTEGER, DIMENSION(:), POINTER :: IRN, JCN
        DOUBLE PRECISION, DIMENSION(:), POINTER :: COLSCA, ROWSCA, pad0
!
!       ------------------------------------
!       Case of distributed assembled matrix
!       matrix on entry:
!       ------------------------------------
        INTEGER NZ_loc, pad1
        INTEGER, DIMENSION(:), POINTER :: IRN_loc, JCN_loc
        COMPLEX*16, DIMENSION(:), POINTER :: A_loc, pad2
!
!    ----------------------------------------
!    Unassembled input matrix: User interface
!    ----------------------------------------
        INTEGER NELT, pad3
        INTEGER, DIMENSION(:), POINTER :: ELTPTR
        INTEGER, DIMENSION(:), POINTER :: ELTVAR
        COMPLEX*16, DIMENSION(:), POINTER :: A_ELT, pad4
!
!    ---------------------------------------------
!    Symmetric permutation :
!               PERM_IN if given by user (optional)
!    ---------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: PERM_IN
!
!
! ******************
! INPUT/OUTPUT data
! ******************
!    --------------------------------------------------------
!    RHS / SOL_LOC
!    -------------
!       right-hand side and solution
!    -------------------------------------------------------
        COMPLEX*16, DIMENSION(:), POINTER :: RHS, REDRHS
        COMPLEX*16, DIMENSION(:), POINTER :: RHS_SPARSE
        COMPLEX*16, DIMENSION(:), POINTER :: SOL_LOC
        INTEGER, DIMENSION(:), POINTER :: IRHS_SPARSE
        INTEGER, DIMENSION(:), POINTER :: IRHS_PTR
        INTEGER, DIMENSION(:), POINTER :: ISOL_LOC
        INTEGER LRHS, NRHS, NZ_RHS, LSOL_LOC, LREDRHS
        INTEGER pad5
!    ----------------------------
!    Control parameters,
!    statistics and output data
!    ---------------------------
        INTEGER ICNTL(40)
        INTEGER INFO(40)
        INTEGER INFOG(40)
        DOUBLE PRECISION COST_SUBTREES
        DOUBLE PRECISION CNTL(15)
        DOUBLE PRECISION RINFO(20)
        DOUBLE PRECISION RINFOG(20)
!    ---------------------------------------------------------
!    Permutations computed during analysis:
!       SYM_PERM: Symmetric permutation
!       UNS_PERM: Column permutations (optionnal)
!    ---------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: SYM_PERM, UNS_PERM
!
!    -------------------------------------
!    Case of distributed matrix on entry:
!    ZMUMPS potentially provides mapping
!    -------------------------------------
        INTEGER, DIMENSION(:), POINTER :: MAPPING
!
!    -------------------------------
!    Deficiency and null space basis
!    -------------------------------
        COMPLEX*16, DIMENSION(:,:), POINTER :: NULL_SPACE
        INTEGER Deficiency, pad6
!    -----
!    Schur
!    -----
        INTEGER NPROW, NPCOL, MBLOCK, NBLOCK
        INTEGER SCHUR_MLOC, SCHUR_NLOC, SCHUR_LLD
        INTEGER SIZE_SCHUR
        INTEGER, DIMENSION(:), POINTER :: LISTVAR_SCHUR
        COMPLEX*16, DIMENSION(:), POINTER :: SCHUR
        COMPLEX*16, DIMENSION(:), POINTER :: SCHUR_CINTERFACE
!    --------------
!    Version number
!    --------------
        CHARACTER(LEN=16) VERSION_NUMBER
!    -----------
!    Out-of-core
!    -----------
        CHARACTER(LEN=256) :: OOC_TMPDIR
        CHARACTER(LEN=64) :: OOC_PREFIX
!    ------------------------------------------
!    To save the matrix in matrix market format
!    ------------------------------------------
        CHARACTER(LEN=256) WRITE_PROBLEM
!
!
! **********************
! INTERNAL Working data
! *********************
        INTEGER INST_Number
!       For MPI
        INTEGER COMM_NODES, MYID_NODES, COMM_LOAD
        INTEGER  MYID, NPROCS, NSLAVES
        INTEGER ASS_IRECV
        INTEGER, DIMENSION(:), POINTER :: POIDS
        INTEGER LBUFR
        INTEGER LBUFR_BYTES
        INTEGER, DIMENSION(:), POINTER ::  BUFR
!       For analysis/facto/solve phases
        INTEGER MAXIS1, pad7
        INTEGER KEEP(500)
        INTEGER*8 KEEP8(150)
!       IS is used for the factors + workspace for contrib. blocks
        INTEGER, DIMENSION(:), POINTER :: IS
!       is1 (maxis1) contains working arrays computed
!       and used only during analysis
        INTEGER, DIMENSION(:), POINTER :: IS1
!       The following data/arrays are computed during the analysis
!       phase and used during the factorization and solve phases.
        INTEGER LNA
        INTEGER NBSA
        INTEGER,POINTER,DIMENSION(:)::STEP, NE_STEPS, ND_STEPS
!  Info for pruning tree
        INTEGER,POINTER,DIMENSION(:)::Step2node
!  ---------------------
        INTEGER,POINTER,DIMENSION(:)::FRERE_STEPS, DAD_STEPS
        INTEGER,POINTER,DIMENSION(:)::FILS, PTRAR, FRTPTR, FRTELT
        INTEGER,POINTER,DIMENSION(:)::NA, PROCNODE_STEPS
!       The two pointer arrays computed in facto and used by the solve
!          (except the factors) are PTLUST_S and PTRFAC.
        INTEGER, DIMENSION(:), POINTER :: PTLUST_S
        INTEGER(8), DIMENSION(:), POINTER :: PTRFAC
!       main real working arrays for factorization/solve phases
        COMPLEX*16, DIMENSION(:), POINTER :: S
!       Information on mapping
        INTEGER, DIMENSION(:), POINTER :: PROCNODE
!       Input matrix ready for numerical assembly
!           -arrowhead format in case of assembled matrix
!           -element format otherwise
        INTEGER, DIMENSION(:), POINTER :: INTARR
        COMPLEX*16, DIMENSION(:), POINTER :: DBLARR
!       Element entry: internal data
        INTEGER NELT_LOC, LELTVAR, NA_ELT, pad8
        INTEGER, DIMENSION(:), POINTER :: ELTPROC
!       Candidates and node partitionning
        INTEGER, DIMENSION(:,:), POINTER :: CANDIDATES
        INTEGER, DIMENSION(:),   POINTER :: ISTEP_TO_INIV2
        INTEGER, DIMENSION(:),   POINTER :: FUTURE_NIV2
        INTEGER, DIMENSION(:,:), POINTER :: TAB_POS_IN_PERE
        LOGICAL, DIMENSION(:), POINTER :: I_AM_CAND
!       For heterogeneous architecture
        INTEGER, DIMENSION(:), POINTER :: MEM_DIST
!       Compressed RHS
        INTEGER, DIMENSION(:),   POINTER :: POSINRHSCOMP
        COMPLEX*16, DIMENSION(:), POINTER :: RHSCOMP
!       For C interface
!   Info on the subtrees to be used during factorization
        DOUBLE PRECISION, DIMENSION(:),   POINTER :: MEM_SUBTREE
        INTEGER, DIMENSION(:),   POINTER :: MY_ROOT_SBTR
        INTEGER, DIMENSION(:),   POINTER :: MY_FIRST_LEAF
        INTEGER, DIMENSION(:),   POINTER :: MY_NB_LEAF
        INTEGER, DIMENSION(:),   POINTER :: DEPTH_FIRST
        DOUBLE PRECISION, DIMENSION(:),   POINTER :: COST_TRAV
        INTEGER NBSA_LOCAL, zwave
        INTEGER(8) :: MAX_SURF_MASTER
        INTEGER :: LWK_USER, zozo
        COMPLEX*16, DIMENSION(:), POINTER :: WK_USER
!    For simulating parallel out-of-core stack.
        DOUBLE PRECISION, DIMENSION(:),POINTER ::CB_SON_SIZE
!   Instance number used/managed by the C/F77 interface
        INTEGER INSTANCE_NUMBER
!    OOC management data that must persist from factorization to solve.
        INTEGER OOC_MAX_NB_NODES_FOR_ZONE
        INTEGER, DIMENSION(:,:),   POINTER :: OOC_INODE_SEQUENCE
        INTEGER(8),DIMENSION(:,:), POINTER :: OOC_SIZE_OF_BLOCK
        INTEGER*8, DIMENSION(:,:),   POINTER :: OOC_VADDR
        INTEGER,DIMENSION(:), POINTER :: OOC_TOTAL_NB_NODES
        INTEGER,DIMENSION(:), POINTER :: OOC_NB_FILES
        CHARACTER,DIMENSION(:,:), POINTER :: OOC_FILE_NAMES
        INTEGER,DIMENSION(:), POINTER :: OOC_FILE_NAME_LENGTH
!    Indices of nul pivots
        INTEGER,DIMENSION(:), POINTER :: PIVNUL_LIST
!    Internal control array
        DOUBLE PRECISION DKEEP(30)
!    Array needed to manage additionnal candidate processor
        INTEGER, DIMENSION(:,:), POINTER :: SUP_PROC
!   ------------------------
!   Root structure(internal)
!   ------------------------
        TYPE (ZMUMPS_ROOT_STRUC) :: root
      END TYPE ZMUMPS_STRUC
    TYPE(ZMUMPS_STRUC), SAVE :: mumps_par
    LOGICAL, PARAMETER :: mklpardiso=.FALSE.
    COMPLEX(kind(0.0d0)), ALLOCATABLE, TARGET, SAVE :: coak(:)
    INTEGER :: i, j, k, nnzi, idum(1), kzi
    COMPLEX(kind(0.0d0)) :: td, dum(1)
    REAL(kind(0.0d0)) :: absmax, tt
    COMPLEX(kind(0.0d0)), external :: ZDOTCM
    REAL(kind(0.0d0)), external :: DZNRM2
    INTEGER , parameter :: IONE=1
    MAXFCT=1
    MNUM=1
    MTYPE=13
    NRHS=1
    MSGLVL=0
    IF (n.EQ.1) THEN
        IF (ijob.EQ.1) THEN
           IF ( nlev.LE.2 ) THEN
              aaa=cmplx(1,0,kind(0.0d0))/a(1)
           ELSE
             absmax=MAXVAL(ABS(dt(nlev-1)%a(                                &
                         dt(nlev-1)%ia(1) : dt(nlev-1)%ia(nn(nlev-1)+1)-1 )))
             IF (ABS(a(1)).GT.100*epsmach*absmax) THEN
                aaa=cmplx(1,0,kind(0.0d0))/a(1)
             ELSE
                aaa=cmplx(0,0,kind(0.0d0))
             END IF
           END IF
        ELSE IF (ijob.EQ.2) THEN
           f(1)=aaa*f(1)
        END IF
        RETURN
    END IF
    IF (ijob.EQ.-2) THEN
       IF (mklpardiso) THEN
          PHASE=-1
          CALL PARDISO( PT, MAXFCT, MNUM, MTYPE, PHASE, N      &
               , dum, idum, idum                               &
               , idum, NRHS, IPARM, MSGLVL, dum, dum, ERROR )
          IF (ERROR.NE.0) THEN
             PRINT 999, ERROR, PHASE
             STOP
          END IF
          DEALLOCATE(PT,X)
       ELSE
          mumps_par%JOB = -2
          CALL ZAGMG_MUMPS(mumps_par)
          IF (transint) THEN
             DEALLOCATE(mumps_par%JCN_loc)
             IF (nlev.EQ.1) DEALLOCATE(mumps_par%IRN_loc,mumps_par%A_loc)
          ELSE
             DEALLOCATE(mumps_par%IRN_loc)
             IF (nlev.EQ.1) DEALLOCATE(mumps_par%JCN_loc,mumps_par%A_loc)
          END IF
       END IF
       IF (coasing) DEALLOCATE(coak)
    ELSE IF (ijob.EQ.1) THEN
       absmax=MAXVAL(ABS(a(ia(1):ia(n+1)-1)))
       IF (mklpardiso) THEN
          ALLOCATE(PT(64))
          do i = 1, 64
            PT(i)%DUMMY = 0
          end do
          IPARM(1:64)=0
          IPARM(1)=1
          IPARM(2)=2
          IPARM(6)=1
          IPARM(10)=13
          IPARM(11)=1
          IF (transint) THEN
             IPARM(12)=2
          ELSE
             IPARM(12)=0
          END IF
          IPARM(18)=-1
          IPARM(19)=-1
          IPARM(20)=0
          ERROR  = 0
          PHASE=12
          CALL zagmg_sortrows(n,a,ja,ia)
          CALL PARDISO( PT, MAXFCT, MNUM, MTYPE, PHASE, N      &
               , a, ia, ja                                     &
               , idum, NRHS, IPARM, MSGLVL, dum, dum, ERROR )
          IF (ERROR.NE.0) THEN
             PRINT 999, ERROR, PHASE
             STOP
          END IF
          flop=1e6*dble(IPARM(19))
          iflop=2*dble(IPARM(18))-n
          ALLOCATE(X(n))
       ELSE
          mumps_par%COMM = 0
          mumps_par%JOB = -1
          mumps_par%SYM = 0
          mumps_par%PAR = 1
          CALL ZAGMG_MUMPS(mumps_par)
          mumps_par%ICNTL(2)=-1
          mumps_par%ICNTL(3)=-1
          mumps_par%ICNTL(4)=0
          mumps_par%ICNTL(14)=80
          mumps_par%ICNTL(18)=3
          mumps_par%ICNTL(24)=1
          nnzi=ia(n+1)-ia(1)
          mumps_par%NZ_loc=nnzi
          IF (transint) THEN
             ALLOCATE( mumps_par%JCN_loc(nnzi) )
             do i=1,n
                do j=ia(i),ia(i+1)-1
                   mumps_par%JCN_loc(j)=i
                end do
             end do
             IF (nlev.GT.1) THEN
               mumps_par%IRN_loc => dt(nlev)%ja(1:nnzi)
               mumps_par%A_loc => dt(nlev)%a(1:nnzi)
             ELSE
               ALLOCATE( mumps_par%IRN_loc(nnzi) , mumps_par%A_loc(nnzi) )
               mumps_par%IRN_loc = ja(1:nnzi)
               mumps_par%A_loc = a(1:nnzi)
             END IF
          ELSE
             ALLOCATE( mumps_par%IRN_loc(nnzi) )
             do i=1,n
                do j=ia(i),ia(i+1)-1
                   mumps_par%IRN_loc(j)=i
                 end do
             end do
             IF (nlev.GT.1) THEN
               mumps_par%JCN_loc => dt(nlev)%ja(1:mumps_par%NZ_loc)
               mumps_par%A_loc => dt(nlev)%a(1:mumps_par%NZ_loc)
             ELSE
               ALLOCATE( mumps_par%JCN_loc(nnzi) , mumps_par%A_loc(nnzi) )
               mumps_par%JCN_loc = ja(1:nnzi)
               mumps_par%A_loc = a(1:nnzi)
             END IF
          END IF
          mumps_par%N=n
          mumps_par%JOB = 4
          CALL ZAGMG_MUMPS(mumps_par)
          flop=mumps_par%RINFO(3)
          i=mumps_par%INFO(27)
          IF (i.GT.0) THEN
             iflop=2*dble(i)-n
          ELSE
             iflop=-2.0d6*dble(i)-n
          END IF
       END IF
       ALLOCATE(coak(n))
       coak(1:n)=cmplx(1,0,kind(0.0d0))
       IF (mklpardiso) THEN
          PHASE=33
          CALL PARDISO( PT, MAXFCT, MNUM, MTYPE, PHASE, N             &
               , a, ia, ja                                            &
               , idum, NRHS, IPARM, MSGLVL, coak, X, ERROR )
          IF (ERROR.NE.0) THEN
             PRINT 999, ERROR, PHASE
             STOP
          END IF
       ELSE
          mumps_par%JOB = 3
          mumps_par%RHS => coak(1:n)
          CALL ZAGMG_MUMPS(mumps_par)
       END IF
       tt=(DZNRM2(N,coak,IONE))
       IF (  SQRT(dble(N))/(ABS(tt)*absmax).LE.100*epsmach  ) THEN
          coasing=.TRUE.
          coak(1:n)=(1.0d0/tt)*coak(1:n)
       ELSE
          coasing=.FALSE.
          DEALLOCATE(coak)
       END IF
    ELSE IF (ijob.EQ.2) THEN
       IF (mklpardiso) THEN
          PHASE=33
          IF (PRESENT(A)) THEN
             CALL PARDISO( PT, MAXFCT, MNUM, MTYPE, PHASE, N      &
                  , a, ia, ja                                     &
                  , idum, NRHS, IPARM, MSGLVL, f, X, ERROR )
          ELSE
             CALL PARDISO( PT, MAXFCT, MNUM, MTYPE, PHASE, N      &
                  , dt(nlev)%a, dt(nlev)%ia, dt(nlev)%ja          &
                  , idum, NRHS, IPARM, MSGLVL, f, X, ERROR )
          END IF
          IF (ERROR.NE.0) THEN
             PRINT 999, ERROR, PHASE
             STOP
          END IF
       ELSE
          mumps_par%JOB = 3
          mumps_par%RHS => f(1:n)
          CALL ZAGMG_MUMPS(mumps_par)
       END IF
       flop=flop+iflop
       IF (coasing) THEN
          td=ZDOTCM(N,coak,IONE,f,IONE)
          CALL ZAXPY(N,-td,coak,IONE,f,IONE)
       END IF
    END IF
    RETURN
999 FORMAT('Error detected in PARDISO; error code:', i2,' ; Phase:',i2)
  END SUBROUTINE zagmg_DIRseq
  SUBROUTINE zagmg_sortrows(n,a,ja,ia)
    IMPLICIT NONE
    INTEGER :: n, ia(n+1), ja(*), i, k, jj
    COMPLEX(kind(0.0d0)) :: a(*), aa
    LOGICAL exc
    DO i=1,n
       exc=.TRUE.
       DO WHILE(exc)
          exc=.FALSE.
          DO k=ia(i)+1,ia(i+1)-1
             IF(ja(k).LT.ja(k-1)) THEN
                jj=ja(k-1)
                aa=a(k-1)
                ja(k-1)=ja(k)
                a(k-1)=a(k)
                ja(k)=jj
                a(k)=aa
                exc=.TRUE.
             END IF
          END DO
       END DO
    END DO
    RETURN
  END SUBROUTINE zagmg_sortrows
END MODULE zagmg_ALLROUTINES
!!!!!!!!!!!!!!!!!!! MAIN DRIVER
  SUBROUTINE zagmg( n,a,ja,ia,f,x,ijob,iprint,nrest,iter,tol )
    USE zagmg_mem
    USE zagmg_ALLROUTINES
    IMPLICIT NONE
    INTEGER    :: n,ia(n+1),ijob,iprint,nrest,iter
    INTEGER  :: ja(*)
    COMPLEX (kind(0.0d0)) :: a(*),f(n),x(n)
    REAL (kind(0.0d0)) :: tol
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Arguments
!  =========
!
!  N       (input) INTEGER.
!          The dimension of the matrix.
!
!  A       (input/output) COMPLEX (kind(0.0d0)). Numerical values of the matrix.
!  IA      (input/output) INTEGER. Pointers for every row.
!  JA      (input/output) INTEGER. Column indices.
!
!          Significant only if IJOB==0,1,2,10,12,100,101,102,110,112
!
!          Detailed description of the matrix format
!
!              On input, IA(I), I=1,...,N, refers to the physical start
!              of row I. That is, the entries of row I are located
!              in A(K), where K=IA(I),...,IA(I+1)-1. JA(K) carries the
!              associated column indices. IA(N+1) must be defined in such
!              a way that the above rule also works for I=N (that is,
!              the last valid entry in arrays A,JA should correspond to
!              index K=IA(N+1)-1).
!
!          AGMG ASSUMES THAT ALL DIAGONAL ENTRIES ARE POSITIVE
!
!              That is, AGMG assumes that, for K=IA(I),...,IA(I+1)-1,
!              some of JA(K) is equal to I with corresponding A(K)>0
!              (value of the diagonal element, which must be positive).
!
!              A,IA,JA are "output" parameters because on exit the
!              entries of each row may occur in a different order
!              (The matrix is mathematically the same, but stored in
!              different way).
!
!
!  F       (input/output) COMPLEX (kind(0.0d0)).
!          On input, the right hand side vector f.
!          Overwritten on output.
!          Significant only if IJOB is not equal to 1 or 101.
!
!  X       (input/output) COMPLEX (kind(0.0d0)).
!          On input and if IJOB== 10, 12, 110, 112, 212: initial guess
!          On output, the computed solution
!              (IJOB==3: result of the application of the preconditioner)
!
! IJOB     (input) INTEGER. Tells AGMG what has to be done.
!          0: performs setup + solve + memory release, no initial guess
!         10: performs setup + solve + memory release, initial guess in x(1:n)
!          1: performs setup only
!             (preprocessing: prepares all parameters for subsequent solves)
!          2: solves only with the preconditioner based on previous setup,
!             using the system matrix provided in A, JA, IA; no initial guess
!         12: solves only with the preconditioner based on previous setup,
!             using the system matrix provided in A, JA, IA;
!             initial guess in x(1:n)
!        202: solves only with the preconditioner based on previous setup and
!             the same matrix as the one provided for set up;  no initial guess
!             using the system matrix provided during previous setup
!        212: solves only with the preconditioner based on previous setup and
!             the same matrix as the one provided for set up;
!             initial guess in x(1:n)
!          3: the vector returned in x(1:n) is not the solution of the linear
!                 system, but the result of the action of the multigrid
!                 preconditioner to the right hand side in f(1:n)
!         -1: erases the setup and releases internal memory
!
!   IJOB == 100,110,101,102,112: same as, respectively, IJOB==0,10,1,2,12
!       but, use the TRANSPOSE of the input matrix in A, JA, IA.
!
!   !!! IJOB==2,3,12,102,112,202,212 require that one has previously called
!       AGMG with IJOB==1 or IJOB==101
!
!   !!! (change with respect to versions 2.x) !!!
!       The preconditioner defined when calling AGMG
!         with IJOB==1 or IJOB==101 is entirely kept in internal memory.
!       Hence the arrays A, JA and IA are not accessed upon subsequent calls
!         with IJOB==3 and IJOB==202,212
!       Upon subsequent calls with IJOB==2,12,102,112: a matrix needs to
!            be supplied in arrays A, JA, IA. It will be used to
!            perform matrix vector product within the main iterative
!            loop (and only for this).
!            Hence the system is solved with this matrix which
!            may differ from the matrix in A, JA, IA that was supplied
!            upon the previous call with IJOB==1 or IJOB==101;
!            then AGMG attempts to solve a linear system with the "new"
!            matrix (supplied when IJOB==2,12,102 or 112) using the
!            preconditioner set up for the "old" one (supplied when
!            IJOB==1 or 101).
!         The same remarks apply to IJOB >= 100 or not: the value IJOB==1
!            or 101 determines whether the preconditioner set up and stored
!            in internal memory is based on the matrix or its transpose;
!            the value IJOB==2,12 or 102,112 is used to determine whether
!            the linear system to be solved is with the matrix or its
!            transpose, independently of the set up.
!            Hence one may set up a preconditioner for a matrix and use it
!            for its transpose.
!       These functionalities (set up a preconditioner and use it for another
!            matrix) are provided for the sake of generality but should be
!            used with care; in general, set up is fast with AGMG and hence
!            it is recommended to rerun it even if the matrix changes only
!            slightly.
!       In the more standard case where one needs successive solves with a
!       same matrix, it recommended to use IJOB==202 or 212 (in which case
!       the arrays A, JA, IA are not accessed) instead of IJOB==2 or 12
!       (in which case the input matrix has to be kept equal to the one
!       provided for set up). This  slightly faster and helps to save memory
!       since the arrays that store the input matrix may be deallocated after
!       set up.
!
! IPRINT   (input) INTEGER.
!              Indicates the unit number where information is to be printed
!              (N.B.: 5 is converted to 6). If nonpositive, only error
!              messages are printed on standard output. Warning messages about
!              insufficient convergence (in the allowed  maximum number of
!              iterations) are further suppressed supplying a negative number.
!
! NREST    (input) INTEGER.
!             Restart parameter for GCR (an implementation of GMRES)
!             Nonpositive values are converted to NREST=10 (default)
!
! !!  If NREST==1, Flexible CG is used instead of GCR (when IJOB==0,10,2,
!             12,100,110,102,112) and also (IJOB==0,1,100,101) performs some
!             simplification during the set up based on the assumption
!             that the matrix supplied in A, JA, IA is symmetric (there is
!             then no more difference between IJOB==1 and IJOB==101).
!
! !!!  NREST=1 Should be used if and only if the matrix is really SYMMETRIC
! !!!         (and positive definite).
!
!  ITER    (input/output) INTEGER.
!          On input, the maximum number of iterations. Should be positive.
!          On output, actual number of iterations.
!            If this number of iteration was insufficient to meet convergence
!            criterion, ITER will be returned negative and equal to the
!            opposite of the number of iterations performed.
!          Significant only if IJOB==0, 2, 10, 12, 100, 102, 110, 112, 202, 212
!
!  TOL     (input) REAL (kind(0.0d0))
!          The tolerance on residual norm. Iterations are stopped whenever
!               || A*x-f || <= TOL* || f ||
!          Should be positive and less than 1.0
!          Significant only if IJOB==0, 2, 10, 12, 100, 102, 110, 112, 202, 212
!
!
!!!!! Remark !!!! Except insufficient number of iterations to achieve
!                 convergence (characterized by a negative value returned
!                 in ITER), all other detected errors are fatal and lead
!                 to a STOP statement.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER, SAVE :: nza
    INTEGER :: i,j,ijb,k,llev,nzext,itef,ierr,errc
    REAL(kind(0.0d0)) :: cputmp
    REAL(kind(0.0d0)) :: resid,tmp
    COMPLEX(kind(0.0d0)), POINTER, DIMENSION(:) :: ascaled
    COMPLEX(kind(0.0d0)) :: adum(1)
    CHARACTER*40 :: valenv
    wfo=.TRUE.
    iout=iprint
    IF (iprint.LE.0) THEN
       iout=6
       wfo=.FALSE.
    ELSE IF (iprint.EQ.5) THEN
       iout=6
    END IF
    ijb=MOD(ijob,100)
    itef=iter
    IF (ijob.NE.-1 .AND. &
        ijob.NE. 0 .AND. &
        ijob.NE. 1 .AND. &
        ijob.NE. 2 .AND. &
        ijob.NE. 3 .AND. &
        ijob.NE. 4 .AND. &
        ijob.NE. 5 .AND. &
        ijob.NE. 6 .AND. &
        ijob.NE. 7 .AND. &
        ijob.NE. 8 .AND. &
        ijob.NE. 9 .AND. &
        ijob.NE.10 .AND. &
        ijob.NE.12 .AND. &
        ijob.NE.100 .AND. &
        ijob.NE.101 .AND. &
        ijob.NE.102 .AND. &
        ijob.NE.110 .AND. &
        ijob.NE.112 .AND. &
        ijob.NE.202 .AND. &
        ijob.NE.212  &
        ) THEN
       WRITE (6,1000) IRANK,ijob
       STOP
    END IF
1000 FORMAT(i3,'*',' FATAL ERROR: ijob=',i5, &
         ' is not allowed')
    IF (MOD(ijob,10).GE.2 .AND. .NOT.preprocessed) THEN
       WRITE (6,1001) IRANK,ijob
       STOP
    END IF
1001 FORMAT(i3,'*',' FATAL ERROR: setup not done: ijob=',i3, &
         ' is not allowed')
    IF (iprint.EQ.-978 .AND. ijob.GE.0) THEN
       valenv(1:10)='AGMG_dump_'
       valenv(19:19)='_'
       CALL DATE_AND_TIME(valenv(11:18),valenv(20:29))
       valenv(28:31)='.mat'
       valenv(33:38)=' '
       OPEN(UNIT=91825,FILE=valenv,FORM='FORMATTED')
       IF (ijob.LE.200) THEN
          WRITE(91825,'(2i18,1x,A4,i1,2i5,i8,ES18.10)')                       &
               n,ia(n+1)-1,'Dump',0,ijob,nrest,iter,tol
         DO i=1,n
            DO j=ia(i),ia(i+1)-1
               WRITE(91825,'(2i18,ES18.10)') i,ja(j),a(j)
            END DO
         END DO
       ELSE
          WRITE(91825,*) n,0,'DUMP',ijob,nrest,iter,tol
       END IF
       CLOSE (91825)
       IF (MOD(ijob,10).NE.1) THEN
          valenv(29:31)='rhs'
          OPEN (UNIT=91825,FILE=valenv,FORM='FORMATTED')
          DO i=1,n
             WRITE(91825,'(ES18.10)') f(i)
          END DO
          CLOSE(91825)
          valenv(29:31)='ini'
          OPEN (UNIT=91825,FILE=valenv,FORM='FORMATTED')
          DO i=1,n
             WRITE(91825,'(ES18.10)') x(i)
          END DO
          CLOSE(91825)
       END IF
       valenv(29:31)='out'
       OPEN (UNIT=91825,FILE=valenv,FORM='FORMATTED')
       iout=91825
       WRITE(iout,910) IRANK
       wfo=.TRUE.
    END IF
    IF (ijob.LT.0 .AND. .NOT.preprocessed) RETURN
    IF (iprint.GT.0 .AND. ijob.GE.0 .AND. MOD(ijob,10).LE.1) THEN
       WRITE(iout,910) IRANK
    END IF
 910 FORMAT(i3,'*ENTERING AGMG **********************************',&
         '***************************')
    wff=wfo.AND.(IRANK.LE.0)
    woo=IRANK.LE.0
    IF (iprint.LT.0 .AND. iprint.NE.-978) woo=.FALSE.
    IF (ijb.LT.0) GOTO 500
    IF (ijob.LE.201) trans=ijob.GE.100
    IF (MOD(ijb,10).GE.2) THEN
    ELSE
      IF (preprocessed) THEN
       CALL zagmg_relmem(-1)
       solve=.FALSE.
       eltm=0.0d0
       cputm=0.0d0
       flop=0.0d0
       kstat=0
      END IF
      CALL zagmg_mestime(-1,0.0d0,0.0d0)
      spd=nrest.EQ.1
      spdcoarse=spd
      transint=trans
      nlev=0
      imult=1
      ALLOCATE(dt(1)%idiag(n+1))
      CALL zagmg_partroword(n,a,ja,ia,dt(1)%idiag)
      nza=ia(n+1)-ia(1)
      ALLOCATE(scald(n),ascaled(ia(n+1)-ia(1)))
      CALL zagmg_diagscal(n,a,ja,ia,dt(1)%idiag,ascaled)
      CALL zagmg_setupL1(n,ascaled,ja,ia)
      DEALLOCATE(ascaled)
      CALL zagmg_smoothsetup
      preprocessed=.TRUE.
      CALL zagmg_mestime(1,cputmp,eltmp)
      IF(wfo)THEN
   !CPU_TIME: next line may be uncommented if implemented
   !          (see subroutine mestime)
   !              WRITE(iout,996) IRANK,cputmp
       WRITE(iout,997) IRANK,eltmp
       WRITE(iout,'()')
      END IF
      IF (MOD(ijb,10).EQ.1) THEN
        IF (nlev.EQ.1) THEN
          ALLOCATE(dt(1)%a(ia(n+1)-ia(1)),dt(1)%ja(ia(n+1)-ia(1)),dt(1)%ia(n+1))
          dt(1)%a=a(ia(1):ia(n+1)-1)
          dt(1)%ja=ja(ia(1):ia(n+1)-1)
          dt(1)%ia=ia(1:n+1)
        END IF
        IF (iprint.EQ.-978) CLOSE(91825)
        RETURN
      END IF
    END IF
    CALL zagmg_mestime(-2,0.0d0,0.0d0)
    resid=tol
    IF (MOD(ijb,10).GE.3) THEN
       IF (wfo) THEN
          WRITE(iout,901) IRANK
       END IF
       CALL zagmg_applyprec(n,f,x,MOD(ijb,10),a,ja,ia)
       CALL zagmg_mestime(2,cputmp,eltmp)
       cputm=cputm+cputmp
       eltm=eltm+eltmp
       solve=.TRUE.
       IF (iprint.EQ.-978) CLOSE(91825)
       RETURN
    ELSE IF (nrest.NE.1) THEN
       nrst=nrest
       IF (nrst.LE.0) nrst=10
       CALL zagmg_GCR(n,f,x,itef,resid,a,ja,ia,ijob)
    ELSE
       CALL zagmg_FlexCG(n,f,x,itef,resid,a,ja,ia,ijob)
    END IF
    iter=itef
    solve=.FALSE.
    print *,cputm,eltm
    CALL zagmg_mestime(2,cputm,eltm)
    IF (wfo) THEN
       IF (wff .AND. iter.NE.0) THEN
          llev=nlev-1
          DO i=2,llev
             WRITE(iout,955) min(i,nlev),kstat(2,i-1),kstat(2,i),       &
                  dble(kstat(2,i))/dble(kstat(2,i-1)),kstat(1,i)
          END DO
       END IF
       WRITE(iout,'()')
       IF (RESID.GT.0 .AND. iter.NE.0) THEN
         tmp=(flop*log(1.0d0/10)/log(RESID))
         WRITE(iout,952) IRANK,tmp/dble(2*nza)
       END IF
   !CPU_TIME: next line may be uncommented if implemented
   !          (see subroutine mestime)
   !            WRITE(iout,998) IRANK,cputm
       WRITE(iout,999) IRANK,eltm
       WRITE(iout,'()')
       IF (RESID.GT.0 .AND. iter.NE.0) THEN
          WRITE (iout,902) IRANK
       END IF
    END IF
    IF (MOD(ijb,10).GT.0) THEN
       solve=.FALSE.
       eltm=0.0d0
       cputm=0.0d0
       flop=0.0d0
       kstat=0
       IF (iprint.EQ.-978) CLOSE(91825)
       RETURN
    END IF
500 CONTINUE
    IF (solve .AND. wfo) THEN
       WRITE(iout,'()')
       WRITE(iout,990) IRANK
       WRITE(iout,'()')
       IF (wff) THEN
          llev=nlev-1
          DO i=2,llev
             WRITE(iout,955) min(i,nlev),kstat(2,i-1),kstat(2,i),       &
                  dble(kstat(2,i))/dble(kstat(2,i-1)),kstat(1,i)
          END DO
       END IF
       WRITE(iout,'()')
       WRITE(iout,953) IRANK,flop/dble(2*nza)
       !CPU_TIME: next line may be uncommented if implemented
       !          (see subroutine mestime)
       !            WRITE(iout,998) IRANK,cputm
       WRITE(iout,999) IRANK,eltm
       WRITE(iout,'()')
    END IF
    CALL zagmg_relmem(ijob)
    preprocessed=.FALSE.
    solve=.FALSE.
    eltm=0.0d0
    cputm=0.0d0
    flop=0.0d0
    kstat=0
    IF (wfo) THEN
       WRITE (iout,903) IRANK
    END IF
    IF (iprint.EQ.-978) CLOSE(91825)
    RETURN
901 FORMAT(i3,'*ONE APPLICATION OF AGMG PRECONDITIONER')
902 FORMAT(i3,            &
  ' (*) 1 work unit represents the cost of 1 (fine grid) residual evaluation ')
903 FORMAT(i3,'*LEAVING AGMG * (MEMORY RELEASED) ***************',&
         '***************************')
952 FORMAT(i3,'*','       Number of work units:',f9.2,             &
              ' per digit of accuracy (*)')
953 FORMAT(i3,'*',' Total number of work units:',f9.2, '(*)')
955 FORMAT('****     level',i2,'   #call=',i6,'   #cycle=',i6,    &
         '   mean=',f7.2,'    max=',i3)
990 FORMAT(i3,'*GLOBAL STATISTICS for preconditioner application:')
996 FORMAT(i3,'*','           Setup time (CPU):   ',ES10.2,     &
         ' seconds')
997 FORMAT(i3,'*','       Setup time (Elapsed):   ',ES10.2,     &
         ' seconds')
998 FORMAT(i3,'*','        Solution time (CPU):   ',ES10.2,     &
         ' seconds')
999 FORMAT(i3,'*','    Solution time (Elapsed):   ',ES10.2,     &
         ' seconds')
  END SUBROUTINE zagmg
!!!!!!!!!!!!!!!!!!! END of MAIN DRIVER
        COMPLEX(kind(0.0d0)) FUNCTION ZDOTCM(N,ZX,INCX,ZY,INCY)
          INTEGER :: INCX,INCY,N
          COMPLEX(kind(0.0d0)) :: ZX(*),ZY(*)
          COMPLEX(kind(0.0d0)) :: ZTEMP
          INTEGER I,IX,IY
          ZTEMP = cmplx(0,0,kind(0.0d0))
          ZDOTCM = cmplx(0,0,kind(0.0d0))
          IF (N.LE.0) RETURN
          IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
          IX = 1
          IY = 1
          IF (INCX.LT.0) IX = (-N+1)*INCX + 1
          IF (INCY.LT.0) IY = (-N+1)*INCY + 1
          DO I = 1,N
             ZTEMP = ZTEMP + conjg(ZX(IX))*ZY(IY)
             IX = IX + INCX
             IY = IY + INCY
          END DO
          ZDOTCM = ZTEMP
          RETURN
20        CONTINUE
          DO I = 1,N
             ZTEMP = ZTEMP + conjg(ZX(I))*ZY(I)
          END DO
          ZDOTCM = ZTEMP
          RETURN
        END FUNCTION ZDOTCM
