module global_parameter
    
    implicit none
    
    logical,parameter::output_div=.false.
    logical,parameter::output_ccv=.true.
    logical,parameter::output_tipper=.false.
   logical,parameter::div_cor=.true. ! open as default
   integer,parameter ::corr_intv = 200
   real*8,parameter ::corr_tol = 1.d-3  ! dose not have to be strict    
     real(8), parameter:: pi=3.141592653589793238462643d0
      real*8,parameter ::sigma_air =1.d-9, sigma_hhs=1.d-2
    !real(8), parameter:: f = 1.d-3 !! resistivity
    real(8), parameter:: mu = 4*pi*10.d0**(-7) !!permeability of vacuum
    !!gridsize growing factor in vertical and horizontal directions(ajusted by users)
    real(8), parameter:: lambda_v = 1.15d0,lambda_h = 1.08d0 
    integer,parameter :: bgf_type = 1 ! type of calculating the background field: 1--analytical(recommend); 2--FEM 
    character,parameter ::pre_type='S' ! I--ILU(0); J--Jacobi; S--SSOR, S is recommended 
    integer,parameter :: fwd_acc = 0  ! using EXCMG or not(0--not using)


end module
