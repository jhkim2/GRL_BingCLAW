! Module containing parameters for Bingham model
module bing_module
	
	real(kind=8) :: rho_a            ! Density of ambient fluid
	real(kind=8) :: rho_s            ! Density of slide

        ! Herschel-Bulkley model parameters
	real(kind=8) n_param             ! Bing rheology parameter
	real(kind=8) alpha_1             ! Shape factor of integration
	real(kind=8) :: alpha_2          ! Shape factor of integration
	real(kind=8) :: beta             ! Shape factor of integration
	real(kind=8) :: visc_dyn         ! Dynamic viscosity
        real(kind=8) :: visc_kin         ! Kinematic viscosity
	real(kind=8) :: gamma_r  

        ! Added Mass parameters
	real(kind=8) :: c_mass           ! Added Mass
	real(kind=8) :: cm_coeff         ! Added Mass Constant

        ! Hydrodrag parameters
        logical      :: hydrodrag
	real(kind=8) :: cF_hyd 
	real(kind=8) :: cP_hyd 

        ! Remolding parameters
        logical      :: remolding 
        real(kind=8) :: tauy_i            ! initial yield strength
        real(kind=8) :: tauy_r            ! residual yield strength
        real(kind=8) :: remold_coeff      ! remolding coefficient

        ! Viscous friction 
        logical, parameter :: visc_fric = .False.
        real(kind=8), parameter :: a_vis= 0.d0 ! Artificial Viscosity

end
