subroutine setprob()

    use regions_module, only: set_regions
    use gauges_module, only: set_gauges
    use geoclaw_module, only: set_geo
    use topo_module, only: read_topo_settings, read_dtopo_settings
    use qinit_module, only: set_qinit, qinit_style
    use fixedgrids_module, only: set_fixed_grids
    use refinement_module, only: set_refinement
    use storm_module, only: set_storm
    use friction_module, only: setup_variable_friction
    use bing_module

    implicit none

    integer :: iunit
    character(len=25) fname

    iunit = 7
    fname = 'setprob.data'
!   # open the unit with new routine from Clawpack 4.4 to skip over
!   # comment lines starting with #:
    call opendatafile(iunit, fname)

    read(7,*) rho_a
    read(7,*) rho_s
    read(7,*) n_param
    read(7,*) gamma_r
    read(7,*) c_mass
    read(7,*) hydrodrag
    read(7,*) cF_hyd 
    read(7,*) cP_hyd    
    read(7,*) remolding 
    read(7,*) tauy_i 
    read(7,*) tauy_r 
    read(7,*) remold_coeff 
    read(7,*) qinit_style

    visc_kin = visc_dyn/rho_s
    alpha_1 = (1.d0/n_param+1.d0)/(1.d0/n_param+2.d0) 
    alpha_2 = 1.d0 - 2.d0/(1.d0/n_param + 2.d0) + 1.d0/(2.d0/n_param+3.d0) 
    beta = (1.d0+1.d0/n_param)**n_param 
    cm_coeff = 1.d0+c_mass*rho_a/rho_s 

    call set_geo()                    !# sets basic parameters g and coord system
    call set_refinement()             !# sets refinement control parameters
    call read_dtopo_settings()        !# specifies file with dtopo from earthquake
    call read_topo_settings()         !# specifies topography (bathymetry) files
    call set_qinit()                  !# specifies file with dh if this used instead
    call set_fixed_grids()            !# Fixed grid settings
    call set_storm()                  ! Set storm parameters
    call setup_variable_friction()    ! Set variable friction parameters

end subroutine setprob
