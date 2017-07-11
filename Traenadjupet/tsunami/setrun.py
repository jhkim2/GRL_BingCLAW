"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

import os
import numpy as np

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")

bouss = False                # Dispersive terms ?
B_param = 1.0 / 15.0         # Parameter for the Boussinesq eqns
use_bous_sw_thresh = False    # Use the switching threshold
bous_sw_thresh = 0.8         # Threshold for the transition to SWE

#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)


    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------
    
    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    probdata.add_param('bouss'             , bouss       ,'Include dispersive terms?')
    probdata.add_param('B_param'           , B_param     ,'Parameter for the Boussinesq eq')
    probdata.add_param('use_bous_sw_thresh', use_bous_sw_thresh,'Use the switching threshold')
    probdata.add_param('bous_sw_thresh'    , bous_sw_thresh      ,'Eta to depth ratio')

    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------
    rundata = setgeo(rundata)

    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    #------------------------------------------------------------------
    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # Set single grid parameters first.
    # See below for AMR parameters.


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:
    clawdata.lower[0] = 0. 
    clawdata.upper[0] = 25. 

    clawdata.lower[1] = 60. 
    clawdata.upper[1] = 72. 


    # Number of grid cells: Coarsest grid
    clawdata.num_cells[0] = 25*30
    clawdata.num_cells[1] = 12*30

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 4

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2

    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0


    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False              # True to restart from prior results
    clawdata.restart_file = '../_checkpoint/fort.chk00251'  # File to use for restart data

    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 1

    if clawdata.output_style==1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.num_output_times = 60
        clawdata.tfinal = 60.*60.*10
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        clawdata.output_times = [0.5, 1.0]

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 3
        clawdata.output_t0 = True
        

    clawdata.output_format = 'ascii'      # 'ascii' or 'netcdf' 

    clawdata.output_q_components = 'all'   # need all
    clawdata.output_aux_components = 'all'  # eta=h+B is in q
    clawdata.output_aux_onlyonce = False    # output aux arrays each frame



    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 1

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.016

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.8
    clawdata.cfl_max = 0.9

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 50000

    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'
    
    # For unsplit method, transverse_waves can be 
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2

    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    clawdata.limiter = ['mc', 'mc', 'mc']

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms
    
    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 'godunov'


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'

    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'



    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 2

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif clawdata.checkpt_style == 1:
        # Checkpoint only at tfinal.
        pass

    elif clawdata.checkpt_style == 2:
        # Specify a list of checkpoint times.  
        clawdata.checkpt_times = [10.]

    elif clawdata.checkpt_style == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5


    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 1

    # List of refinement ratios at each level (length at least mxnest-1)
    amrdata.refinement_ratios_x = [5,2,2]
    amrdata.refinement_ratios_y = [5,2,2]
    amrdata.refinement_ratios_t = [5,2,2]


    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    amrdata.aux_type = ['center','capacity','yleft','center','center']


    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag2refine = True

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.700000

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0  

    #  ----- For developers ----- 
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = True       # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting
    
    # More AMR parameters can be set -- see the defaults in pyclaw/data.py

    # ---------------
    # Regions:
    # ---------------
    rundata.regiondata.regions = []
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    rundata.regiondata.regions.append([1, 3, 0., 1.e10, -180.,180.,-90,90.])
    rundata.regiondata.regions.append([3, 3, 0., 3600.*4., 6, 11, 67, 70])
    #rundata.regiondata.regions.append([2, 4, 0., 1.e10, -7.,-5.,61,63])
    #rundata.regiondata.regions.append([2, 4, 0., 1.e10, -3.,0.,60,62])
    #rundata.regiondata.regions.append([2, 4, 0., 1.e10, -5.,-3.,58,60])
    #rundata.regiondata.regions.append([2, 4, 0., 1.e10, 4.,8.,60,63])
    #rundata.regiondata.regions.append([2, 4, 0., 1.e10, 8.,11.,63,65])

    # ---------------
    # Gauges:
    # ---------------
    rundata.gaugedata.gauges = []

    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    rundata.gaugedata.gauges.append([14, 9.78 , 63.91, 0., 1.e10])  #Bjurn
    rundata.gaugedata.gauges.append([15, 10.76, 63.55, 0., 1.e10])  #Frost
    rundata.gaugedata.gauges.append([17,  5.85, 62.58, 0., 1.e10])  #More
    rundata.gaugedata.gauges.append([18,  6.08, 62.48, 0., 1.e10])  #Sula
    rundata.gaugedata.gauges.append([19, 5.563, 62.41, 0., 1.e10])  #Bergsoy
    rundata.gaugedata.gauges.append([20, 6.596, 62.15, 0., 1.e10])  #Sykkylven
    rundata.gaugedata.gauges.append([21, 8.53, 62.7, 0., 1.e10])   #Sunndalsora
    rundata.gaugedata.gauges.append([22, 7.996,62.76, 0., 1.e10])  #Eresfjord
    rundata.gaugedata.gauges.append([23, 4.78, 60.81, 0., 1.e10])  #Austreim1
    rundata.gaugedata.gauges.append([24, 4.7, 60.67, 0., 1.e10])   #Austreim2

    rundata.gaugedata.gauges.append([26, 10.6 ,   64.99, 0., 1.e10])   #Vikna
    rundata.gaugedata.gauges.append([27, 11.9 ,   65.3 , 0., 1.e10])   #Hommelsto
    rundata.gaugedata.gauges.append([28, 12.65,   66.46, 0., 1.e10])   #Luroy
    rundata.gaugedata.gauges.append([29, 11.97,   66.57, 0., 1.e10])   #Traena
    rundata.gaugedata.gauges.append([30, 14.09,   67.3 , 0., 1.e10])   #Bodo
    rundata.gaugedata.gauges.append([31, 11.76,   67.43, 0., 1.e10])   #Rost
    rundata.gaugedata.gauges.append([32, 12.93,   68.1 , 0., 1.e10])   #Moskenesoy
    rundata.gaugedata.gauges.append([33, 14.08,   68.13, 0., 1.e10])   #Vestvaagoy
    rundata.gaugedata.gauges.append([34, 14.73,   67.91, 0., 1.e10])   #Steigen
    rundata.gaugedata.gauges.append([35, 15.62,   69.24, 0., 1.e10])   #Andoya
    rundata.gaugedata.gauges.append([36, 17.29,   69.67, 0., 1.e10])    #Senja
    rundata.gaugedata.gauges.append([37, 23.82,   71., 0., 1.e10])  #Rolvsoya

    """
    rundata.gaugedata.gauges.append([1, -6.55,	61.47, 0., 1.e10]) #FaroeSouth
    rundata.gaugedata.gauges.append([2, -1.254,	60.56, 0., 1.e10]) #ShetlandNorthFjord
    rundata.gaugedata.gauges.append([3, -1.254, 60.73, 0., 1.e10]) #ShetlandNorth
    rundata.gaugedata.gauges.append([4, -0.6539, 60.8, 0., 1.e10]) #ShetlandNorthEast
    rundata.gaugedata.gauges.append([5, -0.7706, 60.28, 0., 1.e10]) #ShetlandEast
    rundata.gaugedata.gauges.append([6, -3.887, 58.76, 0., 1.e10]) #ScotlandNorth
    rundata.gaugedata.gauges.append([7, -2.821,58.46, 0., 1.e10]) #ScotlandNortheast
    rundata.gaugedata.gauges.append([8, -3.821,57.95, 0., 1.e10]) #ScotlandInverness1
    rundata.gaugedata.gauges.append([9, -3.654,57.68, 0., 1.e10]) #ScotlandInverness2
    rundata.gaugedata.gauges.append([10, -1.967,57.95, 0., 1.e10]) #ScotlandEast1
    rundata.gaugedata.gauges.append([11, -2.287,56.7, 0., 1.e10]) #ScotlandEast2
    rundata.gaugedata.gauges.append([12, -2.437,56.15, 0., 1.e10]) #ScotlandEast3
    rundata.gaugedata.gauges.append([13, -1.721,55.91, 0., 1.e10]) #ScotlandEast4
    rundata.gaugedata.gauges.append([14, 9.78 , 63.91, 0., 1.e10])  #Bjurn
    rundata.gaugedata.gauges.append([15, 10.76, 63.55, 0., 1.e10])  #Frost
    rundata.gaugedata.gauges.append([16, 10.98, 63.91, 0., 1.e10])  #Verrabotn
    rundata.gaugedata.gauges.append([17,  5.85, 62.58, 0., 1.e10])  #More
    rundata.gaugedata.gauges.append([18,  6.08, 62.48, 0., 1.e10])  #Sula
    rundata.gaugedata.gauges.append([19, 5.563, 62.41, 0., 1.e10])  #Bergsoy
    rundata.gaugedata.gauges.append([20, 6.596, 62.15, 0., 1.e10])  #Sykkylven
    rundata.gaugedata.gauges.append([21, 8.53, 62.7, 0., 1.e10])   #Sunndalsora
    rundata.gaugedata.gauges.append([22, 7.996,62.76, 0., 1.e10])  #Eresfjord
    rundata.gaugedata.gauges.append([23, 4.78, 60.81, 0., 1.e10])  #Austreim1
    rundata.gaugedata.gauges.append([24, 4.7, 60.67, 0., 1.e10])   #Austreim2
    #rundata.gaugedata.gauges.append([25, 5.03, 59.81, 0., 1.e10])  #Bomlo
    """

    """
-6.55	61.47	FaroesSouth
-1.254	60.56	ShetlandNorthFjord
-1.254	60.73	ShetlandNorth
-0.6539	60.8	ShetlandNortheast
-0.7706	60.28	ShetlandEast
-3.887	58.76	ScotlandNorth
-2.821	58.46	ScotlandNortheast
-3.821	57.95	ScotlandInverness1
-3.654	57.68	ScotlandInverness2
-1.967	57.95	ScotlandEast1
-2.287	56.7	ScotlandEast2
-2.437	56.15	ScotlandEast3
-1.721	55.91	ScotlandEast4
9.78	63.91	Bjugn
10.76	63.55	Frosta
10.98	63.91	Verrabotn
5.85	62.58	More
6.08	62.48	Sula
5.563	62.41	Bergsoy
6.596	62.15	Sykkylven
8.53	62.7	Sunndalsora
7.996	62.76	Eresfjord
4.78	60.81	Austreim1
4.7	60.67	Austreim2
5.03	59.81	Bomlo

9.78	63.91	Bjugn
5.85	62.58	More
6.08	62.48	Sula
5.563	62.41	Bergsoy
10.6    64.99   Vikna
11.9    65.3    Hommelsto
12.65   66.46   Luroy
11.97   66.57   Traena
14.09   67.3    Bodo
11.76   67.43   Rost
12.93   68.1    Moskenesoy
14.08   68.13   Vestvaagoy
14.73   67.91   Steigen
15.62   69.24   Andoya
17.29   69.67   Senja
    """

    return rundata
    # end of function setrun
    # ----------------------


#-------------------
def setgeo(rundata):
#-------------------
    """
    Set GeoClaw specific runtime parameters.
    For documentation see ....
    """

    try:
        geo_data = rundata.geo_data
    except:
        print "*** Error, this rundata has no geo_data attribute"
        raise AttributeError("Missing geo_data attribute")
       
    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2
    geo_data.earth_radius = 6367.5e3

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0.
    geo_data.dry_tolerance = 1.e-3
    geo_data.friction_forcing = True
    geo_data.manning_coefficient =.025
    geo_data.friction_depth = 1e2

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 1.e-1
    refinement_data.deep_depth = 1e2
    refinement_data.max_level_deep = 3

    # == settopo.data values ==
    topo_data = rundata.topo_data
    # for topography, append lines of the form
    #    [topotype, minlevel, maxlevel, t1, t2, fname]
    # == settopo.data values ==
    rundata.topo_data.topofiles = []
    topofiles = rundata.topo_data.topofiles
    # for topography, append lines of the form
    #    [topotype, minlevel, maxlevel, t1, t2, fname]
    #fname = '/media/jihwan/EXT_NGI/Topo/Storegga/storegga_etopo1.asc'
    #fname = '/media/jihwan/EXT_NGI/Topo/Storegga/paleoNorthAtlantic.tt3'
    fname = '../landslides/vol900/c1/slide_eta.tt3'
    topo_data.topofiles.append([3, 1, 3, 0., 1.e10, fname])
    fname = '/media/jihwan/jik/Topo/Storegga/PaleoNorthAtlantic.tt3'
    topofiles.append([-3, 1, 1, 0., 1.e10, fname])
    fname = '/media/jihwan/jik/Topo/etopo1/W30E30S50N75.asc'
    topofiles.append([3, 1, 1, 0., 1.e10, fname])

    # == setdtopo.data values ==
    dtopo_data = rundata.dtopo_data
    # for moving topography, append lines of the form :   (<= 1 allowed for now!)
    #   [topotype, minlevel,maxlevel,fname]
    dtopo_path = '/media/jihwan/jik/work/ngi/Numerical_models/VP2HD_Case_Studies/Traenadjupet/topo/vol900/c1/dtopo_12000to2500_gamma0.005.tt3'
    dtopo_data.dtopofiles.append([3,3,3,dtopo_path])
    #dtopo_data.dt_max_dtopo = 0.25


    # == setqinit.data values ==
    # == setqinit.data values ==
    rundata.qinit_data.qinit_type = 0
    rundata.qinit_data.qinitfiles = []
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]
    #rundata.qinit_data.qinitfiles.append([1, 3, fname])

    # == setfixedgrids.data values ==
    fgmax_files = rundata.fgmax_data.fgmax_files
    # for fixed grids append to this list names of any fgmax input files
    #fgmax_files.append('fgmax_grid.txt')
    #rundata.fgmax_data.num_fgmax_val = 2 # Save depth only

    return rundata
    # end of function setgeo
    # ----------------------



if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()

