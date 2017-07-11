
import os

import numpy as np
import shutil
from setrun import setgeo,setrun
from setplot import setplot
from clawpack.clawutil import runclaw,data


tauy_i = [3500., 7000., 10000., 15000.]
tauy_r = [2000., 3500.]
remold_coeff = [5e-3, 5e-4, 5e-2]

remolding = True
tauy_i = [7000., 10000., 12000., 15000.]
tauy_r = [1000.,2000.,3000.,3500.]
remold_coeff = [5e-2,5e-3,5e-4,5e-5]

remolding = True
tauy_i = [15000.,12000.,10000.,7000.]
tauy_r = [1000.,2000.,3000.,3500.]
remold_coeff = [5e-2,5e-3,5e-4,5e-5]

remolding = True
tauy_i = [11000.,11500.,12000.,12500.,13000.]
tauy_r = [3000.]
remold_coeff = [7e-3,1e-4,3e-4,5e-4,7e-4]

remolding = True
tauy_i = [11000.,11500.,12000.,12500.,13000.]
tauy_r = [3000.]
remold_coeff = [7e-3,1e-4,3e-4,5e-4,7e-4]

remolding = True
tauy_i = [11000.,11500.,12000.,12500.]
tauy_r = [2500.]
remold_coeff = [7e-3,3e-4,5e-4,7e-4]

for tyi in tauy_i:
    for tyr in tauy_r:
        for rc in remold_coeff:
                rundata = setrun('geoclaw')
                rundata.dtopo_data.dtopofiles = []
                dtopo_path = '/media/jihwan/jik/work/ngi/Numerical_models/VP2HD_Case_Studies/Storegga/topo/test19/dtopo_%dto%d_gamma%s.tt3' %(tyi,tyr,rc)
                rundata.dtopo_data.dtopofiles.append([3,3,3,dtopo_path])

                rundata.write()

                outdir = '/media/jihwan/jik/work/ngi/Numerical_models/VP2HD_Case_Studies/Storegga/results/test19/_output_%d%d%s' %(tyi,tyr,rc) 
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                runclaw.runclaw('xgeoclaw', outdir, True, False,'.', False, False, None)
"""

tauy_i = [100., 500.]
cF_hyd = [0.01, 0.025] 

rundata = setrun('geoclaw')

for tyi in tauy_i:
    for cfh in cF_hyd:
        rundata = setrun('geoclaw')
        rundata.dtopo_data.dtopofiles = []
        dtopo_path = '/media/jihwan/jik/work/ngi/Bing_FVM/Bing_Case_Studies/Storegga/topo/test3/dtopo_%d_cfh%s.tt3' %(tyi,cfh)
        rundata.dtopo_data.dtopofiles.append([3,3,3,dtopo_path])

        rundata.write()
        outdir = '/media/jihwan/jik/work/ngi/Bing_FVM/Bing_Case_Studies/Storegga/results/test3/_output_%d_%s' %(tyi,cfh) 
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        runclaw.runclaw('xgeoclaw', outdir, True, False,'.', False, False, None)
               
"""

