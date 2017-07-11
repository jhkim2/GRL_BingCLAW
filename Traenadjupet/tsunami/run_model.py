
import os

import numpy as np
import shutil
from setrun import setgeo,setrun
from setplot import setplot
from clawpack.clawutil import runclaw,data

remolding = True
tauy_i = [12000]
tauy_r = [2500.]
remold_coeff = [3e-4]
cF_hyd = [0.01]
cP_hyd = [1.0]

for tyi in tauy_i:
    for tyr in tauy_r:
        for rc in remold_coeff:

                rundata = setrun('geoclaw')
                rundata.dtopo_data.dtopofiles = []
                dtopo_path = '/media/jihwan/jik/work/ngi/Numerical_models/VP2HD_Case_Studies/Traenadjupet/topo/vol900/c2/dtopo_%dto%d_gamma%s.tt3' %(tyi,tyr,rc)
                rundata.dtopo_data.dtopofiles.append([3,3,3,dtopo_path])

                rundata.write()

                outdir = '/media/jihwan/jik/work/ngi/Numerical_models/VP2HD_Case_Studies/Traenadjupet/results/vol900/c2/_output_%d%d%s' %(tyi,tyr,rc) 
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                runclaw.runclaw('xgeoclaw', outdir, True, False,'.', False, False, None)

