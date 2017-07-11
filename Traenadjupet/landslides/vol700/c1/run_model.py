
import os

import numpy as np
import shutil
from setrun import setgeo,setrun
from setplot import setplot
from clawpack.clawutil import runclaw,data
from makedtopo import *

remolding = True
n_p = [0.5]
tauy_i = [12000]
tauy_r = [3000.]
remold_coeff = [5e-3]
cF_hyd = [0.01]
cP_hyd = [1.0]
c_mass = [0.]

rundata = setrun('geoclaw')

for tyi in tauy_i:
    for tyr in tauy_r:
        for rc in remold_coeff:
                rundata.probdata.add_param('tauy_i'      ,  tyi      ,'Initial yield strength')
                rundata.probdata.add_param('tauy_r'      ,  tyr      ,'Residual yield strength')
                rundata.probdata.add_param('remold_coeff',  rc       ,'Remolding param. gamma')

                rundata.write()
                outdir = '/media/jihwan/jik/work/ngi/Numerical_models/VP2HD_Case_Studies/Traenadjupet/results/vol700/c1/slide/_%dto%d_gamma%s' %(tyi,tyr,rc)
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                runclaw.runclaw('xgeoclaw', outdir, True, False,'.', False, False, None)

                fg_outdir = '/media/jihwan/jik/work/ngi/Numerical_models/VP2HD_Case_Studies/Traenadjupet/topo/vol700/c1'
                if not os.path.exists(fg_outdir):
                    os.makedirs(fg_outdir)
                outfile = fg_outdir+'/dtopo_%dto%d_gamma%s.tt3' %(tyi,tyr,rc)
                dtopo_tt3_v2(outdir,outfile)

