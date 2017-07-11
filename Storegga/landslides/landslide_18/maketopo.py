
"""
Module to create topo and qinit data files for this example.
"""

from clawpack.geoclaw import topotools
from numpy import *

slope_angle = -2.*pi/180

def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints = 601
    nypoints = 201
    xlower = -1000.e0
    xupper = 11000.e0
    ylower = -2000.e0
    yupper =  2000.e0
    outfile= "slope.tt3"     
    topotools.topo3writer(outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def make_q1():
    """
    Output topography file for the entire domain
    """
    nxpoints = 201
    nypoints = 201
    xlower = -1000.e0
    xupper =  3000.e0
    ylower = -2000.e0
    yupper =  2000.e0
    outfile= "q1.tt3"     
    topotools.topo3writer(outfile,q1,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def topo(x,y):

    x1 =  0.
    z1 = -100.
    z2 = -1500.
    x2 = (z2-z1)/slope_angle + x1
    from numpy import where
    z = where(x<x2,slope_angle*(x-x1)+z1,z2)
    return z

def q1(x,y):

    from numpy import where
    # 2d test
    z = where((x-500)**2+y**2<=500.0**2, 25.0*(1.0-((x-500.0)/500.0)**2-(y/500.0)**2),0.)
    # 1d test
    #z = where((x-750)**2<=250.0**2, 24.0*(1.0-((x-750.0)/250.0)**2),0.) 
    #z = where(abs(x-450.)<450., 10.,0.) 
    return z

if __name__=='__main__':
    maketopo()
    make_q1()

