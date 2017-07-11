
def setgeo(rundata):
    """
    Set GeoClaw specific runtime parameters.
    For documentation see ....
    """

    try:
        geodata = rundata.geodata
    except:
        print "*** Error, this rundata has no geodata attribute"
        raise AttributeError("Missing geodata attribute")

    # == setgeo.data values ==
    geodata.igravity = 1
    geodata.gravity = 9.81
    geodata.icoordsys = 1

    # == settsunami.data values ==
    geodata.sealevel = 0.
    geodata.drytolerance = 1.e-3
    geodata.wavetolerance = 1.e-2
    geodata.depthdeep = 1.e2
    geodata.maxleveldeep = 3
    geodata.ifriction = 0
    geodata.coeffmanning = 0.
    geodata.frictiondepth = 20.

    # == settopo.data values ==
    geodata.topofiles = []
    geodata.topofiles.append([2, 1, 1, 0., 1.e10, 'bowl.topotype2'])

    # == setdtopo.data values ==
    geodata.dtopofiles = []

    # == setqinit.data values ==
    geodata.qinitfiles = []
    # for qinit perturbations append lines of the form  
    #   [minlev, maxlev, fname]
    geodata.qinitfiles.append([1, 2, 'hump.xyz'])

    # == setregions.data values ==
    geodata.regions = []
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    geodata.regions.append([1, 2, 0.e0, 1.e10, -100.,100., -100.,100.])

    # == setgauges.data values ==
    geodata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y]
    geodata.gauges.append([1, 50., 50.])

    # == setfixedgrids.data values ==
    geodata.fixedgrids = []
    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,\
    #  ioutarrivaltimes,ioutsurfacemax]

    return rundata


