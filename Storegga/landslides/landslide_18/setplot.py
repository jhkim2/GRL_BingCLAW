""" 
Set up the plot figures, axes, and items to be done for each frame.
This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters. 
""" 

import os

import numpy as np

jik_outdir = '../_output'
tau_i = 12000.
tau_r = 2000.
param_shear = 5e-4

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps, geoplot
    import matplotlib.pyplot as plt

    plotdata.clearfigures()  # clear any old figures,axes,items data


    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from pyclaw.plotters import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)
    

    #-----------------------------------------
    # Figure for pcolor plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='pcolor', figno=0)
    #plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotfigure.kwargs = {'figsize': (10,10.)}
    plotaxes.scaled = False
    plotaxes.xlimits = [ -2. , 7.] #'auto'
    plotaxes.ylimits = [ 62. , 67.] #'auto'
    #plotaxes.afteraxes = addgauges
    
    def topo_q(current_data):
        q=current_data.q
        b=current_data.aux
        var = q[0,:,:]+b[0,:,:]
        return var
        
    def topo(current_data):
        b=current_data.aux
        var = b[0,:,:]
        return var

    drytol = 1e-4

    def land0(current_data):
        b=current_data.aux
        var = b[0,:,:]
        var1 = np.where(var>0., var , np.nan)
        var = np.ma.array (var1, mask=np.isnan(var1))
	return var

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = land0 #geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 2000.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0
    
    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.outdir = jik_outdir
    def depth(current_data):
        from numpy import ma, where
        q = current_data.q
        h = q[0,:,:]

        var1 = np.where(h>1., h , np.nan)
        var = np.ma.array (var1, mask=np.isnan(var1))
	return var
	
    #plotitem.plot_var = surface_or_depth
    plotitem.plot_var = depth
    #plotitem.pcolor_cmap = 'Reds'
    colors = plt.cm.get_cmap('jet')
    colors = plt.cm.get_cmap('copper_r')
    colors = plt.cm.get_cmap('gist_earth_r')
    plotitem.pcolor_cmap = colors
    plotitem.pcolor_cmin =  1.
    plotitem.pcolor_cmax =  100.
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1

    # Add contour lines of bathymetry:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    #plotitem.outdir = jik_outdir
    plotitem.show = True
    plotitem.plot_var = topo_q #geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels =[-4000,-3000,-2500,-2000,-1500,-1000,-900,-800,-700,-600,-500,-400,-300,-200,-100] 
    plotitem.amr_contour_colors = ['y']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [1]  # show contours only on finest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0

    def fixup(current_data):
        import pylab
        #addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        #pylab.plot((2,6),(65.2,63.),'k-',linewidth=2)
        pylab.title('Height at %4.2f hours' % t, fontsize=20)
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)
    plotaxes.afteraxes = fixup


    #-----------------------------------------
    # Figure for yield stress (tau) plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='tau', figno=1)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'yield stress'
    plotfigure.kwargs = {'figsize': (7,7)}
    plotaxes.scaled = False
    plotaxes.xlimits = [-3, 7.0] #'auto'
    plotaxes.ylimits = [62.,69.] #'auto'
    #plotaxes.afteraxes = addgauges
    
    def tau_y(current_data):
        q=current_data.q
        h = q[0,:,:]
        gamma = q[5,:,:]
        var1 = np.where(h>0.,tau_r+(tau_i-tau_r)*np.exp(-gamma*param_shear), np.nan)
        var = np.ma.array (var1, mask=np.isnan(var1))
        return var

    def fixup(current_data):
        import pylab
        #addgauges(current_data)
        t = current_data.t
        t = t / 60.  # mins
        pylab.title('Yield stress at %4.2f minutes' % t, fontsize=20)
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)
    plotaxes.afteraxes = fixup

    drytol = 1e-4
    
    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = tau_y
    colors = plt.cm.get_cmap('jet_r')
    plotitem.pcolor_cmap = colors
    plotitem.pcolor_cmin =  tau_r
    plotitem.pcolor_cmax =  tau_i
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1

    def topo_q(current_data):
        q=current_data.q
        b=current_data.aux
        var = q[0,:,:]+b[0,:,:]
        return var
        
    def topo(current_data):
        b=current_data.aux
        var = b[0,:,:]
        return var

    # Add contour lines of bathymetry:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = topo_q #geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = [-4000,-3000,-2500,-2000,-1500,-1000,-900,-800,-700,-600,-500,-400,-300,-200,-100] 
    plotitem.amr_contour_colors = ['y']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [0]  # show contours only on finest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    def land0(current_data):
        b=current_data.aux
        var = b[0,:,:]
        var1 = np.where(var>0., var , np.nan)
        var = np.ma.array (var1, mask=np.isnan(var1))
	return var

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = land0 #geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 2000.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0
    

    #-----------------------------------------
    # Figure for pcolor plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='vector field', figno=2)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'vp'
    #plotfigure.kwargs = {'figsize': (18.5,10.5)}
    plotaxes.scaled = False
    plotaxes.xlimits = 'auto' # [-5.,7.]
    plotaxes.ylimits = 'auto' #[62.,69.] #'auto'
    #plotaxes.afteraxes = addgauges

    def fixup(current_data):
        import matplotlib.pyplot as plt
        import pylab
        from matplotlib.pyplot import cm
        t = current_data.t
        t = t / 3600.  # hours
        q=current_data.q
        X = current_data.x
        Y = current_data.y
        h = q[0,:,:]
        u = np.where(h>0., q[3,:,:]/h, np.nan)
        u = np.ma.array (u, mask=np.isnan(u))
        v = np.where(h>0., q[4,:,:]/h, np.nan)
        v = np.ma.array (v, mask=np.isnan(v))
        speed = np.sqrt(u**2 + v**2)
        UN = u/speed
        VN = v/speed        
        plt.quiver(X, Y, UN, VN, speed, cmap=cm.seismic, headlength=7)
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)
        pylab.colorbar()
    plotaxes.afteraxes = fixup

    #-----------------------------------------
    # Figure for Cross section plot
    #-----------------------------------------

    plotfigure = plotdata.new_plotfigure(name='Slice', figno=3)
    #plotfigure.kwargs = {'figsize': (20,8)}
    plotfigure.show = False

    def x_slice(current_data):
        x = current_data.x[:,0]
        y = current_data.y[0,:]
        q = current_data.q
        var = np.zeros(len(x))
        for i in range(len(x)):
            if (abs(x[i])<6):
                j = int(-.55*(i-420)+90)
                var[i] = np.where(q[0,i,j]>0.,q[6,i,j],np.NaN)
                var[i] = q[6,i,j]
        return x,var

    def y_slice(current_data):
        y = current_data.y[0,:]
        q = current_data.q
        eta = q[6,1,:]
        return y,eta

    def B_slice(current_data):
        x = current_data.x[:,0]
        y = current_data.y[0,:]
        q = current_data.q
        var = np.zeros(len(x))
        for i in range(len(x)):
            if (abs(x[i])<6):
                j = int(-.55*(i-420)+90)
                var[i] = q[6,i,j]-q[0,i,j]
        return x,var

    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('water surface')
    plotaxes.title = 'eta'
    plotaxes.xlimits = [2,6]
    plotaxes.ylimits = [-2500.,0.]
    
    # Water Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.outdir = jik_outdir
    plotitem.map_2d_to_1d = x_slice
    plotitem.color = 'r'
    plotitem.plotstyle = '-'
    plotitem.kwargs = {'linewidth':2}

    # Water Surface
    #plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    #plotitem.outdir = '../bous_cfl9_dx05'
    #plotitem.map_2d_to_1d = x_slice
    #plotitem.color = 'r'
    #plotitem.plotstyle = '--'
    #plotitem.kwargs = {'linewidth':2}


    # Topography
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.outdir = jik_outdir
    plotitem.map_2d_to_1d = B_slice
    plotitem.color = 'k'
    plotitem.plotstyle = '-'
    plotitem.kwargs = {'linewidth':2}

    def fixup(current_data):
        import pylab
        import numpy as np
        t = current_data.t
        t = t / 3600.  # hours
        pylab.title('Cross section at %4.2f hours' % t, fontsize=20)
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)
    plotaxes.afteraxes = fixup


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=300, \
                    type='each_gauge')
                    
    plotfigure.show = True
    plotfigure.kwargs = {'figsize': (12,2)}

    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-0.01, .2]
    plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    #plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.outdir='../_output'
    #plotitem.plot_var = 3
    #plotitem.plotstyle = 'b-'
    #plotitem.kwargs = {'linewidth':1}
    
    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'
    plotitem.kwargs = {'linewidth':1}

    # Plot topo as green curve:
    #plotitem = plotaxes.new_plotitem(plot_type='1d_plot')

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo
        
    #plotitem.plot_var = gaugetopo
    #plotitem.plotstyle = 'g-'
    
    def add_zeroline(current_data):
        from pylab import plot, legend
        t = current_data.t
        plot(t, 0*t, 'k')


    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = []             # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.format = 'ascii'                # Format of output
    # plotdata.format = 'netcdf'             

    return plotdata

    
