from pylab import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from ALSutils import *
"""
For day, base_name and scan_number, if only one exist just enter once,
otherwise it will have to be repeated. For instance, for 3 files in one day,
just enter that day, but if there are 2 in one day, and one in another,
it requires: day = [day1, day1, day2]
"""
day = ['2015 02 08']
       
base_name = ['LSNO_SLSL_727']

scan_number =['2611']

setup, fname = build_setup(day, base_name, scan_number)

setup = set_bkg(setup)

epoints = get_epoints(setup)

print setup

process_setup(setup, fname, epoints, offset = 0, save_files = False, to_plot = [69,70], plot = 2,
              lim = [20,62], dx = -1, statistic = 'mean', 
              base_fileout = base_name[0] + '/' + base_name[0] + '_', clean = 999999,
              curvature = [1.63894559e-03 , -3.41817820e-01 ,  8.73668115e+02], 
              cal = [873, 5.719e-3])

#savetxt("setup.txt", setup[:,:], fmt = '%0.1lf')

for i in range (0,len(epoints[:])):   
    inname = base_name[0] + '/' + base_name[0] + '_' + '%0.1lf' %(epoints[i]) + 'eV.txt'
    data = loadtxt(inname)
    if i == 0:
        M = np.zeros((len(data[:,0]),len(epoints[:])))
        data[:,1] = data[:,1]*3/2
        
    M[:,i] = data[:,1]

eloss = []
eloss[:] = data[:,0]
title = r'LaSrNiO$_4$ SLSL on LSAO, #727, 10 fu, $\theta$ = 60$^{\circ}$'
outname = base_name[0] + '/' + base_name[0] + '_' + '3dmap.png' 

plot_3d(M, epoints, eloss, limits = [180,340], title = title, save_fig = False, outname = outname)
