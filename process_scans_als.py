from pylab import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from ALSutils import *
from datetime import datetime

"""
 For day, base_name and scan_number, if only one exist just enter once,
 otherwise it will have to be repeated. For instance, for 3 files in one day,
 just enter that day, but if there are 2 in one day, and one in another,
 it requires: day = [day1, day1, day2]
"""

t0 =  datetime.now()

day = ['2015 02 08']
       
base_name = ['LSNO_SLSL_727']

scan_number =['2611']
              

dopol = False        
polcorr = [4.60165401e-01, -1.25130458e+03, 3.88839849e+07]
process_data = True
plot3d = False

setup, fname = build_setup(day, base_name, scan_number)

setup = set_bkg(setup)

epoints = get_epoints(setup)

#epoints = [75.0]

#print setup

if dopol == True:
    polcorr = calc_bkg_corr(setup,fname, scan_time = 10, lim = [20,62], clean = [6,5,3], 
                  curvature = [1.63894559e-03 , -3.41817820e-01 ,  8.73668115e+02],
                  order = 2)
    print 'pol. correction parameters = ', polcorr 
    
if process_data == True:
    process_setup(setup, fname, epoints, offset = 50, save_files = False, to_plot = [], plot = 0,
              lim = [20,62], dx = 0.025, statistic = 'mean', 
              base_fileout = base_name[0] +  '/' + base_name[0] + '_', clean = [6,5,3],
              curvature = [1.63894559e-03 , -3.41817820e-01 ,  8.73668115e+02], 
              cal = [873, 5.719e-3], pol_corr = polcorr)

if plot3d == True:           
    eloss = []
    for i in range (0,len(epoints[:])):   
        inname = base_name[0] + '/' + base_name[0] + '_' + '%0.1lf' %(epoints[i]) + 'eV.txt'
        data = loadtxt(inname)
        if i == 0:
            M = np.zeros((len(data[:,0]),len(epoints[:])))
            eloss[:] = data[:,0]
            data[:,1] = 3/2*data[:,1]
            
        f = intp.interp1d(data[:,0], data[:,1], kind='linear', 
                                bounds_error=False, fill_value=0.0)
            
        M[:,i] = f(eloss)
    
    title = r'LaSrNiO$_4$ SLSL on LSAO, #727, 10 fu, $\theta$ = 60$^{\circ}$'
    outname = base_name[0] + '/' + base_name[0] + '_' + '3dmap_aligned.png' 
    
    plot_3d(M, epoints, eloss, limits = [150,350], xlim = [], title = title, save_fig = True, outname = outname)
    
tf =  datetime.now()

print tf-t0
