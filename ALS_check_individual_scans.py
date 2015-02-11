from pylab import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from processCCD_image import *
from ALSutils import *

day = []
scan_number = []
base_name = []
bkg = []
for i in range(4):
    day.append('2015 02 08')
    scan_number.append('2611')
    base_name.append('LSNO_SLSL_727')
    bkg.append(2)

filenos = [13,37,55]       

#bkg = [2,25]
e_inc = 65

setup, fname, epoints = set_setup(day, base_name, scan_number, filenos, e_inc, bkg = bkg)

#print setup
#print epoints
#print fname

process_setup(setup, fname, epoints, offset = 0, save_files = False, plot = 2,
              lim = [20,62], dx = -1, statistic = 'mean', 
              base_fileout = 'RIXS_data_', clean = 70,
              curvature = [1.63894559e-03 , -3.41817820e-01 ,  8.73668115e+02], 
              cal = [873, 5.719e-3])
