from pylab import *
from pyspec import fit, fitfuncs
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os, pdb, glob
from processCCD_image import *
from datetime import datetime

def get_frame_list(file_info, ai_info):
    """
    Take output from get_scan_info(folder)
    and returns a nested list of frames
    Frames are ordered by energy and grouped if they are at the same energy    
    """
    
    # remove backgrounds
    not_background = ai_info[:,3] == 1.0
    file_info = file_info[not_background]
    ai_info = ai_info[not_background]

    # sort by beamline energy
    bl_energy_order = np.argsort(ai_info[:,1])
    file_info = file_info[bl_energy_order]
    ai_info = ai_info[bl_energy_order]
    
    frame_list = list(ai_info[:,0])
    bl_energy_list = list(ai_info[:,1])
    
    frames = []
    while len(frame_list) > 0:
        frame_set = [frame_list.pop(0)]
        currE = bl_energy_list.pop(0)
        for j in range(len(frame_list)):
            try:
                if np.round(bl_energy_list[j], decimals=1) == np.round(currE, decimals=1):
                    frame_set.append(frame_list.pop(j))
                    bl_energy_list.pop(j)
            except IndexError:
                pass
        frames.append(frame_set)   
    return frames

def get_scans_info(folder_list):
    file_info, ai_info = get_scan_info(folder_list[0])
    for folder in folder_list[1:]:
        file_info_, ai_info_ = get_scan_info(folder)
        file_info = np.vstack((file_info, file_info_))
        ai_info = np.vstack((ai_info, ai_info_))
    return file_info, ai_info

def get_scan_info(folder):
    """ Finds all fits files in a given folder and returns:
    file_info - array with pathnames and filenames
    ai_info - array with frames, bl_energy, spec_energy, andor
    """
    #pdb.set_trace()
    searchterm = os.path.join(folder, '*.fits')
    path_list = []
    for path in glob.glob(searchterm):
        if 'Captured' not in path:
            path_list.append(path)
    path_list = np.array(path_list)
    
    
    file_list = []
    for path in path_list:
        file_list.append(os.path.split(path)[1])

    file_list = np.array(file_list)
    
    aifilename = folder + file_list[0][0:-10] + 'AI.txt'  
    ai_info = np.loadtxt(aifilename, usecols = (0,1,2,3), skiprows = 10, dtype=float)
    if len(ai_info.shape) == 1:
        fake = ai_info + np.NaN 
        ai_info = np.vstack((ai_info, fake))
    #pdb.set_trace()
    
    frames = []
    bl_energy = []
    spec_energy = []
    andor = []   # 0 for background 1 for data
    scans = []
    #pdb.set_trace()
    for filename in file_list:
        fileno  = int(filename[-10:-5])
        #try:
        row_index = ai_info[:,0] == fileno
        #except IndexError:
        #row_index = ai_info[0] == fileno
        frames.append(ai_info[row_index, 0][0])
        bl_energy.append(ai_info[row_index, 1][0])
        spec_energy.append(ai_info[row_index, 2][0])
        andor.append(ai_info[row_index, 3][0])
        #print filename[-15:-11]
        scans.append(int(filename[-15:-11]))
    
    #pdb.set_trace()
    file_info = np.vstack((path_list, file_list)).transpose()
    ai_info = np.vstack((frames, bl_energy, spec_energy, andor, scans)).transpose()
    return file_info, ai_info

def get_ai(ai_info, frame):
    frames = np.array(ai_info[:, 0])
    row_index = frames == frame
    return ai_info[row_index, :][0]  # hacK?

def get_file(ai_info, file_info, frame):
    row_index = ai_info[:,0] == frame
    return file_info[row_index, 0][0] # hacK?

def get_a_background(ai_info, file_info):
    row_index = ai_info[:,3] == 0
    return  file_info[row_index, 0][0]# hacK?

def plot_3d(M, epoints, eloss, limits = [], xlim = [], ylim = [], 
            title = '', save_fig = False, outname = '3dmap.png'):
    """
    Plots 3d incident energy x energy loss graph. M is a matrix with the data
    as M[incident energy, energy loss]. epoints is the incident energies, and
    eloss the energy loss.
    """

    plt.figure()
    if len(limits[:]) == 0:
        fig = plt.contourf(epoints, eloss, M, 100)
        plt.colorbar()
    else:
        tks = []
        for i in range (6):
            tks.append(limits[0] + i*(limits[1]-limits[0])/5)
            
        z = np.linspace(limits[0],limits[1], 100, endpoint=True)
        fig = plt.contourf(epoints, eloss, M, z)
        plt.colorbar(ticks=tks)

    plt.xlabel('Incident Energy (eV)', fontsize=15)
    plt.ylabel('Energy Loss (eV)', fontsize=15)
    if len(xlim) != 0:    
        plt.xlim(xlim[0],xlim[1])
    if len(ylim) != 0:
        plt.ylim(ylim[0], ylim[1])
    plt.title(title)
    if save_fig is True:
        plt.savefig(outname, format='png', dpi=1000)
    
def calc_xas(infile, ewindow):
    """
    Calculate the XAS from summing the RIXS data along eloss. "infile" contains
    the path to the eloss x intensity scans for each incident energy. "ewindow"
    (=[e0,ef]) sets the the limits of energy loss to be integrated (e0 to ef).
    
    INCLUDE ERROR CALCULATION!
    """
    
    for i in range (len(infile[:])):
        
        data = loadtxt(infile[i])
        
        if i == 0:
            xas = np.zeros((len(infile[:])))
        
        for j in range (len(data[:,0])):
            if data[j,0] > ewindow[0] and data[j,0] < ewindow[1]:
                xas[i] = xas[i] + data[j,1]
        
    return xas
    
def calc_eloss(infile, energy, ewindow):
    """
    Calculate the incident energy x intensity for a specific energy loss. 
    "infile" contains the path to the eloss x intensity scans for each incident 
    energy. "energy" is the energy loss to be used, and "ewindow" sets the range
    in energy loss to be binned. For a given incident energy, the energy loss
    from "energy"-"ewindow to "energy"+"ewindow" is averaged.
    """
    
    eloss = []
    error = []
    for j in range (len(infile)):
        data = loadtxt(infile[j])
        eloss.append(0.0)
        error.append(0.0)
        w = 0
        for i in range(len(data)):
            if np.abs(data[i,0]-energy) < ewindow:
                eloss[j] = eloss[j] + data[i,1]
                error[j] = error[j] + data[i,2]**2
                w = w + 1
        eloss[j] = eloss[j]/w + 0.0
        error[j] = error[j]**0.5/w
        
    eloss = np.array(eloss)
    error = np.array(error)
        
    return eloss, error
                

def build_setup(day, base_name, scan_number):
    """
    The CCD images are processed based on the "setup". This file contains
    all scans for a given sample and temperature which should be the same (apart
    from different incident energies). It consists of 5 colums:"incident energy",
    "file number", "background position", "shutter status" (0-closed, 1-open), 
    "scan order" (first in "scan_number" is 0).
    
    "day", "base_name" and "scan_number" are vectors with info about the scans.
    IT ASSUMES THAT THE SCANS ARE SAVED IN: data/CCD Scan "scan_number"/, AND
    THAT THEY HAVE THE SAME NAME FORMAT AS THEY ARE SAVED AT ALS.
    
    Outputs "setup" and "fname", the later has the path to all scans in "setup".
    """
    
    for w in range(len(scan_number)):
        
        if len(day[:]) == 1:
            begin = 'data/' + day[0]
        else:
            begin = 'data/' + day[w]
            
        if len(scan_number[:]) == 1:
            mid = '/CCD Scan ' + scan_number[0]
            end_file = scan_number[0]
            end_ai = scan_number[0] + '-AI.txt'
        else:
            mid = '/CCD Scan ' + scan_number[w]
            end_file = scan_number[w]
            end_ai = scan_number[w] + '-AI.txt'
            
        if len(base_name[:]) == 1:
            mid2 = '/' + base_name[0] + '_'
            ai = begin + mid + '/' + base_name[0] + '_' + end_ai
        else:
            mid2 = '/' + base_name[w] + '_'
            ai = begin + mid + '/' + base_name[w] + '_' + end_ai
        
        ai_info = loadtxt(ai, usecols = (0,1,3), skiprows = 10)  
    
        scan = np.zeros((len(ai_info[:,0]),1))
        
        scan[:] = w
        
        ai_info = np.hstack((ai_info,scan))
             
        if w == 0:
            fname = []
            for i in range (len(ai_info[:,0])):
                fname.append(begin + mid + mid2 + end_file + '-%05d.fits' %(ai_info[i,0]))
            
            setup = np.empty((len(ai_info[:,0]), 4))
            setup[:,:] = ai_info[:,:]
        if w != 0:
            for i in range (len(ai_info[:,0])):
                fname.append(begin + mid + mid2 + end_file + '-%05d.fits' %(ai_info[i,0]))
            
            setup = np.vstack((setup,ai_info))
    
    tmp_bkg = np.zeros((len(setup[:,0]),1))
    
    setup = np.hstack((setup,tmp_bkg))
    
    setup[:,[0,1]] = setup[:,[1,0]]
    setup[:,[2,3]] = setup[:,[3,2]]
    setup[:,[2,4]] = setup[:,[4,2]]
    
    return setup, fname
    
def set_setup(day, base_name, scan_number, filenos, e_inc, bkg = []):
    """
    Similar to the "build_setup" (read the one above!). This is used to set
    the "setup" using specific scans (in filenos). Besides "setup and "fname",
    it exports "epoints" that have the incident energies.
    """
    
    setup = np.zeros((len(filenos[:]),5))
    fname = []
    epoints = []
    for i in range(len(filenos[:])):
        setup[i,0] = i*0.1+e_inc
        setup[i,1] = filenos[i]
        
        if len(bkg[:]) != 0:
            setup[i,2] = i + len(filenos[:])
            fname.append('data/' + day[i] + '/CCD Scan ' + scan_number[i] + '/' +
        base_name[i] + '_' + scan_number[i] + '-%05d.fits' %(filenos[i]))
 
        else:
            setup[i,2] = -1
            fname = bkg
        
        setup[i,3] = 1
        setup[i,4] = 0
        epoints.append(i*0.1+e_inc)
    
    if len(bkg[:]) != 0:
        for i in range(len(filenos[:])):
            fname.append('data/' + day[i] + '/CCD Scan ' + scan_number[i] + '/' +
            base_name[i] + '_' + scan_number[i] + '-%05d.fits' %(bkg[i]))
            
    
    return setup, fname, epoints

def set_bkg(setup, bkg_method = []):
    """
    Takes the "setup" and finds the closest background for every scan. It can
    do it in bkg_method: overall_closest, and same_file. In the former, for each
    scan, the closest background scan as displayed in "setup" will be assigned,
    while in the later, the closest background within the same "scan_number" will
    be assigned. A different bkg_method can be assigned for each "scan_number".
    The "same_file" method is used by default.
    
    THERE MUST BE AN EASIER WAY TO DO THIS!
    """
    repeat = 0
    k = 0
    for i in range (len(setup[:,0])):
        if setup[i,3] == 1.0 and repeat == 0: #finds repeat number and sets first bkg
            k = i+1
            repeat = 1
            while k > 0 and k < len(setup[:,0]) and setup[k,3] != 0.0:
                if np.abs(setup[i,0] - setup[k,0]) < 0.005:
                    k = k + 1
                    repeat = repeat + 1
                else:
                    k = -1
    
    #print repeat
    num_scans = int(np.max(setup[:,4]) + 1)
    
    if len(bkg_method) == 0:
        for w in range (num_scans):
            bkg_method.append('same file') #if no bkg method is given will use 'same file'                   
                        
    for w in range (len(bkg_method[:])):                    
        if 'overall closest' == bkg_method[w]:
            for i in range (len(setup[:,0])/repeat):
                i = repeat*i
                if setup[i,3] == 1.0 and setup[i,4] == w:
                    diff = 9999        
                    for j in range(len(setup[:,0])/repeat): #find closest bkg
                        j = j*repeat
                        if setup[j,3] == 0.0 and np.abs(i-j) < diff and np.abs(i-j) != 0:
                            diff = i - j
                    
                    for j in range(repeat):
                        setup[i+j,2] = i - diff + j
        
        if 'same file' == bkg_method[w]:
            for i in range (len(setup[:,0])/repeat):
                i = repeat*i
                if setup[i,3] == 1.0 and setup[i,4] == w:
                    diff = 9999        
                    for j in range(len(setup[:,0])/repeat): #find closest bkg
                        j = j*repeat
                        if setup[j,3] == 0.0 and np.abs(i-j) < diff and np.abs(i-j) != 0 and setup[j,4] == setup[i,4]:
                            diff = i - j
                    
                    for j in range(repeat):
                        setup[i+j,2] = i - diff + j

    return setup
    
def get_epoints(setup, delta = 0.005):
    """
    Collects all different incident energies for each scan, and return them into
    "epoints". All energies within +- delta are taken to be the same.
    """
    
    aux = np.empty((len(setup[:,0]),len(setup[0,:])))    
    
    aux[:,:] = setup[:,:]
    
    s = "i8"
    for i in range(len(aux[0,:])-1):
        s+=',i8'
        
    aux.view(s).sort(order=['f0'], axis=0)
    
    epoints = []
    epoints.append(aux[0,0])
    for i in range(1,len(aux[:,0])):
        if np.abs(aux[i,0] - aux[i-1,0]) > delta:
            epoints.append(aux[i,0])
    
    return epoints
            
def process_setup(setup, fname, epoints = [], offset = 100, save_files = True,
                    to_plot = [], plot = 0, lim = [1,9999], dx = -1, statistic = 'mean',
                    base_fileout = 'RIXS_data_', clean = [], curvature = [],
                    cal = [], cal_func = 'slit', scan_time = 10, pol_corr = []):

    """
    Process images defined on setup. 
    
    Setup file consists of 5 colums: Energy (or an indice),
    file number, background scan (position of it in fname), shutter (0 closed, 1 open), and
    scan number (its position in the scan_number list, not relevant for this part).
    
    fname - list with path for the images.
    
    epoints - Energies (or indices) to be used.
    
    offset - offset for plotting.
    
    save_files = True/False
    
    to_plot = [] - Select energies/indices to plot
    
    plot = 0/1/2 - 0 = processed spectrum, 1 = unprocessed spectrum, 2 = background
    
    lim = [] - Limits of the image to be used.
    
    dx, statistic - variables for binning.
    
    base_file_out - Path to saved files, the energy/indice in epoints will be added
    to the end of the path.
    
    clean = 0 - Cosmic rays removal parameter.
    
    curvature= [] - Parameters for curvature correction.
    
    cal = [] - Parameters for pixel -> energy calibration
    
    cal_func = 'slit' or 'gauss' - function to be fitted to the elastic peak
    in order to align it to zero energy loss.
    
    scan_time = 10 - Time (in minutes) collecting each scan. This is important
    for the polinomial correction time calculation.
    
    pol_corr = [] - Polinomial correction for the background. This is needed
    whenever the CDD warms up since the background keeps changing over an 
    extended period.
    """
    plt.figure()  
    
    #Images are processed according to their incident energy
    for i in range (len(epoints)):
        
        print 'Started energy %0.1lf eV...' %(epoints[i])
        #t0 = datetime.now()
        
        fname_list = []
        fname_list_BG = []
        time_BG = []
        #Find all images with same incident energy and their backgrounds
        for j in range (len(setup[:,0])):
            if setup[j,3] == 1.0 and np.abs(setup[j,0] - epoints[i]) < 0.005:
                fname_list.append(fname[j])
                fname_list_BG.append(fname[int(setup[j,2])])
                time_BG.append(scan_time*j)
        
        print 'Processing %d scans...' %(len(fname_list))
        
        ims = CCD(fname_list, photon_E=epoints[i], poly_order=2, binpix=2,
                  fname_list_BG = fname_list_BG, exclude=lim)
                  
        if len(curvature[:]) == 3:
            ims.curvature = curvature
                 
       # t0_clean =  datetime.now()
        if len(clean) != 0:
            ims.clean_std_new(clean) # This process is taking ~ 0.3 sec/scan!
        else:
            ims.clean_std(clean[0])
        #tf_clean = datetime.now()
        
        """
        If pol_corr will be done, it renormalize the background images based on 
        the integrated background.
        """
        if len(pol_corr) != 0:        
            ims.get_specs()
    
            w = 0
            for spec in ims.BGspecs:
                sum_bkg = 0
                for j in range (len(spec[1][:])):
                    sum_bkg = sum_bkg + spec[1][j]
                pol = 0
                for j in range (len(pol_corr)):
                    pol = pol_corr[j]*time_BG[w]**(len(pol_corr)-j-1) + pol
                #print pol/sum_bkg
                ims.BGimages[w] = ims.BGimages[w]*pol/sum_bkg
                w = w + 1
                
        ims.sub_backgrounds()
        ims.get_specs()
        cal[0] = ims.shift_e(func = cal_func)
        ims.sum_specs()
        
        if len(cal[:]) == 2:
            ims.calibrate(cal[0],cal[1])
            
        if dx <= 0.0:
            ims.bin_points(dx, statistic = statistic)
        
        if len(to_plot[:]) == 0:
            if plot == 0:
                ims.plot_spectrum(offset = i*offset)
            if plot == 1:
                ims.plot_specs()
            if plot == 2:
                ims.plot_BGimages()
        else:
            for w in range (len(to_plot)):
                if np.abs(epoints[i] - to_plot[w]) < 0.005:
                    if plot == 0:
                        ims.plot_spectrum(offset = i*offset)
                    if plot == 1:
                        ims.plot_specs()
                    if plot == 2:
                        ims.plot_BGimages()

        if save_files is True:
            fileout = base_fileout + '%0.1lfeV.txt' %(epoints[i])
            ims.fileout = fileout
            ims.write_file()
            print 'saved file: ' + fileout
        
        
        #tf =  datetime.now()

        #print 'Clean_std_new time = ', tf_clean - t0_clean
        #print 'Total time = ', tf-t0
        print 'Done!'

def calc_bkg_corr(setup,fname, scan_time = 10, lim = [1,9999], clean = 1E5,
                  curvature = [], order = 0):
                      
    """
    Calculates the necessary background correction based on the integrated 
    intensity of the background scans over time. It fits a "order" polynomial
    to the integrated background x time, and retuns it as a p[len(order)] vector.
    """

    sum_bkg = []
    time = []
    j = 0
    for i in range (len(fname)):
        
        fname_list = []
        fname_list_BG = []
        if setup[i,3] == 0.0:
            fname_list.append(fname[i])
            
            ims = CCD(fname_list, photon_E=70, poly_order=2, binpix=2,
                      fname_list_BG = fname_list_BG, exclude=lim)
            
            if len(curvature[:]) == 3:
                ims.curvature = curvature
            
            ims.clean_std_new(clean)
            
            ims.get_specs()
            ims.sum_specs()
            
            sum_bkg.append(0.0)
            time.append(scan_time*i*1.0)
            for w in range (len(ims.spectrum[1])):
                sum_bkg[j] = sum_bkg[j] + ims.spectrum[1][w]
                
            j = j + 1
    
    p = np.polyfit(time,sum_bkg,order)    
    
    xnew = np.linspace(time[0], time[len(time)-1], 100)
    ynew = []
    for i in range (len(xnew)):
        y = 0
        for j in range (len(p)):
            y = p[j]*xnew[i]**(len(p)-j-1) + y
        ynew.append(y)
    
    plt.plot(time,sum_bkg, 'bs')
    plt.plot(xnew,ynew)
    plt.xlabel('Time')
    plt.ylabel('Bkg')
    time_lim = 0.1*np.abs(time[len(time)-1]-time[0])
    bkg_lim = 0.1*(np.max(sum_bkg) - np.min(sum_bkg))
    plt.xlim(time[0]-time_lim,time[len(time)-1]+time_lim)
    plt.ylim(np.min(sum_bkg) - bkg_lim, np.max(sum_bkg) + bkg_lim)

    return p      
        
        
        