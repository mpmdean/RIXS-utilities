from pylab import *
from pyspec import fit, fitfuncs
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os, pdb, glob
from processCCD_image import *

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
    Needs improvement!
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
    #plt.colorbar()
    plt.xlabel('Incident Energy (eV)', fontsize=15)
    plt.ylabel('Energy Loss (eV)', fontsize=15)
    if len(xlim) != 0:    
        plt.xlim(xlim[0],xlim[1])
    if len(ylim) != 0:
        plt.ylim(ylim[0], ylim[1])
    plt.title(title)
    if save_fig is True:
        plt.savefig(outname, format='png', dpi=1000)
    #plt.show()
    
def calc_xas(infile, outfile, epoints, ewindow):
    
    for i in range (len(infile[:])):
        
        data = loadtxt(infile[i])
        
        if i == 0:
            xas = np.zeros((len(infile[:])))
        
        for j in range (len(data[:,0])):
            if data[j,0] > ewindow[0] and data[j,0] < ewindow[1]:
                xas[i] = xas[i] + data[j,1]
        
    return xas

def build_setup(day, base_name, scan_number):
    
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
    THERE MUST BE AN EASIER WAY TO DO THIS!
    """
    dark_number = 0
    for i in range(len(setup[:,0])):
        if setup[i,3] == 0.0:
            dark_number = dark_number + 1
            
    dark = np.zeros((dark_number,2))
    for i in range(len(setup[:,0])):
        if setup[i,3] == 0.0:
            dark[dark_number-1,0] = setup[i,1]
            dark[dark_number-1,1] = setup[i,4]
            dark_number = dark_number - 1
    
    repeat = 0
    k = 0
    for i in range (len(setup[:,0])):
        if setup[i,3] == 1.0 and repeat == 0: #finds repeat number and sets first bkg
            k = i+1
            repeat = 1
            while k > 0 and k < len(setup[:,0]):
                if np.abs(setup[i,0] - setup[k,0]) < 0.005:
                    k = k + 1
                    repeat = repeat + 1
                else:
                    k = -1
    
#    print repeat
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
                    base_fileout = 'RIXS_data_', clean = 0, curvature = [],
                    cal = [], scan_time = 10, pol_corr = []):

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
    """
    plt.figure()  
    
    for i in range (len(epoints[:])):
        
        fname_list = []
        fname_list_BG = []
        time_BG = []
        for j in range (len(setup[:,0])):
            if setup[j,3] == 1.0 and np.abs(setup[j,0] - epoints[i]) < 0.005:
                fname_list.append(fname[j])
                fname_list_BG.append(fname[int(setup[j,2])])
                time_BG.append(scan_time*(setup[j,1]-1))
            
        ims = CCD(fname_list, photon_E=epoints[i], poly_order=2, binpix=2,
                  fname_list_BG = fname_list_BG, exclude=lim)

        if len(curvature[:]) == 3:
            ims.curvature = curvature
                  
        ims.clean(clean) #ims.clean_std(clean) ???? 
        
        if len(pol_corr) != 0:        
            ims.get_specs()
    
            w = 0
            for spec in ims.BGspecs:
                sum_bkg = 0
                for j in range (len(spec[1][:])):
                    sum_bkg = sum_bkg + spec[1][j]
                pol = pol_corr[0] + pol_corr[1]*time_BG[w] + pol_corr[2]*time_BG[w]**2
                #print pol, sum_bkg
                ims.BGimages[w] = ims.BGimages[w]*pol/sum_bkg
                w = w + 1
            
        ims.sub_backgrounds()
        ims.get_specs()
        ims.sum_specs()
        
        if len(cal[:]) == 2:
            ims.calibrate(cal[0],cal[1])
            
        if dx <= 0.0:
            dx = ims.spectrum[0][0]-ims.spectrum[0][1]
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
            
           
            
def fit_pol(x,p,mode='eval'):
   if mode == 'eval':
      out = p[0] + p[1]*x + p[2]*x**2
   elif mode == 'params':
      out = ['offset','first order', 'second order']
   elif mode == 'name':
      out = "Exponential"
   elif mode == 'guess':
      g = peakguess(x, p)
      out = g[4:6]
   else:
      out = []
   return out

def calc_bkg_corr(setup,fname, scan_time = 10, lim = [1,9999], clean = 1E5,
                  curvature = []):

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
            
            ims.clean(clean)
            
            ims.get_specs()
            ims.sum_specs()
            
            sum_bkg.append(0.0)
            time.append(scan_time*i*1.0)
            for w in range (len(ims.spectrum[1])):
                sum_bkg[j] = sum_bkg[j] + ims.spectrum[1][w]
                
            j = j + 1


    err = []
    err[:] = np.sqrt(sum_bkg[:])

    M = np.zeros((len(time),3))
    
    M[:,0] = time[:]
    M[:,1] = sum_bkg[:]
    M[:,2] = err[:]    

    funcs = [fit_pol]
    
    p0 = sum_bkg[0]
    p1 = (sum_bkg[len(sum_bkg)-1] - sum_bkg[0])/(time[len(time)-1]-time[0])
    p2 = -p1
    
    guess = array([p0,p1,p2])
    
    f = fit.fit(x= M[:,0], y=M[:,1], e = M[:,2], funcs=funcs, guess=guess, ifix= [0,0,0])  
    f.go()
    clf()
    errorbar(f.datax, f.datay, f.datae, fmt='k.')
    fx,fy = f.evalfitfunc(nxpts = 200)
    plt.plot(fx, fy, 'r-')       
    
    plt.plot(M[:,0],M[:,1], 'bs')
    plt.xlabel('Time')
    plt.ylabel('Bkg')
    time_lim = 0.1*np.abs(time[len(time)-1]-time[0])
    bkg_lim = 0.1*(np.max(sum_bkg) - np.min(sum_bkg))
    plt.xlim(time[0]-time_lim,time[len(time)-1]+time_lim)
    plt.ylim(np.min(sum_bkg) - bkg_lim, np.max(sum_bkg) + bkg_lim)

    return f.result          
        
        
        