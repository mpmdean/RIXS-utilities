import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import processCCD_image

"""
Assumes that it does 1 LH then 1 LV
"""

def getParam(filename, param):
    fid = open(filename)
    alltxt = fid.read().lower()  # everything in lowercase
    param = param.lower()
    for line in alltxt.splitlines():
        if line.find(param) != -1:
            variable, valuestr = line.split(' = ')
    try:
        return float(valuestr)
    except NameError:
        print "{} not found in file {}".format(param, filename)


base_name = '../RIXS/Ni_LTNAO_'
scan_list = []
for scan_LH, scan_LV in zip(range(589, 603, 2), range(590,603, 2)): 
    scan_list.append(scan_LV)
    scan_list.append(scan_LH)

clean = [10,6,4]
previous_Ein = 852.8
elastic_guess = 550 - 189.1
curvature = [2.01257084e-06, -1.09139974e-02, 3.70021820e+02]

filebkg = '{}{:04d}'.format(base_name, 474)

plt.figure()
for scanno in scan_list:
    print '\n'
    filename = '{}{:04d}'.format(base_name, scanno)
    try:
        ims = processCCD_image.CCD([filename], photon_E=851, poly_order=2, binpix=1,
                  fname_list_BG = [filebkg], exclude=[1,9999])
        print "found {}".format(filename)
    except IOError:
        print "Didn't find {}".format(filename)
        continue
        
    ims.curvature = curvature
    ims.clean_std_new(clean)
    ims.sub_backgrounds()
    ims.get_specs()
    ims.sum_specs()
    
    Ein = ims.file_dictionary['beamline_monoenergy']
    print "Energy changed by {:.1f}".format(Ein-previous_Ein)
    elastic_guess = elastic_guess - (Ein-previous_Ein)/0.0771
    
    if ims.file_dictionary['beamline_idpolar'] == 1.0:
        elastic_list = ims.fit_elastic(cen=None, sigma=2, x_window = [elastic_guess-10, elastic_guess+10])
        if np.abs(elastic_list - elastic_guess) < 5.0:
            print "LV Sucessful fit: guess = {:.1f}, value = {:.1f}".format(elastic_guess, elastic_list)
            elastic = elastic_list
        else:
            print "LV FIT FAILED!!!!!!. Using guess {}".format(elastic_guess)
            elastic = elastic_guess
    elif ims.file_dictionary['beamline_idpolar'] == 0.0:
        print "LH using guess {}".format(elastic_guess)
        elastic = elastic_guess

    ims.calibrate(elastic, 0.0767)
    ims.plot_spectrum(label=filename[-4:])
    
    previous_Ein = Ein
    
    ims.fileout = "../RIXS/LTNAO_L3_L2/LTNAO_{:04d}.dat".format(scanno)
    ims.write_file()

plt.legend()
plt.show()