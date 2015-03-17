import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.interpolate

plt.ion()

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

def getM(scan_list, base_name, eloss, ein, IDPolar_val):
    epoints = []
    for scanno in scan_list:
        filename = '{}{:04d}.dat'.format(base_name, scanno)
        if getParam(filename, 'idpolar') == IDPolar_val:
            epoints.append(getParam(filename, 'monoenergy'))
    
    epoints = np.array(epoints)
    print epoints
    M_at_epoints = np.zeros((len(eloss),len(epoints)))
    
    i = 0
    for scanno in scan_list:
        filename = '{}{:04d}.dat'.format(base_name, scanno)
        if getParam(filename, 'idpolar') == IDPolar_val:
            print filename
            data = np.loadtxt(filename)
            plt.plot(-1*data[:,0], data[:,1], label=filename[-8:-4])
            f = scipy.interpolate.interp1d(-1*data[:,0], data[:,1], kind='linear', bounds_error=True, fill_value=np.nan) # last two added MPMD
            M_at_epoints[:,i] = f(eloss)
            i += 1
    
    f = scipy.interpolate.interp2d(epoints, eloss, M_at_epoints, kind='linear', bounds_error=False, fill_value=np.nan)
    M = f(ein, eloss)
    return M
    #return M, M_at_epoints, epoints

# vectors to interpolate on
eloss = np.linspace(-1, 8, 900)
ein = np.arange(851, 857.3, 0.1)

#  LTNAO
#scan_list = range(412,467+1)
scan_list = range(430, 435)

base_name = '../RIXS/LTNAO_L3_l2/LTNAO_'


IDPolar_val = 0.0 # LH
print "############ RA_LH ###############"
plt.clf()
M_LTNAO_LH = getM(scan_list, base_name, eloss, ein, IDPolar_val)
plt.show()
print "############ RA_LV ###############"
IDPolar_val = 1.0 # LV
M_LTNAO_LV = getM(scan_list, base_name, eloss, ein, IDPolar_val)
plt.legend()
