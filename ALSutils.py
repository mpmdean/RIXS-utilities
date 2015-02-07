import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import os, pdb, glob

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

def get_scan_info(folder):
    """ Finds all fits files in a given fold and returns:
    file_info - array with pathnames and filenames
    ai_info - array with frames, bl_energy, spec_energy, andor
    """
    searchterm = os.path.join(folder, '*.fits')
    path_list = np.array(glob.glob(searchterm))
    
    file_list = []
    for path in path_list:
        file_list.append(os.path.split(path)[1])
    
    file_list = np.array(file_list)
    
    aifilename = folder + file_list[0][0:-10] + 'AI.txt'  
    ai_info = np.loadtxt(aifilename, usecols = (0,1,2,3), skiprows = 10, dtype=float)
    
    frames = []
    bl_energy = []
    spec_energy = []
    andor = []   # 0 for background 1 for data
    for filename in file_list:
        fileno  = int(filename[-10:-5])
        row_index = ai_info[:,0] == fileno
        frames.append(ai_info[row_index, 0][0])
        bl_energy.append(ai_info[row_index, 1][0])
        spec_energy.append(ai_info[row_index, 2][0])
        andor.append(ai_info[row_index, 3][0])
    
    file_info = np.vstack((path_list, file_list)).transpose()
    ai_info = np.vstack((frames, bl_energy, spec_energy, andor)).transpose()
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