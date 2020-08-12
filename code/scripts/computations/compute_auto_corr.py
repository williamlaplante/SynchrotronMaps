'''
Exports correlation values and corresponding thetas in an npz file.
'''
import sys
sys.path.append('../../')
from corrfunc import compute_corr, C
import numpy as np
import pandas as pd
from processing import read_ref_map, read_dust_map
import concurrent.futures
import healpy as hp
import time

    
#Fixing parameters
nsides = [128,256,512,1024]

resol = np.degrees(hp.nside2resol(128)) #degrees 
x = np.arange(resol,6,resol) #degrees
thicknesses = np.append((x[1:]-x[:-1]), (x[-1]-x[-2]))

#wrapper function for multiprocessing
def compute_corr_wrap(p):
    return compute_corr(*p)


for zmin, zmax in zip([0.1,0.3,0.5,1.2], [0.2,0.4,0.6,1.3]):
    col_name = str(zmin) +'z'+ str(zmax)
    
    for nside in nsides:
        ref_map = read_ref_map(nside, zmin, zmax)
        with concurrent.futures.ProcessPoolExecutor() as executor:
            args = ((ref_map, ref_map, R, dr) for R, dr in zip(x, thicknesses))
            out = list(executor.map(compute_corr_wrap, args))
            y,err = zip(*out)
            filename = "./out/auto_corr/" + str(zmin) +'z'+ str(zmax) + '.' + str(nside) + '.%d.npz' %time.time()
            np.savez(filename, x=x, y=y, err=err)
            

