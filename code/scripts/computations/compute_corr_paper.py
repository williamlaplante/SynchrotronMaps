'''
Exports correlation values and corresponding thetas in an npz file.
'''
import sys
sys.path.append('../../')
from corrfunc import compute_corr
import numpy as np
from processing import read_dust_map, read_ref_map
import concurrent.futures
import healpy as hp
import time

 t1 = time.perf_counter()
 
print('Fixing parameters...')
#Fixing parameters
nside = 4096
resol = np.degrees(hp.nside2resol(nside)) #degrees 
x = np.exp(np.linspace(np.log(resol/1.2), np.log(2.5), np.log2(nside)+1)) #log scale 
thicknesses = np.append((x[1:]-x[:-1]), (x[-1]-x[-2]))

#wrapper function for multiprocessing
def compute_corr_wrap(p):
	return compute_corr(*p)

print('Reading in dust map with resolution %d ...' %nside)
dust_map = read_dust_map(nside)

print('Entering zrange loop...')

for zmin, zmax in zip([0.1,0.3,0.5,1.2], [0.2,0.4,0.6,1.3]):
    ref_map = read_ref_map(nside, zmin, zmax)
    
    print('Starting multiprocessing...')

    with concurrent.futures.ProcessPoolExecutor() as executor:
        args = ((ref_map, dust_map, R, dr) for R, dr in zip(x, thicknesses))
        out = list(executor.map(compute_corr_wrap, args))

    print('Ending multiprocessing...')
    
	y,err = zip(*out)
	filename = "./out/cross_corr_paper/" + str(zmin) +'z'+ str(zmax) + '.' + str(nside) + '.%d.npz' %time.time()
	np.savez(filename, x=x, y=list(y), err=list(err))

t2 = time.perf_counter()

print('Program terminated in %f second(s)' %(t2-t1))
