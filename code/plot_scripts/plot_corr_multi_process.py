import sys
sys.path.append('../')
from corrfunc import compute_corr
import matplotlib.pyplot as plt
import numpy as np
from processing import generate_ref_maps, generate_dust_maps
import concurrent.futures
import healpy as hp
import time

start = time.perf_counter()
    
#Fixing parameters
nside = 1024	
resol = np.degrees(hp.nside2resol(nside)) #degrees 

#retrieving maps
ref_dict = generate_ref_maps(nside)
dust_map = generate_dust_maps(nside)[0]

#fixing plotting parameters
plt.figure(figsize=(10,8))
plt.xlabel('Degrees (arcmin)')
plt.ylabel('Correlation')
colors = ['red', 'pink', 'cyan', 'blue']

x = np.exp(np.linspace(np.log(resol/1.2), np.log(2.5), np.log2(nside)+1)) #log scale 
thicknesses = np.append((x[1:]-x[:-1]), (x[-1]-x[-2]))

def compute_corr_wrap(p):
	return compute_corr(*p)

#computing the correlation values for all reference maps against 1998 dust map
for zrange,color in zip(ref_dict.keys(), colors):
	with concurrent.futures.ProcessPoolExecutor() as executor:
	    args = ((ref_dict[zrange], dust_map, r, thick) for r, thick in zip(x,thicknesses))
	    out = list(executor.map(compute_corr_wrap, args))
	    y,err = zip(*out)
	    plt.errorbar(x*60, list(y), yerr=list(err), color=color, label=zrange)


plt.legend()
plt.title('Correlation plot for dust maps with NSIDE='+str(nside))
plt.xscale('log')

filename = '../temp/corr_dust_maps_'+str(nside)+'.%d.jpg' % time.time()
plt.savefig(filename) 

end = time.perf_counter()

print('Finished in '+str(end-start)+' second(s)')
