import sys
sys.path.append('../')
from corrfunc import compute_corr
import matplotlib.pyplot as plt
import numpy as np
from processing import generate_ref_maps, generate_dust_maps
import healpy as hp
import time

start = time.perf_counter()

for nside in [128]:
    
    #Fixing parameters
    resol = np.degrees(hp.nside2resol(nside)) #degrees 
    dr = 2*resol
    linear = False

    #retrieving maps
    ref_dict = generate_ref_maps(nside)
    dust_map = generate_dust_maps(nside)[0]

    #fixing plotting parameters
    plt.figure(figsize=(10,8))
    plt.xlabel('Degrees (arcmin)')
    plt.ylabel('Correlation')
    colors = ['red', 'pink', 'cyan', 'blue']

    if linear :
        x = np.arange(resol, 3.5, dr) #linear scale
    else:
        x = np.exp(np.linspace(np.log(resol/1.2), np.log(2.5), np.log2(nside)+1)) #log scale 
        thicknesses = np.append((x[1:]-x[:-1]), (x[-1]-x[-2]))


    #computing the correlation values for all reference maps against 1998 dust map
    for zrange,color in zip(ref_dict,colors):
            out = [compute_corr(ref_dict[zrange], dust_map, r, thick) for r, thick in zip(x,thicknesses)]
            y,err = zip(*out)
            plt.errorbar(x*60,list(y),yerr=list(err), color=color, label=zrange)

    
    plt.legend()
    plt.title('Correlation plot for dust maps with NSIDE='+str(nside))

    if linear:
        plt.xscale('linear')
    else:
        plt.xscale('log')

    filename = '../temp/corr_dust_maps_'+str(nside)+'.%d.jpg' % time.time()
    plt.savefig(filename) 
    plt.clf()

end = time.perf_counter()

print('Finished in '+str(end-start)+' second(s)')
