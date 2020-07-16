from corrfunc import compute_corr
import matplotlib.pyplot as plt
import numpy as np
from processing import generate_ref_maps, generate_dust_maps
import healpy as hp

#Fixing parameters
nside = 1024
resol = np.degrees(hp.nside2resol(nside)) #degrees
dr = 4*resol

#retrieving maps
ref_dict = generate_ref_maps(nside)
dust_map = generate_dust_maps(nside)[0]

#fixing plotting parameters
plt.figure(figsize=(10,8))
plt.xlabel('Degrees (arcmin)')
plt.ylabel('Correlation')
colors = ['red', 'pink', 'cyan', 'blue']
x = np.arange(resol, 3.5, dr)

#computing the correlation values for all reference maps against 1998 dust map
for zrange,color in zip(ref_dict,colors):
    out = [compute_corr(ref_dict[zrange], dust_map, r, dr) for r in x]
    y,err = zip(*out)
    plt.errorbar(x*60,list(y),yerr=list(err), color=color, label=zrange)
    
    
plt.legend()
plt.title('Correlation plot for dust maps with NSIDE='+str(nside))

plt.savefig('temp/corr_dust_maps_1024.jpg')
