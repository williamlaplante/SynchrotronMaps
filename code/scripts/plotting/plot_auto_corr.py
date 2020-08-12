import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import time
import os

t1 = time.perf_counter()

print('Program starts')
print('--------------')

path = "../computations/out/auto_corr/"

print('Setting up plotting parameters...')

figure, axes = plt.subplots(nrows=2, ncols=4, figsize=(20,8), sharex=True)
figure.suptitle('Auto-correlation of LSS maps for different z-ranges.')
axes[0,0].set_ylabel('Correlation')


print('Plotting...')

for file in os.listdir(path):
    category = file[:3]
    resolution = file.split('.')[3]
    
    if category == "0.1":
        out = np.load(path+file)
        axes[0,0].errorbar(out['x'], out['y'], yerr = out['err'], label=resolution)
        axes[0,0].set_title('0.1 < z < 0.2')
        axes[0,0].set_xlabel('Angular separation (deg)')
        axes[0,0].legend()
        
    elif category == "0.3":
        out = np.load(path+file)
        axes[0,1].errorbar(out['x'], out['y'], yerr = out['err'], label=resolution)
        axes[0,1].set_title('0.3 < z < 0.4')
        axes[0,1].set_xlabel('Angular separation (deg)')
        axes[0,1].legend()

    elif category == "0.5":
        out = np.load(path+file)
        axes[0,2].errorbar(out['x'], out['y'], yerr = out['err'], label=resolution)
        axes[0,2].set_title('0.5 < z < 0.6')
        axes[0,2].set_xlabel('Angular separation (deg)')
        axes[0,2].legend()

    elif category == "1.2":
        out = np.load(path+file)
        axes[0,3].errorbar(out['x'], out['y'], yerr = out['err'], label=resolution)
        axes[0,3].set_title('1.2 < z < 1.3')
        axes[0,3].set_xlabel('Angular separation (deg)')
        axes[0,3].legend()

    else:
        raise Exception('Houston, we have a problem!')

figure.tight_layout()

print('Saving plot...')

figure.savefig('../results/auto_corr_by_z.%d.jpg' % time.time(), format='jpeg')

t2 = time.perf_counter()

print('Program finished in %f second(s).' % (t2-t1))
