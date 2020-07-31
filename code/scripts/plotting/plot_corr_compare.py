import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import time

t1 = time.perf_counter()
print('Program starts')
print('--------------')

out_1024 = np.load('../computations/out/arr_1024.1596072265.npz')
out_2048 = np.load('../computations/out/arr_2048.1596152119.npz')

print('Setting up plotting parameters...')
figure, axes = plt.subplots(nrows=2, ncols=4, figsize=(20,8), sharex=True)
x = out_2048['x_values']*60 #out_1024 has the same x_values

print('Plotting...')
for j in range(4):
    axes[0,j].errorbar(x, out_1024['y_values_z'+str(j)], yerr=out_1024['y_error_z'+str(j)], label='nside=1024')
    axes[0,j].errorbar(x, out_2048['y_values_z'+str(j)], yerr=out_2048['y_error_z'+str(j)], label='nside=2048')
    axes[0,j].set_title('z%d' %j)
    axes[0,j].set_xscale('log')
    axes[0,j].set_xlabel('Degrees (arcmin)')
    axes[0,j].legend()

axes[1,0].errorbar(x, out_1024['y_values_gen_hc'], yerr=out_1024['y_error_gen_hc'],linestyle='--' , marker='o', label='nside=1024 ; bruteforce')
axes[1,0].plot(x, out_1024['y_values_gen_legendre']/9.8, label='legendre (times factor of 1/9.8)')
axes[1,0].legend()
axes[1,0].set_xlabel('Degrees (arcmin)')
axes[1,1].errorbar(x, out_2048['y_values_gen_hc'], yerr=out_2048['y_error_gen_hc'],linestyle='--', marker='o', label='nside=2048 ; bruteforce')
axes[1,1].plot(x, out_2048['y_values_gen_legendre']/9.35, label='legendre (times factor of 1/9.35)')
axes[1,1].legend()
axes[1,1].set_xlabel('Degrees (arcmin)')

for i in range(2) : axes[i,0].set_ylabel('Correlation')

figure.tight_layout()
figure.suptitle('Correlation of dust maps with LSS for resolutions nside=1024, 2048 for different z-ranges')

print('Saving plot...')
figure.savefig('../results/comp_corr_by_z.%d.jpg' % time.time(), format='jpeg')

t2 = time.perf_counter()

print('Program finished in %f second(s).' % (t2-t1))