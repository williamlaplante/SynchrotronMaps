'''
Exports correlation values and corresponding thetas in an npz file.
'''
import sys
sys.path.append('../../')
from corrfunc import compute_corr, C
import numpy as np
import pandas as pd
from processing import generate_ref_maps, generate_dust_maps
import concurrent.futures
import healpy as hp
import time

    
#Fixing parameters
nsides = [1024,2048]	
resol = np.degrees(hp.nside2resol(nsides[1])) #degrees 
x = np.exp(np.linspace(np.log(resol/1.2), np.log(2.5), np.log2(nsides[1])+1)) #log scale 
thicknesses = np.append((x[1:]-x[:-1]), (x[-1]-x[-2]))

#wrapper function for multiprocessing
def compute_corr_wrap(p):
	return compute_corr(*p)

for nside in nsides:
	#retrieving maps
	ref_dict = generate_ref_maps(nside)
	dust_map = generate_dust_maps(nside)[0]

	Y = pd.DataFrame()
	Y_err = pd.DataFrame()
	#computing the correlation values for all reference maps against 1998 dust map
	for zrange,i in zip(ref_dict.keys(), range(len(ref_dict.keys()))):
		with concurrent.futures.ProcessPoolExecutor() as executor:
			args = ((ref_dict[zrange], dust_map, r, thick) for r, thick in zip(x,thicknesses))
			out = list(executor.map(compute_corr_wrap, args))
			y,err = zip(*out)
			col_name_Y = 'y_values_z'+str(i)
			col_name_Y_err = 'y_error_z'+str(i)
			Y[col_name_Y] = list(y)
			Y_err[col_name_Y_err] = list(err)

			if i==0:
				C_l = [abs(np.random.normal())/(10*n) for n in range(1,21)]
				m = hp.synfast(C_l, nside)
				args =((m,m,r,thick) for r,thick in zip(x,thicknesses))
				out = list(executor.map(compute_corr_wrap, args))
				y,err = zip(*out)
				y_values_gen_legendre = [C(r, C_l) for r in x]
				y_values_gen_hc = list(y)
				y_error_gen_hc = list(err)

	
	filename = './out/arr_'+str(nside)+'.%d.npz' % time.time()
	np.savez(filename,x_values=x, y_values_gen_legendre = y_values_gen_legendre, y_values_gen_hc = y_values_gen_hc, y_error_gen_hc = y_error_gen_hc, **{col_name: np.array(col_vals) for col_name, col_vals in Y.iteritems()}, **{col_name: np.array(col_vals) for col_name, col_vals in Y_err.iteritems()})

'''
The arrays are composed as follows (not following same ordering, but same labeling):

for nside in [1024,2048]:
	arr_compute_corr_nside= [x_values, y_values_z0, y_error_z0,
							y_values_z1, y_error_z1,
							y_values_z2, y_error_z2,
							y_values_z3, y_error_z3,
							y_values_gen_hc, y_error_gen_hc, y_values_gen_legendre]


x_values : array of theta values constructed from nside=2048. always kept constant.

y_values_z(i) : array of correlation values associated with the x_values for the z(i)th redshift range.
				The redshift ranges are z0 : 0.1<z<0.2, z1 : 0.3<z<0.4, z2 : 0.5<z<0.6, z3 : 1.2<z<1.3.

y_values_gen_hc/y_error_gen_hc : correlation values associated with x_values obtained by pre-setting Cl's, generating a map
				  from the Cl's, and using compute_corr to get correlation values and error.

y_values_gen_legendre : correlation values obtained using legendre polynomials for theta in corresponding x_values.
'''

