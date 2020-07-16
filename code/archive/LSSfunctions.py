#Functions related to the NYU-VAGC Large-Scale Structure data

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

#The region in the sky we will be looking at
RArange = [100,270]
DECrange = [-10,75]

#Data related functions

def getLSSdata(sample_type, ccformat=False, zmin=0, zmax=100, div=100):
    
    '''
    Gets the data from the bright, full, void or safe catalog. 
    
    There is the option of getting the data in a cross-correlation format matching the
    formatting of GSM's 'gettemp' function. This latter option returns a 2d density array. 
    One can also specify the range of z values wanted for the ccformat option.
    
    '''
    if sample_type.lower() not in ['bright','full','void','safe']:
        return
    else: 
        path = '/Users/williiamlaplante/Research/SynchrotronMaps/data/'+sample_type+'/'
        filename = 'post_catalog.dr72'+sample_type+'0.fits'
        hdul = fits.open(path+filename)
        data = hdul[1].data
        hdul.close()
        z = data['Z']
        index = np.where((zmin<z)&(z<zmax))
        
        if ccformat:
            RA = data['RA']
            DEC = data['DEC']
            arr = np.flip(np.histogram2d(DEC[index], RA[index],
                         bins=div,
                         range=[DECrange, RArange])[0],axis=0)
            return arr
        
        return data[index]
    
    
def cclss(sample_type, zmin, zmax, div, centered=False):
    '''
    Same as cctemp but for the lss samples.
    
    '''
    lssdata = getLSSdata(sample_type=sample_type, ccformat=True, zmin=zmin, zmax=zmax, div=div)
    
    if centered:
        lssdata = lssdata - lssdata.mean()
    
    return lssdata
    
# Plotting related functions

def plotLSS(sample_type,zmin=0,zmax=2):
    '''
    Plots the large scale structure data, that is galaxy position in a RA / DEC grid,
    but only for a specific range of z values.
    
    '''
    data = getLSSdata(sample_type, zmin=zmin, zmax=zmax)
    plt.figure(figsize=(10,10))
    plt.plot(data['RA'],data['DEC'],'.',markersize=0.5,color='black')
    plt.xlabel('R.A.')
    plt.ylabel('D.E.C.')
    plt.xlim(RArange[0], RArange[1])
    plt.ylim(DECrange[0], DECrange[1])
    return


def plotdNdz(sample_type, numpoints=100):

    z = getLSSdata(sample_type)['Z']
    zval = np.linspace(z.min(), z.max(), numpoints)
    
    x = np.array([(zval[i+1]+zval[i])/2 for i in np.arange(len(zval)-1)])
    y = np.array([((zval[i]<z) & (z<=zval[i+1])).sum() for i in np.arange(0,len(zval)-1)])
    
    plt.plot(x,y,color='black',label=sample_type)
    dz = (z.max()-z.min())/numpoints
    dzstring = str(float("{:.0e}".format(dz)))
    plt.ylabel(r'$\frac{dN}{dz} \times \Delta z$ = '+dzstring)
    plt.xlabel('z')
    plt.legend()
    
    return 
