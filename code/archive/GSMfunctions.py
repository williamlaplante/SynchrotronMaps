#Functions related to the Global Sky Model

import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from astropy.coordinates import SkyCoord
from astropy import units as u
import itertools
from astropy.io import fits
from pygsm import GlobalSkyModel2016

#The region in the sky we will be looking at
RArange = [100,270]
DECrange = [-10,75]

#NSIDE parameter should be kept at 1024 for GSM maps
NSIDE = 1024

#Data related functions

def gettemp(log=True, ccformat=False, freq = 500, div=100):
    
    '''
    Gets the temperatures of the global sky model.
    
    Returns an array of temperatures. Note that NSIDE=1024 in the model.
    
    There is the option of getting a 1d array in cross-correlation format(ccformat).
    The format of this array is a 2d array of the temperatures, which are exactly 
    the same temperatures illustrated by the plotting function 'plotGSMregion'. That is,
    it is exactly the 2d array yielding the colormap of the figure.

    '''
    
    gsm = GlobalSkyModel2016(freq_unit='MHz')

    if log:
        temp = np.log(gsm.generate(freq)) #Generate the temperatures
    else:
        temp = gsm.generate(freq)
    
    
    if ccformat:
        pixel_range = getindex(div)
        temp = temp[pixel_range]
        return np.flip(np.reshape(temp, (div,div), 'F'),axis=0)
    
    return temp


def cctemp(div, log=False, freq=500, centered=False):
    '''
    returns the temperature data in the form of a matrix. This function 
    is to combine transformations made on the data sample, such as 
    normalization and centering.
    
    '''
    tempdata = gettemp(log=log, ccformat=True,freq=freq,div=div)
    
    if centered:
        tempdata = tempdata - tempdata.mean()
        
    return tempdata


def _getrange(div=100):
    '''
    Returns a range of ra and dec values.
    '''
    ra_values = np.linspace(RArange[0],RArange[1],div)
    dec_values = np.linspace(DECrange[0],DECrange[1],div)
    return ra_values,dec_values



def _getgrid(div=100):
    '''
    Returns the cartesian product of the ra and dec values
    '''
    ra_values,dec_values = _getrange(div)
    return map(list,zip(*itertools.product(ra_values,dec_values)))



def _convertgrid2Gal(div=100): 
    #Convert each set of points in the grid to galactic coordinates
    '''
    Converts each grid point (ra,dec) to galactic coordinates.
    
    Note: This part takes some time to run, perhaps later you could try 
    other solutions to improve it.
    '''
    ra_values, dec_values = _getgrid(div)
    return SkyCoord(ra=ra_values, dec=dec_values, frame='fk5', unit='deg').galactic


def getindex(div=100, nside=1024): 
    '''
    Converts each galactic coordinate to a pixel number, i.e. we retrieve
    the index of the wanted temperatures.
    Useful to use as an index to retrieve the temperatures sought 
    for data analysis or other stuff.
    '''
    coord = _convertgrid2Gal(div)
    return hp.ang2pix(nside, coord.l.value, coord.b.value, lonlat=True)


#Plotting functions

def plotGSMregion(div=100, pixel_size=20,log=False, freq=500):
    temp = gettemp(log=log,ccformat=False,freq=freq,div=div)
    pixel_range = getindex(div)
    x,y = _getgrid(div)
    colormap = temp[pixel_range]
    plt.figure(figsize=(10,8))
    plt.scatter(x,y,c=colormap,s=pixel_size,marker='s')
    plt.colorbar()
    plt.xlabel('RA')
    plt.ylabel('DEC')
    return


def illustrate():
    NSIDE = 1024
    temp = DataFunctions.gettemp()
    center = SkyCoord(ra=190*u.deg, dec=30*u.deg).galactic
    vec = hp.ang2vec(center.l.value,center.b.value,lonlat=True)
    ipix_disc = hp.query_disc(nside=NSIDE, vec=vec, radius=np.radians(10))
    temp[ipix_disc] = temp.max()
    hp.mollview(temp)
    return 
