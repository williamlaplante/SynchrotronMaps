import numpy as np
import pandas as pd
from astropy.io import fits
from pygsm import GlobalSkyModel2016
import PlotFunctions

def getLSSdata(num, ccformat=False, div=100, zmin=0, zmax=100):
    '''
    gets the data from a fits file in the data/num directory, where 
    num is a number from 0 to 49. One can specifiy the wanted z-range.
    There is also the option of getting the data in a 1d array, formatted
    properly for cross-correlation.
    
    '''
    if not isinstance(num,int):
        raise Exception('num should be an integer')
    elif num<0 or num>49:
        raise Exception('num variable should be between 0 and 49 inclusively')
    else:
        path = '/Users/williiamlaplante/Research/SynchrotronMaps/data/' + str(num) 
        filename = '/post_catalog.dr72all'+str(num)+'.fits'
        hdul = fits.open(path+filename)
        data = hdul[1].data
        hdul.close()
        z = data['Z']        
    
        if ccformat:
            xmin,xmax = (100,270)
            ymin,ymax = (-10,75)
            ra = data['RA']
            dec = data['DEC']
            arr = np.flip(np.histogram2d(data['RA'][np.where((zmin<z)&(z<zmax))],
                                         data['DEC'][np.where((zmin<z)&(z<zmax))],
                                         bins=div,
                                         range=[[xmin, xmax], [ymin, ymax]])[0], axis=1)

            llsdata = np.ndarray.flatten(np.flip(np.flip(arr,axis=1)),'F')
            return llsdata - np.mean(llsdata)
        else:
            data = data[np.where((zmin<z)&(z<zmax))]
            
    return data



def getinfo(num):
    '''
    Gets the info on a fits file, i.e. extracts its header. 
    '''
    if not isinstance(num,int):
        raise Exception('num should be an integer')
    elif num<0 or num>49:
        raise Exception('num variable should be between 0 and 49 inclusively')
        return []
    else:
        path = '/Users/williiamlaplante/Research/SynchrotronMaps/data/' + str(num) 
        filename = '/post_catalog.dr72all'+str(num)+'.fits'
        hdul = fits.open(path+filename)
        info = hdul[0].header
        hdul.close()
    return info



def gettemp(log=True, ccformat=False, freq = 500, div=100):
    
    '''
    Gets the temperatures of the global sky model.
    
    Returns an array of temperatures. Note that NSIDE=1024 in the model.
    
    There is the option of getting a 1d array of the temperature data only
    in the wanted RA, DEC range, that is RA = (100,270), DEC = (-10,75). The
    returned array gives the temperature of the 2d grid from right to left,
    down to up.

    '''
    
    gsm = GlobalSkyModel2016(freq_unit='MHz')

    if log:
        temp = np.log(gsm.generate(freq)) #Generate the temperatures
    else:
        temp = gsm.generate(freq)
    
    if ccformat:
        pixel_range = PlotFunctions.getindex(div)
        temp = temp[pixel_range]
        temp = temp - np.mean(temp)
    
    return temp


def makesomenoise(mean=0, std=0.1, div=100):
    '''
    Makes a noisy array of size divxdiv.
    '''
    
    return np.ndarray.flatten(np.random.normal(loc=mean, scale=std, size=(div,div)))


def getzranges():
    
    '''
    Returns the zranges for each data file.
    
    '''
    
    z = []
    for i in range(50):
        a = getinfo(i)
        z.append((float(a['ZMIN']),float(a['ZMAX'])))
    return z


def combinedata(num_zvalues=10, div=100):
    '''
    Combines the data of all the data files in ccformat, ordered in zrange values. That is, for a
    specific range of z values, all data files have some data in that range. We extract that data and
    put it into the same rows, and each columns are associated with a data file.
    
    '''
    dframe = pd.DataFrame([])
    zmin, zmax = (0.0010009038, 0.4969387) #hardcoded zmin and zmax for speed
    zrange = np.linspace(zmin,zmax,num_zvalues)
    z = ['zmid' + str(i) for i in range(1,len(zrange))]
    indexvals = np.arange(div**2)
    col = np.array([])
    for i in range(50):
        for j in range(len(zrange)-1):
            col = np.append(col, getLSSdata(i, ccformat=True, div=100, zmin=zrange[j], zmax=zrange[j+1]))
        
        dframe[str(i)] = pd.Series(col)
        col = np.array([])
    
    index = pd.MultiIndex.from_product([z,indexvals], names=['zrange', 'data'])
    dframe.index = index
    return dframe.fillna(0)
    




