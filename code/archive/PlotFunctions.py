import DataFunctions
import StatsFunctions
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from astropy.coordinates import SkyCoord
from astropy import units as u
import itertools

#The region in the sky we will be looking at
RArange = [100,270]
DECrange = [-10,75]


def plotLSS(sample_type,zmin=0,zmax=2):
    '''
    Plots the large scale structure data, that is galaxy position in a RA / DEC grid,
    but only for a specific range of z values.
    
    '''
    alldata = DataFunctions.getLSSdata(sample_type)
    data = alldata[np.where((zmin<=alldata['Z'])&(alldata['Z']<=zmax))]
    plt.figure(figsize=(10,10))
    plt.plot(data['RA'],data['DEC'],'.',markersize=0.5,color='black')
    plt.xlabel('R.A.')
    plt.ylabel('D.E.C.')
    plt.xlim(RArange[1], RArange[0])
    plt.ylim(DECrange[0], DECrange[1])
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



def getindex(div=100): 
    '''
    Converts each galactic coordinate to a pixel number, i.e. we retrieve
    the index of the wanted temperatures.
    Useful to use as an index to retrieve the temperatures sought 
    for data analysis or other stuff.
    '''
    coord = _convertgrid2Gal(div)
    NSIDE = 1024
    return hp.ang2pix(NSIDE, coord.l.value, coord.b.value, lonlat=True)



def plotGSMregion(div=100, pixel_size=20, ra_minmax = [100,270], dec_minmax = [-10,75]):
    temp = DataFunctions.gettemp(False)
    pixel_range = getindex(div)
    x,y = _getgrid(div)
    colormap = temp[pixel_range]
    plt.figure(figsize=(10,8))
    plt.scatter(x,y,c=colormap,s=pixel_size,marker='s')
    plt.xlim(270,100)
    plt.colorbar()
    plt.xlabel('RA')
    plt.ylabel('DEC')
    return



def _getangles():
    '''
    Retrieves the angle of each pixel in the mollweide projection for NSIDE=1024 and transforms them into fk5 coordinates.
    
    Returns a SkyCoord object with RA and DEC angles of each pixel in the mollweide projection.
    '''
    NSIDE=1024
    temp = DataFunctions.gettemp()
    index = np.arange(len(temp))
    gal_angles = hp.pix2ang(NSIDE,index,lonlat=True)
    return SkyCoord(l=gal_angles[0]*u.degree, b=gal_angles[1]*u.degree, frame='galactic').fk5



def _retrieveindex():
    '''
    Returns an array of indices, where each index is associated with an (RA,DEC) coordinate
    respecting predefined boundaries. Note that each each (RA,DEC) is associated with a temperature.
    
    '''
    angles = _getangles()
    ra_index = np.where((angles.ra.value>100)&(angles.ra.value<270))
    dec_index = np.where((angles.dec.value>-10)&(angles.dec.value<70))
    return np.intersect1d(ra_index,dec_index)



def plotGSMregion_2():
    '''
    Putting everything together, we can plot temperatures at a specific region in the sky in RA and DEC coordinates.
    '''
    temp = DataFunctions.gettemp(False)
    angles = _getangles()
    index = _retrieveindex()
    index2 = []
    for i in np.arange(0,len(index),5): #Skips each n-values. To make the code faster when doing trial runs.
        index2.append(index[i])
    
    plt.figure(figsize=(10,8))
    plt.scatter(angles.ra.value[index2], angles.dec.value[index2],s=1.5, c=temp[index2])
    plt.colorbar()
    plt.xlim(270,100)
    plt.xlabel('RA')
    plt.ylabel('DEC')
    return



def plotNz(num, numpoints=100):
    data = DataFunctions.getLSSdata(num)['RA']
    z = DataFunctions.getLSSdata(num)['Z']
    zval = np.linspace(z.min(), z.max(), numpoints)
    x = np.array([(zval[i+1]+zval[i])/2 for i in np.arange(len(zval)-1)])
    y = np.array([len(data[np.where((zval[i]<z) & (z<=zval[i+1]))]) for i in np.arange(0,len(zval)-1)])
    
    plt.plot(x,y,color='black')
    dz = (z.max()-z.min())/numpoints
    dzstring = str(float("{:.0e}".format(dz)))
    plt.ylabel(r'$\frac{dN}{dz} \times \Delta z$ = '+dzstring)
    plt.xlabel('z')
    return

def plotcorr(num, numpoints=15):
    '''
    plots correlation between temperature and galaxy positions
    at increasing intervals of z for the data in the directory data/num.

    '''
    temperatures = DataFunctions.gettemp(log=True, ccformat=True, freq=500, div=100)
    zvalues = DataFunctions.getLSSdata(num)['Z']
    
    zrange = np.linspace(zvalues.min(),zvalues.max(),numpoints)
    zmid = [(zrange[i]+zrange[i+1])/2 for i in np.arange(len(zrange)-1)]
    
    y = np.empty(len(zrange)-1)
    
    for i in range(len(zrange)-1):
        data = DataFunctions.getLSSdata(num=num, ccformat=True, div=100, zmin=zrange[i], zmax=zrange[i+1])
        y[i] = StatsFunctions.crosscorrelation(data,temperatures)
        
    plt.figure(figsize=(10,8))
    plt.plot(zmid, y,lw=0.4)
    ymin = [-0.0015]*(numpoints-1)
    ymax = [0.0015]*(numpoints-1)
    plt.fill_between(zmid,ymin,ymax, color='green', alpha=0.2)
    plt.xlabel('z')
    plt.ylabel('Correlation')
    
    return 


def plotall(numpoints=15):
    
    zmin, zmax = (0.0010009038, 0.4969387) #hardcoded zmin and zmax for speed
    zrange = np.linspace(zmin,zmax,numpoints)
    zmid = [(zrange[i]+zrange[i+1])/2 for i in np.arange(len(zrange)-1)]
    temperatures = DataFunctions.gettemp(log=True, ccformat=True, freq=500, div=100)

    plt.figure(figsize=(12,8))
    plt.xlabel('z')
    plt.ylabel('Cross-Correlation')
    
    for num in range(50):
        y = np.empty(len(zrange)-1)

        for i in range(len(zrange)-1):
            data = DataFunctions.getLSSdata(num=num, ccformat=True, div=100, zmin=zrange[i], zmax=zrange[i+1])
            y[i] = StatsFunctions.crosscorrelation(data,temperatures)
        plt.plot(zmid, y, lw=0.4)
     
    return
    