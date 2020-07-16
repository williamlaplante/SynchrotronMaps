#Functions related to grid operations

import numpy as np
import itertools
import healpy as hp
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.visualization import ZScaleInterval, PercentileInterval

class Grid2d:
    
    def __init__(self, xrange, yrange, shape):
        
        if isinstance(xrange, list) and isinstance(yrange, list):
            
            if len(xrange)!=2 or len(yrange)!=2:
                raise Exception("Range length is not equal to 2.")
            else:
                self.xrange = xrange
                self.yrange = yrange
        else:
            raise TypeError('Input range should be of type tuple.')
        
        if isinstance(shape, tuple) and len(shape)==2 and all(np.array(shape)>0):
            self.shape = shape
        else:
            raise Exception('Please input a correct 2d shape.')
        
        self.frame = None
            
            
            
    def set_frame(self, frame):
        '''
        Sets the coordinate system of the grid.
        
        '''
        self.frame = frame
    
        return 
    
    def getrange(self):
        
        '''
        Returns a range of values in the x and y dimensions with equal spacing. 

        Returns:
        --------------
        Tuple of the form (x,y) where x and y are numpy arrays divided according to dimension size.
        '''
        return np.linspace(self.xrange[0],self.xrange[1],self.shape[0]), np.linspace(self.yrange[0],self.yrange[1],self.shape[1])



    def getgrid(self):
        '''
        Returns unzipped cartesian product of x,y values. That is, if xrange=[0,1]=yrange, then it returns the cartesian product [(0,0), (0,1), (1,0), (1,1)].
        '''
        xrange,yrange = self.getrange()
        return list(itertools.product(xrange,yrange))



    def transform(self, frame1, frame2, unit): 
        '''
        Converts grid points of the form (x,y) from coordinate system of frame1 to frame2.
        
        Parameters:
        ---------------
        units : units of the coordinate system created at instantiation of class. 
        frame1 : frame of the coordinate system created at instantiation of class. 
        frame2 : frame to which frame1 will be transformed.
        Returns:
        -------------
        A SkyCoord object with the of the transformed coordinates.
      
        '''
        frames = ['galactic', 'fk5', 'icrs']
        units  = ['deg', 'rad']
        
        if unit not in units: raise Expection('Units incorrect.')
        if frame1 not in frames or frame2 not in frames: raise Exception('Frames incorrect.')
            
        coord = self.getgrid()
        x, y = zip(*coord)

        return SkyCoord(list(x), list(y), frame=frame1, unit=unit).transform_to(frame2)


        
class HealpyObj():
    
    def __init__(self, data, nest, coord, lonlat=True):
        self.NSIDE = hp.npix2nside(len(data))
        self.nest = nest
        self.coord = coord
        self.data = data
        self.lonlat = lonlat
        self.scale = None
        
        
    def set_scale(self, scale):
        self.scale = scale
        
    
    def scale_data(self, data):
        if self.scale == None:
            return data
        
        if self.scale.lower()=='log':
            return np.log(data)
            
        elif self.scale.lower()=='zscale':
            zscale = ZScaleInterval()
            return zscale(data)
            
        elif self.scale.lower()=='percentile':
            interval = PercentileInterval(90)
            return interval(data)
        else:
            return data
            
        
    def getindex(self, Grid2d): 
        '''
        For each coordinate point on a 2d grid, it retrieves the index (pixel) of that coordinate on a healpix map.

        Parameters:
        -----------------
        Grid2d : Grid2d object

        '''
    
        if self.coord == 'G':
            coord = Grid2d.transform(Grid2d.frame, 'galactic', unit='deg')
            return hp.ang2pix(self.NSIDE, coord.l.value, coord.b.value,nest=self.nest, lonlat=self.lonlat) 
        else:
            return -1 #Need to implement here for non-galactic transformations.

    
    def mollview(self, title):
        '''
        Plots a mollview projection of a healpix map. 
        
        Parameters:
        ------------
        title: title of the plot
        '''
        hp.mollview(self.scale_data(self.data), nest=self.nest, cmap='viridis', coord=self.coord, title=title) 
        
        return
    
    
    def getregion(self, Grid2d, arr2d=False):
        '''
        For each coordinate pair on a 2d grid, it retrieves the corresponding value at that coordinate on a healpix map. If arr2d is true, returns it in the form of a 2d array matching the conventions of the original grid. 
        
        That is, if we take the cartesian product of [0,1] and [0,1], we will get [(0,0), (0,1), (1,0), (1,1)] and thus we reshape it so that the matrix follows the following convetion: origin at bottom left, x,y dimensions going left to right and down to up respectively.
        
        '''
        if arr2d:
            return np.flip(np.reshape(self.data[self.getindex(Grid2d)], Grid2d.shape, 'F'),axis=0)
        
        else:
            return self.data[self.getindex(Grid2d)]
        
        
    def plotimage(self, Grid2d):
        '''
        Plots the region of an healpix map on a 2d grid as an image. The region is specified with a Grid2d object.
        
        Note: This function is almost identical to plotregion(), except it uses a matrix as its source to plot the image. This is mainly to confirm that the matrix manipulations are correct in getregion() and not to deal with markersize stuff from plotregion().
        '''
        arr = self.getregion(Grid2d, arr2d=True)
        plt.figure(figsize=(8,8))
        plt.imshow(self.scale_data(arr), cmap='viridis')
        plt.axis('off')
        plt.colorbar()
        
        return 
    
    
    def plotregion(self, Grid2d, markersize): #Try to automate markersize!
        '''
        Plots a region of an healpix map on a 2d grid. The region is specified with a Grid2d object.
        
        Parameters:
        -------------
        Grid2d : Grid2d object specifying the necessary information for the grid
        
        markersize : size of each pixel 
        '''
        
        
        colormap = (self.data[self.getindex(Grid2d)]).copy()
        
        x,y = zip(*Grid2d.getgrid())
        plt.figure(figsize=(10,8))
        plt.scatter(x, y, c=self.scale_data(colormap), s=markersize, marker='s')
        plt.colorbar()
        plt.xlabel('x')
        plt.ylabel('y')
        
        return
    
    
    def illustrate(self, Grid2d):
        
        '''
        Shows on a mollweide projection where a grid is located. 
        
        '''
        
        
        '''
        BROKEN FUNCTION
        Need to work on this function. There is a problem for degree vs radian,
        and there is another problem for the coordinate system used (the center and radius are based on RA,DEC and not on galactic.
        
        '''
        data = (self.data).copy()
        
        centerx =  np.radians((max(Grid2d.xrange) + min(Grid2d.xrange))/2)
        centery =  np.radians((max(Grid2d.yrange) + min(Grid2d.yrange))/2)
        
        radius = min(max(Grid2d.xrange) - centerx, max(Grid2d.yrange)-centery)
                    
        center_gal = SkyCoord(centerx, centery,frame='fk5',unit='deg').galactic
        
        
        vec = hp.ang2vec(center_gal.l.value,center_gal.b.value,lonlat=self.lonlat) 
        
        ipix_disc = hp.query_disc(nside=self.NSIDE, vec=vec, radius=radius)
        data[ipix_disc] = max(data)
        
        hp.mollview(data, nest=self.nest, cmap='jet', coord=self.coord) 


        
        
        
class DataGrid(Grid2d):
    
    def __init__(self, xrange, yrange, shape, xdata, ydata):
        
        self.xdata = xdata
        self.ydata = ydata
        self.scale = None
        
        Grid2d.__init__(self, xrange, yrange, shape) 
    
    
    def scale_data(self, data):
   
        if self.scale == None:
            return data
        
        if self.scale.lower()=='log':
            return np.log(data)
            
        elif self.scale.lower()=='zscale':
            zscale = ZScaleInterval()
            return zscale(data)
            
        elif self.scale.lower()=='percentile':
            interval = PercentileInterval(90)
            return interval(data)
        else:
            return data
        
        
    def overdensity_field(self):
        '''
        Gives the overdensity field of the data as a 2d array in the convention followed by the usual x,y plots. That is, the origin is in the bottom left, with x and y dimensions having increasing values from left to right and down to up respectively.
        
        '''
       
        arr = np.histogram2d(self.ydata, self.xdata, bins=self.shape, 
                                     range=[self.yrange, self.xrange])[0]
        arr = np.flip(arr,axis=0)
                              
        return (arr - arr.mean())/arr.mean()
        
    def plot_image(self):
        '''
        Plot an image of the overdensity field.  
        
        '''
        
        arr = self.overdensity_field()
        
        plt.figure(figsize=(12,10))
        plt.imshow(self.scale_data(arr), cmap='viridis')
        plt.axis('off')
        plt.colorbar()
     
    def plot_scatter(self):
        '''
        Scatter plot of the x,y data.
        '''
        plt.figure(figsize=(12,10))
        plt.plot(self.xdata, self.ydata, '.', markersize=0.5,color='black')
        plt.xlim(self.xrange[0], self.xrange[1])
        plt.ylim(self.yrange[0], self.yrange[1])
        
        return
        
        
        
        
        