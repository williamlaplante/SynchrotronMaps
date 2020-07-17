import numpy as np
import healpy as hp

def query_annulus(m, nside, vec, R, dr, nest=False, idx=False):
    '''
    queries an annulus centered at vec with radius R < r < R + dr on a map m with resolution nside
    and ordering nest (true for nest, false for ring).

    Parameters:
    ---------------
    m : The input healpix map (if idx=False, the map is not considered)
    nside : The resolution of the map.
    vec : A vector [x,y,z] representing the center of the annulus.
    R : The radius of the inner circle of the annulus in degrees.
    dr : The thickness of the annulus in degrees.
    nest : Ordering of the map. If true, nested, else ring.
    idx : When true, returns the index of the queried pixels. Else, returns queried pixels from map. 

    Notes:
    -  When idx=True, the parameters taken into account are not extracted from the map's properties.
       For example, if there is an input map while idx=True, nside won't be determined with
       hp.get_nside(m).

    - If the map contains hp.UNSEEN pixels, this function ignores them. This is a good way to avoid
      querying some unwanted pixels.

    '''
    
    disk1 = hp.query_disc(nside, vec, np.radians(R+dr), nest=nest)
    disk2 = hp.query_disc(nside, vec, np.radians(R), nest=nest)
    index = np.setdiff1d(disk1, disk2, assume_unique=True)
    if index.size == 0:
        raise Exception('Empty annulus. Please increase thickness of annulus or map resolution.')
    if idx:
        return index
    else:
        circle = m[index]
        return circle[circle!=hp.UNSEEN]





def get_annuli_len(nside, R, dr, step=1, nest=False):
    '''
    Gets the average size of the returned array when calling query_annulus(). 
    
    Parameters:
    ------------
    nside, R, dr, nest : Identical to the parameters in the query_annulus() function.
    step : The number of pixels to skip when hopping from one center point to another.
    '''
    npix = hp.nside2npix(nside)
    X,Y,Z = hp.pix2vec(nside, np.arange(0,npix,step), nest=nest)
    l = 0
    for vec in zip(X,Y,Z):
        l+=hp.query_disc(nside, vec, np.radians(R+dr),nest=nest).size - hp.query_disc(nside,vec,np.radians(R),nest=nest).size

    return l/X.size



def compute_corr(map1, map2, R, dr, nest=False):
    '''
    Computes two-point or auto correlation function w(|r2-r1|) = w(R).  

    Parameters:
    ------------
    map1: first healpix map
    
    map2: second healpix map
    
    theta: distance |r1-r2| in degrees
    
    dr: thickness of annuli in degrees
    
    nest: ordering of the healpix maps
    
    '''
    if hp.get_nside(map1)!=hp.get_nside(map2):
        raise Exception('Maps are of two different resolutions.')
        
    #parameters
    nside = hp.get_nside(map1)
    npix = hp.nside2npix(nside) 

    if dr < np.degrees(hp.nside2resol(nside)):
        dr = np.degrees(hp.nside2resol(nside))


    if R==0: 
        return (map1*map2).mean(), (map1*map2).std()/np.sqrt(npix)
    else: 
        ring_len = get_annuli_len(nside, R, dr, step=4, nest=nest)

    
    #Get the desired pixels and the vectors associated with them.
    pix1, pix2 = np.arange(npix)[map1!=hp.UNSEEN], np.arange(npix)[map2!=hp.UNSEEN]
    pix_inter =  np.intersect1d(pix1,pix2)
    X,Y,Z = hp.pix2vec(nside, pix_inter, nest=nest)
    

    it = (query_annulus(map2, nside, vec, R, dr, nest).sum()/ring_len for vec in zip(X,Y,Z))
    
    result = map1[pix_inter]* np.fromiter(it, dtype=np.float)
    
    return result.mean(), result.std()/np.sqrt(npix)
