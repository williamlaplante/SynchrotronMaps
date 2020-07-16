import healpy as hp
import numpy as np
import Statsfunctions as statsfunc
from astropy.io import fits

full_path_ref = '/Users/williiamlaplante/Research/SynchrotronMaps/data/ReferenceSamples/'
full_path_dust = '/Users/williiamlaplante/Research/SynchrotronMaps/data/DustMapsData/'


def _remove_local(m,nest=False, rad=1):
    '''
    Removes local averages from field. 
    
    Parameters:
    --------------
    m : healpy map on which we remove local averages
    
    nest : ordering of the map
    
    rad : radius of the local region
    
    '''
    nside = hp.get_nside(m)
    npix = hp.nside2npix(nside)
    X,Y,Z = hp.pix2vec(nside,np.arange(npix),nest=nest)
    
    return np.array([m[hp.vec2pix(nside,x,y,z, nest=nest)] - m[hp.query_disc(nside, [x,y,z], np.radians(rad))].mean() for x,y,z in zip(X,Y,Z)])


def generate_ref_maps(NSIDE):
    '''
    gets the full reference sample used in the dust map paper, and divides it into four
    sub-reference maps, which are divided by z intervals. Then it returns overdensity fields
    of these subsets of data in the form of healpix maps.
    '''
    refsample = fits.open(full_path_ref + 'full_ref_sample.fits')[1].data
    ra = refsample['RA']
    dec = refsample['DEC']
    z = refsample['Z']
    
    ra1, dec1 = ra[np.where((0.1<z) & (z<0.2))], dec[np.where((0.1<z) & (z<0.2))] 
    ra2, dec2 = ra[np.where((0.3<z) & (z<0.4))], dec[np.where((0.3<z) & (z<0.4))]
    ra3, dec3 = ra[np.where((0.5<z) & (z<0.6))], dec[np.where((0.5<z) & (z<0.6))]
    ra4, dec4 = ra[np.where((1.2<z) & (z<1.3))], dec[np.where((1.2<z) & (z<1.3))] 
    ref1 = statsfunc.coord2hp(ra1,dec1, NSIDE)
    ref2 = statsfunc.coord2hp(ra2,dec2, NSIDE)
    ref3 = statsfunc.coord2hp(ra3,dec3, NSIDE)
    ref4 = statsfunc.coord2hp(ra4,dec4, NSIDE)
    
    return {'0.1<z<0.2' : ref1, '0.3<z<0.4' : ref2, '0.5<z<0.6' : ref3, '1.2<z<1.3' : ref4}


def generate_dust_maps(NSIDE):
    '''
    Gets the re
    
    '''
    ebv1998 = np.load(full_path_dust + 'EBV_1998.npy')
    map1 = hp.reorder(ebv1998, n2r=True) #convert ordering to ring
    map1 = hp.ud_grade(map1, NSIDE) #degrade or upgrade resolution to fit other maps
    map1 = _remove_local(map1) #remove local average
    
    ebv2015 = np.load(full_path_dust + 'EBV_2015.npy')
    map2 = hp.reorder(ebv2015, n2r=True) #convert ordering to ring
    map2 = hp.ud_grade(map2, NSIDE) #degrade resolution to nside of 128
    map2 = _remove_local(map2) #remove local average
    
    return [map1, map2]
