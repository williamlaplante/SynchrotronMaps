import healpy as hp
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from array import array

full_path_ref = '/Users/williiamlaplante/Research/SynchrotronMaps/data/ReferenceSamples/'
full_path_dust = '/Users/williiamlaplante/Research/SynchrotronMaps/data/DustMapsData/'


def change_coord(m, coord):
    """ Change coordinates of a HEALPIX map

    Parameters
    ----------
    m : map or array of maps
      map(s) to be rotated
    coord : sequence of two character
      First character is the coordinate system of m, second character
      is the coordinate system of the output map. As in HEALPIX, allowed
      coordinate systems are 'G' (galactic), 'E' (ecliptic) or 'C' (equatorial)

    Example
    -------
    The following rotate m from galactic to equatorial coordinates.
    Notice that m can contain both temperature and polarization.
    >>>> change_coord(m, ['G', 'C'])
    """
    # Basic HEALPix parameters
    npix = m.shape[-1]
    nside = hp.npix2nside(npix)
    ang = hp.pix2ang(nside, np.arange(npix))

    # Select the coordinate transformation
    rot = hp.Rotator(coord=reversed(coord))

    # Convert the coordinates
    new_ang = rot(*ang)
    new_pix = hp.ang2pix(nside, *new_ang)

    return m[..., new_pix]


def coord2hp(ra, dec, nside, nest=False, lonlat=True, remove_zeros=True, galactic=True):
    '''
    Transformation from a set of ra, dec coordinates in degrees to an overdensity field in healpix format.
    
    Parameters:
    -----------------
    
    ra: ra coordinates as numpy array
    
    dec: dec coordinates as numpy array
    
    nside: nside parameter of healpix map
    
    nest: Ordering of healpix map. If true, nested, else, ring.
    
    lonlat: Degrees or radians, and latitude vs co-latitude. See healpy documentation for specifics.
    
    remove_zeros: If true, pixels holding 0 count will hold the value of hp.UNSEEN. Useful to exclude 0's from further calculations.
    
    galactic: If true, transforms the map from equatorial to galactic coordinates. Else it returns it in equatorial coordinates.
    
    
    '''
    
    arr = hp.ang2pix(nside,ra,dec,nest=nest,lonlat=lonlat)
    count_arr = np.histogram(arr, bins=hp.nside2npix(nside), range=[0,hp.nside2npix(nside)])[0].astype('float')
    
    if count_arr.mean()==0:
        raise Exception('Empty set of coordinates.')
    
    mean = count_arr[count_arr!=0].mean()
    count_arr[count_arr==0] = hp.UNSEEN
    
    if galactic:
        count_arr = change_coord(count_arr, ['C', 'G'])
        
    count_arr[count_arr!=hp.UNSEEN] = count_arr[count_arr!=hp.UNSEEN]/mean - 1
    return count_arr


def remove_local(m,nest=False, rad=1):
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
    '''
    Note: 
    - histogram of z values in full_ref_sample.fits reveals identical shape to fig 5 in paper.
    - sample size of full_ref_sample.fits is ~1.13 mil, which is exactly as in the paper.   
    - full_ref_sample.fits has the property that its galactic coordinates are all above 50 deg, and its limited on the zrange.
    '''
    refsample = fits.open(full_path_ref + 'full_ref_sample.fits')[1].data
    ra = refsample['RA']
    dec = refsample['DEC']
    z = refsample['Z']

    
    #We isolate the ranges as they do in the paper to perform plotting
    ra1, dec1 = ra[np.where((0.1<z) & (z<0.2))], dec[np.where((0.1<z) & (z<0.2))] 
    ra2, dec2 = ra[np.where((0.3<z) & (z<0.4))], dec[np.where((0.3<z) & (z<0.4))]
    ra3, dec3 = ra[np.where((0.5<z) & (z<0.6))], dec[np.where((0.5<z) & (z<0.6))]
    ra4, dec4 = ra[np.where((1.2<z) & (z<1.3))], dec[np.where((1.2<z) & (z<1.3))]
    
    #we convert to healpix format with histogram binning
    ref1 = coord2hp(ra1,dec1, NSIDE)
    ref2 = coord2hp(ra2,dec2, NSIDE)
    ref3 = coord2hp(ra3,dec3, NSIDE)
    ref4 = coord2hp(ra4,dec4, NSIDE)
    
    return {'0.1<z<0.2' : ref1, '0.3<z<0.4' : ref2, '0.5<z<0.6' : ref3, '1.2<z<1.3' : ref4}




def generate_ref_maps_v2(nside):
    
    refsample = fits.open(full_path_ref + 'full_ref_sample.fits')[1].data
    ra,dec,z = refsample['RA'],refsample['DEC'],refsample['Z']
    correct_idx = map_cut_idx(nside)
    wrong_idx = np.setdiff1d(np.arange(hp.nside2npix(nside)), correct_idx)
    ra_ranges = [ra[np.where((z1<z) & (z<z2))] for z1,z2 in zip([0.1,0.3,0.5,1.2],[0.2,0.4,0.6,1.3])]
    dec_ranges  = [dec[np.where((z1<z) & (z<z2))] for z1,z2 in zip([0.1,0.3,0.5,1.2],[0.2,0.4,0.6,1.3])]

    maps = []
    for i in range(4):
        arr = hp.ang2pix(nside, ra_ranges[i], dec_ranges[i], nest=False,lonlat=True)
        count_arr = np.histogram(arr, bins=hp.nside2npix(nside), range=[0,hp.nside2npix(nside)])[0].astype('float')
        count_arr[correct_idx] = count_arr[correct_idx]/count_arr[correct_idx].mean() - 1
        count_arr[wrong_idx] = hp.UNSEEN
        count_arr = change_coord(count_arr, ['C', 'G'])
        maps.append(count_arr)
        
    return maps
    
    
def generate_dust_maps(NSIDE):
    '''
    Gets the re
    
    '''
    ebv1998 = np.load(full_path_dust + 'EBV_1998.npy')
    map1 = hp.reorder(ebv1998, n2r=True) #convert ordering to ring
    map1 = hp.ud_grade(map1, NSIDE) #degrade or upgrade resolution to fit other maps
    map1 = remove_local(map1) #remove local average
    
    ebv2015 = np.load(full_path_dust + 'EBV_2015.npy')
    map2 = hp.reorder(ebv2015, n2r=True) #convert ordering to ring
    map2 = hp.ud_grade(map2, NSIDE) #degrade resolution to nside of 128
    map2 = remove_local(map2) #remove local average
    
    return [map1, map2]


def generate_dust_map_1998(nside):
    ebv1998 = np.load(full_path_dust + 'EBV_1998.npy')
    map1 = hp.reorder(ebv1998, n2r=True) #convert ordering to ring
    map1 = hp.ud_grade(map1, nside) #degrade or upgrade resolution to fit other maps
    map1 = remove_local(map1) #remove local average
    return map1



def map_cut_idx(nside, gal_cut=50):
    '''
    returns the index of all values in a healpix map of resolution nside respecting a galactic cut of gal_cut.
    
    '''
    theta, phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    idx = SkyCoord(theta,phi,unit='deg',frame='fk5').galactic.b.value >= gal_cut
    return np.arange(hp.nside2npix(nside))[idx]


def read_ref_map(nside, zmin, zmax):
    if zmin not in [0.1,0.3,0.5,1.2] or zmax not in [0.2,0.4,0.6,1.3]:
        raise Exception('Please input a correct zmin or zmax.')
    filename = str(zmin)+'z'+str(zmax)+'_'+str(nside)+'_'+'overdensity_field.bin'
    filepath = path_ref = "/Users/williiamlaplante/Research/SynchrotronMaps/data/processed_maps/ref_maps/" + filename
        
    with open(filepath, "rb") as f:
        float_array = array('d')
        float_array.frombytes(f.read())
        
    return np.frombuffer(float_array, dtype=float)

def read_dust_map(nside):
    filename = "ebv_" + str(nside) + "_1998.bin"
    filepath = "/Users/williiamlaplante/Research/SynchrotronMaps/data/processed_maps/dust_maps/"+filename
            
    with open(filepath, "rb") as f:
        float_array = array('d')
        float_array.frombytes(f.read())
        
    return np.frombuffer(float_array, dtype=float)