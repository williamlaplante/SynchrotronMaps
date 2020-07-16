def get_ring_len(nside, theta, step=1, nest=False):
    '''
    Gets the average ring length for some radius theta and resolution nside.
    
    '''
    dr = hp.nside2resol(nside)
    npix = hp.nside2npix(nside)
    X,Y,Z = hp.pix2vec(nside, np.arange(0,npix,step), nest=nest)
    l = 0
    for vec in zip(X,Y,Z):
        l+=hp.query_disc(nside, vec, np.radians(theta)+dr,nest=nest).size - hp.query_disc(nside,vec,np.radians(theta),nest=nest).size

    return l/len(X)

    
    
def compute_corr(map1, map2, theta, nest=False):
    '''
    
    Parameters:
    ------------
    map1: first healpix map
    
    map2: second healpix map
    
    theta: distance |r1-r2| in degrees
    
    thickness: thickness of annuli in degrees
    
    nest: ordering of the healpix maps
    
    '''
    if hp.get_nside(map1)!=hp.get_nside(map2):
        raise Exception('Maps are of two different resolutions.')
        
    #set healpix maps parameters
    nside = hp.get_nside(map1)
    npix = hp.nside2npix(nside)
    
    #Get the vectors associated with pixels. We exclude the hp.UNSEEN points from our results.
    pix1 = np.arange(npix)[map1!=hp.UNSEEN] 
    pix2 = np.arange(npix)[map2!=hp.UNSEEN]
    pix_inter =  np.intersect1d(pix1,pix2)
    X,Y,Z = hp.pix2vec(nside, pix_inter, nest=nest)
    
    #Fix the ring parameters
    ring_len = get_ring_len(nside, theta, step=4, nest=nest)
    dr = 6*hp.nside2resol(nside)
    
    def query_annulus(vec):
        disk1 = hp.query_disc(nside, vec, np.radians(theta)+dr, nest=nest)
        disk2 = hp.query_disc(nside, vec, np.radians(theta), nest=nest)
        index = np.setdiff1d(disk1, disk2, assume_unique=True)
        circle = map2[index]
        if index.size==0:
            err = 'queried annulus of null size with radius '+str(np.radians(theta))+'rad and thickness '+str(dr)+'rad'
            raise Exception(err)
        return circle[circle!=hp.UNSEEN]
    
    it = (query_annulus(vec).sum()/ring_len for vec in zip(X,Y,Z))
    
    result = map1[pix_inter]* np.fromiter(it, dtype=np.float)
    
    return result.mean(), result.std()/np.sqrt(npix)


def plot_corr(dust_map, *ref_maps):

    x = np.exp(np.linspace(-4.2,0.8,15))
    
    colors = ['red', 'pink', 'cyan', 'blue']
    
    plt.figure(figsize=(8,6))
    for ref_map, color in tqdm(zip(ref_maps, colors),total=len(ref_maps)):
        
        out=[compute_corr(ref_map, dust_map, theta, 6*hp.nside2resol(256)) for theta in x] 
        
        y, err = zip(*out)
        plt.errorbar(x*60,list(y), yerr=list(err), color=color)
        
        
    plt.xscale('log')
    plt.xlabel('degrees (arcmin)')
    plt.ylabel('correlation')
    return 



#Making the maps:
NSIDE = 256

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

ebv1998 = np.load(full_path_dust + 'EBV_1998.npy')

map1 = hp.reorder(ebv1998, n2r=True) #convert ordering to ring

map1 = hp.ud_grade(map1, NSIDE) #degrade or upgrade resolution to fit other maps

map1 = remove_local(map1) #remove local average


plot_corr(map1, ref1 ,ref2, ref3, ref4)