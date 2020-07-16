#Functions related to statistics applied on the data.

import numpy as np
import healpy as hp
from progress.bar import Bar
from tqdm.notebook import tqdm
from scipy.special import legendre
from itertools import product


#Dot product function (normalized)
def dot_product(X,Y):
    normx = np.linalg.norm(X)
    normy = np.linalg.norm(Y)
    
    if normx == 0 or normy == 0:
        return 0
    else:
        return np.multiply(X/normx, Y/normy).sum()
        

#Correlation function
def correlation(X,Y):
    '''
    Correlation function. Assumes the input x and y are two numpy matrices of identical 
    shape. The formula is 
    
    corr = COV(x,y)/(STD(x)*STD(y)) = COV(x,y)/sqrt(VAR(x)*VAR(y))
    
    Note that 
    
    COV(x,y) = E(x.*y) - E(x)E(y),
    
    where .* is pointwise multiplication between the two matrices.
    '''
    
    
    cov = np.multiply(X,Y).mean() - X.mean()*Y.mean()
    
    if X.std()==0 or Y.std()==0:
        return 0
    else:
        return cov/(X.std()*Y.std())


#Algorithm dividing into subarrays two matrices and computing a function that takes submatrices as arguments.
def sub_apply(a, b, func, subarr_size = 10, filter_zeros=True): 
    '''
    Algorithm:
    Takes two matrices of size NxN. Then it divides each matrix into submatrices of size 'subarr_size' x 'subarr_size',
    and computes the given function on submatrices from one matrix and the other element-wise (elements being the 
    submatrices). From this an array of scalar values is obtained. The mean and std are computed and returned.
    
    Parameters :
    ----------------------
    
    a: First matrix of size (N,N)
    
    b: Second matrix of size (N,N)
    
    func: function that will be applied on subarrays.
    
    subarr_size: The wanted size of the subarrays. Must be a multiple of N.
    
    filter_zeros: If true, removes all correlation values that are exactly equal to 0 from the array.
    
    '''
    
    if a.shape != b.shape: #We want same dimensions for both matrices
        raise Exception('The matrices are of different shape.')
        return 0
    
    if a.shape[0] != a.shape[1]: #We want square matrices
        raise Exception('The matrices are not squares. They have shape: {}'.format(a.shape))
        return 0
    
    N,_ = a.shape
    vals = []
    
    if N % subarr_size != 0: #We want the size of the subarray to be a multiple of the size of the whole array
        return 0
    
    step = subarr_size #This is solely to express understanding of how subarray size translates to algorithm.
 
    for i in range(0,N,step):
        for j in range(0,N,step):
            vals.append(func(a[i:i+step, j:j+step],b[i:i+step, j:j+step]))
            
    if filter_zeros:
        vals = list(filter(lambda x: x!=0, vals)) #Remove all 0 values.
    
    return np.mean(vals), np.std(vals)


#Same as sub_apply, except we also compute function between "not-element-wise" submatrices.
def cross_sub_apply(arr_fixed, arr_iter, func, subarr_size=10, filter_zeros=True, weights=False): 
    
    '''
    Similar in concept to sub_apply. The difference is that instead of applying the function subarray-wise only,
    that is, first subarray of array #1 taken with first subarray of array #2, second with second, etc...
    we have that the first subarray of arr_fixed is taken with all of the subarrays of arr_iter, and we average over
    these values to get our first scalar value. We repeat for the second subarray of fixed_arr,the third, etc...
    and this way we get (N/subarr_size)^2 scalar values. We then return the mean and std of these values.
    
    
    Parameters:
    --------------
    arr_fixed: The array that is fixed. 
    
    arr_iter: The array over which we iterate.
    
    func: The function applied on the subarrays.
    
    subarr_size: If subarray size is N^2, then subarr_size is N.
    
    weights: If true, we compute a weighted average on values gotten from a single iteration. The weights are given by 1/(1+distance).
    
    '''
    if arr_fixed.shape != arr_iter.shape: #We want same dimensions for both matrices
        raise Exception('The matrices are of different shape.')
 
    if arr_fixed.shape[0] != arr_fixed.shape[1]: #We want square matrices
        raise Exception('The matrices are not squares. They have shape: {}'.format(a.shape))
    
    N,_ = arr_fixed.shape
    vals = []
    
    if N % subarr_size != 0: #We want the size of the subarray to be a multiple of the size of the whole array
        raise Exception('The subarray size must be a multiple of the array size')
    
    step = subarr_size #This is solely to express understanding of how subarray size translates to algorithm.
 
    for i in range(0,N,step):
        for j in range(0,N,step):
            subvals = []
            distance = []
            
            for n in range(0,N,step):
                for m in range(0,N,step):
                    subvals.append(func(arr_fixed[i:i+step, j:j+step], arr_iter[n:n+step, m:m+step]))
                    distance.append(max([abs(n-i)/step, abs(m-j)/step]))
            if filter_zeros:
                subvals = np.array(subvals)
                distance = np.array(distance)[subvals!=0]
                subvals = subvals[subvals!=0]
            
            if len(subvals)!=0:
                if weights:
                    subvals_avg = np.average(subvals, weights=1/(1+np.array(distance)))
                else:
                    subvals_avg = np.mean(subvals)
                    
                vals.append(subvals_avg)
    
    return np.mean(vals), np.std(vals)


#Applies dot product function above on elements lying on radius
def circular_product(X,Y, filter_zeros=True):
    '''
    Algorithm: We pick an element at position i,j from both matrices X,Y. We pick the boundary elements of a square of size 2*theta+1
    (basically, the side of the square is a distance theta from the point i,j) centered at i,j from both X and Y. We compute the
    normalized dot product between the two boundaries gotten from X and Y. we then increase theta and repeat the latter action. This 
    gives us a set of subvalues associated with position i,j. We average those subvalues, giving us a single value associated with
    position i,j. We then repeat for each position in the matrix. Finally, we have a total of NxN values. We return the mean and std.

    Note that boundary conditions are not circular. They are flattening. That is, we return 0 when going in negative indices, and we
    return N-1 when going over the array size. (For keen observers, we see that the function actually returns N instead of N-1 when 
    slice=True. This is because when slicing numpy array like X[a:b], a <= index < b. 

    Parameters:
    -------------
    X: The first matrix of size NxN.

    Y: The second matrix of size NxN.

    filter_zeros: If true, we remove the 0's from the sub_value array. (array containing values for fixed position, changing theta)


    '''
        
    if X.shape != Y.shape: 
        raise Exception('The matrices are of different shape.')
        return 0
    
    if X.shape[0] != X.shape[1]: 
        raise Exception('The matrices are not squares. They have shape: {}'.format(a.shape))
        return 0
    
    N,_ = X.shape
    
    def b(index, slice=False):
    
        if index<=0:
            return 0

        if slice:
            if index>N: return N
            else : return index
        else:
            if index>=N : return N-1
            else : return index
    
    vals = []
    
    for i in range(N):
        for j in range(N):
            subvals = []
            for theta in range(1,N):
                upper_border = X[b(i-theta), b(j-theta, True) : b(j+theta, True)]
                lower_border = X[b(i+theta), b(j-theta, True) : b(j+theta, True)]
                left_border = X[b(i-theta, True) : b(i+theta, True), b(j-theta)]
                right_border = X[b(i-theta, True) : b(i+theta, True), b(j+theta)]
                
                X_circle = np.concatenate((upper_border, lower_border, left_border, right_border), axis=0)
                
                upper_border = Y[b(i-theta), b(j-theta, True) : b(j+theta, True)]
                lower_border = Y[b(i+theta), b(j-theta, True) : b(j+theta, True)]
                left_border = Y[b(i-theta, True) : b(i+theta, True), b(j-theta)]
                right_border = Y[b(i-theta, True) : b(i+theta, True), b(j+theta)]
                
                Y_circle = np.concatenate((upper_border, lower_border, left_border, right_border), axis=0)
                
                subvals.append(dot_product(X_circle, Y_circle))
                
            if filter_zeros:
                subvals = list(filter(lambda x: x!=0, subvals))
            
            if subvals:    
                vals.append(np.mean(subvals))
    
    return np.mean(vals), np.std(subvals)
                

#Computes arr[i,j] - <arr_i,j> where arr_i,j is a subarray of size kxk centered at i,j. i.e., localized centering.
def filter_arr(arr, k):
    '''
    Takes an array of size NxM, and for each point in the array, it centers it locally. i.e., removes the mean from the point.
    
    Conditions are as follows: 
        -We pick a radius of k for all points except for border stripes of size k. For the latter case its a bit different.
        -Boundary conditions are not cyclic. 
        
    Parameters:
    ---------------
    arr : The array we wish to apply a function on locally.
    k : The radius of the subarrays. i.e. if k=1, we look at 3x3 sub-arrays, k=2 means 4x4, etc...
    
    '''
    
    N,M = arr.shape
    final_arr = np.empty([N,M])
      
    #Center of array
    
    for i in range(k,N-k):
        for j in range(k, M-k):
            final_arr[i,j] = arr[i,j] - arr[i-k:i+k+1, j-k:j+k+1].mean()
        

    #Four sides, excluding corners
    
    #Vertical lefmost stripe
    for i in range(k,N-k):
        for j in range(0, k):
            final_arr[i,j] = arr[i,j] - arr[i-k:i+k+1, 0:k].mean()
    
    #Vertical rightmost stripe
    for i in range(k,N-k):
        for j in range(M-k, M):
            final_arr[i,j] = arr[i,j] - arr[i-k:i+k+1, M-k:M].mean()
 
    #Horizontal uppermost stripe
    for i in range(0,k):
        for j in range(k, M-k):
            final_arr[i,j] = arr[i,j] - arr[0:k, j-k:j+k+1].mean()
            
    #Horizontal lowermost stripe
    for i in range(N-k,N):
        for j in range(k, M-k):
            final_arr[i,j] = arr[i,j] - arr[N-k:N, j-k:j+k+1].mean()
    
    #Upper-left corner
    for i in range(0,k):
        for j in range(0,k):
            final_arr[i,j] = arr[i,j] - arr[0:k, 0:k].mean()
    
    #Upper-right corner
    for i in range(0,k):
        for j in range(M-k,M):
            final_arr[i,j] = arr[i,j] - arr[0:k, M-k:M].mean()
    
    #Lower-left corner
    for i in range(N-k, N):
        for j in range(0,k):
            final_arr[i,j] = arr[i,j] - arr[N-k:N, 0:k].mean()
    
    #Lower-right corner
    for i in range(N-k, N):
        for j in range(M-k,M):
            final_arr[i,j] = arr[i,j] - arr[N-k:N, M-k:M].mean()
        
    return final_arr

def get_circle(mat, center_point, radius, full=False):
    '''
    Used to retrieve points on an approximated circle on a 2d array (matrix). Center point 
    is an index and radius is in term of indices.
    
    Parameters:
    ------------
    mat: A numpy matrix of size (N,M)
    
    center_point: The center of the circle, given as an index tuple (i,j)
    
    radius: radius of approximated circle
    
    '''
    if full:
        f = lambda x:x
    else:
        f = abs
        
    j, i = np.meshgrid(range(mat.shape[1]), range(mat.shape[0]))
        
    return mat[f((i-center_point[0])**2 + (j-center_point[1])**2 - radius**2) < radius]

    

def get_circle_coord(x, y, center_point, radius, full=True):
    '''
    Used to retrieve all points of an x-y grid found on the boudary of a circle.
    
    Parameters:
    ------------
    x: x coordinate
    
    y: y coordinate
    
    center_point: The center of the circle, given in same units as x,y coordinates
    
    radius: radius of circle in units of x,y coordinates
    
    full: If false, returns only coordinates on the boundary.
    
    '''
    if full:
        mask = (x-center_point[0])**2 + (y-center_point[1])**2 <= radius**2
    else:
        epsilon=radius
        mask = abs((x-center_point[0])**2 + (y-center_point[1])**2 - radius**2) < epsilon
    
    return x[mask], y[mask]




def w(iter_map, circle_map, theta):
    '''
    Two point correlation function approximated. 
    
    parameters:
    ------------
    iter_map: map from which we draw the scalar values at each index.
    
    circle_map: map from which we draw the circle values centered at index of scalar value of iter_map.
    
    theta: radius of circle drawn / distance of two point correlation function.
    
    returns w(theta) and the error on w(theta).
    
    The error comes from the variance of the dot product over position in the sky. Hence, for a fixed theta, it quantifies how much
    variation there is over the sky for the "circular" dot product.
    '''
    
    if iter_map.shape != circle_map.shape:
        raise Exception('matrices are of different shape.')
    
    
    it = (get_circle(circle_map, index, theta).mean() for index in product(range(circle_map.shape[0]), range(circle_map.shape[1])))
    arr = np.fromiter(it, dtype=np.float).reshape(iter_map.shape)
    res = iter_map*arr
    
    return res.mean(), res.std()/np.sqrt(res.size)



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
    
    galactic: If true, transforms the map to galactic coordinates. Else it returns it in equatorial coordinates.
    
    
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




def field2healpix(ra, dec, nside, prog=True):
    '''
    Converts a density field to a healpix map. Assumes density field is in degrees,
    and has coordinates system ra (x-axis) and dec (y-axis).
    
    
    Parameters:
    --------------
    nside: resolution of wanted healpix map
    
    ra: ra coordinates forming the field
    
    dec: dec coordinates forming the field
    

        
    '''
    resol = np.degrees(hp.nside2resol(nside, arcmin=False))
    radius = resol/2
    N = len(ra)
    
    index = np.arange(hp.nside2npix(nside))
    theta, phi = hp.pix2ang(nside, index, nest=False, lonlat=True)
    
    if prog:
        it = (((ra-cp[0])**2 + (dec-cp[1])**2 <= radius**2).sum() for cp in tqdm(zip(theta,phi), total=len(theta)))
    else:
        it = (((ra-cp[0])**2 + (dec-cp[1])**2 <= radius**2).sum() for cp in zip(theta,phi))
        
    return np.fromiter(it, dtype=np.float)


def C(theta, C_l):
    '''
    Two point correlation function defined by its fourier transform (basically) for some fixed theta
    
    Parameters:
    -------------
    theta: separation angle.
        
    C_l: array of C_l values computed with anafast. 
    
    '''
    it = ( ((2*l+1)/4*np.pi) * legendre(l)(np.cos(theta)) for l in range(len(C_l)) )
    return np.dot(C_l, np.fromiter(it, dtype=np.float))



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












#CORRELATION FUNCTION








#The database-independent two-point correlation function.

def compute_corr_v1(map1, map2, theta, nest=False):
    '''
    version 1 of w(theta). This version gets the whole list of possible m(r1)m(r2) using the ring method 
    to then compute <m(r1)m(r2)>. It is very explicit and thus slow. On the other hand, it is clear what the
    function does.
    
    The function is symmetric and linear. (experimentally tested, and it is also logical)
    
    Parameters:
    ------------
    map1: first healpix map
    
    map2: second healpix map
    
    theta: distance |r1-r2| in radians
    
    nest: ordering of the healpix maps
    
    '''
    if hp.get_nside(map1)!=hp.get_nside(map2):
        raise Exception('Maps are of two different resolutions.')
        
    #set healpix maps parameters
    nside = hp.get_nside(map1)
    npix = hp.nside2npix(nside)
    
    #Get the vectors associated with pixels
    X,Y,Z = hp.pix2vec(nside, np.arange(npix), nest=False)
    
    #Function used to retrieve mean of circle with center vec and radius theta in map2
    def get_circle(vec):
        dr = 10*hp.nside2resol(nside)  
        
        disk1 = hp.query_disc(nside, vec, theta+dr, nest=nest)
        disk2 =  hp.query_disc(nside, vec, theta, nest=nest)
        index = np.setdiff1d(disk1,disk2, assume_unique=True)
        
        circle = map2[index].astype(np.float)
        return circle[circle!=hp.UNSEEN]
               
    #Function used to retrieve scalar at vector position x,y,z in map1
    def get_scalar(x,y,z):
        scalar = map1[hp.vec2pix(nside, x, y, z, nest=nest)]
        if scalar == hp.UNSEEN:
            return 0
        else:
            return scalar
   
    vals = []
    for x,y,z in tqdm(zip(X,Y,Z), total=X.size):
        vals.extend(get_scalar(x,y,z)*get_circle([x,y,z]))
    
    ans = np.array(vals)
    
    return ans[ans!=0].mean(), ans[ans!=0].std()







#Function used to store the index ring matrices in the CorrelatonFunctionDB file.

def make_ring_matrices(nside, thetas, nest=False, save=True):
    
    for theta in tqdm(thetas, position=0):
        npix = hp.nside2npix(nside)
        dr = 2*hp.nside2resol(nside)
        R = np.radians(theta)
        
        X,Y,Z = hp.pix2vec(nside, np.arange(npix), nest=nest)
        
        l=0

        #We look for the smallest size of the indices first. This is to avoid indices mismatch between rows.
        for x,y,z in tqdm(zip(X,Y,Z),total=npix,position=1):
            disk1_size = hp.query_disc(nside, [x,y,z], R+dr, nest=nest).size
            disk2_size =  hp.query_disc(nside, [x,y,z], R, nest=nest).size
            if (disk1_size-disk2_size)>l: l=(disk1_size-disk2_size)

        #We initialize our matrix of indices.
        A = np.empty((npix,l), np.int)

        #We loop over all pixels in the map, and retrieve the indices of the values
        #on the circle centered on that pixel and of radius theta.
        #We thus fill up a matrix of indices, where each row is the indices of the values on that circle.

        for x,y,z,i in tqdm(zip(X,Y,Z,range(npix)),total=npix, position=2):
            disk1 = hp.query_disc(nside, [x,y,z], R+dr, nest=nest)
            disk2 =  hp.query_disc(nside, [x,y,z], R, nest=nest)
            idx = np.setdiff1d(disk1,disk2, assume_unique=True)
            
            if idx.size<l:
                #pad with nan values up to size l
                num = l-idx.size
                idx = np.append(idx, [npix]*num)
            
            A[i,:] = idx

        if save:
            full_path = '/Users/williiamlaplante/Research/SynchrotronMaps/data/CorrelationFunctionDB/'
            if nest: n = '_nest'
            else: n = '_ring'
            np.save(full_path + str(nside) + '_' + str(theta) + n + '.npy', A)
    
    return 







#Function used to compute the correlation function values beteen two maps. It assumes the index array is in the DB.
def compute_corr_v2(m1,m2, theta, nest=False):

    
    if hp.get_nside(m1)!=hp.get_nside(m2):
        raise Exception('Maps are of different resolution.')
    
        
    if theta==0:
        return (m1*m2).mean(), (m1*m2).std()
    
    nside = hp.get_nside(m1)
    full_path = '/Users/williiamlaplante/Research/SynchrotronMaps/data/CorrelationFunctionDB/'
    if nest:
        n = '_nest'
    else:
        n = '_ring'
     
    try:
        theta_list = list(str(theta))
        while len(theta_list)<6:
            theta_list.append('0')
            
        A = np.load(full_path + str(nside) + '_' + ''.join(theta_list)[:6] + n + '.npy')
    except FileNotFoundError:
        raise Exception('There is no index matrix for this theta and resolution, you must create one first.')
        return
    
    arr1 = np.append(m1, hp.UNSEEN)[A]
    arr1[arr1==hp.UNSEEN] = np.nan
    m2[m2==hp.UNSEEN] = np.nan
    arr1 = arr1*m2.reshape(-1,1)
    return np.nanmean(arr1, dtype=np.float), np.nanstd(arr1, dtype=np.float)