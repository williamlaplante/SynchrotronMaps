import unittest
from corrfunc import compute_corr
import numpy as np
import healpy as hp

class TestCorrfunc(unittest.TestCase):
    
    def test_compute_corr(self):
        
        #Make two maps of same resolution
        nside = 32
        a = np.random.normal(size=hp.nside2npix(nside))
        b = np.random.normal(size=hp.nside2npix(nside))
        c = np.random.normal(size=hp.nside2npix(nside))
        
        #Set parameters for corrfunc
        R = 7 #degrees
        dr = 2*np.degrees(hp.nside2resol(nside)) #degrees
        
        #test for symmetry
        mean1, err1 = compute_corr(a,b,R,dr)
        mean2, err2 = compute_corr(b,a,R,dr)
        self.assertAlmostEqual(mean1, mean2, places=5)
        self.assertAlmostEqual(err1, err2, places=3)
        
        #tests for linearity
        mean1, err1 = compute_corr(5*a,4*b,R,dr)
        mean2, err2 = compute_corr(a,b,R,dr)
        self.assertAlmostEqual(mean1, 20*mean2, places=5)
        
        mean1, err1 = compute_corr(a+b,c, R, dr)
        mean2, err2 = compute_corr(a,c,R,dr)
        mean3, err3 = compute_corr(b,c,R,dr)
        self.assertAlmostEqual(mean1, (mean2+mean3))
                                  
                                  
        #test for field of ones
        a = np.ones(hp.nside2npix(nside))
        b = a
        
        self.assertAlmostEqual(compute_corr(a,b,R,dr)[0], 1, places=4)
        
        #test for field of zeros
        a = np.zeros(hp.nside2npix(nside))
        self.assertAlmostEqual(compute_corr(a,b,R,dr)[0], 0)
        
    def test_query_annulus(self):
        #resolution vs dr test
        
        #r<resolution -> full circle
        
        #combination of r1-r2 and r2-r3 same as r1-r3
        
        #how it deals hp.UNSEEN ; different tests for maps with boundaries
        pass
        
        
    def test_get_annuli_len(self):
        #tests for consistency with query_annulus
        pass

if __name__ == '__main__':
    unittest.main()
        
        
