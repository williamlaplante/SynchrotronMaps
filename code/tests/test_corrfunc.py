import sys
sys.path.append('../')
import unittest
from corrfunc import compute_corr, query_annulus
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
        

        #test for annuli with thickness that is too small for map resolution
        a = np.random.normal(size=hp.nside2npix(nside))
        b = np.random.normal(size=hp.nside2npix(nside))
        thick = np.degrees(hp.nside2resol(nside))

        self.assertEqual(compute_corr(a,b,R,thick)[0], compute_corr(a,b,R,thick/10)[0])

        #test for when R=0
        self.assertEqual(compute_corr(a,b,0,thick)[0], (a*b).mean())
        
        #test for when radius is smaller than resolution
        self.assertEqual(compute_corr(a,b,thick/6,thick)[0], compute_corr(a,b,thick/4,thick)[0])

    def test_query_annulus(self):
        nside=32
        resol = np.degrees(hp.nside2resol(nside))
        center = hp.pix2vec(nside,0)
        m = np.random.normal(size=hp.nside2npix(nside))

        #resolution vs dr test
        resol = np.degrees(hp.nside2resol(nside))
        self.assertRaises(Exception, query_annulus, m, nside, center, 5, resol/100, idx=True)

        #when r << resolution -> same as full circle. (i.e. when the inner radius is below the lower limit which is a func of resol)
        a1 = np.append(np.zeros(1),query_annulus(m,nside,center,resol/100,10,idx=True)).astype(int)
        a2 = hp.query_disc(nside,center,np.radians(10))
        self.assertTrue(np.allclose(a1,a2))

        #combination of r1-r2 and r2-r3 same as r1-r3
        r1_r2 = query_annulus(m,nside,center,5,1,idx=True)
        r2_r3 = query_annulus(m,nside,center,6,1,idx=True)
        r1_r3 = query_annulus(m,nside,center,5,2,idx=True)
        arr1 = np.sort(np.concatenate((r1_r2,r2_r3)))
        arr2 = np.sort(r1_r3)
        self.assertEqual(arr1.size, arr2.size)
        self.assertTrue(np.allclose(arr1,arr2))
   
        #how it deals hp.UNSEEN ; different tests for maps with boundaries

        pass
        
        
    def test_get_annuli_len(self):
        #tests for consistency with query_annulus
        pass

if __name__ == '__main__':
    unittest.main()
        
        
