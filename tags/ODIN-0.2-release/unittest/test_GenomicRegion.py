from __future__ import print_function
import unittest
from rgt.GenomicRegion import GenomicRegion


class TestGenomicRegion(unittest.TestCase):
    def test_overlap(self):
        r = GenomicRegion(chrom=1, initial=10, final=15)
        
        #usual cases
        r2 = GenomicRegion(chrom=1, initial=20, final=25)
        self.assertFalse(r.overlap(r2))
        
        r2 = GenomicRegion(chrom=1, initial=0, final=5)
        self.assertFalse(r.overlap(r2))
        
        r2 = GenomicRegion(chrom=1, initial=7, final=12)
        self.assertTrue(r.overlap(r2))
        
        r2 = GenomicRegion(chrom=1, initial=12, final=18)
        self.assertTrue(r.overlap(r2))
        
        r2 = GenomicRegion(chrom=1, initial=12, final=14)
        self.assertTrue(r.overlap(r2))
        
        #r2 within r
        r2 = GenomicRegion(chrom=1, initial=11, final=13)
        self.assertTrue(r.overlap(r2))
        
        #border cases
        #GenomicRegions touch, but do not overlap
        r2 = GenomicRegion(chrom=1, initial=5, final=10)
        self.assertFalse(r.overlap(r2))
        
        #here, they overlap
        r2 = GenomicRegion(chrom=1, initial=5, final=11)
        self.assertTrue(r.overlap(r2))
        
        #they touch, do not overlap
        r2 = GenomicRegion(chrom=1, initial=15, final=20)
        self.assertFalse(r.overlap(r2))
        
        #they overlap in 1 bp (14th)
        r2 = GenomicRegion(chrom=1, initial=14, final=20)
        self.assertTrue(r.overlap(r2))
        
        #they have zero length
        r = GenomicRegion(chrom=1, initial=10, final=10)
        r2 = GenomicRegion(chrom=1, initial=10, final=10)
        self.assertFalse(r.overlap(r2))
        
        #they have zero length
        r = GenomicRegion(chrom=1, initial=10, final=10)
        r2 = GenomicRegion(chrom=1, initial=11, final=11)
        self.assertTrue(r.overlap(r2))
        
        #they have zero length
        r = GenomicRegion(chrom=1, initial=10, final=10)
        r2 = GenomicRegion(chrom=1, initial=5, final=10)
        self.assertTrue(r.overlap(r2))
        
    def test_extend(self):
        #normal extend
        r = GenomicRegion(chrom=1, initial=10, final=20)
        r.extend(5,15)
        self.assertEqual(r.initial, 5)
        self.assertEqual(r.final, 35)
        
        #use negative values to extend
        r2 = GenomicRegion(chrom=1, initial=10, final=20)
        r2.extend(-5,-1)
        self.assertEqual(r2.initial, 15)
        self.assertEqual(r2.final, 19)
        
        #extend to under zero
        r3 = GenomicRegion(chrom=1, initial=10, final=20)
        r3.extend(15,0)
        self.assertEqual(r3.initial, 0)
        
        #extend so that inital and final coordinate change
        r4 = GenomicRegion(chrom=1, initial=10, final=20)
        r4.extend(-50,-50)
        self.assertEqual(r4.initial, 0)
        self.assertEqual(r4.final, 60)
        
    def test_len(self):
        r = GenomicRegion(chrom=1, initial=10, final=20)
        self.assertEqual(len(r), 10)
        
    def test_cmp(self):
        r = GenomicRegion(chrom=1, initial=10, final=20)
        
        r2 = GenomicRegion(chrom=1, initial=12, final=22)
        self.assertTrue(r < r2)
        
        r2 = GenomicRegion(chrom=1, initial=8, final=18)
        self.assertTrue(r > r2)
        
        r2 = GenomicRegion(chrom=1, initial=10, final=12)
        self.assertTrue(r > r2)
        
        r2 = GenomicRegion(chrom=1, initial=12, final=14)
        self.assertTrue(r < r2)
        
        r2 = GenomicRegion(chrom='X', initial=4, final=8)
        self.assertTrue(r < r2)
        
        r2 = GenomicRegion(chrom=1, initial=10, final=18)
        self.assertTrue(r >= r2)
    
if __name__ == '__main__':
    unittest.main()
    
    