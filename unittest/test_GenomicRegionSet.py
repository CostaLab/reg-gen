from __future__ import print_function
from __future__ import division
import sys
import unittest
from rgt.GenomicRegion import *
from rgt.GenomicRegionSet import *
import os
from rgt.Util import GenomeData
from rgt.Util import OverlapType


"""Unit Test"""

class TestGenomicRegionSet(unittest.TestCase):
    
    def region_sets(self,listA,listB):
        """ Setting two GenomicRegionSets as self.setA and self.setB for each case test. """
        self.setA = GenomicRegionSet('for Unit Test')
        for i in range(len(listA)):
            self.setA.add(GenomicRegion(chrom=listA[i][0], initial=listA[i][1], final=listA[i][2]))
        
        self.setB = GenomicRegionSet('for Unit Test')
        for i in range(len(listB)):
            self.setB.add(GenomicRegion(chrom=listB[i][0], initial=listB[i][1], final=listB[i][2]))
    
    def test_extend(self):
        """
        Two empty sets
        A : none 
        R : none
        """
        self.region_sets([],
                         [])
        self.setA.extend(100,100)
        self.assertEqual(len(self.setA.sequences), 0)
        """
        One region
        A :   -----
        R : ---------
        """
        self.region_sets([['chr1',5,10]],
                         [])
        result = self.setA
        result.extend(4,4)
        self.assertEqual(len(result.sequences), 1)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 14)
        """
        Many region
        A :   -----   ------         -----    -----
        R : --------=---------     ------------------
        """
        self.region_sets([['chr1',5,10],['chr1',15,20],['chr1',40,50],['chr1',65,75]],
                         [])
        result = self.setA
        result.extend(5,5)
        self.assertEqual(len(result.sequences), 4)
        self.assertEqual(result.sequences[0].initial, 0)
        self.assertEqual(result.sequences[0].final, 15)
        self.assertEqual(result.sequences[1].initial, 10)
        self.assertEqual(result.sequences[1].final, 25)
        self.assertEqual(result.sequences[2].initial, 35)
        self.assertEqual(result.sequences[2].final, 55)
        self.assertEqual(result.sequences[3].initial, 60)
        self.assertEqual(result.sequences[3].final, 80)
        """
        Many region in different chromosome
        A :   -----   ------         -----    -----
        R : none
        """
        self.region_sets([['chr1',5,10],['chr2',15,20],['chr3',40,50],['chr4',65,75]],
                         [])
        result = self.setA
        result.extend(5,5)
        self.assertEqual(len(result.sequences), 4)
        self.assertEqual(result.sequences[0].initial, 0)
        self.assertEqual(result.sequences[0].final, 15)
        self.assertEqual(result.sequences[0].chrom, 'chr1')
        self.assertEqual(result.sequences[1].initial, 10)
        self.assertEqual(result.sequences[1].final, 25)
        self.assertEqual(result.sequences[1].chrom, 'chr2')
        self.assertEqual(result.sequences[2].initial, 35)
        self.assertEqual(result.sequences[2].final, 55)
        self.assertEqual(result.sequences[2].chrom, 'chr3')
        self.assertEqual(result.sequences[3].initial, 60)
        self.assertEqual(result.sequences[3].final, 80)
        self.assertEqual(result.sequences[3].chrom, 'chr4')
        """
        One region
        A :   -----
        R : ---------
        """
        self.region_sets([['chr1',100,200]],
                         [])
        result = self.setA
        result.extend(10,10,percentage=True)
        self.assertEqual(len(result.sequences), 1)
        self.assertEqual(result.sequences[0].initial, 90)
        self.assertEqual(result.sequences[0].final, 210)
        
    def test_sort(self):
        self.region_sets([['chr1',15,20],['chr1',40,50],['chr1',65,75],['chr1',5,10]],
                         [])
        self.setA.sort()
    
    def test_intersect(self):
        """
        Two empty sets
        A : none 
        B : none
        R : none
        """
        self.region_sets([],
                         [])
        result = self.setA.intersect(self.setB)
        self.assertEqual(len(result.sequences), 0)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.ORIGINAL)
        self.assertEqual(len(result.sequences), 0)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.COMP_INCL)
        self.assertEqual(len(result.sequences), 0)
        """
        One empty set
        A :   -----
        B : none
        R : none
        """
        self.region_sets([['chr1',5,10]],
                         [])
        result = self.setA.intersect(self.setB)
        self.assertEqual(len(result.sequences), 0)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.ORIGINAL)
        self.assertEqual(len(result.sequences), 0)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.COMP_INCL)
        self.assertEqual(len(result.sequences), 0)
        """
        A : none
        B :   -----
        R : none
        """
        self.region_sets([],
                         [['chr1',5,10]])
        result = self.setA.intersect(self.setB)
        self.assertEqual(len(result.sequences), 0)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.ORIGINAL)
        self.assertEqual(len(result.sequences), 0)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.COMP_INCL)
        self.assertEqual(len(result.sequences), 0)
        """
        No overlapping
        A : ------      ---------               ------- 
        B :        ----          ------  ------   
        R : none
        """
        self.region_sets([['chr1',1,5],['chr1',11,20],['chr1',33,38]],
                         [['chr1',7,9],['chr1',20,25],['chr1',26,31]])
        result = self.setA.intersect(self.setB)
        self.assertEqual(len(result.sequences), 0)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.ORIGINAL)
        self.assertEqual(len(result.sequences), 0)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.COMP_INCL)
        self.assertEqual(len(result.sequences), 0)
        """
        End-to-end attach
        A : ------      ------
        B :       ------
        R : none
        """
        self.region_sets([['chr1',1,5],['chr1',11,20]],
                         [['chr1',5,11]])
        result = self.setA.intersect(self.setB)
        self.assertEqual(len(result.sequences), 0)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.ORIGINAL)
        self.assertEqual(len(result.sequences), 0)

        result = self.setA.intersect(self.setB, mode=OverlapType.COMP_INCL)
        self.assertEqual(len(result.sequences), 0)
        """
        No length attach
        A : .      .
        B :    .   .
        R : none
        """
        self.region_sets([['chr1',2,2],['chr1',20,20]],
                         [['chr1',5,5],['chr1',20,20]])
        result = self.setA.intersect(self.setB)
        self.assertEqual(len(result.sequences), 1)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.ORIGINAL)
        self.assertEqual(len(result.sequences), 1)

        result = self.setA.intersect(self.setB, mode=OverlapType.COMP_INCL)
        self.assertEqual(len(result.sequences), 1)
        
        """
        Perfect overlapping
        A : ------
        B : ------
        R : ------
        """
        self.region_sets([['chr1',1,10],['chr1',500,550],['chr1',600,650],['chr1',700,750],['chr1',725,800]],
                         [['chr1',1,10],['chr1',500,550],['chr1',600,650],['chr1',700,750],['chr1',725,800]])
        result = self.setA.intersect(self.setB, mode=OverlapType.OVERLAP, rm_duplicates=True)
        
        self.assertEqual(len(result.sequences), 4)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 10)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.ORIGINAL)
        self.assertEqual(len(result.sequences), 5)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 10)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.COMP_INCL)
        self.assertEqual(len(result.sequences), 5)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 10)
        """
        One overlapping region
        A : ------
        B :     --------
        R1:     --       (overlap)
        R2: ------       (original)
        R3:              (comp_incl)
        """
        self.region_sets([['chr1',1,10]],
                         [['chr1',7,20]])
        result = self.setA.intersect(self.setB)
        self.assertEqual(len(result.sequences), 1)
        self.assertEqual(result.sequences[0].initial, 7)
        self.assertEqual(result.sequences[0].final, 10)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.ORIGINAL)
        self.assertEqual(len(result.sequences), 1)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 10)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.COMP_INCL)
        self.assertEqual(len(result.sequences), 0)
        """
        Two simple overlapping regions
        A : -------      --------
        B :     -------------
        R1:     ---      ----     (overlap)
        R2: -------      -------- (original)
        R3:                       (comp_incl)
        """
        self.region_sets([['chr1',1,10],['chr1',26,35]],
                         [['chr1',7,30]])
        result = self.setA.intersect(self.setB)
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 7)
        self.assertEqual(result.sequences[0].final, 10)
        self.assertEqual(result.sequences[1].initial, 26)
        self.assertEqual(result.sequences[1].final, 30)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.ORIGINAL)
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 10)
        self.assertEqual(result.sequences[1].initial, 26)
        self.assertEqual(result.sequences[1].final, 35)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.COMP_INCL)
        self.assertEqual(len(result.sequences), 0)
        """
        Two separately overlapping regions 
        A : -------      --------
        B :     -----        --------
        R1:     ---          ----     (overlap)
        R2: -------      --------     (original)
        R3:                           (comp_incl)
        """
        self.region_sets([['chr1',1,10],['chr1',26,35]],
                         [['chr1',7,15],['chr1',30,40]])
        result = self.setA.intersect(self.setB)
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 7)
        self.assertEqual(result.sequences[0].final, 10)
        self.assertEqual(result.sequences[1].initial, 30)
        self.assertEqual(result.sequences[1].final, 35)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.ORIGINAL)
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 10)
        self.assertEqual(result.sequences[1].initial, 26)
        self.assertEqual(result.sequences[1].final, 35)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.COMP_INCL)
        self.assertEqual(len(result.sequences), 0)
        """
        Many various overlapping (mixed)
        A :   ------------------            --------   ---------
        B : ----   -------    ------            ----------      
        R1:   --   -------    --                ----   ---       (overlap)
        R2:   ------------------            --------   --------- (original)
        R3:                                                      (comp_incl) 
        """
        self.region_sets([['chr1',3,30],['chr1',50,60],['chr1',70,85]],
                         [['chr1',1,5],['chr1',10,19],['chr1',27,35],['chr1',55,75]])
        
        result = self.setA.intersect(self.setB)
        self.assertEqual(len(result.sequences), 5)
        self.assertEqual(result.sequences[0].initial, 3)
        self.assertEqual(result.sequences[0].final, 5)
        self.assertEqual(result.sequences[1].initial, 10)
        self.assertEqual(result.sequences[1].final, 19)
        self.assertEqual(result.sequences[2].initial, 27)
        self.assertEqual(result.sequences[2].final, 30)
        self.assertEqual(result.sequences[3].initial, 55)
        self.assertEqual(result.sequences[3].final, 60)
        self.assertEqual(result.sequences[4].initial, 70)
        self.assertEqual(result.sequences[4].final, 75)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.ORIGINAL)
        self.assertEqual(len(result.sequences), 3)
        self.assertEqual(result.sequences[0].initial, 3)
        self.assertEqual(result.sequences[0].final, 30)
        self.assertEqual(result.sequences[1].initial, 50)
        self.assertEqual(result.sequences[1].final, 60)
        self.assertEqual(result.sequences[2].initial, 70)
        self.assertEqual(result.sequences[2].final, 85)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.COMP_INCL)
        self.assertEqual(len(result.sequences), 0)
        """
        Different chromosomes
        A : chr1  -------
        B : chr2  -------
        R : none
        """
        self.region_sets([['chr1',1,10]],
                         [['chr2',1,10]])
        result = self.setA.intersect(self.setB)
        self.assertEqual(len(result.sequences), 0)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.ORIGINAL)
        self.assertEqual(len(result.sequences), 0)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.COMP_INCL)
        self.assertEqual(len(result.sequences), 0)
        """
        Completely included overlapping
        A : ---------------------------
        B : ----    ------       -----------
        R1: ----    ------       ------      (overlap)
        R2: ---------------------------      (original)
        R3:                                  (comp_incl)
        """
        self.region_sets([['chr1',1,50]],
                         [['chr1',1,5],['chr1',10,19],['chr1',45,60]])
        result = self.setA.intersect(self.setB)
        self.assertEqual(len(result.sequences), 3)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 5)
        self.assertEqual(result.sequences[1].initial, 10)
        self.assertEqual(result.sequences[1].final, 19)
        self.assertEqual(result.sequences[2].initial, 45)
        self.assertEqual(result.sequences[2].final, 50)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.ORIGINAL)
        self.assertEqual(len(result.sequences), 1)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 50)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.COMP_INCL)
        self.assertEqual(len(result.sequences), 0)
        """
        A : ----    ------       -----------
        B : ---------------------------
        R1: ----    ------       ------      (overlap)
        R2: ----    ------       ----------- (original)
        R3: ----    ------                   (comp_incl)
        """
        self.region_sets([['chr1',1,5],['chr1',10,19],['chr1',45,60]],
                         [['chr1',1,50]])
        result = self.setA.intersect(self.setB)
        self.assertEqual(len(result.sequences), 3)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 5)
        self.assertEqual(result.sequences[1].initial, 10)
        self.assertEqual(result.sequences[1].final, 19)
        self.assertEqual(result.sequences[2].initial, 45)
        self.assertEqual(result.sequences[2].final, 50)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.ORIGINAL)
        self.assertEqual(len(result.sequences), 3)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 5)
        self.assertEqual(result.sequences[1].initial, 10)
        self.assertEqual(result.sequences[1].final, 19)
        self.assertEqual(result.sequences[2].initial, 45)
        self.assertEqual(result.sequences[2].final, 60)
        
        
        result = self.setA.intersect(self.setB, mode=OverlapType.COMP_INCL)
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 5)
        self.assertEqual(result.sequences[1].initial, 10)
        self.assertEqual(result.sequences[1].final, 19)

        """
        A : --------------         -------
                ------
        B :       -----          ----------------
        R1:       -----            -------      (overlap)
                  ----
        R2: --------------         -------      (original)
                ------
        R3:                        -------      (comp_incl)
        """
        self.region_sets([['chr1',1,50],['chr1',20,40],['chr1',70,80]],
                         [['chr1',25,45],['chr1',65,95]])
        result = self.setA.intersect(self.setB)
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 25)
        self.assertEqual(result.sequences[0].final, 45)
        self.assertEqual(result.sequences[1].initial, 70)
        self.assertEqual(result.sequences[1].final, 80)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.ORIGINAL)
        self.assertEqual(len(result.sequences), 3)
        self.assertEqual(result.sequences[1].initial, 1)
        self.assertEqual(result.sequences[1].final, 50)
        self.assertEqual(result.sequences[0].initial, 20)
        self.assertEqual(result.sequences[0].final, 40)
        self.assertEqual(result.sequences[2].initial, 70)
        self.assertEqual(result.sequences[2].final, 80)
        
        result = self.setA.intersect(self.setB, mode=OverlapType.COMP_INCL)
        self.assertEqual(len(result.sequences), 1)
        self.assertEqual(result.sequences[0].initial, 70)
        self.assertEqual(result.sequences[0].final, 80)

    def test_closest(self):
        """
        Two empty sets
        A : none 
        B : none
        R : none
        """
        self.region_sets([],
                         [])
        result = self.setA.closest(self.setB)
        self.assertEqual(len(result.sequences), 0)
        """
        One empty set
        A :   -----
        B : none
        R : none
        """
        self.region_sets([['chr1',5,10]],
                         [])
        result = self.setA.closest(self.setB)
        self.assertEqual(len(result.sequences), 0)
        """
        A : none
        B :   -----
        R : none
        """
        self.region_sets([],
                         [['chr1',5,10]])
        result = self.setA.closest(self.setB)
        self.assertEqual(len(result.sequences), 0)
        """
        Overlapping within set
        A : -----====-----
        B :      ----
        R : none
        """
        self.region_sets([['chr1',1,10],['chr1',6,15]],
                         [['chr1',6,10]])
        result = self.setA.closest(self.setB)
        self.assertEqual(result, False)
        """
        A :      ----
        B : -----====-----
        R : none
        """
        self.region_sets([['chr1',6,10]],
                         [['chr1',1,10],['chr1',6,15]])
        result = self.setA.closest(self.setB)
        self.assertEqual(result, False)
        """
        No overlapping
        A : ------      ---------               ------- 
        B :        ----          ------  ------   
        R :                      ------
        """
        self.region_sets([['chr1',1,5],['chr1',11,20],['chr1',33,38]],
                         [['chr1',7,9],['chr1',20,25],['chr1',26,31]])
        result = self.setA.closest(self.setB)
        self.assertEqual(len(result.sequences), 1)
        self.assertEqual(result.sequences[0].initial, 20)
        self.assertEqual(result.sequences[0].final, 25)
        """
        End-to-end attach
        A : ------      ------
        B :       ------
        R :       ------
        """
        self.region_sets([['chr1',1,5],['chr1',11,20]],
                         [['chr1',5,11]])
        result = self.setA.closest(self.setB)
        self.assertEqual(len(result.sequences), 1)
        self.assertEqual(result.sequences[0].initial, 5)
        self.assertEqual(result.sequences[0].final, 11)
        """
        Perfect overlapping
        A : ------
        B : ------
        R : none
        """
        self.region_sets([['chr1',1,10]],
                         [['chr1',1,10]])
        result = self.setA.closest(self.setB)
        self.assertEqual(result, False)
        """
        One overlapping region
        A : ------
        B :     --------
        R : none
        """
        self.region_sets([['chr1',1,10]],
                         [['chr1',7,20]])
        result = self.setA.closest(self.setB)
        self.assertEqual(result, False)
        """
        Two simple overlapping regions
        A : -------      --------
        B :     -------------
        R : none
        """
        self.region_sets([['chr1',1,10],['chr1',26,35]],
                         [['chr1',7,30]])
        result = self.setA.closest(self.setB)
        self.assertEqual(result, False)
        """
        Two separately overlapping regions 
        A : -------      --------
        B :     -----        --------
        R : none
        """
        self.region_sets([['chr1',1,10],['chr1',26,35]],
                         [['chr1',7,15],['chr1',30,40]])
        result = self.setA.closest(self.setB)
        self.assertEqual(result, False)
        """
        Many various overlapping (mixed)
        A :   ------------------            --------   ---------
        B : ----   -------    ------            ----------      
        R : none
        """
        self.region_sets([['chr1',3,30],['chr1',50,60],['chr1',70,85]],
                         [['chr1',1,5],['chr1',10,19],['chr1',27,35],['chr1',55,75]])
        result = self.setA.closest(self.setB)
        self.assertEqual(result, False)
        """
        Different chromosomes
        A : chr1  -------
        B : chr2  -------
        R : chr2  -------
        
        """
        self.region_sets([['chr1',1,10]],
                         [['chr2',1,10]])
        result = self.setA.closest(self.setB)
        self.assertEqual(len(result.sequences), 0)
        """
        Completely included overlapping
        A : ---------------------------
        B : ----    ------       -----------
        R : none
        """
        self.region_sets([['chr1',1,50]],
                         [['chr1',1,5],['chr1',10,19],['chr1',45,60]])
        result = self.setA.closest(self.setB)
        self.assertEqual(result, False)
        """
        A : ----    ------       -----------
        B : ---------------------------
        R : none
        """
        self.region_sets([['chr1',1,5],['chr1',10,19],['chr1',45,60]],
                         [['chr1',1,50]])
        result = self.setA.closest(self.setB)
        self.assertEqual(result, False)
        """
        A : ----         ------                  ---
        B :        ---              -----
        R :        ---
        """
        self.region_sets([['chr1',1,5],['chr1',27,45],['chr1',85,95]],
                         [['chr1',15,20],['chr1',55,65]])
        result = self.setA.closest(self.setB)
        self.assertEqual(len(result.sequences), 1)
        self.assertEqual(result.sequences[0].initial, 15)
        self.assertEqual(result.sequences[0].final, 20)
    
    def test_remove_duplicates(self):
        """
        A : ===== -----
        R : ----- -----
        """
        self.region_sets([['chr1',1,10],['chr1',1,10],['chr1',15,25]],
                         [])
        self.setA.remove_duplicates()
        result = self.setA
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 10)
        self.assertEqual(result.sequences[1].initial, 15)
        self.assertEqual(result.sequences[1].final, 25)
        """
        A : =====--- -----
        R : =====--- -----
        """
        self.region_sets([['chr1',1,10],['chr1',1,15],['chr1',20,25]],
                         [])
        self.setA.remove_duplicates()
        result = self.setA
        self.assertEqual(len(result.sequences), 3)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 10)
        self.assertEqual(result.sequences[1].initial, 1)
        self.assertEqual(result.sequences[1].final, 15)
        self.assertEqual(result.sequences[2].initial, 20)
        self.assertEqual(result.sequences[2].final, 25)
        """
        A : ===== ----- ------  ====
        R : ----- ----- ------  ----
        """
        self.region_sets([['chr1',1,10],['chr1',1,10],['chr1',15,25],['chr1',30,35],['chr1',40,45],['chr1',40,45]],
                         [])
        self.setA.remove_duplicates()
        result = self.setA
        self.assertEqual(len(result.sequences), 4)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 10)
        self.assertEqual(result.sequences[1].initial, 15)
        self.assertEqual(result.sequences[1].final, 25)
        self.assertEqual(result.sequences[2].initial, 30)
        self.assertEqual(result.sequences[2].final, 35)
        self.assertEqual(result.sequences[3].initial, 40)
        self.assertEqual(result.sequences[3].final, 45)

    def test_window(self):
        """
        A :             -------
        B : ------[ 99 ]       [   199   ]---
        window = 100
        R :       -                           only one base overlaps with extending A
        """   
        self.region_sets([['chr1',200,300]],
                         [['chr1',1,101],['chr1',499,550]])
        result = self.setA.window(self.setB,adding_length=100)
        self.assertEqual(len(result.sequences), 1)
        self.assertEqual(result.sequences[0].initial, 100)
        self.assertEqual(result.sequences[0].final, 101)
        """
        A :             -------
        B : ------[ 99 ]       [   199   ]---
        window = 200
        R : ------                        -   
        left-hand side is covered, and the right-hand side is only one base overlapped
        """   
        self.region_sets([['chr1',200,300]],
                         [['chr1',1,101],['chr1',499,550]])
        result = self.setA.window(self.setB,adding_length=200)
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 1)  # GenomicRegion.extend will choose 1 rather than 0
        self.assertEqual(result.sequences[0].final, 101)
        self.assertEqual(result.sequences[1].initial, 499)
        self.assertEqual(result.sequences[1].final, 500)
        """
        A :                         ----    ----
        B :             --------                    ----
        window = 1000 (default)
        R :                 ----                    ----
        """   
        self.region_sets([['chr1',3000,3500],['chr1',4000,4500]],
                         [['chr1',1500,2500],['chr1',5000,5500]])
        result = self.setA.window(self.setB)
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 2000)
        self.assertEqual(result.sequences[0].final, 2500)
        self.assertEqual(result.sequences[1].initial, 5000)
        self.assertEqual(result.sequences[1].final, 5500)
        """
        A :                         ----    ----
        B :             --------                    ----
        window = 2000
        R :             --------                    ----
                            ----                    ----
        window = 100
        R : none
        """   
        self.region_sets([['chr1',3000,3500],['chr1',4000,4500]],
                         [['chr1',1500,2500],['chr1',5000,5500]])
        result = self.setA.window(self.setB,adding_length=2000)
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 1500)
        self.assertEqual(result.sequences[0].final, 2500)
        self.assertEqual(result.sequences[1].initial, 5000)
        self.assertEqual(result.sequences[1].final, 5500)
        result = self.setA.window(self.setB,adding_length=100)
        self.assertEqual(len(result.sequences), 0)
        
    def test_subtract(self):
        """
        A : none
        B :    ------
        R : none
        """
        self.region_sets([],
                         [['chr1',6,15]])
        result = self.setA.subtract(self.setB)
        self.assertEqual(len(result.sequences), 0)
        """
        A :    ------
        B : none
        R :    ------
        """
        self.region_sets([['chr1',6,15]],
                         [])
        result = self.setA.subtract(self.setB)
        self.assertEqual(len(result.sequences), 1)
        self.assertEqual(result.sequences[0].initial, 6)
        self.assertEqual(result.sequences[0].final, 15)
        """
        A : ------
        B :    ------
        R : ---
        """
        self.region_sets([['chr1',1,10]],
                         [['chr1',6,15]])
        result = self.setA.subtract(self.setB)
        self.assertEqual(len(result.sequences), 1)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 6)
        """
        A :    ------
        B : ------
        R :       ---
        """
        self.region_sets([['chr1',6,15]],
                         [['chr1',1,10]])
        result = self.setA.subtract(self.setB)
        self.assertEqual(len(result.sequences), 1)
        self.assertEqual(result.sequences[0].initial, 10)
        self.assertEqual(result.sequences[0].final, 15)
        """
        A :    ---
        B : ---------
        R : none
        """
        self.region_sets([['chr1',6,10]],
                         [['chr1',1,15]])
        result = self.setA.subtract(self.setB)
        self.assertEqual(len(result.sequences), 0)
        """
        A : ---------
        B :    ---
        R : ---   ---
        """
        self.region_sets([['chr1',1,15]],
                         [['chr1',6,10]])
        result = self.setA.subtract(self.setB)
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 6)
        self.assertEqual(result.sequences[1].initial, 10)
        self.assertEqual(result.sequences[1].final, 15)
        """
        A :    ------
        B :    ------
        R : none
        """
        self.region_sets([['chr1',6,15]],
                         [['chr1',6,15]])
        result = self.setA.subtract(self.setB)
        self.assertEqual(len(result.sequences), 0)
        """
        A :   ----------              ------
        B :          ----------                    ----
        R :   -------                 ------
        """
        self.region_sets([['chr1',5,30],['chr1',70,85]],
                         [['chr1',20,50],['chr1',100,110]])
        result = self.setA.subtract(self.setB)
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 5)
        self.assertEqual(result.sequences[0].final, 20)
        self.assertEqual(result.sequences[1].initial, 70)
        self.assertEqual(result.sequences[1].final, 85)
        """
        A :        ------   -----
        B :    ------
        R :          ----   -----
        """
        self.region_sets([['chr1',20,30],['chr1',35,55]],
                         [['chr1',10,23],['chr1',100,110]])
        result = self.setA.subtract(self.setB)
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 23)
        self.assertEqual(result.sequences[0].final, 30)
        self.assertEqual(result.sequences[1].initial, 35)
        self.assertEqual(result.sequences[1].final, 55)
        """
        A :   ch1     ---------------------
              ch2     -------------------------
        B :   ch1             ------
              ch2                        ------
        R :   ch1     --------      -------
              ch2     -------------------
        """
        self.region_sets([['chr1',0,30000],['chr2',0,35000]],
                         [['chr1',20000,23000],['chr2',31000,35000]])
        result = self.setA.subtract(self.setB)
        self.assertEqual(len(result.sequences), 3)
        self.assertEqual(result.sequences[0].initial, 0)
        self.assertEqual(result.sequences[0].final, 20000)
        self.assertEqual(result.sequences[1].initial, 23000)
        self.assertEqual(result.sequences[1].final, 30000)
        self.assertEqual(result.sequences[2].initial, 0)
        self.assertEqual(result.sequences[2].final, 31000)
        """
        A :   -----------------------------------------------------------
        B :    ---    ---------         ----           ----
        R :   -   ----         ---------    -----------    --------------
        """
        self.region_sets([['chr1',5,1000]],
                         [['chr1',10,15],['chr1',30,70],['chr1',120,140],['chr1',200,240]])
        result = self.setA.subtract(self.setB)
        self.assertEqual(len(result.sequences), 5)
        
        """
        A :   -----------------------              ------
                   -----     -----  -----------
        B :    ---    ---------         ----           ----
        R :   -   ----         ------              ----
                   ---         ---  ----    ---
        """
        self.region_sets([['chr1',5,100],['chr1',20,40],['chr1',60,80],['chr1',95,150],['chr1',180,220]],
                         [['chr1',10,15],['chr1',30,70],['chr1',120,140],['chr1',200,240]])
        result = self.setA.subtract(self.setB)
        #print(result.sequences)
        self.assertEqual(len(result.sequences), 8)
        self.assertEqual(result.sequences[0].initial, 5)
        """
        A :   -----------------------------------------------------------
        B :    ---    ---------         ----           ----
        R :   -   ----         ---------    -----------    --------------
        """
        self.region_sets([['chr1',5,1000],['chr2',5,1000],['chr4',5,1000]],
                         [['chr1',10,15],['chr1',30,70],['chr1',120,140],['chr1',200,240],
                          ['chr2',10,15],['chr2',30,70],['chr2',120,140],['chr2',200,240],
                          ['chr4',10,15],['chr4',30,70],['chr4',120,140],['chr4',200,240]])
        result = self.setA.subtract(self.setB)
        self.assertEqual(len(result.sequences), 15)
        
        
    def test_merge(self):
        """
        A : none
        R : none
        """
        self.region_sets([],
                         [])
        self.setA.merge()
        result = self.setA
        self.assertEqual(len(result.sequences), 0)
        """
        A : -----  -----
        R : -----  -----
        """
        self.region_sets([['chr1',1,10],['chr1',15,25]],
                         [])
        self.setA.merge()
        result = self.setA
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 10)
        self.assertEqual(result.sequences[1].initial, 15)
        self.assertEqual(result.sequences[1].final, 25)
        """
        A1: ------------   ----
        A2:    -----
        R : ------------   ----
        """
        self.region_sets([['chr1',1,30],['chr1',11,20],['chr1',40,50]],
                         [])
        self.setA.merge()
        result = self.setA
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 30)
        self.assertEqual(result.sequences[1].initial, 40)
        self.assertEqual(result.sequences[1].final, 50)
        """
        A1: --------       ----
        A2:    ---------
        R : ------------   ----
        """
        self.region_sets([['chr1',1,30],['chr1',20,40],['chr1',50,60]],
                         [])
        self.setA.merge()
        result = self.setA
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 40)
        self.assertEqual(result.sequences[1].initial, 50)
        self.assertEqual(result.sequences[1].final, 60)
        """
        A : =======
        R : -------
        """
        self.region_sets([['chr1',1,30],['chr1',1,30]],
                         [])
        self.setA.merge()
        result = self.setA
        self.assertEqual(len(result.sequences), 1)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 30)
        
    def test_cluster(self):
        """
        Empty sets
        A : none 
        R : none
        """
        self.region_sets([],
                         [])
        result = self.setA.cluster(10)
        self.assertEqual(len(result.sequences), 0)
        """
        A :  ------- 
        R :  -------
        """
        self.region_sets([['chr1',1,10]],
                         [])
        result = self.setA.cluster(10)
        self.assertEqual(len(result.sequences), 1)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 10)
        """
        A :  -----
                  ------
        R :  -----------
        """
        self.region_sets([['chr1',1,10],['chr1',10,20]],
                         [])
        result = self.setA.cluster(10)
        self.assertEqual(len(result.sequences), 1)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 20)
        """
        A :  -----  -----
        R1:  -----  -----
        R2:  ------------
        """
        self.region_sets([['chr1',1,10],['chr1',15,25]],
                         [])
        result = self.setA.cluster(1)
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 10)
        self.assertEqual(result.sequences[1].initial, 15)
        self.assertEqual(result.sequences[1].final, 25)
        result = self.setA.cluster(5)
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 10)
        self.assertEqual(result.sequences[1].initial, 15)
        self.assertEqual(result.sequences[1].final, 25)
        result = self.setA.cluster(6)
        self.assertEqual(len(result.sequences), 1)
        self.assertEqual(result.sequences[0].initial, 1)
        self.assertEqual(result.sequences[0].final, 25)
        """
        A :  ---- ----  ----   ----    ----
        R1:  ---------  ----   ----    ----
        R2:  ---------------   ----    ----
        R3:  ----------------------    ----
        R4:  ------------------------------
        R5:  ------------------------------
        """
        self.region_sets([['chr1',1,10],['chr1',15,25],['chr1',35,45],
                          ['chr1',60,70],['chr1',90,100]],
                         [])
        result = self.setA.cluster(6)
        self.assertEqual(len(result.sequences), 4)
        result = self.setA.cluster(11)
        self.assertEqual(len(result.sequences), 3)
        result = self.setA.cluster(16)
        self.assertEqual(len(result.sequences), 2)
        result = self.setA.cluster(21)
        self.assertEqual(len(result.sequences), 1)
        result = self.setA.cluster(26)
        self.assertEqual(len(result.sequences), 1)
        
    def test_flank(self):
        """
        A :        -----
        R1:     ---     ---
        """
        self.region_sets([['chr1',60,75]],
                         [])
        result = self.setA.flank(10)
        self.assertEqual(len(result.sequences), 2)
        self.assertEqual(result.sequences[0].initial, 50)
        self.assertEqual(result.sequences[0].final, 60)
        self.assertEqual(result.sequences[1].initial, 75)
        self.assertEqual(result.sequences[1].final, 85)
        """
        A :        -----     ----
        R1:   -----     =====    ----
        """
        self.region_sets([['chr1',60,75],['chr1',90,100]],
                         [])
        result = self.setA.flank(15)
        self.assertEqual(len(result.sequences), 4)
        self.assertEqual(result.sequences[0].initial, 45)
        self.assertEqual(result.sequences[0].final, 60)
        self.assertEqual(result.sequences[1].initial, 75)
        self.assertEqual(result.sequences[1].final, 90)
        self.assertEqual(result.sequences[2].initial, 75)
        self.assertEqual(result.sequences[2].final, 90)
        self.assertEqual(result.sequences[3].initial, 100)
        self.assertEqual(result.sequences[3].final, 115)
        
    def test_jaccard(self):
        """
        self           --8--      ---10---      -4-
        y         ---10---             ---10---
        intersect      -5-             -4-    
        similarity:   ( 5 + 4 )/[(8 + 10 + 4) + (10 +10) - (5 + 4 )]
                      = 9/33
        """
        self.region_sets([['chr1',50,58],['chr1',70,80],['chr1',90,94]],
                         [['chr1',45,55],['chr1',76,86]])
        result = self.setA.jaccard(self.setB)
        self.assertEqual(result, 9/33)
    
    def test_get_genome_data(self):
        """hg19"""
        result = GenomicRegionSet("hg19")
        result.get_genome_data(organism="hg19")
        self.assertEqual(len(result.sequences), 23)
        """hg19, with Mitochondria chromosome"""
        result = GenomicRegionSet("hg19")
        result.get_genome_data(organism="hg19",chrom_M=True)
        self.assertEqual(len(result.sequences), 24)
        
    def test_random_regions(self):
        
        self.region_sets([['chr1',0,10000],['chr2',0,20000],['chrX',0,30000]],
                         [])
        result = self.setA.random_regions(organism="mm9", 
                                          total_size=100, 
                                          overlap_result=False, 
                                          overlap_input=False)
        result.sort()
        #print("-"*80)
        #print("The result random regions are: ")
        #for s in result.sequences:
        #    print("\t%s\t%10d\t%10d%10d" % (s.chrom,s.initial,s.final,s.__len__()))
        #print("Overlaps within result: ",result.within_overlap())
        
        
        self.region_sets([['chr1',0,10000],['chr2',0,20000],['chrX',0,30000]],
                         [])
        result = self.setA.random_regions(organism="mm9", 
                                          total_size=100, 
                                          overlap_result=True, 
                                          overlap_input=False)
        result.sort()
        #print("-"*80)
        #print("The result random regions are: ")
        #for s in result.sequences:
        #    print("\t%s\t%10d\t%10d%10d" % (s.chrom,s.initial,s.final,s.__len__()))
        #print("Overlaps within result: ",result.within_overlap())
        
        self.region_sets([['chr1',0,10000],['chr2',0,20000],['chrX',0,30000]],
                         [])
        result = self.setA.random_regions(organism="mm9", 
                                          total_size=100, 
                                          overlap_result=False, 
                                          overlap_input=True)
        result.sort()
        #print("-"*80)
        #print("The result random regions are: ")
        #for s in result.sequences:
        #    print("\t%s\t%10d\t%10d%10d" % (s.chrom,s.initial,s.final,s.__len__()))
        #print("Overlaps within result: ",result.within_overlap())
        
        self.region_sets([['chr1',0,10000],['chr2',0,20000],['chrX',0,30000]],
                         [])
        result = self.setA.random_regions(organism="mm9", 
                                          total_size=100, 
                                          overlap_result=True, 
                                          overlap_input=True)
        result.sort()
        #print("-"*80)
        #print("The result random regions are: ")
        #for s in result.sequences:
        #    print("\t%s\t%10d\t%10d%10d" % (s.chrom,s.initial,s.final,s.__len__()))
        #print("Overlaps within result: ",result.within_overlap())
        
        self.region_sets([['chr1',0,1000],['chr2',0,2000],['chrX',0,3000]],
                         [])
        result = self.setA.random_regions(organism="mm9", 
                                          multiply_factor=100, 
                                          overlap_result=False, 
                                          overlap_input=False)
        result.sort()
        #print("-"*80)
        #print("The result random regions are: ")
        #for s in result.sequences:
        #    print("\t%s\t%10d\t%10d%10d" % (s.chrom,s.initial,s.final,s.__len__()))
        #print("Overlaps within result: ",result.within_overlap())
        
        self.region_sets([['chr1',0,1000],['chr2',0,2000],['chrX',0,3000]],
                         [])
        result = self.setA.random_regions(organism="mm9", 
                                          multiply_factor=100, 
                                          overlap_result=False, 
                                          overlap_input=False,
                                          chrom_M=True)
        result.sort()
        #print("-"*80)
        #print("The result random regions are: ")
        #for s in result.sequences:
        #    print("\t%s\t%10d\t%10d%10d" % (s.chrom,s.initial,s.final,s.__len__()))
        #print("Overlaps within result: ",result.within_overlap())

"""
    
    def test_projection_test(self):
        
        self.region_sets([['chr1',0,10000000],['chr2',0,20000000],['chr3',0,30000000]],
                         [['chr1',10,11],['chr2',10,11],['chr3',10,11]])
        result = self.setA.projection_test(self.setB)
        #print(result)
        #self.assertEqual(result, 11/31)
"""     
        
if __name__ == "__main__":

    unittest.main()