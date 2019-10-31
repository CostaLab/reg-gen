
import os.path
from ctypes import *
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.GenomicRegion import GenomicRegion

me = os.path.abspath(os.path.dirname(__file__))
lib = cdll.LoadLibrary(os.path.join(me, "librgt.so"))

# Determine whether GenomicRegions overlap
overlap_c = lib.overlap
overlap_c.argtypes = [POINTER(c_char), c_int, c_int, POINTER(c_char), c_int, c_int]
overlap_c.restype = c_bool

# Compare GenomicRegions
compare_c = lib.compareGenomicRegions
compare_c.argtypes = [POINTER(c_char), c_int, c_int, POINTER(c_char), c_int, c_int]
compare_c.restype = c_int

# Compute jaccard index
jaccardC = lib.jaccard
jaccardC.argtypes = [POINTER(c_char_p), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_char_p), POINTER(c_int), POINTER(c_int), c_int]
jaccardC.restype = c_double

# Intersect genomic regions
intersect_overlap_c = lib.intersectGenomicRegionSetsOverlap
intersect_overlap_c.argtypes = [POINTER(c_char_p), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_char_p), POINTER(c_int), POINTER(c_int), c_int, POINTER(POINTER(c_int)), POINTER(POINTER(c_int)), POINTER(POINTER(c_int)), POINTER(c_int)]
intersect_overlap_c.restype = None

intersect_original_c = lib.intersectGenomicRegionSetsOriginal
intersect_original_c.argtypes = [POINTER(c_char_p), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_char_p), POINTER(c_int), POINTER(c_int), c_int, POINTER(POINTER(c_int)), POINTER(POINTER(c_int)), POINTER(POINTER(c_int)), POINTER(c_int)]
intersect_original_c.restype = None

intersect_completely_included_c = lib.intersectGenomicRegionSetsCompletelyIncluded
intersect_completely_included_c.argtypes = [POINTER(c_char_p), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_char_p), POINTER(c_int), POINTER(c_int), c_int, POINTER(POINTER(c_int)), POINTER(POINTER(c_int)), POINTER(POINTER(c_int)), POINTER(c_int)]
intersect_completely_included_c.restype = None

chromA = c_char_p("chr1")
initialA = c_int(0)
finalA = c_int(2)

chromB = c_char_p("chr1")
initialB = c_int(5)
finalB = c_int(7)

chromC = "chr1"
initialC = c_int(0)
finalC = c_int(4)

chromD = "chr1"
initialD = c_int(5)
finalD = c_int(6)

overlapping = overlap_c(chromA, initialA, finalA, chromB, initialB, finalB)
print("Overlapping? ", overlapping)

compare = compare_c(chromA, initialA, finalA, chromB, initialB, finalB)
print("Compare = ", compare)


def jaccardIndex(gnrsA, gnrsB):
    # Convert to ctypes
    chroms = [gr.chrom for gr in gnrsA.sequences]
    chromsA = (c_char_p * len(chroms))(*chroms)

    chroms = [gr.chrom for gr in gnrsB.sequences]
    chromsB = (c_char_p * len(chroms))(*chroms)

    ints = [gr.initial for gr in gnrsA.sequences]
    initialsA = (c_int * len(ints))(*ints)

    ints = [gr.initial for gr in gnrsB.sequences]
    initialsB = (c_int * len(ints))(*ints)

    ints = [gr.final for gr in gnrsA.sequences]
    finalsA = (c_int * len(ints))(*ints)

    ints = [gr.final for gr in gnrsB.sequences]
    finalsB = (c_int * len(ints))(*ints)

    # Call C-function
    return jaccardC(chromsA, initialsA, finalsA, len(gnrsA), chromsB, initialsB, finalsB, len(gnrsB))

set1 = GenomicRegionSet("A")
set1.add(GenomicRegion("chr1", 0, 10))
set1.add(GenomicRegion("chr1", 15, 20))
set1.add(GenomicRegion("chr1", 30, 45))
print(set1.sequences)
set2 = GenomicRegionSet("B")
set2.add(GenomicRegion("chr1", 0, 5))
set2.add(GenomicRegion("chr1", 10, 25))
set2.add(GenomicRegion("chr1", 35, 45))
print(set2.sequences)

jaccard2 = jaccardIndex(set1, set2)
print("jaccard2", jaccard2)


def intersect(gnrsA, gnrsB, overlap_type):
    # Convert to ctypes
    lenA = len(gnrsA)
    lenB = len(gnrsB)
    lenR = min(lenA, lenB)

    chromsA_python = [gr.chrom for gr in gnrsA.sequences]
    chromsA_c = (c_char_p * lenA)(*chromsA_python)

    chromsB_python = [gr.chrom for gr in gnrsB.sequences]
    chromsB_c = (c_char_p * lenB)(*chromsB_python)

    initialsA_python = [gr.initial for gr in gnrsA.sequences]
    initialsA_c = (c_int * lenA)(*initialsA_python)

    initialsB_python = [gr.initial for gr in gnrsB.sequences]
    initialsB_c = (c_int * lenB)(*initialsB_python)

    finalsA_python = [gr.final for gr in gnrsA.sequences]
    finalsA_c = (c_int * lenA)(*finalsA_python)

    finalsB_python = [gr.final for gr in gnrsB.sequences]
    finalsB_c = (c_int * lenB)(*finalsB_python)

    indices_c = POINTER(c_int)((c_int * lenR)())
    initialsR_c = POINTER(c_int)((c_int * lenR)())
    finalsR_c = POINTER(c_int)((c_int * lenR)())
    sizeR_c = c_int()


    # Call C-function
    if overlap_type == 0:
        intersect_overlap_c(chromsA_c, initialsA_c, finalsA_c, lenA, chromsB_c, initialsB_c, finalsB_c, lenB,
                            pointer(indices_c), pointer(initialsR_c), pointer(finalsR_c), byref(sizeR_c))
    elif overlap_type == 1:
        intersect_original_c(chromsA_c, initialsA_c, finalsA_c, lenA, chromsB_c, initialsB_c, finalsB_c, lenB,
                             pointer(indices_c), pointer(initialsR_c), pointer(finalsR_c), byref(sizeR_c))
    elif overlap_type == 2:
        intersect_completely_included_c(chromsA_c, initialsA_c, finalsA_c, lenA, chromsB_c, initialsB_c, finalsB_c, lenB,
                                        pointer(indices_c), pointer(initialsR_c), pointer(finalsR_c), byref(sizeR_c))

    result = GenomicRegionSet(gnrsA.name)
    for i in range(sizeR_c.value):
        result.add(GenomicRegion(chromsA_python[indices_c[i]], initialsR_c[i], finalsR_c[i]))

    return result

set1 = GenomicRegionSet("A")
set1.add(GenomicRegion("chr1", 0, 10))
set1.add(GenomicRegion("chr1", 15, 20))
set1.add(GenomicRegion("chr1", 30, 45))
print(set1.sequences)
set2 = GenomicRegionSet("B")
set2.add(GenomicRegion("chr1", 0, 5))
set2.add(GenomicRegion("chr1", 10, 25))
set2.add(GenomicRegion("chr1", 35, 45))
print(set2.sequences)

for mode in (0, 1, 2):
    inter = intersect(set1, set2, mode)
    print("intersect, mode =", mode, ": ", str(inter.sequences))

'''
http://johnstowers.co.nz/2011/07/15/interfacing-python-c-opencv-via-ctypes/
FROM EXAMPLE CODE:

func = lib.test_get_data_nulls
func.restype = POINTER(c_char)
func.argtypes = [POINTER(c_int)]

l = c_int()
data = func(byref(l))

print(data, l, data.contents)

lib.test_data_print(data,l)

func_out = lib.test_get_data_nulls_out
func_out.argtypes = [POINTER(POINTER(c_char)), POINTER(c_int)]
func.restype = None

l2 = c_int()
data2 = POINTER(c_char)()
func_out(byref(data2), byref(l2))

print(data2, l2, data2.contents)

lib.test_data_print(data2, l2)

print("equal ", data[0] == data2[0], data[1] == data2[1], data[2] == data2[2], data[3] == data2[3], data[4] == data2[4])

func = lib.test_get_fixed_array_size_2
func.argtypes = [POINTER(c_double)]
func.restype = None

data = (c_double * 2)()
func(data)
x,y = data
print("array ", x, y)
'''
