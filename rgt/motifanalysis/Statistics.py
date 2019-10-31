
###################################################################################################
# Libraries
###################################################################################################

# Python 3 compatibility


import time

# Internal
from ..GeneSet import GeneSet
from ..Util import OverlapType
from ..GenomicRegionSet import GenomicRegionSet

# External
from numpy import asarray, argsort, sum, arange, nonzero, minimum

###################################################################################################
# Functions
###################################################################################################


def ecdf(x):
    """
    Auxiliary function for multiple_test_correction
    """
    nobs = len(x)
    return arange(1, nobs+1) / float(nobs)


def multiple_test_correction(pvals, alpha=0.05, method='indep'):
    """ 
    p-value correction for false discovery rate.
   
    This covers Benjamini/Hochberg for independent or positively correlated and
    Benjamini/Yekutieli for general or negatively correlated tests. Both are
    available in the function multipletests, as method=`fdr_bh`, resp. `fdr_by`.

    If there is prior information on the fraction of true hypothesis, then alpha
    should be set to alpha * m/m_0 where m is the number of tests,
    given by the p-values, and m_0 is an estimate of the true hypothesis.
    (see Benjamini, Krieger and Yekuteli)

    The two-step method of Benjamini, Krieger and Yekutiel that estimates the number
    of false hypotheses will be available (soon).

    Method names can be abbreviated to first letter, 'i' or 'p' for fdr_bh and 'n' for
    fdr_by.

    Author: Josef Pktd, H Raja and Vincent Davis (scikits.statsmodels.sandbox.stats.multicomp)

    Keyword arguments:
    pvals -- List of p-values from the individual tests.
    alpha -- Error rate (float). (default 0.05)
    method -- {'indep', 'negcorr')
        
    Return:
    rejected -- List of booleans. True if a hypothesis is rejected, False otherwise.
    pvalue_corrected -- A list with the p-values adjusted for multiple hypothesis testing to limit FDR.
    """

    pvals = asarray(pvals)

    pvals_sortind = argsort(pvals)
    pvals_sorted = pvals[pvals_sortind]
    sortrevind = pvals_sortind.argsort()

    if method in ['i', 'indep', 'p', 'poscorr']:
        ecdffactor = ecdf(pvals_sorted)
    elif method in ['n', 'negcorr']:
        cm = sum(1./arange(1, len(pvals_sorted)+1))
        ecdffactor = ecdf(pvals_sorted) / cm
    else:
        raise ValueError('only indep and negcorr implemented')
    reject = pvals_sorted < ecdffactor*alpha
    if reject.any():
        rejectmax = max(nonzero(reject)[0])
    else:
        rejectmax = 0
    reject[:rejectmax] = True

    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    pvals_corrected[pvals_corrected > 1] = 1
    return reject[sortrevind], pvals_corrected[sortrevind]


def fisher_table(motif_name, regions, mpbs, gene_set=False, mpbs_set=False):
    """
    TODO

    Keyword arguments:
    motif_name -- TODO
    regions -- TODO
    mpbs -- TODO
    gene_set -- TODO
    mpbs_set -- TODO

    Return:
    a -- number of input regions intersecting with mpbs
    b -- number of input regions not intersecting with mpbs
    gene_set -- TODO
    mpbs_set -- GenomicRegionSet of mpbs regions intersecting with input regions
    """

    # Fetching motif
    mpbs_motif = GenomicRegionSet(name="mpbs_motif")
    for region in mpbs:
        if motif_name in region.name:
            mpbs_motif.add(region)

    # Performing intersections
    if len(mpbs_motif) > 0:
        # regions which are overlapping with mpbs_motif
        intersect_original = regions.intersect(mpbs_motif, mode=OverlapType.ORIGINAL, rm_duplicates=True)

        l_intersect = len(intersect_original)

        # Fetching genes
        if gene_set:
            gene_set_res = GeneSet(motif_name)
            for genomic_region in intersect_original:
                if genomic_region.name:
                    gene_list = [e if e[0] != "." else e[1:] for e in genomic_region.name.split(":")]
                    for g in gene_list:
                        gene_set_res.genes.append(g)
            gene_set_res.genes = list(set(gene_set_res.genes))  # Keep only unique genes
        else:
            gene_set_res = None

        # Fetching mpbs
        if mpbs_set:
            mpbs_set_res = mpbs_motif.intersect(regions, mode=OverlapType.ORIGINAL, rm_duplicates=True)
        else:
            mpbs_set_res = None

        return l_intersect, len(regions) - l_intersect, gene_set_res, mpbs_set_res

    else:
        gene_set_res = GeneSet(motif_name) if gene_set else None
        mpbs_set_res = GenomicRegionSet(mpbs_motif.name) if mpbs_set else None

        return 0, len(regions), gene_set_res, mpbs_set_res


def get_fisher_dict(motif_names, regions, mpbs, gene_set=False, mpbs_set=False):
    """
    TODO

    Keyword arguments:
    motif_names -- TODO
    regions -- TODO
    mpbs -- TODO
    gene_set -- TODO
    mpbs_set -- TODO

    Return:
    res1_dict -- dictionary containing number of input regions intersecting with mpbs regions. Keys: motif names
    res2_dict --  dictionary containing number of input regions not intersecting with mpbs regions. Keys: motif names
    geneset_dict -- TODO
    mpbs_dict -- GenomicRegionSet of mpbs regions intersecting with input regions. Keys: motif names
    """

    # Calculating statistics for EV
    res1_dict = dict()
    res2_dict = dict()

    geneset_dict = dict()
    mpbs_dict = dict()

    for motif in motif_names:
        table = fisher_table(motif, regions, mpbs, gene_set=gene_set, mpbs_set=mpbs_set)

        # number of input regions intersecting mpbs regions
        res1_dict[motif] = table[0]
        # number of input regions NOT intersecting mpbs regions
        res2_dict[motif] = table[1]

        if gene_set:
            geneset_dict[motif] = table[2]

        # GenomicRegionSet of mpbs regions intersecting input regions
        if mpbs_set:
            mpbs_dict[motif] = table[3]

    # Return
    return res1_dict, res2_dict, geneset_dict, mpbs_dict
