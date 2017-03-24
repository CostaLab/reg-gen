
###################################################################################################
# Libraries
###################################################################################################

# Python
from __future__ import print_function

# Internal
from rgt.GeneSet import GeneSet
from rgt.Util import OverlapType
from rgt.GenomicRegionSet import GenomicRegionSet

# External
from numpy import asarray, argsort, sum, arange, nonzero, minimum, maximum, int64, any, nan, inf, abs

###################################################################################################
# Functions
###################################################################################################


def ecdf(x):
    """
    Auxiliary function for multiple_test_correction
    """
    nobs = len(x)
    return arange(1,nobs+1)/float(nobs)


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
        raise ValueError('only indep and necorr implemented')
    reject = pvals_sorted < ecdffactor*alpha
    if reject.any():
        rejectmax = max(nonzero(reject)[0])
    else:
        rejectmax = 0
    reject[:rejectmax] = True

    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    pvals_corrected[pvals_corrected>1] = 1
    return reject[sortrevind], pvals_corrected[sortrevind]


def fisher_table((motif_name, regions, mpbs, return_geneset, output_mpbs_file)):
    """
    TODO

    Keyword arguments:
    m -- TODO
    region_file_name -- TODO
    mpbs_file_name -- TODO
    return_geneset -- TODO
    output_mpbs_file -- TODO

    Return:
    a -- TODO
    b -- TODO
    gene_set -- TODO
    mpbs_list -- TODO
    """

    # Initialization
    return_vec = []

    # Fetching motif
    mpbs_motif = GenomicRegionSet(name="grep_motif")
    for region in mpbs.sequences:
        if motif_name in region.name:
            mpbs_motif.add(region)

    # Performing intersections
    if len(mpbs_motif) > 0:
        # regions which are overlapping with mpbs_motif
        intersect_original = regions.intersect(mpbs_motif, mode=OverlapType.ORIGINAL, rm_duplicates=True)
        # regions which are not overlapping with regions from mpbs_motif
        subtract_overlap = regions.subtract(mpbs_motif, whole_region=True)

        # Counting the number of regions
        return_vec.append(len(intersect_original))
        return_vec.append(len(subtract_overlap))

        # Fetching genes
        if return_geneset:
            gene_set = GeneSet(motif_name)
            for genomic_region in intersect_original.sequences:
                if genomic_region.name:
                    gene_list = [e if e[0] != "." else e[1:] for e in genomic_region.name.split(":")]
                    for g in gene_list:
                        gene_set.genes.append(g)
            gene_set.genes = list(set(gene_set.genes))  # Keep only unique genes
            return_vec.append(gene_set)

        # Fetching mpbs
        if output_mpbs_file:
            mpbs_list = []
            fetch_mpbs = mpbs_motif.intersect(regions, mode=OverlapType.ORIGINAL, rm_duplicates=True)
            for region in fetch_mpbs.sequences:
                s = []
                s.extend([region.chrom, str(region.initial), str(region.final), region.name, str(region.data), region.orientation])
                mpbs_list.append(s)
            return_vec.append(mpbs_list)

    else:
        return_vec.append(0)
        return_vec.append(len(regions))
        gene_set = GeneSet(motif_name)
        return_vec.append(gene_set)
        mpbs_list = []
        return_vec.append(mpbs_list)

    # Return
    return return_vec


def get_fisher_dict(motif_names, region_file_name, mpbs_file_name, return_geneset=False,
                    output_mpbs_file=None, color="0,130,0"):
    """
    TODO

    Keyword arguments:
    motif_names -- TODO
    region_file_name -- TODO
    mpbs_file_name -- TODO
    temp_file_path -- TODO
    return_geneset -- TODO
    output_mpbs_file -- TODO

    Return:
    res1_dict -- TODO
    res2_dict -- TODO
    geneset_dict -- TODO
    """

    # Initialization
    regions = GenomicRegionSet(name="regions")
    regions.read_bed(region_file_name)
    mpbs = GenomicRegionSet(name="mpbs")
    mpbs.read_bed(mpbs_file_name)

    # Sort region and mpbs bed files
    regions.sort()
    mpbs.sort()

    # Calculating statistics for EV
    res1_dict = dict()
    res2_dict = dict()

    if return_geneset:
        geneset_dict = dict()

    for mpbs_name_group in motif_names:

        # Creating data input
        curr_data_input = [[m, regions, mpbs, return_geneset, output_mpbs_file] for m in mpbs_name_group]

        # Evaluating fisher table with multiprocessing
        curr_res = [fisher_table(x) for x in curr_data_input]
        for i in range(0, len(mpbs_name_group)):
            res1_dict[mpbs_name_group[i]] = curr_res[i][0]
            res2_dict[mpbs_name_group[i]] = curr_res[i][1]
            if return_geneset:
                geneset_dict[mpbs_name_group[i]] = curr_res[i][2]
            if output_mpbs_file:
                for vec in curr_res[i][3]:
                    output_mpbs_file.write("\t".join(vec + [vec[1], vec[2], color]) + "\n")

    # Return
    if return_geneset:
        return res1_dict, res2_dict, geneset_dict

    return res1_dict, res2_dict

