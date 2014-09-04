
###################################################################################################
# Libraries
###################################################################################################

# Python
from multiprocessing import Pool

# Internal
from .. GeneSet import GeneSet
from .. Util import ErrorHandler, OverlapType

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

def fisher_table((grs1, grs2)):
    """ 
    TODO

    Keyword arguments:
    grs1 -- TODO
    grs2 -- TODO
        
    Return:
    len(intersect_grs) -- TODO
    len(subtract_grs) -- TODO
    gene_set -- TODO
    """

    # Performing intersections
    intersect_grs = grs1.intersect(grs2, mode=OverlapType.ORIGINAL)
    subtract_grs = grs1.subtract(intersect_grs)

    # Fetching genes
    gene_set = GeneSet(grs2.name)
    for gr in intersect_grs:
        if(gr.name):
            gene_list = [e if e[0]!="." else e[1:] for e in gr.name.split(":")]
            for g in gene_list: gene_set.genes.append(g)

    # Keep only unique genes
    gene_set.genes = list(set(gene_set.genes))

    # Fetching mpbs
    intersect_grs2 = grs2.intersect(grs1, mode=OverlapType.ORIGINAL)

    return [len(intersect_grs), len(subtract_grs), gene_set, intersect_grs2]

def get_fisher_dict(grouped_mpbs_dict_keys, region_set, mpbs_dict):
    """ 
    TODO

    Keyword arguments:
    grouped_mpbs_dict_keys -- TODO
    region_set -- TODO
    mpbs_set -- TODO
        
    Return:
    res1_dict -- TODO
    res2_dict -- TODO
    geneset_dict -- TODO
    """

    # Calculating statistics for EV
    res1_dict = dict()
    res2_dict = dict()
    geneset_dict = dict()
    res_mpbs_dict = dict()
    for mpbs_name_group in grouped_mpbs_dict_keys:

        # Creating data input
        curr_data_input = [[region_set, mpbs_dict[m]] for m in mpbs_name_group if m]
        curr_proc_nb = len(curr_data_input)

        # Evaluating randomic c and d with multiprocessing
        pool = Pool(curr_proc_nb)
        curr_res = pool.map(fisher_table,curr_data_input)
        pool.close()
        pool.join()
        for i in range(0,len(mpbs_name_group)):
            res1_dict[mpbs_name_group[i]] = curr_res[i][0]
            res2_dict[mpbs_name_group[i]] = curr_res[i][1]
            geneset_dict[mpbs_name_group[i]] = curr_res[i][2]
            res_mpbs_dict[mpbs_name_group[i]] = curr_res[i][3]

    return res1_dict, res2_dict, geneset_dict, res_mpbs_dict


