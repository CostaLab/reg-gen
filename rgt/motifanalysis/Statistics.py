
###################################################################################################
# Libraries
###################################################################################################

# Python
from os import waitpid
from os.path import basename, join, dirname
from subprocess import Popen, check_output

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

def fisher_table((motif_name,region_file_name,mpbs_file_name,return_geneset,output_mpbs_file)):
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
    to_remove = []
    return_vec = []

    # Fetching motif
    grep_file_name = mpbs_file_name+motif_name+"_grep.bed"; to_remove.append(grep_file_name)
    p1 = Popen("grep \"\t\""+motif_name+"\"\t\" "+mpbs_file_name+" > "+grep_file_name, shell=True)
    waitpid(p1.pid, 0)

    # Performing intersections
    a_file_name = mpbs_file_name+motif_name+"_A.bed"
    b_file_name = mpbs_file_name+motif_name+"_B.bed"
    n_lines_grep = int(check_output(['wc', '-l', grep_file_name]).split()[0])
    if(n_lines_grep > 0):
        p2 = Popen("intersectBed -a "+region_file_name+" -b "+grep_file_name+" -wa -u > "+a_file_name, shell=True)
        waitpid(p2.pid, 0)
        p3 = Popen("intersectBed -a "+region_file_name+" -b "+grep_file_name+" -wa -v > "+b_file_name, shell=True)
        waitpid(p3.pid, 0)
        to_remove.append(a_file_name); to_remove.append(b_file_name)

        # Counting the number of lines
        a = int(check_output(['wc', '-l', a_file_name]).split()[0])
        b = int(check_output(['wc', '-l', b_file_name]).split()[0])
        return_vec.append(a); return_vec.append(b)

        # Fetching genes
        if(return_geneset):
            gene_set = GeneSet(motif_name)
            a_file = open(a_file_name,"r")
            for line in a_file:
                ll = line.strip().split("\t")
                if(ll[3]):
                    gene_list = [e if e[0]!="." else e[1:] for e in ll[3].split(":")]
                    for g in gene_list: gene_set.genes.append(g)
            a_file.close()
            gene_set.genes = list(set(gene_set.genes)) # Keep only unique genes
            return_vec.append(gene_set)

        # Fetching mpbs
        if(output_mpbs_file):
            mpbs_list = []
            mpbs_temp_file_name = mpbs_file_name+motif_name+"_mpbstemp.bed"; to_remove.append(mpbs_temp_file_name)
            p4 = Popen("intersectBed -a "+grep_file_name+" -b "+region_file_name+" -wa -u > "+mpbs_temp_file_name, shell=True)
            waitpid(p4.pid, 0)
            mpbs_temp_file = open(mpbs_temp_file_name,"r")
            for line in mpbs_temp_file: mpbs_list.append(line.strip().split("\t"))
            mpbs_temp_file.close()
            return_vec.append(mpbs_list)

    else:
        b = int(check_output(['wc', '-l', region_file_name]).split()[0])
        return_vec.append(0); return_vec.append(b)
        gene_set = GeneSet(motif_name)
        return_vec.append(gene_set)
        mpbs_list = []
        return_vec.append(mpbs_list)

    # Remove all files
    for e in to_remove:
        p5 = Popen("rm "+e, shell=True)
        waitpid(p5.pid, 0)

    # Return
    return return_vec

def get_fisher_dict(motif_names, region_file_name, mpbs_file_name, temp_file_path, return_geneset=False, output_mpbs_file=None, color="0,130,0"):
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
    to_remove = []
    region_name = ".".join(basename(region_file_name).split(".")[:-1])
    mpbs_name = ".".join(basename(mpbs_file_name).split(".")[:-1])

    # Sort region and mpbs bed files
    region_file_name_sort = join(temp_file_path,region_name+"_sort.bed"); to_remove.append(region_file_name_sort)
    mpbs_file_name_sort = join(temp_file_path,mpbs_name+"_sort.bed"); to_remove.append(mpbs_file_name_sort)
    p1 = Popen("sort -k1,1 -k2,2n "+region_file_name+" > "+region_file_name_sort, shell=True)
    waitpid(p1.pid, 0)
    p2 = Popen("sort -k1,1 -k2,2n "+mpbs_file_name+" > "+mpbs_file_name_sort, shell=True)
    waitpid(p2.pid, 0)

    # Calculating statistics for EV
    res1_dict = dict()
    res2_dict = dict()
    if(return_geneset): geneset_dict = dict()
    for mpbs_name_group in motif_names:

        # Creating data input
        curr_data_input = [[m,region_file_name_sort,mpbs_file_name_sort,return_geneset,output_mpbs_file] for m in mpbs_name_group]
        curr_proc_nb = len(curr_data_input)

        # Evaluating fisher table with multiprocessing
        curr_res = [fisher_table(xx) for xx in curr_data_input]
        for i in range(0,len(mpbs_name_group)):
            res1_dict[mpbs_name_group[i]] = curr_res[i][0]
            res2_dict[mpbs_name_group[i]] = curr_res[i][1]
            if(return_geneset): geneset_dict[mpbs_name_group[i]] = curr_res[i][2]
            if(output_mpbs_file):
                for vec in curr_res[i][3]: output_mpbs_file.write("\t".join(vec+[vec[1],vec[2],color])+"\n")

    # Remove all files
    for e in to_remove:
        p5 = Popen("rm "+e, shell=True)
        waitpid(p5.pid, 0)

    # Return
    if(return_geneset): return res1_dict, res2_dict, geneset_dict
    return res1_dict, res2_dict


