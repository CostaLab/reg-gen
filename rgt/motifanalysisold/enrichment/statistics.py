#################################################################################################
# Functions to evaluate statistics on data
#################################################################################################

#################################################################################################
##### LIBRARIES #################################################################################
#################################################################################################

# Local Libraries
import html

# Distal Libraries
from .. util import *

# External Libraries
from fisher import pvalue
import statsmodels.sandbox.stats.multicomp as sm

#################################################################################################
##### FUNCTIONS #################################################################################
#################################################################################################

def fisherMultiple(thresholdPvalue,combinationList,multipleAlpha,realDict,randDict,geneDict,enrichedOnly=True):
    """Evaluate statistics and prints on file.

    Keyword arguments:
    thresholdPvalue -- Threshold p-value to determine enriched bindings.
    combinationList -- List of combination numbers used.
    multipleAlpha -- Alpha number for multiple testing correction.
    realDict -- Evidence-based statistics dictionary.
    randDict -- Random statistics dictionary.
    geneDict -- Dictionary of gene lists.
    enrichedOnly -- Use only enriched motifs on cobinding. (default True)

    Returns:
    resultTableList -- List of resultTables.
    """

    # Separate realDict, randDict and geneDict based on combination number
    realDictSep = dict([(c,dict()) for c in combinationList])
    randDictSep = dict([(c,dict()) for c in combinationList])
    geneDictSep = dict([(c,dict()) for c in combinationList])
    for comb in combinationList:
        for k in realDict.keys():
            if(len(k.split(",")) == comb):
                realDictSep[comb][k] = realDict[k]
                randDictSep[comb][k] = randDict[k]
                geneDictSep[comb][k] = geneDict[k]

    # Iterating on combination numbers
    resultTableList = []
    enrichedList = []
    for comb in combinationList:

        # Calculating statistics
        resultsTable = []
        for k in realDictSep[comb].keys():
            a = float(realDictSep[comb][k][0])
            b = float(realDictSep[comb][k][1])
            c = float(randDictSep[comb][k][0])
            d = float(randDictSep[comb][k][1])
            try:
                pValue = pvalue(a,b,c,d)
            except:
                pValue.right_tail = 1
            if(a+b > 0): per = a/(a+b)
            else: per = 0.0
            if(c+d > 0): bper = c/(c+d)
            else: bper = 0.0
            newGeneList = []
            for g in geneDict[k]:
                if(g[0] != "."): newGeneList.append(g)
            resultsTable.append(k.split(",")+[pValue.right_tail,a,b,c,d,round(per*100.0,2),round(bper*100.0,2),newGeneList])

        # Performing multiple test and adding corrected p-values to result table
        try:
            [h,pc,a,b] = sm.multipletests([e[comb] for e in resultsTable], alpha=multipleAlpha, returnsorted=False)
        except ValueError:
            pc = [e[comb] for e in resultsTable]
        for i in range(0,len(resultsTable)): resultsTable[i].insert(comb+1,pc[i])    

        # Sorting results tables
        resultsTable = sort.sortTableRowsByCol(resultsTable, col=comb+1, order="asc")

        # ITER 1 - Evaluating enriched list | ITER > 1 - Removing cobindings that contains motifs that are not in enriched list
        if(enrichedOnly):
            if(comb == 1): enrichedList = [e[0] for e in resultsTable if e[2] <= thresholdPvalue]
            else:
                newResultsTable = []
                for vec in resultsTable:
                    containsNonEnriched = False
                    for e in vec[:comb]:
                        if(e not in enrichedList):
                            containsNonEnriched = True
                            break
                    if(not containsNonEnriched): newResultsTable.append(vec)
                resultsTable = newResultsTable

        # Append resultTable to resultTableLists
        resultTableList.append(resultsTable)

    return resultTableList


