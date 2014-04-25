#################################################################################################
# Creates graphical outputs.
#################################################################################################

#################################################################################################
##### LIBRARIES #################################################################################
#################################################################################################

# External Libraries
import numpy as np
import scipy.stats as st
import matplotlib as mpl
mpl.use('Agg') # As the server does not have the $DISPLAY env variable, this fix had to be added.
import matplotlib.pyplot as plt
import pylab

#################################################################################################
##### FUNCTIONS #################################################################################
#################################################################################################

def mpbsScoreHistogram(mpbsList,colorList,alphaList,labelList,outputLocation,bins=100,outExt="png",useLegend=True):
    """Creates a histogram of MPBS motif matching scores.

    Keyword arguments:
    mpbsList -- List of mpbs dictionaries.
    colorList -- List of colors for each histogram.
    alphaList -- List of alphas for each histogram.
    labelList -- List of labels for each histogram.
    outputLocation -- Location of the output graphs.
    bins -- Number of bins for histogram. (default 100)
    outExt -- Output extension for the images. (default 'png')
    useLegend -- Whether to use legend or not. (default True)

    Returns:
    outputLocation<factorName>.<outExt> -- Histogram images for each factor.
    """

    # Fetching data
    dataVec = [dict() for e in range(0,len(mpbsList))]
    for i in range(0,len(dataVec)):
        for factor in mpbsList[i].keys():
            for chrName in mpbsList[i][factor].keys():
                for e in mpbsList[i][factor][chrName]:
                    if(factor in dataVec[i].keys()): dataVec[i][factor].append(e[3])
                    else: dataVec[i][factor] = [e[3]]

    # Creating graphs
    for factor in dataVec[0].keys():

        # Creating figure
        fig = plt.figure(figsize=(8,5), facecolor='w', edgecolor='k')
        ax = fig.add_subplot(111)
        
        # Iterating on the different mpbs dictionaries
        histList = []
        realLabelList = []
        for i in range(0,len(dataVec)):

            # Verifying if factor exists
            if(not (factor in dataVec[i].keys())): continue
            realLabelList.append(labelList[i])
            
            # Creating hist list
            histVec = [0] * (bins+1)
            for e in dataVec[i][factor]: histVec[int((float(e)/1000.0)*bins)] += 1
            histList.append(histVec)

            # Plotting
            ax.hist(dataVec[i][factor], bins, facecolor=colorList[i], alpha=alphaList[i], label=labelList[i])

        # Evaluating correlations
        corrString = ""; pValueString = ""
        for i in range(0,len(histList)):
            for j in range(i+1,len(histList)):
                corrVec = st.pearsonr(histList[i],histList[j])
                corrString += realLabelList[i]+"+"+realLabelList[j]+" = "+str(corrVec[0])+"  "
                pValueString += realLabelList[i]+"+"+realLabelList[j]+" = "+str(corrVec[1])+"  "

        # Saving figure
        ax.set_title(factor)
        if(corrString != ""): ax.set_xlabel(corrString[:-2]+"\n"+pValueString[:-2])
        if(useLegend): ax.legend()
        fig.savefig(outputLocation+factor+"."+outExt, format=outExt, dpi=300, bbox_inches='tight')


