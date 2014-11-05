#################################################################################################
# HTML output functions.
#################################################################################################

#################################################################################################
##### LIBRARIES #################################################################################
#################################################################################################

# External Libraries
import HTML

#################################################################################################
##### FUNCTIONS #################################################################################
#################################################################################################

def printHTML(logoIndexList,genesIndex,headerList,logoLocation,resultsMatrix,outputFileName):
    """Print results to a html file.

    Keyword arguments:
    logoIndexList -- List of indexes of the results Matrix that contains the factor names.
    genesIndex -- Index of the results Matrix that contains the gene list.
    headerList -- List of header titles.
    logoLocation -- Location of the logo graphs in png.
    resultsMatrix -- Matrix containing the results to be printed.
    outputFileName -- Path + name of the output html file.

    Returns:
    outputFileName -- File containing the html results.
    """

    # Initial table
    table = [headerList]

    # Adding data to table 
    for vec in resultsMatrix:
        resVec = []
        for i in range(0,len(vec)):
            if(i in logoIndexList):
                resVec.append(str(vec[i]))
                resVec.append("<img src='"+logoLocation+str(vec[i])+".png' width=200 >")
            elif(i == genesIndex):
                geneL = "+".join(vec[i])
                resVec.append("<a href='http://biit.cs.ut.ee/gprofiler/index.cgi?significant=1&sort_by_structure=1&ordered_query=0&organism=mmusculus&query="+geneL+"' > GO Enrichment </a>")
            else: resVec.append(str(vec[i]))
        table.append(resVec)

    # Printing table
    htmlcode = HTML.table(table)
    outputFile = open(outputFileName,"w")
    for line in htmlcode: outputFile.write(line)
    outputFile.close()
    return 0


