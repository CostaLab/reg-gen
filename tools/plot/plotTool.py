# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os
lib_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(lib_path)
import numpy
from scipy.stats import mstats, wilcoxon, mannwhitneyu, rankdata
import time, datetime
import argparse
import HTML
from collections import *
import statsmodels.sandbox.stats.multicomp as sm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import getpass
import fnmatch

# Local Libraries
# Distal Libraries
from rgt.GenomicRegion import *
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
from rgt.Util import *
from rgt.CoverageSet import *

# Local test
dir = os.getcwd()

"""
Plot tool for generating analytic figures.

Author: Joseph Kuo

"""

def print2(parameter, string):
    """ Show the message on the console and also save in a list for future backup. """
    print(string)
    parameter.append(string)

def bedCoverage(bed,reads):
    """ Return coverage matrix of multiple reads on one bed. 
    bed --> GenomicRegionSet
    
    """
    c=[]
    for rp in reads:
        r = os.path.abspath(rp)   # Here change the relative path into absolute path
        cov = CoverageSet(r,bed)
        cov.coverage_from_genomicset(r)
        #cov.normRPM()
        c.append(cov.coverage)
        print("    processing: "+rp)
    return numpy.transpose(c)  

def printTable(namesCol,namesLines,table,fileName):
    """ Write a txt file from the given table. """
    d = os.path.dirname(fileName)
    try:
        os.stat(d)
    except:
        os.mkdir(d)
             
    f = open(fileName,"w")
    f.write("\t"+("\t".join(namesCol))+"\n")
    for i,line in enumerate(table):
        f.write(namesLines[i]+"\t"+("\t".join([str(j) for j in line]))+"\n")
    f.close()

def output(f, directory, folder, filename, extra=None):
    """Output the file in the defined folder """
    pd = os.path.join(dir,directory,folder)
    try:
        os.stat(os.path.dirname(pd))
    except:
        os.mkdir(os.path.dirname(pd))
    try:
        os.stat(pd)
    except:
        os.mkdir(pd)    
    
    # Saving 
    if extra == None:
        f.savefig(os.path.join(pd,filename), bbox_inches='tight',dpi=300)
    else:
        f.savefig(os.path.join(pd,filename), bbox_extra_artists=(extra), bbox_inches='tight',dpi=300)
    
    if args.pdf:
        pp = PdfPages(os.path.join(pd,filename) + '.pdf')
        pp.savefig(f, bbox_extra_artists=(extra),bbox_inches='tight') 
    
        # Set the file's metadata via the PdfPages object:
        d = pp.infodict()
        d['Title'] = args.t
        d['CreationDate'] = datetime.datetime.today()
        pp.close()
    
    if args.show:
        plt.show()

def output_parameters(parameter, directory, folder, filename):
    pd = os.path.join(dir,directory,folder)
    try:
        os.stat(os.path.dirname(pd))
    except:
        os.mkdir(os.path.dirname(pd))
    try:
        os.stat(pd)
    except:
        os.mkdir(pd)    
    if parameter:
        logp = open(os.path.join(pd,"parameters.log"),'w')
        for s in parameter:
            logp.write("{}\n".format(s))
        logp.close()
    
def output_array(array, directory, folder, filename):
    """ Write a txt file from the given array. """
    pd = os.path.join(dir,directory,folder)
    try:
        os.stat(os.path.dirname(pd))
    except:
        os.mkdir(os.path.dirname(pd))
    try:
        os.stat(pd)
    except:
        os.mkdir(pd)   
             
    f = open(os.path.join(pd,filename),"w")
    for i,line in enumerate(array):
        f.write(("\t".join(j for j in line))+"\n")
    f.close()

def integrate_html(directory):
    try:
        os.stat(os.path.join(dir,directory))
    except:
        print("No such directory.")
        return
    f = open(os.path.join(dir,directory,"index.html"),'w')
    table = []
    # Header
    table.append(['<font size="7">' + directory + "</font>"])
    # Each row is a plot with its data
    for root, dirnames, filenames in os.walk(os.path.join(dir,directory)):

        for filename in fnmatch.filter(filenames, '*.html'):
            if filename != 'index.html':
                table.append(["<a href='"+os.path.join(root, filename)+"'>"+filename.split('.')[0]+"</a>"])
    htmlcode = HTML.table(table)
    for line in htmlcode: f.write(line)
    f.close()

#################################################################################################
##### PARAMETERS ################################################################################
#################################################################################################

parser = argparse.ArgumentParser(description='Provides various plotting tools.\nAuthor: Joseph Kuo, Ivan Gesteira Costa Filho')

helpinput = 'The file name of the input Experimental Matrix file. Recommended to add more columns for more information for ploting. For example, cell type or factors.'
helpoutput = 'The directory name for the output files. For example, project name.'
helpt = 'The title shown on the top of the plot and also the folder name...'


subparsers = parser.add_subparsers(help='sub-command help',dest='mode')

parser_boxplot = subparsers.add_parser('boxplot',help='Boxplot based on the bam and bed files for gene association analysis.')
parser_boxplot.add_argument('input',help=helpinput)
parser_boxplot.add_argument('output', help=helpoutput)
parser_boxplot.add_argument('-t', default='Boxplot', help=helpt)
parser_boxplot.add_argument('-g',choices=['read','region','col4','col5','col6'], default='read', help="The way to group data into different subplots.(Default:read)")
parser_boxplot.add_argument('-s',choices=['read','region','col4','col5','col6'], default='col4', help="The way to sort data which shared the same color.(Default:col4)")
parser_boxplot.add_argument('-c',choices=['read','region','col4','col5','col6'], default='region', help="The way to color data within the same subplot.(Default:region)")
parser_boxplot.add_argument('-nqn', action="store_true", help='No quantile normalization in calculation.')
parser_boxplot.add_argument('-pdf', action="store_true", help='Save the figure in pdf format.')
parser_boxplot.add_argument('-html', action="store_true", help='Save the figure in html format.')
parser_boxplot.add_argument('-sl', type=float, default=0.01, help='Define the significance level for multiple test. Default: 0.01')
parser_boxplot.add_argument('-show', action="store_true", help='Show the figure in the screen.')

###########################################################################
parser_lineplot = subparsers.add_parser('lineplot', help='Generate lineplot with various modes.')

choice_method = ['midpoint','leftend','rightend','bothends'] # Be consist as the arguments of GenomicRegionSet.relocate_regions
choice_line = ['regions-read','reads-region']
choice_group=['col4','col5','col6']
#choice_direction = ['left', 'right', 'both']

parser_lineplot.add_argument('input', help=helpinput)
parser_lineplot.add_argument('output', default='lineplot', help=helpoutput)
parser_lineplot.add_argument('-t', default='Lineplot', help=helpt)
parser_lineplot.add_argument('-c', choices=choice_method, default='midpoint', 
                             help='Define the center to calculate coverage on the regions. Options are: '+', '.join(choice_method) + 
                             '.(Default:midpoint) The bothend mode will flap the right end region for calculation.')
parser_lineplot.add_argument('-l', choices=choice_line, default='regions-read', 
                             help="Define the lines' content of each single plot. Options are: "+', '.join(choice_line) + '(Default:regions-read)')
parser_lineplot.add_argument('-g', choices=choice_group, default=None, 
                             help='Define the grouping way for reads and regions which will be shown as the labels of rows. For example, choosing the column where cell type information is, all the read(s) and region(s) are grouped together to form the plot. Options are: '+', '.join(choice_group) + '.(Default:col4)')
parser_lineplot.add_argument('-f', choices=choice_group, default=None,
                             help='Define the factor column in Experimental Matrix which will be shown as the labels of columns and legends.')
parser_lineplot.add_argument('-e', type=int, default=2000, help='Define the extend length of interested region for plotting.(Default:2000)')
parser_lineplot.add_argument('-r', type=int, default=200, help='Define the readsize for calculating coverage.(Default:200)')
parser_lineplot.add_argument('-s', type=int, default=50, help='Define the stepsize for calculating coverage.(Default:50)')
parser_lineplot.add_argument('-b', type=int, default=100, help='Define the binsize for calculating coverage.(Default:100)')
parser_lineplot.add_argument('-pdf', action="store_true", help='Save the figure in pdf format.')
parser_lineplot.add_argument('-html', action="store_true", help='Save the figure in html format.')
parser_lineplot.add_argument('-show', action="store_true", help='Show the figure in the screen.')
parser_lineplot.add_argument('-table', action="store_true", help='Store the tables of the figure in text format.')

###########################################################################
parser_heatmap = subparsers.add_parser('heatmap', help='Generate heatmap with various modes.')

choice_method = ['midpoint','leftend','rightend','bothends'] # Be consist as the arguments of GenomicRegionSet.relocate_regions
choice_line = ['regions-read','reads-region']
choice_group=['col4','col5','col6']
parser_heatmap.add_argument('input', help=helpinput)
parser_heatmap.add_argument('output', default='lineplot', help=helpoutput)
parser_heatmap.add_argument('-t', default='Heatmap', help=helpt)
parser_heatmap.add_argument('-sort',type=int, default=None, help='Define the way to sort the signals. \
Default is no sorting at all, the signals arrange in the order of their position; \
"0" is sorting by the average ranking of all signals; \
"1" is sorting by the ranking of 1st column; "2" is 2nd and so on... ')
parser_heatmap.add_argument('-c', choices=choice_method, default='midpoint', 
                             help='Define the center to calculate coverage on the regions. Options are: '+', '.join(choice_method) + 
                             '.(Default:midpoint) The bothend mode will flap the right end region for calculation.')
parser_heatmap.add_argument('-g', choices=choice_group, default=None, 
                             help='Define the grouping way for reads and regions. For example, choosing the column where \
                             cell type information is, all the read(s) and region(s) are grouped together to form the plot. \
                             Options are: '+', '.join(choice_group) + '.(Default:col4)')
parser_heatmap.add_argument('-f', choices=choice_group, default=None,
                             help='Define the factor column in Experimental Matrix.')
parser_heatmap.add_argument('-e', type=int, default=2000, help='Define the extend length of interested region for plotting.(Default:2000)')
parser_heatmap.add_argument('-r', type=int, default=200, help='Define the readsize for calculating coverage.(Default:200)')
parser_heatmap.add_argument('-s', type=int, default=50, help='Define the stepsize for calculating coverage.(Default:50)')
parser_heatmap.add_argument('-b', type=int, default=100, help='Define the binsize for calculating coverage.(Default:100)')
parser_heatmap.add_argument('-log', action="store_true", help='Set the signal in logarithm scale.')
parser_heatmap.add_argument('-pdf', action="store_true", help='Save the figure in pdf format.')
parser_heatmap.add_argument('-html', action="store_true", help='Save the figure in html format.')
parser_heatmap.add_argument('-show', action="store_true", help='Show the figure in the screen.')


parser_integrate = subparsers.add_parser('integrate', help='Offer some tools to integrate all the plots.')
parser_integrate.add_argument('-ihtml', help='Integrate all the html files within the given directory and generate index.html for all plots.')

args = parser.parse_args()


################################################################################################

if args.mode =='integrate':
    integrate_html(args.ihtml)
    sys.exit("Integration finished.")




#################################################################################################
##### INPUT #####################################################################################
#################################################################################################

t0 = time.time()
# Input parameters dictionary
parameter = [] 

parameter.append("Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
parameter.append("User: " + getpass.getuser())
parameter.append("Project: " + args.output)
parameter.append("\nCommand:\n   $ python " + " ".join(sys.argv))
parameter.append("")



# Read the Experimental Matrix
exps = ExperimentalMatrix()
exps.read(args.input)
beds = exps.get_regionsets() # A list of GenomicRegionSets
bednames = exps.get_regionsnames()
reads = exps.get_readsfiles()
readsnames = exps.get_readsnames()
fieldsDict = exps.fieldsDict



################################ boxplot  #######################################################

if args.input and args.mode=='boxplot':
    
    def quantile_normalization(matrix):
        """ Return the np.array which contains the normalized values
        """
        
        rank_matrix = []
        for c in range(matrix.shape[1]):
            col = matrix[:,c]
            rank_col = mstats.rankdata(col)
            rank_matrix.append(rank_col)
    
        ranks = numpy.array(rank_matrix)
        trans_rank = numpy.transpose(ranks)
        
        # Calculate for means of ranks
        print("    Calculating for the mean of ranked data...")
        sort_matrix = numpy.sort(matrix,axis=0)
        means = []
        for r in range(matrix.shape[0]):
            row = [x for x in sort_matrix[r,:]]
            means.append(numpy.mean(row))
    
        # Replace the value by new means
        print("    Replacing the data value by normalized mean...")
        normalized_table = numpy.around(trans_rank)
        for i, v  in enumerate(means):
            normalized_table[normalized_table == i+1] = v
        #print(rounded_rank)
        return normalized_table
    
    def tables_for_plot(norm_table,all_bed,beds):
        """ Return a Dict which stores all tables for each bed with bedname as its key. """
        tableDict = {} # Storage all tables for each bed with bedname as the key
        conList = []   # Store containers of beds
        iterList = []
        
        for i,bed in enumerate(beds):
            tableDict[bed.name] = []
            bed.sort()
            conList.append(bed.__iter__())
            iterList.append(conList[-1].next())
            
        for i, r in enumerate(all_bed.sequences):
            for j in range(len(beds)):
                while r > iterList[j]:
                    try:
                        iterList[j] = conList[j].next()
                    except:
                        break
                if r == iterList[j]:
                    tableDict[beds[j].name].append(norm_table[i])
                elif r < iterList[j]:
                    continue
        return tableDict
    
    def unique(a):
        seen = set()
        return [seen.add(x) or x for x in a if x not in seen]
    
    def group_tags(exps, groupby, colorby, sortby):
        """Generate the tags for the grouping of plot
        
        Parameters:
            groupby = 'bam','bed','cell',or 'factor'
            colorby = 'bam','bed','cell',or 'factor'
            sortby = 'bam','bed','cell',or 'factor'
        """
        
        if groupby == "read":
            l = [exps.get_type(i,"factor") for i in exps.get_readsnames()]
            group_tags = unique(l)
        elif groupby == "region":
            l = [exps.get_type(i,"factor") for i in exps.get_regionsnames()]
            group_tags = unique(l)
        else:
            group_tags = exps.fieldsDict[exps.fields[int(groupby[-1])-1]].keys()
            
        if colorby == "read":
            l = [exps.get_type(i,"factor") for i in exps.get_readsnames()]
            color_tags = unique(l)
        elif colorby == "region":
            l = [exps.get_type(i,"factor") for i in exps.get_regionsnames()]
            color_tags = unique(l)
        else:
            color_tags = exps.fieldsDict[exps.fields[int(colorby[-1])-1]].keys()
        
        if sortby == "read":
            l = [exps.get_type(i,"factor") for i in exps.get_readsnames()]
            sort_tags = unique(l)
        elif sortby == "region":
            l = [exps.get_type(i,"factor") for i in exps.get_regionsnames()]
            sort_tags = unique(l)
        else:
            sort_tags = exps.fieldsDict[str(exps.fields[int(sortby[-1])-1])].keys()
        return group_tags, color_tags, sort_tags
    
    def group_data(tables, exps, group_tags, color_tags, sort_tags):
        
        plotDict = {}  # Extracting the data from different bed_bams file
        cues = {}   # Storing the cues for back tracking
        for bedname in tables.keys():
            plotDict[bedname] = {}
            mt = numpy.array(tables[bedname])
            for i,readname in enumerate(exps.get_readsnames()):
                plotDict[bedname][readname] = mt[:,i]
                #print(plotDict[bedname][readname])
                x = tuple(exps.get_types(readname) + exps.get_types(bedname))
                cues[x] = [bedname, readname]
        #print(cues.keys())
        
        sortDict = {}  # Storing the data by sorting tags
        for g in group_tags:
            #print("    "+g)
            sortDict[g] = {}
            for a in sort_tags:
                #print("        "+c)
                sortDict[g][a] = {}
                for c in color_tags:
                    #print("            "+a)
                    sortDict[g][a][c] = []
                    for k in cues.keys():
                        if set([g,a,c]) == set(k):
                            sortDict[g][a][c] = plotDict[cues[k][0]][cues[k][1]]
                            
        return sortDict 
    
    def box_plot(sortDict, group_tags, color_tags, sort_tags):
        """ Return boxplot from the given tables.
        
        """
        ## Initialize the figure
        #plt.title("Boxplot of Gene Association analysis")
        
        #colors = plt.cm.Accent(numpy.linspace(0, 1, 12))
        colors = [ 'lightgreen', 'pink', 'cyan', 'lightblue', 'tan']
        #colors = [(0, 35/255, 138/255),(132/255, 29/255, 20/255)]
        # Subplt by group_tags
        f, axarr = plt.subplots(1, len(group_tags), dpi=300)
        canvas = FigureCanvas(f)
        canvas.set_window_title(args.t)
        #f.suptitle(args.t, fontsize=20)
        try: axarr = axarr.reshape(-1)
        except: axarr = [axarr]
        plt.subplots_adjust(bottom=0.3)
        
        axarr[0].set_ylabel("Count number (log)")
        for i, g in enumerate(group_tags):
            axarr[i].set_title(g, y=0.94)
            axarr[i].set_yscale('log')
            axarr[i].tick_params(axis='y', direction='out')
            axarr[i].yaxis.tick_left()
            axarr[i].yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.7, zorder=1)
            d = []  # Store data within group
            color_t = []  # Store tag for coloring boxes
            x_ticklabels = []  # Store ticklabels
            for k, a in enumerate(sort_tags):
                for j, c in enumerate(color_tags):
                    if sortDict[g][a][c] == []:  # When there is no matching data, skip it
                        continue
                    else:
                        d.append([x+1 for x in sortDict[g][a][c]])
                        color_t.append(colors[j])
                        x_ticklabels.append(a)  #  + "." + c
            
            # Fine tuning boxplot
            bp = axarr[i].boxplot(d, notch=False, sym='o', vert=True, whis=1.5, positions=None, widths=None, 
                                  patch_artist=True, bootstrap=None)
            z = 10 # zorder for bosplot
            plt.setp(bp['whiskers'], color='black',linestyle='-',linewidth=0.8,zorder=z)
            plt.setp(bp['fliers'], markerfacecolor='gray',color='none',alpha=0.3,markersize=1.8,zorder=z)
            plt.setp(bp['caps'],color='none',zorder=z)
            plt.setp(bp['medians'], color='black', linewidth=1.5,zorder=z+1)
            legends = []
            for patch, color in zip(bp['boxes'], color_t):
                patch.set_facecolor(color) # When missing the data, the color patch will exceeds
                patch.set_zorder(z)
                legends.append(patch)
                
            # Fine tuning subplot
            axarr[i].set_xticks([len(color_tags)*n + 1 + (len(color_tags)-1)/2 for n,s in enumerate(sort_tags)])
            #plt.xticks(xlocations, sort_tags, rotation=90, fontsize=10)
            axarr[i].set_xticklabels(sort_tags, rotation=0, fontsize=10)
            
            axarr[i].set_ylim(bottom=0.95)
            for spine in ['top', 'right', 'left', 'bottom']:
                axarr[i].spines[spine].set_visible(False)
            axarr[i].tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
            
            if i > 0:
                plt.setp(axarr[i].get_yticklabels(),visible=False)
                #plt.setp(axarr[i].get_yticks(),visible=False)
                axarr[i].minorticks_off()
                axarr[i].tick_params(axis='y', which='both', left='off', right='off', labelbottom='off')
                    
        plt.setp([a.get_yticklabels() for a in axarr[1:]], visible=False)
        #plt.legend(colors, color_tags, loc=7)
        axarr[-1].legend(legends[0:len(color_tags)], color_tags, loc='center left', handlelength=1, 
                 handletextpad=1, columnspacing=2, borderaxespad=0., prop={'size':10},
                 bbox_to_anchor=(1.05, 0.5))
        f.tight_layout(pad=2, h_pad=None, w_pad=None)
        output(f=f, directory = args.output, folder = args.t, filename="boxplot",extra=plt.gci())
        
        
        ########## HTML ###################
        if args.html:
            pd = os.path.join(dir,args.output,args.t)
            try:
                os.stat(os.path.dirname(pd))
            except:
                os.mkdir(os.path.dirname(pd))
            try:
                os.stat(pd)
            except:
                os.mkdir(pd)    
            f = open(os.path.join(pd,'boxplot.html'),'w')
            table = []
            # Header
            table.append(['<font size="7">' + args.t + "</font>"])
            # Each row is a plot with its data
            table.append(["<img src='boxplot.png' width=800 >"])
            
            #### Calculate p value ####
                
            for g in group_tags:
                #table.append(['<font size="5">' + g + "</font>"])
                indM = 0
                header = []
                data_p = []
                arr = []
                
                for s in sort_tags:
                    for c in color_tags:
                        header.append("{0}.{1}".format(s,c))
                        data_p.append(sortDict[g][s][c])
    
                        for i, d in enumerate(data_p[:indM]):
                            u, p_value = mannwhitneyu(data_p[indM], d)
                            
                            arr.append(p_value)
                        indM = indM + 1
                #print(len(arr))
                [h,pc,a,b] = sm.multipletests(arr, alpha=0.05, returnsorted=False)
                #print(len(h))
                #print(len(pc))
                ar = numpy.chararray([len(color_tags)*len(sort_tags),len(color_tags)*len(sort_tags)], itemsize=10)
                ar[:] = "-"
                k = 0
                for c in color_tags:
                    for s in sort_tags:
                        for i, d in enumerate(header[:k]):
                            ar[k,i] = "{:3.1e}".format(pc[0.5*k*(k-1) + i])
                        k = k + 1
                            
                fig = plt.figure()
                at = plt.subplot(111, frame_on=False) 
                at.set_title(str(g))
                at.set_axis_off()
                
                nrows, ncols = ar.shape
                width, height = 1.0/(ncols+1), 1.0/(nrows+1)
                
                tb = matplotlib.table.Table(at, bbox=[0,0,1,1])
                # Add cells
                for (i,j), val in numpy.ndenumerate(ar):
                    tb.add_cell(i+1, j+1, width, height, text=val, loc='center', facecolor='none')
                    try:
                        if float(val) < args.sl:
                            tb.get_celld()[(i+1,j+1)].get_text().set_color('red')
                    except: pass
                # Add column label
                for i, h in enumerate(header):
                    # Row header
                    tb.add_cell(i+1, 0, width, height, text=h, loc='right', edgecolor='none', facecolor='none')
                    # Column header
                    tb.add_cell(0, i+1, width, height/2, text=h, loc='center', edgecolor='none', facecolor='none') #
                    
                tb.set_fontsize(15)
                at.add_table(tb)
                #plt.subplots_adjust(bottom=0.3)
                fn = os.path.join(dir,args.output,args.t,g+'_p_value.png')
                fig.savefig(fn, bbox_inches='tight',dpi=300)
                table.append(["<img src='"+fn+"' width=800 >"])
            table.append(["<a href='"+os.path.join(dir, args.output,args.t,"parameters.log")+"'>Parameters</a>"])
            htmlcode = HTML.table(table)
            for line in htmlcode: f.write(line)
            f.close()
    ###########################################################################
    #                          BOXPLOT
    ###########################################################################
    # Combine the beds

    
    print2(parameter,"\nStep 1/5: Combining all regions")
    all_bed = GenomicRegionSet("All regions")
    for bed in beds:
        all_bed.combine(bed)
    all_bed.remove_duplicates() #all_bed is sorted!!
    print2(parameter,"    " + str(len(all_bed.sequences)) + " regions from all bed files are combined.")
    t1 = time.time()
    
    print2(parameter,"    --- finished in {0:.3f} secs\n".format(t1-t0))
    
    # Coverage of reads on all_bed
    print2(parameter,"Step 2/5: Calculating coverage of each bam file on all regions")
    all_table = bedCoverage(all_bed,reads) 
    regionnames = [r.chrom+":"+str(r.initial)+"-"+str(r.final) for r in all_bed]
    #printTable(readsnames,regionnames,all_table,os.path.join(dir,"cov/all_bed.txt"))
    t2 = time.time()
    print2(parameter,"    --- finished in {0:.3f} secs\n".format(t2-t1))
    
    # Quantile normalization
    print2(parameter,"Step 3/5: Quantile normalization of all coverage table")
    
    if args.nqn:
        print2(parameter,"    No quantile normalization.")
        norm_table = all_table
    else:
        norm_table = quantile_normalization(all_table)
    #printTable(readsnames,regionnames,norm_table,os.path.join(dir,"cov/norm_bed.txt"))
    t3 = time.time()
    print2(parameter,"    --- finished in {0:.3f} secs\n".format(t3-t2))
    
    # Generate individual table for each bed
    print2(parameter,"Step 4/5: Constructing different tables for box plot")
    tables = tables_for_plot(norm_table,all_bed, beds)
    #print(tables.keys())
    #print(tables)
    
    #for bed in beds:
        #print(len(bed.sequences))
        #regionnames = [r.chrom+":"+str(r.initial)+"-"+str(r.final) for r in bed]
        #printTable(readsnames,regionnames,tables[bed.name],dir+"/cov/"+bed.name+".txt")
    t4 = time.time()
    print2(parameter,"    --- finished in {0:.3f} secs\n".format(t4-t3))
    
    # Plotting
    print2(parameter,"Step 5/5: Plotting")
    group_tags, color_tags, sort_tags = group_tags(exps, groupby=args.g, colorby=args.c, sortby=args.s)
    #print(group_tags, "\n", color_tags, "\n", sort_tags)
    sortDict = group_data(tables, exps, group_tags, color_tags, sort_tags)
    #print(sortDict)
    box_plot(sortDict,group_tags, color_tags, sort_tags)
    t5 = time.time()
    print2(parameter,"    --- finished in {0:.3f} secs\n".format(t5-t4))
    print2(parameter,"Total running time is : " + str(datetime.timedelta(seconds=round(t5-t0))))
    output_parameters(parameter, directory = args.output, folder = args.t, filename="parameters.log")
    
##################################################################################################
################################ lineplot  #######################################################
##################################################################################################    
elif args.input and args.mode=='lineplot':
    print2(parameter, "Parameters:\tExtend length:\t"+str(args.e))
    print2(parameter, "\t\tRead size:\t"+str(args.r))
    print2(parameter, "\t\tBin size:\t"+str(args.b))
    print2(parameter, "\t\tStep size:\t"+str(args.s))
    print2(parameter, "\t\tCenter mode:\t"+str(args.c+"\n"))
    
    
    t0 = time.time()
    # Processing the regions by given parameters
    print2(parameter, "Step 1/3: Processing regions by given parameters")
    processed_beds = []
    processed_bedsF = [] # Processed beds to be flapped
    for bed in beds:
        if args.c == 'bothends':
            newbed = bed.relocate_regions(center='leftend', left_length=args.e, right_length=args.e+int(0.5*args.b))
            processed_beds.append(newbed)
            newbedF = bed.relocate_regions(center='rightend', left_length=args.e+int(0.5*args.b), right_length=args.e)
            processed_bedsF.append(newbedF)
        else:
            newbed = bed.relocate_regions(center=args.c, left_length=args.e, right_length=args.e+int(0.5*args.b))
            processed_beds.append(newbed)
            #for b in newbed.sequences:
            #    print(b.__repr__())
    
    # Grouping files
    if args.g:
        groupedbed = OrderedDict()  # Store all bed names according to their types
        for bed in bednames:
            ty = exps.get_type(bed,exps.fields[int(args.g[3])-1])
            try: groupedbed[ty].append(bed)
            except: groupedbed[ty] =[bed]
        groupedbam = OrderedDict()  # Store all bam names according to their types
        for bam in readsnames:
            ty = exps.get_type(bam,exps.fields[int(args.g[3])-1])
            try: groupedbam[ty].append(bam)
            except: groupedbam[ty] =[bam]
    else:
        groupedbed = OrderedDict()
        groupedbed["All"] = bednames
        groupedbam = OrderedDict()
        groupedbam["All"] = readsnames
    
    t1 = time.time()
    print2(parameter, "    --- finished in {0:.2f} secs".format(t1-t0))
    
    
    # Calculate for coverage
    print2(parameter, "\nStep 2/3: Calculating the coverage to all reads and averaging ")
    len_coverage = int(2* args.e/args.s)
    data = OrderedDict()
    for t in groupedbed.keys():
        print2(parameter, "  For "+t+":")
        data[t] = []
        totn = len(groupedbed[t]) * len(groupedbam[t])
        bi = 0
        for ai, bed in enumerate(groupedbed[t]):
            i = bednames.index(bed)
            data[t].append([]) # New empty list for bed
            for aj, bam in enumerate(groupedbam[t]):
                ts = time.time()
                #print("    calculating the coverage of " + bed + " and " + bam)
                j = readsnames.index(bam)
                c = CoverageSet(bed+"."+bam,processed_beds[i])            
                c.coverage_from_bam(reads[j], read_size = args.r, binsize = args.b, stepsize = args.s)
                # When bothends, consider the fliping end
                if args.c == 'bothends':
                    flap = CoverageSet("for flap", processed_bedsF[i])
                    flap.coverage_from_bam(reads[j], read_size = args.r, binsize = args.b, stepsize = args.s)
                    ffcoverage = numpy.fliplr(flap.coverage)
                    c.coverage = numpy.concatenate((c.coverage, ffcoverage), axis=0)
                # Averaging the coverage of all regions of each bed file
                #avearr = numpy.zeros(len_coverage+1)
                #avearr = [0] * int((2*args.e / args.s + 1))
                #for a in c.coverage:
                #    avearr = avearr + a
                #avearr = avearr/len(c.coverage)
                avearr = numpy.array(c.coverage)
                avearr = numpy.average(avearr, axis=0)
                numpy.transpose(avearr)
                data[t][ai].append(avearr) # Store the array into data list
                bi += 1
                te = time.time()
                print2(parameter, "     Computing ("+ str(bi)+"/"+str(totn)+")\t" + "{0:40}   --{1:<6.1f}secs".format(bed+"."+bam, ts-te))
    t2 = time.time()
    print2(parameter, "    --- finished in {0:.2f} secs".format(t2-t1))
    
    # Plotting
    print2(parameter, "\nStep 3/3: Plotting the lineplots")
    rot = 50
    ticklabelsize = 7
    #color_list = [ 'lightgreen', 'pink', 'cyan', 'lightblue', 'tan']
    color_list = plt.cm.Set2(numpy.linspace(0, 1, 8))
    
    if args.l == 'regions-read':
        f, axs = plt.subplots(len(data.keys()),len(groupedbam.values()[0]), dpi=300) # figsize=(8.27, 11.69)
        if len(data.keys()) == 1 and len(groupedbam.values()) == 1: axs=[axs]
    
        for it, ty in enumerate(data.keys()):
            for i,bam in enumerate(groupedbam[ty]):
                if it == 0: axs[it,i].set_title(bam.split("_")[0],fontsize=11)
                
                #x = range(-args.e, args.e + int(0.5*args.b), args.s)
                # Processing for future output
                if args.table:
                    pArr = numpy.array(["Name","X","Y"]) # Header
                    
                for j, bed in enumerate(groupedbed[ty]): 
                    y = data[ty][j][i]
                    x = numpy.linspace(-args.e, args.e, len(y))
                    #x = range(-args.e + int(0.5*args.b), args.e + int(0.5*args.b)+args.s, args.s)
                    #print("x {}  y {}".format(len(x), len(y)))
                    axs[it, i].plot(x,y, color=color_list[j], lw=1)
                    # Processing for future output
                    if args.table:
                        name = numpy.array([bed]*len(x))
                        xvalue = numpy.array(x)
                        yvalue = numpy.array(y)
                        conArr = numpy.vstack([name,xvalue,yvalue])
                        conArr = numpy.transpose(conArr)
                        pArr = numpy.vstack([pArr, conArr])
                if args.table:
                    output_array(pArr, directory = args.output, folder ="lineplot_tables",filename=ty+"_"+bam)
                
                axs[it,i].set_xlim([-args.e, args.e])
                plt.setp(axs[it, i].get_xticklabels(), fontsize=ticklabelsize, rotation=rot)
                plt.setp(axs[it, i].get_yticklabels(), fontsize=ticklabelsize)
                axs[it, i].locator_params(axis = 'x', nbins = 4)
                axs[it, i].locator_params(axis = 'y', nbins = 4)
                
                for n in groupedbed:
                    legendl = exps.fieldsDict[exps.fields[int(args.f[-1])-1]].keys()
                axs[0,-1].legend(legendl, loc='center left', handlelength=1, handletextpad=1, columnspacing=2, borderaxespad=0., prop={'size':10}, bbox_to_anchor=(1.05, 0.5))
                
        for i,ty in enumerate(data.keys()):
            axs[i,0].set_ylabel("{}".format(ty),fontsize=12, rotation=90)
        f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
        output(f=f, directory = args.output, folder = args.t, filename="lineplot",extra=plt.gci())
    
    if args.html:
        ########## HTML ###################
        pd = os.path.join(dir,args.output,args.t)
        try:
            os.stat(os.path.dirname(pd))
        except:
            os.mkdir(os.path.dirname(pd))
        try:
            os.stat(pd)
        except:
            os.mkdir(pd)    
        f = open(os.path.join(pd,'lineplot.html'),'w')
        table = []
        # Header
        table.append(['<font size="7">' + args.t + "</font>"])
        # Each row is a plot with its data
        table.append(["<img src='lineplot.png' width=800 >"])
        
        table.append(["<a href='"+os.path.join(dir, args.output,args.t,"parameters.log")+"'>Parameters</a>"])
        htmlcode = HTML.table(table)
        for line in htmlcode: f.write(line)
        f.close()
             
    t3 = time.time()
    print2(parameter, "    --- finished in {0:.2f} secs".format(t3-t2))
    print2(parameter, "Total running time is : " + str(datetime.timedelta(seconds=t3-t0)))
    output_parameters(parameter, directory = args.output, folder = args.t, filename="parameters.log")

#################################################################################################    
################################ heatmap  #######################################################
#################################################################################################
elif args.input and args.mode=='heatmap':
    print2(parameter, "Parameters:\tExtend length:\t"+str(args.e))
    print2(parameter, "\t\tRead size:\t"+str(args.r))
    print2(parameter, "\t\tBin size:\t"+str(args.b))
    print2(parameter, "\t\tStep size:\t"+str(args.s))
    print2(parameter, "\t\tCenter mode:\t"+str(args.c+"\n"))
    
    t0 = time.time()
    # Processing the regions by given parameters
    print2(parameter, "Step 1/4: Processing regions by given parameters")
    processed_beds = []
    processed_bedsF = [] # Processed beds to be flapped
    for bed in beds:
        if args.c == 'bothends':
            newbed = bed.relocate_regions(center='leftend', left_length=args.e, right_length=args.e+int(0.5*args.b))
            processed_beds.append(newbed)
            newbedF = bed.relocate_regions(center='rightend', left_length=args.e+int(0.5*args.b), right_length=args.e)
            processed_bedsF.append(newbedF)
        else:
            newbed = bed.relocate_regions(center=args.c, left_length=args.e, right_length=args.e+int(0.5*args.b))
            processed_beds.append(newbed)
            
    # Grouping files
    if args.g:
        groupedbed = OrderedDict()  # Store all bed names according to their types
        for bed in bednames:
            ty = exps.get_type(bed,exps.fields[int(args.g[3])-1])
            try: groupedbed[ty].append(bed)
            except: groupedbed[ty] =[bed]
        groupedbam = OrderedDict()  # Store all bam names according to their types
        for bam in readsnames:
            ty = exps.get_type(bam,exps.fields[int(args.g[3])-1])
            try: groupedbam[ty].append(bam)
            except: groupedbam[ty] =[bam]
    else:
        groupedbed = OrderedDict()
        groupedbed["All"] = bednames
        groupedbam = OrderedDict()
        groupedbam["All"] = readsnames
    
    t1 = time.time()
    print2(parameter, "    --- finished in {0:.2f} secs".format(t1-t0))
    
    
    # Calculate for coverage
    print2(parameter, "\nStep 2/4: Calculating the coverage to all reads and averaging ")
    len_coverage = int(2* args.e/args.s)
    data = OrderedDict()
    for t in groupedbed.keys():
        print2(parameter, "  For "+t+":")
        data[t] = OrderedDict()
        totn = len(groupedbed[t]) * len(groupedbam[t])
        bi = 0
        for ai, bed in enumerate(groupedbed[t]):
            i = bednames.index(bed)
            data[t][bed] = OrderedDict() # New empty list for bed
            for aj, bam in enumerate(groupedbam[t]):
                ts = time.time()
                #print("    calculating the coverage of " + bed + " and " + bam)
                j = readsnames.index(bam)
                c = CoverageSet(bed+"."+bam,processed_beds[i])            
                c.coverage_from_bam(reads[j], read_size = args.r, binsize = args.b, stepsize = args.s)
                # When bothends, consider the fliping end
                if args.c == 'bothends':
                    flap = CoverageSet("for flap", processed_bedsF[i])
                    flap.coverage_from_bam(reads[j], read_size = args.r, binsize = args.b, stepsize = args.s)
                    ffcoverage = numpy.fliplr(flap.coverage)
                    c.coverage = numpy.concatenate((c.coverage, ffcoverage), axis=0)

                if args.log:
                    data[t][bed][bam] = numpy.log10(numpy.vstack(c.coverage)) # Store the array into data list
                else:
                    data[t][bed][bam] = numpy.vstack(c.coverage) # Store the array into data list
                #print(data[t][bed][bam])
                bi += 1
                te = time.time()
                print2(parameter, "     Computing ("+ str(bi)+"/"+str(totn)+")\t" + "{0:40}   --{1:<6.1f}secs".format(bed+"."+bam, ts-te))
    t2 = time.time()
    print2(parameter, "    --- finished in {0:.2f} secs".format(t2-t1))
    
    # Sorting 
    print2(parameter, "\nStep 3/4: Sorting the data for heatmap")
    if args.sort == None:
        pass
    elif args.sort == 0:
        for t in data.keys():
            for i, bed in enumerate(data[t].keys()):
                #print(numpy.sum(data[t][bed].values()[0], axis=1))
                #print(len(numpy.sum(data[t][bed].values()[0], axis=1)))
                
                sumarr = numpy.sum([numpy.sum(d, axis=1) for d in data[t][bed].values()], axis=0)
                #print(sumarr)
                #sumarr = numpy.sum(sumarr, axis=1)
                ind = stats.rankdata(sumarr,method='ordinal') # The index for further sorting
                #numpy.fliplr(ind)
                
                for j, bam in enumerate(data[t][bed].keys()):
                    d = numpy.empty(shape=(data[t][bed][bam].shape))
                    for k, ranki in enumerate(ind):
                        d[-ranki,:] = data[t][bed][bam][k,:]
                    data[t][bed][bam] = d
    else:
        for t in data.keys():
            for i, bed in enumerate(data[t].keys()):
                sumarr = numpy.sum(data[t][bed].values()[args.sort - 1], axis=1)
                #print(sumarr)
                #sumarr = numpy.sum(sumarr, axis=1)
                ind = stats.rankdata(sumarr,method='ordinal') # The index for further sorting
                #list(ind)
                #print(ind)
                for j, bam in enumerate(data[t][bed].keys()):
                    d = numpy.empty(shape=(data[t][bed][bam].shape))
                    for k, ranki in enumerate(ind):
                        d[-ranki,:] = data[t][bed][bam][k,:]
                    data[t][bed][bam] = d
                #print(data[t][bed].values()[0])
    t3 = time.time()
    print2(parameter, "    --- finished in {0:.2f} secs".format(t3-t2))
    # Plotting
    print2(parameter, "\nStep 4/4: Plotting the heatmap")
    #fig = [] # Store all figures  one figure for one type
    tickfontsize = 6
    cm_list = ['Blues', 'Reds', 'Greens', 'Oranges', 'Purples']
    pngfilenames = []
    
    for ti, t in enumerate(data.keys()):
        #fig.append(plt.figure())
        #rows = len(data[t].keys())
        columns = len(data[t].values()[0].keys())
        #fig, axs = plt.subplots(rows,columns, sharey=True, dpi=300)
        #matplotlib.pyplot.subplots_adjust(left=1, right=2, top=2, bottom=1)
        fig = plt.figure()
        plt.suptitle("Heatmap: "+t, y=1.05)
        rows = len(data[t].keys())
        
        ratio = 6
        #gs = gridspec.GridSpec(rows*ratio,columns)
        axs = numpy.empty(shape=(rows+1,columns), dtype=object)

        for bi, bed in enumerate(data[t].keys()):
            for bj, bam in enumerate(data[t][bed].keys()):
                max_value = numpy.amax(data[t][bed][bam])
                
                axs[bi, bj] = plt.subplot2grid(shape=(rows*ratio+1, columns), loc=(bi*ratio, bj), rowspan=ratio)

                if bi == 0: axs[bi, bj].set_title(exps.get_type(bam,'factor'), fontsize=7)
                im = axs[bi, bj].imshow(data[t][bed][bam], extent=[-args.e, args.e, 0,1], aspect='auto',
                                        vmin=0, vmax=max_value, interpolation='nearest', cmap=cm_list[bj])
                axs[bi, bj].set_xlim([-args.e, args.e])
                axs[bi, bj].set_xticks([-args.e, 0, args.e])
                #axs[bi, bj].set_xticklabels([-args.e, 0, args.e]
                plt.setp(axs[bi, bj].get_xticklabels(), fontsize=tickfontsize, rotation=0)
                #plt.setp(axs[bi, bj].get_yticklabels(), fontsize=10)
                #axs[bi, bj].locator_params(axis = 'x', nbins = 2)
                #axs[bi, bj].locator_params(axis = 'y', nbins = 4)
                for spine in ['top', 'right', 'left', 'bottom']:
                    axs[bi, bj].spines[spine].set_visible(False)
                axs[bi, bj].tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
                axs[bi, bj].tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
                #if bj > 0:
                #    plt.setp(axs[bi, bj].get_yticklabels(),visible=False)
                #plt.setp(axarr[i].get_yticks(),visible=False)
                axs[bi, bj].minorticks_off()
                if bj == 0:
                    nregion = len(exps.objectsDict[bed])
                    axs[bi, bj].set_ylabel(exps.get_type(bed,'factor')+" ("+str(nregion) + ")", fontsize=7)
                if bi == rows-1:
                    #divider = make_axes_locatable(axs[bi,bj])
                    #cax = divider.append_axes("bottom", size="5%", pad=0.5)
                    axs[rows,bj] = plt.subplot2grid((rows*ratio+1, columns), (rows*ratio, bj))
                    axs[rows,bj].tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
                    
                    #cbar = grid.cbar_axes[i//2].colorbar(im)
                    cbar = plt.colorbar(im, cax = axs[bi+1,bj], ticks=[0, max_value], orientation='horizontal')
                    cbar.outline.set_linewidth(0.5)
                    cbar.ax.xaxis.set_ticks_position('none')
                    if args.log:
                        cbar.ax.set_xticklabels(['0', '{:1.1f}'.format(max_value)], fontsize=tickfontsize)# horizontal colorbar
                        cbar.set_label('log10', fontsize=tickfontsize)
                    else:
                        cbar.ax.set_xticklabels(['0', int(max_value)], fontsize=tickfontsize)# horizontal colorbar
                    #cbar.set_label('Amplitute of signal')
        fig.tight_layout(pad=1, h_pad=1, w_pad=1)
        output(f=fig, directory = args.output, folder = args.t, filename="heatmap"+ "_" + t,extra=plt.gci())
        pngfilenames.append("heatmap"+ "_" + t)
        
    if args.html:
    ########## HTML ###################
        pd = os.path.join(dir,args.output,args.t)
        try:
            os.stat(os.path.dirname(pd))
        except:
            os.mkdir(os.path.dirname(pd))
        try:
            os.stat(pd)
        except:
            os.mkdir(pd)
        f = open(os.path.join(pd,'heatmap.html'),'w')
        table = []
        # Header
        table.append(['<font size="7">' + args.t + "</font>"])
        # Each row is a plot with its data
        for s in pngfilenames:
            table.append(["<img src='"+ s +".png' width=800 >"])
        
        table.append(["<a href='"+os.path.join(dir, args.output,args.t,"parameters.log")+"'>Parameters</a>"])
        htmlcode = HTML.table(table)
        for line in htmlcode: f.write(line)
        f.close()
        
    t4 = time.time()
    print2(parameter, "    --- finished in {0:.2f} secs".format(t4-t3))
    print2(parameter, "\nTotal running time is : " + str(datetime.timedelta(seconds=t4-t0)))
    
    output_parameters(parameter, directory = args.output, folder = args.t, filename="parameters.log")
    
    
    