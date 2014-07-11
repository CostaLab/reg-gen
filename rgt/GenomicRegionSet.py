from __future__ import print_function
from __future__ import division
from GenomicRegion import *
import random
from copy import deepcopy
import os
import time
from scipy import stats
from Util import GenomeData, OverlapType, AuxiliaryFunctions
from GeneSet import GeneSet

"""
Represent list of GenomicRegions.

Authors: Ivan G. Costa, Manuel Allhoff, Joseph Kuo

Define a list of GenomicRegions. Perform operations of this set.

Methods:

readBed(filename):
Add GenomicRegions defined by row in bed file.

intersect(region):
Return intersection with region.

"""

class GenomicRegionSet:

    def __init__(self, name):
        """Initialize a set of GenomicRegions"""
        self.name = name
        self.sequences = []
        self.sorted = False
        self.fileName=""
        self.genome_path = ""
    
    def get_chrom(self):
        """Return all chromosomes"""
        return [ r.chrom for r in self.sequences ]
    
    def add(self, region):
        """Add GenomicRegion"""
        self.sequences.append(region)
        self.sorted = False
    
    def __len__(self):
        return len(self.sequences)
    
    def __iter__(self):
        return iter(self.sequences)

    def extend(self,left,right):
        """Perform extend step for every element"""
        for s in self.sequences:
            s.extend(left,right)

    def sort(self):
        """Sort Elements by criteria defined by a GenomicRegion"""
        self.sequences.sort(cmp = GenomicRegion.__cmp__)

    def read_bed(self, filename):
        """Read BED file and add every row as a GenomicRegion. 
        Chrom (1), start (2), end (2), name (4) and orientation (6) is used for GenomicRegion. 
        All other columns (5, 7, 8, ...) are put to the data variable of the GenomicRegion.
        The numbers in parentheses are the columns of the BED format.
        See BED format at: http://genome.ucsc.edu/FAQ/FAQformat.html#format1 """
        self.fileName=filename
        with open(filename) as f:
            for line in f:
                try:
                    name, orientation, data = None, None, None
                    line = line.strip("\n")
                    line = line.split("\t")
                    size = len(line)
                    chrom = line[0]
                    start, end = int(line[1]), int(line[2])
                    if start > end:
                        start, end =  end, start
                    if size > 3:
                        name = line[3]
                    if size > 5:
                        orientation = line[5]
                    if size > 5:
                        data = "\t".join( [line[4]] + line[6:] )
                    self.add( GenomicRegion(chrom, start, end, name, orientation, data) )
                except:
                    print("Error at line",line,self.fileName)
            self.sort()
  
    def read_bedgraph(self, filename):
        """Read BEDGRAPH file and add every row as a GenomicRegion. 
        See BEDGRAPH format at: http://genome.ucsc.edu/goldenPath/help/bedgraph.html """
        self.fileName=filename
        with open(filename) as f:
            for line in f:
                try:
                    line = line.strip("\n")
                    line = line.split("\t")
                    assert len(line) == 4

                    chrom, start, end, data = line[0], int(line[1]), int(line[2]), str(line[3])

                    self.add( GenomicRegion(chrom=chrom, initial=start, final=end, data=data) )
                except:
                    print("Error at line", line, self.fileName)
                    
            self.sort()
            
    def random_subregions(self,size):
        """Return a subsampling of the genomic region set with a specific number of regions"""
        z = GenomicRegionSet(self.name + '_random')
        samp = random.sample(range(len(self)),size)
        for i in samp:
            z.add(self.sequences[i])
        z.sort()
        return z                

    def write_bed(self,filename):
        """Write GenomicRegions to BED file"""
        with open(filename, 'w') as f:
            for s in self:
                print(s, file=f)

    def gene_association(self, gene_set = None, organism = "hg19", promoterLength = 1000, threshDist = 50000):
        """Associates coordinates to genes given the following rules:
           1. If the peak is inside gene (promoter+coding) then this peak is associated with that gene.
           2. If a peak is inside overlapping genes, then the peak is annotated with both genes.
           3. If peak is between two genes (not overlapping neither), then both genes are annotated.
           4. If the distance between peak and gene is greater than a threshold distance, then it is not annotated.

        Keyword arguments:
        gene_set -- List of gene names as a GeneSet object. If None, then consider all genes to be enriched. (default None)
        organism -- Organism in order to fetch genomic data. (default hg19)
        promoterLength -- Length of the promoter region. (default 1000)
        threshDist -- Threshold maximum distance for a coordinate to be considered associated with a gene. (default 50000)

        Returns:
        result_grs -- GenomicRegionSet exactely as self, but with the following additional information:
                      1. name: String of genes associated with that coordinate separated by ':'
                      2. data: String of proximity information (if the coordinate matched to the corresponding gene in the
                               previous list in a proximal position (PROX) or distal position (DIST)) separated by ':'
                      The gene will contain a '.' in the beginning of its name if it is not in the gene_set given.
        """
 
        # Initializations
        chromosome_list = list(set(self.get_chrom()))
        genome_data = GenomeData(organism)

        # Separating coordinates by chromosomes
        coord_dict = dict([(chrom,[]) for chrom in chromosome_list])
        for gr in self.sequences:
            coord_dict[gr.chrom].append([gr.initial,gr.final,gr.name,gr.data,gr.orientation])

        # Sorting coord_dict
        for chrom in coord_dict.keys(): coord_dict[chrom] = sorted(coord_dict[chrom], key=lambda x: x[0])

        # Reading assocDict
        gene_regions = GenomicRegionSet("gene_regions")
        gene_regions.read_bed(genome_data.get_association_file())
        assocDict = dict()
        geneFlagDict = dict()
        for gr in gene_regions:
            geneFlagDict[gr.name.upper()] = False
            curr_region = [gr.initial,gr.final,gr.name,0,gr.orientation]
            try: assocDict[gr.chrom].append(curr_region)
            except Exception: assocDict[gr.chrom] = [curr_region]
       
        # Sorting assocDict
        for chrom in assocDict.keys(): assocDict[chrom] = sorted(assocDict[chrom], key=lambda x: x[0])
        
        # Updating assocDict to contain all chromosomes in coord_dict
        for chrom in chromosome_list:
            if(chrom not in assocDict.keys()): assocDict[k] = []

        # Updating geneFlagDict based on gene_set
        if(gene_set):
            gene_list = [g.upper() for g in gene_set.genes]
            for e in gene_list: geneFlagDict[e] = True
        else: 
            for e in geneFlagDict.keys(): geneFlagDict[e] = True

        # Associating coord_dict with assocDict
        aDict = dict() # Results dictionary
        for chrName in chromosome_list:

            aDict[chrName] = [] # Adding new result list
            counter = 0 # Counter to run through gene list

            # Iterating on coordinate list (main list)
            for coord in coord_dict[chrName]: 

                didBreak = False # Check wether it breaked or not.

                # Running linearly through gene list
                while(counter < len(assocDict[chrName])):

                    # Extend the gene coordinates to encompass promoter region based on strandness
                    if(assocDict[chrName][counter][4] == "+"): geneCoord = [assocDict[chrName][counter][0]-promoterLength,assocDict[chrName][counter][1]]
                    else: geneCoord = [assocDict[chrName][counter][0],assocDict[chrName][counter][1]+promoterLength]

                    check = AuxiliaryFunctions.overlap(coord,geneCoord) # Check overlap between coordinate and gene+promoter

                    if(check == 0): # If contain overlap, then check if the coordinate also contains overlap with the next gene

                        # Verify if next gene+promoter also overlap
                        genesList = [assocDict[chrName][counter][2]+"_PROX"] # List of overlapping genes
                        if(counter < len(assocDict[chrName])-1): # If this is not the last gene (there is a 'next') then verify overlap
                            if(assocDict[chrName][counter+1][4] == "+"): geneCoordNext = [assocDict[chrName][counter+1][0]-promoterLength,assocDict[chrName][counter+1][1]]
                            else: geneCoordNext = [assocDict[chrName][counter+1][0],assocDict[chrName][counter+1][1]+promoterLength]
                            checkNext = AuxiliaryFunctions.overlap(coord,geneCoordNext) # Verify overlap between coordinate and next gene
                            if(checkNext == 0): genesList.append(assocDict[chrName][counter+1][2]+"_PROX") # If there is an overlap then add this gene to association list
                        #else: # If this is the last gene (there is no 'next') then dont verify overlap

                        # Verify if genes are enriched
                        for i in range(0,len(genesList)): # Check if genes in genesList are enriched 
                            if(not geneFlagDict[genesList[i][:-5]]): genesList[i] = "."+genesList[i] # If the gene is not enriched, put a dot (.) in front of it (will be used by motif matching)
                            #else: If gene is enriched than let the name without the dot (.) in front of it
                        didBreak = True
                        break

                    elif(check == -1): # If gene is after coordinate, then check if they are within maximum distance for overlap. Also do the same check for the previous gene.

                        # Verify overlap again using maximum distance with current gene+promoter
                        genesList = [] # List of overlapping genes
                        maxGeneCoord = [geneCoord[0]-threshDist,geneCoord[1]+threshDist]
                        maxCheck = AuxiliaryFunctions.overlap(coord,maxGeneCoord)
                        if(maxCheck == 0): genesList.append(assocDict[chrName][counter][2]+"_DIST") # If it overlapped then put current gene in overlap list
                    
                        # Verify overlap again using maximum distance with previous gene+promoter
                        if(counter > 0): # Do this verification only if this is not the first gene
                            if(assocDict[chrName][counter-1][4] == "+"): geneCoordPrev = [assocDict[chrName][counter-1][0]-promoterLength,assocDict[chrName][counter-1][1]]
                            else: geneCoordPrev = [assocDict[chrName][counter-1][0],assocDict[chrName][counter-1][1]+promoterLength]
                            maxGeneCoordPrev = [geneCoordPrev[0]-threshDist,geneCoordPrev[1]+threshDist]
                            maxCheckPrev = AuxiliaryFunctions.overlap(coord,maxGeneCoordPrev)
                            if(maxCheckPrev == 0): genesList.append(assocDict[chrName][counter-1][2]+"_DIST") # If it overlapped then put previous gene in overlap list

                        # Verify if genes are enriched
                        if(len(genesList) == 0): genesList.append(".") # If genesList is empty then put a '.' to represent non-association
                        else: # If genesList is not empty then verify enriched genes
                            for i in range(0,len(genesList)): # Check if genes in genesList are enriched 
                                if(not geneFlagDict[genesList[i][:-5]]): genesList[i] = "."+genesList[i] # If the gene is not enriched, put a dot (.) in front of it (will be used by motif matching)
                                #else: If gene is enriched than let the name without the dot (.) in front of it
                        didBreak = True
                        break

                    #elif(check == -1): If gene is before coordinate, dont do anything! Just move to the next gene.
                    counter += 1

                # If not breaked, then the gene list is over and the coordinate can only overlap the last gene.
                if(not didBreak): 
                    genesList = []
                    if(len(assocDict[chrName]) == 0): # If gene list is empty then dont associate coordinate with any gene
                        genesList.append(".")
                    else: # If gene list is not empty then try to verify overlap between coordinate and last gene
                        if(assocDict[chrName][counter-1][4] == "+"): geneCoordPrev = [assocDict[chrName][counter-1][0]-promoterLength,assocDict[chrName][counter-1][1]]
                        else: geneCoordPrev = [assocDict[chrName][counter-1][0],assocDict[chrName][counter-1][1]+promoterLength]
                        maxGeneCoordPrev = [geneCoordPrev[0]-threshDist,geneCoordPrev[1]+threshDist]
                        maxCheckPrev = AuxiliaryFunctions.overlap(coord,maxGeneCoordPrev)
                        if(maxCheckPrev == 0): # If it overlapped then put previous gene in overlap list and check for enrichment
                            genesList.append(assocDict[chrName][counter-1][2]+"_DIST")
                            if(not geneFlagDict[genesList[0][:-5]]): genesList[0] = "."+genesList[0]
                        else: genesList.append(".") # If does not overlap then put a '.' to represent non-association

                # Write the curent coordinate with its corresponding overlapping genes (enriched or not)
                aDict[chrName].append(coord[:5]) # Writing the raw coordinate until STRAND field only
                aDict[chrName][-1][2] = ":".join(genesList) # Write list of overlapping genes
                aDict[chrName][-1][3] = 0 # Force the SCORE field to be '0'
                aDict[chrName][-1][4] = "." # Force the STRAND field to be '.'

        # Converting aDict to genomic region set
        result_grs = GenomicRegionSet(self.name+"_associated")
        for chrom in aDict.keys():
            for coord in aDict[chrom]:
                curr_genes_list = coord[2].split(":")
                new_genes_list = []
                new_prox_list = []
                for curr_gene in curr_genes_list:
                    new_genes_list.append("_".join(curr_gene.split("_")[:-1]))
                    new_prox_list.append(curr_gene.split("_")[-1])
                new_genes_list = ":".join(new_genes_list)
                new_prox_list = ":".join(new_prox_list)
                result_grs.add(GenomicRegion(chrom, coord[0], coord[1], name=new_genes_list, data=new_prox_list))
             
            return result_grs

    def filter_by_gene_association(self, gene_set = None, organism = "hg19", promoterLength = 1000, threshDist = 50000):
        """Updates self in order to keep only the coordinates associated to genes which are in gene_set.
           It also returns information regarding the mapped genes and proximity information (if mapped in proximal (PROX)
           or distal (DIST) regions). If a region was associated with two genes, one in the gene_set and the other not in
           gene_set, then the updated self still keeps information about both genes. However, the latter won't be reported
           as mapped gene.

           Keyword arguments:
             gene_set -- List of gene names as a GeneSet object. If None, then consider all genes to be enriched. (default None)
             organism -- Organism in order to fetch genomic data. (default hg19)
             promoterLength -- Length of the promoter region. (default 1000)
             threshDist -- Threshold maximum distance for a coordinate to be considered associated with a gene. (default 50000)

           Returns:
             None -- Updates self in order to keep only the coordinates associated to genes which are in gene_set
             all_genes = GeneSet that contains all genes associated with the coordinates
             mapped_genes = GeneSet that contains the genes associated with the coordinates which are in gene_set
             all_proxs = List that contains all the proximity information of genes associated with the coordinates
             mapped_proxs = List that contains all the proximity information of genes associated with the coordinates which are in gene_set
        """

        # Making association with gene_set
        assoc_grs = self.gene_association(gene_set, organism, promoterLength, threshDist)

        # Initializing updated genomic region set
        updated_grs = GenomicRegionSet(assoc_grs.name) # New GRS will contain only coordinates associated to genes which are in gene_set

        # Initializing resulting gene sets
        all_genes = GeneSet("all_genes") # Contains all genes associated with the coordinates
        mapped_genes = GeneSet("mapped_genes") # Contains the genes associated with the coordinates which are in gene_set
        all_proxs = [] # Contains all the proximity information of genes associated with the coordinates
        mapped_proxs = [] # Contains all the proximity information of genes associated with the coordinates which are in gene_set

        # Iterating on associated genomic regions
        for gr in assoc_grs.sequences:

            # If the coordinate wasn't associated with any gene, it is not considered
            if(gr.name == "."): continue

            # Fetching list of genes and proximity information
            curr_genes = gr.name.split(":")
            curr_proxs = gr.data.split(":")

            # Updating all/mapped - genes/proxs
            flag_assoc = False
            for i in range(0,len(curr_genes)):
                g = curr_genes[i]
                p = curr_proxs[i]
                if(g[0] != "."):
                    mapped_genes.genes.append(g)
                    mapped_proxs.append(p)
                    flag_assoc = True
                else: g = g[1:]
                all_genes.genes.append(g)
                all_proxs.append(p)
            
            # Adding current genomic region to updated GRS if there was a gene associated with it in gene_set
            if(flag_assoc): updated_grs.add(deepcopy(gr))

        # Removing duplicates
        all_genes.genes = list(set(all_genes.genes))
        mapped_genes.genes = list(set(mapped_genes.genes))
        all_proxs = list(set(all_proxs))
        mapped_proxs = list(set(mapped_proxs))

        # Updating self
        self.name = updated_grs.name
        self.sequences = updated_grs.sequences
        self.sorted = False

        return all_genes, mapped_genes, all_proxs, mapped_proxs

    """
    def filter_by_gene_association(self,fileName,geneSet,geneAnnotation,genomeSize,promoterLength=1000,threshDist=50000):
        from rgt.motifanalysis.util import bedFunctions, sort
        from rgt.motifanalysis.enrichment.geneAssociation import *
        self.fileName=fileName
        regionsToGenes={}
        coordDict = bedFunctions.createBedDictFromSingleFile(fileName, features=[1,2,3,4,5]) 
        coordDict = sort.sortBedDictionary(coordDict, field=0)
        [dictBed,allBed] = geneAssociationByPromoter(coordDict,geneSet,geneAnnotation,genomeSize,promoterLength,threshDist)  
        #print dictBed
        genes=[]
        totalPeaks=0
        allgenes=[]
        for chr in dictBed.keys():
            for (v1,v2,name,orientation,data) in dictBed[chr]:
                totalPeaks+=1
                names=name.split(":")
                keep=[n for n in names if "." not in n]
                if len(keep) > 0:
                    self.add(GenomicRegion(chr,v1,v2,":".join(keep)))
                genes = genes+keep
                allgenes=allgenes+[n.strip(".") for n in names]
                regionsToGenes[chr+":"+str(v1)+"-"+str(v2)]=[n.strip(".") for n in names]
        #print "Total Peaks", total
        mappedGenes=len(list(set(allgenes)))
        self.sort()
        self.genes=list(set(genes))
        if geneSet == None:
          le=0
        else:
          le=len(geneSet)

        return le, len(self.genes), mappedGenes, totalPeaks,regionsToGenes
    """

    def intersect(self,y,mode=OverlapType.OVERLAP):
        """Return the overlapping regions with three different modes.
        
        (mode = OverlapType.OVERLAP) 
        Return new GenomicRegionSet including only the overlapping regions with y.
        
            Keyword arguments:
            y -- the GenomicRegionSet which to compare with
        
            Return:
            z -- a new GenomicRegionSet including only the overlapping regions
        
            Graphical explanation:
            self       ----------              ------
            y                 ----------                    ----
            Result            ---
            
        (mode = OverlapType.ORIGINAL)
        Return the regions of original GenomicRegionSet which have any intersections with y.
        
            Keyword arguments:
            y -- the GenomicRegionSet which to compare with
            
            Return:
            z -- the regions of original GenomicRegionSet which have any intersections with y
            
            Graphical explanation:
            self       ----------              ------
            y              ----------                    ----
            Result     ----------
            
        (mode = OverlapType.COMP_INCL)
        Return region(s) of the GenomicRegionSet which are 'completely' included by y.
        
            Keyword arguments:
            y -- the GenomicRegionSet which to compare with
            
            Return:
            z -- region(s) of self which are 'completely' included by y
            
            Graphical explanation:
            self        -------------             ------
            y              ----------      ---------------              ----
            Result                                ------
        """
        if len(self) == 0 or len(y) == 0:
            return GenomicRegionSet(self.name + '_' + y.name)
        z = GenomicRegionSet(self.name + '_' + y.name)
        
        # If there is overlap within self or y, they should be merged first. 
        self.sort()
        y.sort()
        # Iteration
        con_self = self.__iter__()
        con_y = y.__iter__()
        s = con_self.next()
        ss = con_y.next()
        cont_loop = True
        pre_z = GenomicRegion(chrom="chr1",initial=0,final=0) # Store the previous region which add in z
        while cont_loop:
            # When the regions overlap
            if s.overlap(ss):
                if mode == OverlapType.OVERLAP:
                    sss = GenomicRegion(chrom=s.chrom, 
                                        initial=max(s.initial, ss.initial), 
                                        final=min(s.final, ss.final))
                    z.add(sss)
                elif mode == OverlapType.ORIGINAL:
                    if s == pre_z:
                        pass
                    else:
                        z.add(s)
                        pre_z = s
                elif mode == OverlapType.COMP_INCL:
                    if s.initial >= ss.initial and s.final <= ss.final:
                        z.add(s)
                    else:
                        pass    
                if s.final >= ss.final:
                    try:
                        ss = con_y.next()
                    except:
                        cont_loop = False
                else:
                    try:
                        s = con_self.next()
                    except:
                        cont_loop = False   
            # When the region have no overlap
            elif s <= ss:
                try:
                    s = con_self.next()
                except:
                    cont_loop = False
            elif s > ss:
                try:      
                    ss = con_y.next()
                except:
                    cont_loop = False
        z.sort()
        return z
    
    def closest(self,y):
        """Return a new GenomicRegionSet including the region(s) of y which is closest to any self region. 
           If there are intersection, return False.
        
        Keyword arguments:
        y -- the GenomicRegionSet which to compare with
        
        Return:
        z -- the region(s) which is nearest to self
        
        """
        if len(self) == 0 or len(y) == 0:
            return GenomicRegionSet('Empty set') 
        elif len(self.intersect(y)) != 0:
            return False
        else:
            z_dict = {}  # For storing the distance and the regions
            self.sort()
            y.sort()
            con_self = self.__iter__()
            con_y = y.__iter__()
            s = con_self.next()
            ss = con_y.next()
            cont_loop = True
            while cont_loop:
                if s.chrom == ss.chrom:
                    # ----
                    #        ------
                    if s < ss:
                        z_dict[ss.initial - s.final] = ss
                        try: s = con_self.next()
                        except: cont_loop = False
                    #          ------
                    #  -----
                    elif s > ss:
                        z_dict[s.initial - ss.final] = ss
                        try: ss = con_y.next()
                        except: cont_loop = False
                elif s.chrom != ss.chrom:
                    if s < ss:
                        try: s = con_self.next()
                        except: cont_loop = False
                    elif s > ss:
                        try: ss = con_y.next()
                        except: cont_loop = False
            
            if len(z_dict.keys()) == 0:
                return GenomicRegionSet('Empty set')
            else:
                minimum_distance = min(z_dict.keys())
                z = GenomicRegionSet('Closest region')
                z.add(z_dict[minimum_distance])
                return z
            
    def remove_duplicates(self):
        """Remove the duplicate regions and remain the unique regions. (No return)"""
        self.sort()
        for i in range(len(self.sequences) - 1):
            try:
                if self.sequences[i] == self.sequences[i+1]:
                    del self.sequences[i+1]
            except:
                continue
            
    def window(self,y,adding_length = 1000):
        """Return the overlapping regions of self and y with adding a specified number 
        (1000, by default) of base pairs upstream and downstream of each region in self. 
        In effect, this allows regions in y that are near regions in self to be detected.
        
        Keyword arguments:
        y -- the GenomicRegionSet which to compare with
        adding_length -- the length of base pairs added to upstream and downstream of self (default 1000)
        
        Return:
        z -- a GenomicRegionSet including the regions of overlapping between extended self and original y.
        
        """
        if len(self) == 0 or len(y) == 0:
            return GenomicRegionSet('None region')
        # Establish an extended GenomicRegionSet
        extended_self = deepcopy(self)
        for i in extended_self:
            i.extend(adding_length,adding_length)
        # Find their intersections
        return extended_self.intersect(y)
    
    def subtract(self,y):
        """Return a GenomicRegionSet excluded the overlapping regions with y.
        
        Keyword arguments:
        y -- the GenomicRegionSet which to subtract by
        
        Return:
        z -- the remaining regions of self after subtraction
        
        Graphical explanation:
        self     ----------              ------
        y               ----------                    ----
        Result   -------                 ------
        
        """
        if len(self) == 0:
            return GenomicRegionSet(self.name + '-' + y.name)
        elif len(y) == 0:
            return self
        else:
            z = GenomicRegionSet(self.name + '-' + y.name)
            self.sort()
            y.sort()
            con_self = self.__iter__()
            con_y = y.__iter__()
            s = con_self.next()
            ss = con_y.next()
            cont_loop = True
            while cont_loop:
                # When the regions overlap
                if s.overlap(ss):
                    if s.initial < ss.initial:
                        z.add(GenomicRegion(chrom=s.chrom, initial=s.initial, final=ss.initial))
                    if s.final > ss.final:
                        s.initial = ss.final
                        try:
                            ss = con_y.next()
                            pass
                        except:
                            z.add(GenomicRegion(chrom=s.chrom, initial=ss.final, final=s.final))
                            try:
                                s = con_self.next()
                                pass
                            except:
                                cont_loop = False
                    else:
                        try:
                            s = con_self.next()
                            pass
                        except:
                            cont_loop = False 
                # When the region have no overlap
                elif s < ss:
                    z.add(GenomicRegion(chrom=s.chrom, initial=s.initial, final=s.final))
                    try:
                        s = con_self.next()
                        pass
                    except:
                        try:
                            ss = con_y.next()
                            pass
                        except:
                            cont_loop = False
                elif s > ss:
                    try:
                        ss = con_y.next()
                        pass
                    except:
                        z.add(GenomicRegion(chrom=s.chrom, initial=s.initial, final=s.final))
                        try:
                            s = con_self.next()
                            pass
                        except:
                            cont_loop = False
            z.sort()
            return z
    
    def subtract_aregion(self,y):
        """Return a GenomicRegionSet excluded the overlapping regions with y.
        
        Keyword arguments:
        y -- the GenomicRegion which to subtract by
        
        Return:
        result -- the remaining regions of self after subtraction
        
        Graphical explanation:
        self     ----------              ------
        y               ----------
        Result   -------                 ------
        
        """
        if len(self) == 0:
            return GenomicRegionSet('None region')

        else:
            z = GenomicRegionSet('Subtracted RegionSet')
            z.add(y)
            result = self.subtract(z)
            return result
    
    def merge(self):
        """Merge the regions within the GenomicRegionSet"""
        self.sort()
        
        if len(self.sequences) in [0, 1]:
            return
        else:
            z = GenomicRegionSet(name = self.name)
            prev_region = self.sequences[0]
            for cur_region in self.sequences[1:]:
                if prev_region.overlap(cur_region):
                    prev_region.initial = min(prev_region.initial, cur_region.initial)
                    prev_region.final = max(prev_region.final, cur_region.final)
                    prev_region.data += '_$_' + cur_region.data #use extra character to distinguish data sets
                else:
                    z.add(prev_region)
                    prev_region = cur_region
            z.add(prev_region)
            self.sequences = z.sequences
    
    def combine(self,region_set):
        """ Adding another GenomicRegionSet without merging the overlapping regions. """
        self.sequences.extend(region_set.sequences)
        #self.name = self.name + "+" + region_set.name
        self.sorted = False
        
    def cluster(self,max_distance):
        """Cluster the regions with a certain distance and return the result as a new GenomicRegionSet.
        
        Keyword arguments:
        max_distance -- the maximum distance between regions within the same cluster
        
        Return:
        z -- a GenomicRegionSet including clusters
        
        Graphical explanation:
        self           ----           ----            ----
                          ----             ----                 ----
        Result(d=1)    -------        ---------       ----      ----
        Result(d=10)   ---------------------------------------------        
        
        """
        self.sort()
        if len(self) == 0:
            return GenomicRegionSet('None region')
        elif len(self) == 1:
            return self
        else:
            z = GenomicRegionSet('Clustered region set')
            previous = self.sequences[0]
            for s in self.sequences[1:]:
                s_ext = deepcopy(s)
                s_ext.extend(max_distance,max_distance)
                if s_ext.overlap(previous):
                    previous.initial = min(previous.initial, s.initial)
                    previous.final = max(previous.final, s.final)
                else:
                    z.add(previous)
                    previous = s
            z.add(previous)
            return z
        
    def flank(self,size):
        """Return two flanking intervals with given size from both ends of each region.
        
        Keyword arguments:
        size -- the length of flanking intervals (default = SAME length as the region)
        
        Return:
        z -- a GenomicRegionSet including all flanking intervals
        
        Graphical explanation:
        self        -----           --            ---
        Result -----     -----    --  --       ---   ---
        
        """
        if len(self) == 0:
            return GenomicRegionSet("Empty")
        else:
            z = GenomicRegionSet("Flanking intervals")
            for s in self:
                s1 = GenomicRegion(name='upstream',chrom=s.chrom,
                                   initial=max(0, s.initial - size),
                                   final=s.initial)
                s2 = GenomicRegion(name='downstream',chrom=s.chrom,
                                   initial=max(0, s.final),
                                   final=s.final + size)  # Adding the limit of chromosome length
                z.add(s1)
                z.add(s2)
            return z
    
    def jaccard(self,query):
        """Return jaccard index, a value of similarity of these two GenomicRegionSet
        
        Keyword arguments:
        query -- the GenomicRegionSet which to compare with
        
        Return:
        similarity -- (Total length of overlapping regions)/(Total length of original regions)
        
        Graphical explanation:
        self           --8--      ---10---    -4-
        query     ---10---             ---10---
        intersect      -5-             -4-    2
        similarity:   ( 5 + 4 + 2)/[(8 + 10 + 4) + (10 +10) - (5 + 4 + 2)]
                      = 11/31
        
        """
        intersects = self.intersect(query)
        #print(intersects.total_coverage(),self.total_coverage(), query.total_coverage(),sep="\t")
        similarity = intersects.total_coverage()/(self.total_coverage() + query.total_coverage() - intersects.total_coverage())
        return similarity
    
    def within_overlap(self):
        """ Check whether there is overlapping within or not. """
        refer_posi = GenomicRegion(name="reference",chrom="chr1",initial=0,final=0)
        self.sort()
        for s in self:
            if s.overlap(refer_posi):
                return True
            refer_posi = s
        return False

    def total_coverage(self):
        """ Return the sum of all lengths of regions. """
        length = 0
        for s in self:
            length = length + len(s)
        return length
    
    def get_genome_data(self,organism, chrom_X=False, chrom_M=False):
        """ Add genome data from database into the GenomicRegionSet. """
        genome = GenomeData(organism)
        chromosome_file = open(genome.get_chromosome_sizes(),'r')
        for line in chromosome_file:
            if chrom_X and "chrX" in line and "random" not in line:
                chrom_region = GenomicRegion(chrom=line.split("\t")[0],
                                             initial=0,
                                             final=int(line.split("\t")[1]))
                self.add(chrom_region)
                continue
            elif chrom_M and "chrM" in line and "random" not in line:
                chrom_region = GenomicRegion(chrom=line.split("\t")[0],
                                             initial=0,
                                             final=int(line.split("\t")[1]))
                self.add(chrom_region)
                continue
            elif "random" not in line and "_" not in line and "chrX" not in line and "chrM" not in line: 
                chrom_region = GenomicRegion(chrom=line.split("\t")[0],
                                             initial=0,
                                             final=int(line.split("\t")[1]))
                self.add(chrom_region)
                continue
            else:
                continue
        chromosome_file.close()
        
    def random_regions(self, organism, total_size=None, multiply_factor=1, 
                       overlap_result=True, overlap_input=True, chrom_X=False, chrom_M=False, 
                       filter_path=None):
        """Return a GenomicRegionSet which contains the random regions generated by given entries and given number on the given organism.
        
        Keyword arguments:
        organism -- which organism's genome to use. (hg19, mm9)
        total_size -- Given the number of result random regions.
        multiply_factor -- this factor multiplies to the number of entries is the number of exporting random regions.
                           ** total_size has higher priority than multiply_factor. **
        overlap_result -- the results whether overlap with each other or not. (True/False)
        overlap_input -- the results whether overlap with input entries or not. (True/False)
        chrom_M -- the result should cover mitochondria chromosome or not. (True/False)
        
        Return:
        z -- a GenomicRegionSet which contains the random regions
        
        """

        # Fetching the chromosome length from data
        chrom_map = GenomicRegionSet("chrom_map")
        chrom_map.get_genome_data(organism, chrom_X=False, chrom_M=False)
        chrom_list = chrom_map.get_chrom()
        if filter_path:
            filter_map = GenomicRegionSet('filter')
            filter_map.read_bed(filter_path)
            chrom_map = chrom_map.subtract(filter_map)
            #print("chrom_map length: ", len(chrom_map.sequences))
            
        # Defining input_map, result_list (containing lengths of all result regions)
        result_list = []
        input_map = self
        input_num = len(self)
        input_list = [len(x) for x in self.sequences] # Contain a list of length of inputs
            
        if type(total_size) == int:
            for i in range(total_size):
                result_list.append(input_list[i % input_num])
        elif multiply_factor > 0:
            for i in range(int(multiply_factor * input_num)):
                result_list.append(input_list[i % input_num])
        
        # Dealing with positions of random regions
        def list_max(chrom_list, result_map):
            """ Generate a list containing maximum length for each chromosome """
            map_max = [] # Store the minimum length corresponding to chrom_list
            for ch in chrom_list:
                map_max.append(max([len(s) for s in result_map.any_chrom(ch)])) 
            return map_max
        
        def randoming(result_map, chrom_list,result_list):
            """ Return the new GenomicRegion as the result of randomization. """
            cont_loop = True
            candidate_chrom = [chrom_list[i] for i,l in enumerate(map_max) if l > length] # screen for the potential chromosomes which has space
            #print(candidate_chrom)
            while cont_loop:
                sample = random.choice(result_map.any_chrom(random.choice(candidate_chrom), len_min=length))
                try:
                    random_posi = random.randint(sample.initial, sample.final - length)
                    #print(sample.chrom)
                    cont_loop = False
                except:
                    pass
            return GenomicRegion(chrom=sample.chrom,initial=random_posi,final=random_posi + length)
        
        ###################################################################
        z = GenomicRegionSet(name="random regions")
        if overlap_input == False and overlap_result == False:
            result_map = chrom_map.subtract(input_map) ### overlap_input == False
            
            map_max = list_max(chrom_list, result_map)
            
            for length in result_list:
                new_region = randoming(result_map, chrom_list, result_list)
                z.add(new_region)
                result_map = result_map.subtract_aregion(new_region) ### overlap_result == False
                map_max[chrom_list.index(new_region.chrom)] = max([len(s) for s in result_map.any_chrom(new_region.chrom)])
                
        elif overlap_input == False and overlap_result == True:
            result_map = chrom_map.subtract(input_map) ### overlap_input == False
            map_max = list_max(chrom_list, result_map)
            
            for length in result_list:
                new_region = randoming(result_map, chrom_list, result_list)
                z.add(new_region)
        
        elif overlap_input == True and overlap_result == False:
            result_map = chrom_map
            map_max = list_max(chrom_list, result_map)
            
            for length in result_list:
                new_region = randoming(result_map, chrom_list, result_list)
                z.add(new_region)
                result_map = result_map.subtract_aregion(new_region) ### overlap_result == False
                map_max[chrom_list.index(new_region.chrom)] = max([len(s) for s in result_map.any_chrom(new_region.chrom)])
                
        elif overlap_result == True and overlap_result == True:
            result_map = chrom_map
            map_max = list_max(chrom_list, result_map)
            
            for length in result_list:
                new_region = randoming(result_map, chrom_list, result_list)
                z.add(new_region)
        return z
    
    #hallo
    def projection_test(self, query, organism, extra=None, background=None):
        """" Return the p value of binomial test. """
        chrom_map = GenomicRegionSet("Genome")
        chrom_map.get_genome_data(organism=organism)
        #print("coverage of reference: ",self.total_coverage(),"\tcoverage of genome: ",chrom_map.total_coverage())
        if background: #backgound should be a GenomicRegionSet
            ss = self.intersect(background, OverlapType.OVERLAP)
            possibility =  ss.total_coverage() / background.total_coverage()
        else:
            possibility = self.total_coverage() / chrom_map.total_coverage() # The average likelihood

        nquery = query.relocate_regions(center='midpoint', left_length=0, right_length=0)
        intersect_regions = self.intersect(nquery,mode=OverlapType.OVERLAP)
        n = len(nquery)
        k = len(intersect_regions)
        #print("intersections: ",k,"\tnumber of query",n,"\tgenetic coverage: ",possibility)
        p = float(stats.binom_test(k, n, possibility))
        if extra:
            try: return possibility, k/n, p  # for the case n = 0
            except: return possibility, 0, p
        else:
            return p

    def any_chrom(self,chrom,len_min=False, len_max=False):
        """ Return a list of regions which belongs to given chromosome. (used in random_regions) """
        if len_min == False and len_max == False:
            return [s for s in self if s.chrom == chrom] 
        elif type(len_min) == int and len_max == False:
            return [s for s in self if s.chrom == chrom and s.final - s.initial >= len_min] 
        elif type(len_max) == int and len_min == False:
            return [s for s in self if s.chrom == chrom and s.final - s.initial <= len_max] 
        else:
            return [s for s in self if s.chrom == chrom and len_min <= s.final - s.initial <= len_max] 
        
    def relocate_regions(self,center='midpoint', left_length=2000, right_length=2000):
        """ Return a new GenomicRegionSet which relocates the regions by given center and extend length.
        
        Parameters:
        center:    
            midpoint  -- locate the new region's center as original region's midpoint
            leftend      -- locate the new region's center as original region's 5' end (if no orientation information, default is left end)
            rightend      -- locate the new region's center as original region's 3' end (if no orientation information, default is right end)
            bothends  -- locate the new region's center as original region's both ends
        left_length:
        right_length:
            The length to extend from the center
        """
        new_regions = GenomicRegionSet('new'+self.name)
        for r in self.sequences:
            # Define the position
            if center == 'midpoint':
                mp = int(0.5*(r.initial + r.final))
                nr = GenomicRegion(chrom=r.chrom, initial=mp, final=mp)
            elif center == 'leftend':
                nr = GenomicRegion(chrom=r.chrom, initial=r.initial, final=r.initial)
            elif center == 'rightend':
                nr = GenomicRegion(chrom=r.chrom, initial=r.final, final=r.final)
            elif center == 'bothends':
                nrL = GenomicRegion(chrom=r.chrom, initial=r.initial, final=r.initial) # newbed on left end
                nrR = GenomicRegion(chrom=r.chrom, initial=r.final, final=r.final) # newbed on right end
            try:new_regions.add(nr)
            except: 
                new_regions.add(nrL)
                new_regions.add(nrR)
                
        # Extend the region
        new_regions.extend(left_length, right_length)
        return new_regions

    def maximum_length(self):
        """Return the length of the maximum region from the GenomicRegionSet. """
        try:
            maxl = len(self.sequences[0])
            for s in self:
                maxl = max(len(s), maxl)
            return maxl
        except:
            print("There is no region in the given GenomicRegionSet, {}".format(self.name))
            return 0 
