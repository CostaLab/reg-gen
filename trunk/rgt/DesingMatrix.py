possible_types=["genes","regions","reads"]

class ExperimentalMatrix:

    def __init__(self):
        self.fields=[]
        self.fieldsDict={}
        self.names=[]
        self.files=[]
        self.types=[]
        self.objects={}
    
    def read(self,file):
        f=open(file)
        header=f.readline()
        header=header.strip("\n")
        header=header.split("\t")
        assert(header[0]=="name")
        assert(header[1]=="type")
        assert(header[2]=="file")
        self.fields=header
        for fi in range(3,len(fields)):
            self.fieldsDict[header[fi]]={}
        self.filedsDict=
        for line in f:
            line=line.strip("\n")
            line=line.split("\t")
            self.names.append(line[0])
            self.files[line[0]]=line[2]
            self.types[line[0]]=line[1]
            for fi in range(3,len(fields)):
                dict=self.fieldsDict[header[fi]]
                try:
                    dict[l[fi]].append(line[0])
                except:
                    dict[l[fi]]=[line[0]]
        self.loadObjects()

    def loadObjects():
        for i,t in enumerate(self.type):
            if t == "regions":
                bed = SetGenomicRegions(self.files[i])
                bed.readBed(self.files[i])
                self.objects[self.names[i]]=bed
            if t == "genes":
                genes= SetGenes()
                genes.read(self.files[i])
                self.objects[self.names[i]]=genes

            
