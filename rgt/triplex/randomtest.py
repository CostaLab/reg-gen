from rgt.RNADNAInteractionSet import *

class RandomTest:
    def __init__(self, txp_path, organism):
    	self.txp = RNADNAInteractionSet(organism=organism, filename=txp_path)
    	
    def target_count(self):
    	"""Count the number of TFFs for each TFO

    	The count number is stored in a dictionary as self.merged_TFO
    	"""
    	self.txp.merged_TFO()
    	self.merged_TFO = OderedDict()
    	for tfo, tffs in iteritems(self.txp.merged):
    		self.merged_TFO[tfo] = len(tffs)
        
    def randomization(repeat):

