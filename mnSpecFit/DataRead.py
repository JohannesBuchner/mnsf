import astropy.io.fits as fits
from astropy.table import Table
import pyqt_fit.kernel_smoothing as smooth
from numpy import logical_and, array, mean, histogram

import matplotlib.pyplot as plt
from spectralTools.step import Step

from spectralTools.binning.tteBinning import tteBinning
from glob import glob


import os
import errno

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

class DataRead:
    
    def __init__(self,dataFile):
        
        self.dataFile = dataFile
        self.data = fits.open(dataFile)
        self.binDict = dict()



    def MakeBins(self):

        for key in self.binDict.keys():
            
            tab = self.binDict[key]
            make_sure_path_exists(self.directory+key)
            self.tab = tab 
            print "Writing:\n\t%s"%self.directory+key+"/"+self.instrument+"_"+self.det+".fits"
            tab.write(self.directory+key+"/"+self.instrument+"_"+self.det+".fits",format="fits",overwrite=True)

            

        
        
        
      
        
        
                
        


        

        
