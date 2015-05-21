import astropy.io.fits as fits
from astropy.table import Table

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
    
    def __init__(self,dataFile,pathExt=""):
        
        self.dataFile = dataFile
        self.pathExt = pathExt
        self.data = fits.open(dataFile)
        self.binDict = dict()



    def MakeBins(self):

        for key in self.binDict.keys():
            
            tab = self.binDict[key]
            make_sure_path_exists(self.pathExt+self.directory+key)
            self.tab = tab 
            print "Writing:\n\t%s"%self.pathExt+self.directory+key+"/"+self.instrument+"_"+self.det+".fits"
            tab.write(self.pathExt+self.directory+key+"/"+self.instrument+"_"+self.det+".fits",format="fits",overwrite=True)

            
            self.InfoTxt(self.pathExt+self.directory+key,tab)
            
    def InfoTxt(self,dir,tab):


        f = open("%s/binInfo.txt"%dir,"w")

        # Write time bins
        f.write("TSTART:\t%.2f\n"%tab.meta["TMIN"])
        f.write("TSTOP:\t%.2f\n"%tab.meta["TMAX"])

        # Write Duration
        f.write("DURATION:\t%.2f\n"%tab.meta["duration"])
        
        # Write Total Counts
        f.write("COUNTS:\t%.1f\n"%(tab["total"]).sum())
        
        # Write Background
        f.write("BKG:\t%.1f\n"%(tab["bkg"]).sum())
        f.close()
        
        
      
        
        
                
        


        

        
