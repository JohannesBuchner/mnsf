from DataRead import DataRead
from PHAMaker import PHAMaker

import astropy.io.fits as fits
from astropy.table import Table
from numpy import logical_and, array, mean, histogram, arange

import matplotlib.pyplot as plt
from binning.step import Step

from binning.tteBinning import tteBinning
from glob import glob


#from phaMake.phaMake import phaMake

import os
import errno
from random import randint


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


detLU = {"NAI_00":"n0",\
         "NAI_01":"n1",\
         "NAI_02":"n2",\
         "NAI_03":"n3",\
         "NAI_04":"n4",\
         "NAI_05":"n5",\
         "NAI_06":"n6",\
         "NAI_07":"n7",\
         "NAI_08":"n8",\
         "NAI_09":"n9",\
         "NAI_10":"na",\
         "NAI_11":"nb",\
         "BGO_00":"b0",\
         "BGO_01":"b1"}


class GBMReader(DataRead):


    
    def ReadData(self,bkgIntervals=[[-20,-.1],[50.,300.]],rspFile=None):
        '''
        Momentarily GBM specific 
        
        '''

        
        self.trigTime=self.data[0].header["TRIGTIME"]
        self.tstart = self.data[0].header["TSTART"]-self.trigTime
        self.tstop = self.data[0].header["TSTOP"]-self.trigTime
        
     
        

        self.instrument="GBM"
        self.det = detLU[self.data[0].header["DETNAM"]]
        


        directory = ""

        for x in self.dataFile.split('/')[:-1]:
            directory+=x+"/"
            
        self.directory = directory

        if rspFile == None:
        
            rsp = glob(directory+"glg_cspec_"+self.det+"*.rsp")[0]
            
        else:

            rsp = rspFile

        print "Found RSP: "+rsp
        self.rsp = rsp
            

            
        self.emin = self.data[1].data["E_MIN"]
        self.emax = self.data[1].data["E_MAX"]

        chans = array(zip(self.emin,self.emax))


        meanChan = array(map(mean,chans))

        self.meanChan = meanChan
        #self.chans = chans


        self.chanWidth = self.emax-self.emin
        
        self.bkgIntervals = bkgIntervals
        
        #phas =self.data[2]["PHA"]
        
        #tt = logical_and(phas>=chans[0],phas<=chans[1])
        self.evts = self.data[2].data["TIME"] - self.trigTime #Filter chans
        
        
        #get the background fit
        
        tb = tteBinning(self.dataFile,int(self.tstart),int(self.tstop),self.bkgIntervals)
        self.dataBinner = tb
        


    def CreateCounts(self,start=0,stop=10,pha=False):
        
        
        self.dataBinner.MakeBackgroundSelectionsForDataBinner()
        self.dataBinner._FitBackground()
        self.bkgMods = self.dataBinner.polynomials

        self.thisStart = start
        self.thisStop  = stop 

        #GO BY TIME BIN

        
        bins = self.dataBinner.bins
        j=0
        for i in range(len(bins)-1):
            
            lob=bins[i]
            hib=bins[i+1]
            
            if lob>=start and hib<=stop:
                
                bkgCounts = []
                bkgError = []
            
                totalCounts = []
            
                for ch in range(128):
                
                    tt = self.data[2].data["PHA"] == ch
                
                
                    ## get evts between times:
                    tt2 = self.evts >= lob
                    tt2 = logical_and(tt2,self.evts <hib)
                
                    tt= logical_and(tt,tt2)
                
                    #Num total counts
                
                    totalCounts.append(len(self.evts[tt]))
                    bkgCounts.append(self.bkgMods[ch].integral(lob,hib))
                    bkgError.append(self.bkgMods[ch].integralError(lob,hib))
                totalCounts =array(totalCounts)/(hib-lob)
                bkgCounts=array(bkgCounts)/(hib-lob)
                bkgError = array(bkgError)/(hib-lob)
                sourceCounts  = totalCounts  - bkgCounts


                directory = "bin%d"%j

                
                
                        
                tab = Table(array(zip(totalCounts,sourceCounts,bkgCounts,bkgError,self.emin,self.emax,self.meanChan)),names=["total","source","bkg","berr","emin","emax","meanChan"])
                tab.meta={"duration":hib-lob,"INST":self.instrument,"DET":self.det,"RSP":self.rsp,"TMIN":lob,"TMAX":hib}
                
                


                tab.meta["binN"]=directory
                tab.meta["fileLoc"] = self.pathExt+self.directory+directory+"/"
                tab.meta["file"]=self.dataFile.split('/')[-1]
                self.binDict[directory] = tab


                if pha:
                    make_sure_path_exists(self.directory+directory)
                    exposure = hib-lob # Need to add deadtime????
                    phaFile = PHAMaker(self.pathExt+self.directory+directory+"/",self.det,exposure,self.rsp,self.emin,self.emax)
                    phaFile.SetCounts(totalCounts*exposure)
                    phaFile.SetBak(bkgCounts,bkgError)
                    phaFile.Write()

                
                j+=1
        self.PlotData()

    def PlotData(self,save = None):

        tBins = []
        for i in range(len(self.dataBinner.bins)-1):
    
            tBins.append([self.dataBinner.bins[i],self.dataBinner.bins[i+1]])
        tBins=array(tBins)


        cnts,_ = histogram(self.dataBinner.evts,bins=self.dataBinner.bins)

        maxCnts = max(cnts/self.dataBinner.binWidth)
        minCnts = min(cnts/self.dataBinner.binWidth) 
        
        fig=plt.figure(randint(1,1E6))
        ax=fig.add_subplot(111)

        Step(ax,tBins,cnts/self.dataBinner.binWidth,"k",.5)


        #Plot the background region

        cnts,_ = histogram(self.dataBinner.filteredEvts,bins=self.dataBinner.bins)
        
        
        Step(ax,tBins,cnts/self.dataBinner.binWidth,"b",.7)


        #Plot the selection region

        
        ax.vlines(self.thisStart,minCnts-50,maxCnts+50,colors='limegreen',linewidth=1.5)
        ax.vlines(self.thisStop,minCnts-50,maxCnts+50,colors='limegreen', linewidth=1.5)
        

        bkg = []
        oneSecBins =arange(self.dataBinner.bins[0],self.dataBinner.bins[-1],1.) 
        for i in range(len(oneSecBins)-1):
    
            b=0
            for j in range(len(self.bkgMods)):
        
                b+= self.bkgMods[j].integral(oneSecBins[i],oneSecBins[i+1])/(oneSecBins[i+1]-oneSecBins[i])
            bkg.append(b)
        meanT = map(mean,zip(oneSecBins[:-1],oneSecBins[1:]))
        ax.plot(meanT,array(bkg),linewidth=2,color="r")

        minX = min(map(lambda x: x[0],self.bkgIntervals))
        maxX = max(map(lambda x: x[1],self.bkgIntervals))


        ax.set_ylim(bottom=minCnts-50.,top=maxCnts+50.)
        ax.set_xlim(left=minX,right=maxX)
        
        if save:
            fig.savefig(save,bbox_inches="tight")
