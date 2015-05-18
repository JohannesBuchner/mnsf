from rsp import rsp
import numpy as np
from joblib import Parallel, delayed

class RSPconvolve:


    def __init__(self,rspFile):


        thisRSP = rsp(rspFile)

        self.drm = thisRSP.drm
        self.channelE = thisRSP.channelE
        self.photonE = thisRSP.photonE


        self._CreatePhotonEvalEnergies()

        del thisRSP

        #energyChans = range(len(self.photonE)) ###Modify!
        

    def _CreatePhotonEvalEnergies(self):

        chanWidth = self.photonE[:,1]-self.photonE[:,0]
        resFrac511 = 0.2 #lifted from RMFIT!
        resExp = -0.15   #lifted from RMFIT!

        binCenter = np.array(map(np.mean,self.photonE))
        

        resFrac = resFrac511*(binCenter/511.)**resExp
        resFWHM = binCenter*resFrac
        numEchans = np.ones(len(binCenter))

        self.lowEval = chanWidth<resFWHM/2.
        
        self.medEval = chanWidth>=resFWHM/2.
        numEchans[self.medEval]=3.
        self.highEval = chanWidth/2.>=resFWHM/3.
        numEchans[self.highEval]=7.

        self.lowEne = binCenter[self.lowEval]
        self.medEne=np.array(map(lambda x,y: [x-0.333333*y,x,x+0.333333*y] ,binCenter[self.medEval],chanWidth[self.medEval]))
        self.highEne=np.array(map(lambda x,y: [x-0.5*y,x-0.333333*y,x-0.16667*y,x,x+0.16667*y,x+0.333333*y,x-0.5*y] ,binCenter[self.highEval],chanWidth[self.highEval]))
        
        self.chanWidth = chanWidth
            

    def SetModel(self,model):
        '''
        Pass a function that is the spectral model for the fit


        '''

        self.model = model

        return



    def SetModelVec(self,cnts):

        self.vec = cnts*self.chanWidth
    
    def CreateModelVector(self):
        '''
        Call the Cython integrator to make a vector that will 
        be convolved with the DRM.


        This function is slow. It has been updated to ignore the 
        full integrator

        '''

        
        self._ConvolveMatrix()


    def _ConvolveMatrix(self):

        self.counts = np.dot(self.drm.transpose(),self.vec)


    def GetCounts(self):

        return self.counts


    
