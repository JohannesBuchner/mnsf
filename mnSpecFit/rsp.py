import astropy.io.fits as pf
from numpy import zeros, array, matrix, abs
import math

class rsp(object):

    def __init__(self,rspFileName):


        rspFile = pf.open(rspFileName)
    
        #angle = rspFile[2].header["DET_ANG"]
#        print "Opening "+rspFileName
#        print "Read RSP file with angle: %lf"%angle
        
        

        self.fileName = rspFileName.split('/')[-1]

        
        

        self.chanData = rspFile[1].data
        try:
            self.numEnergyBins = rspFile[2].header['NUMEBINS']
            self.numDetChans = rspFile[2].header['DETCHANS']
            self.det = rspFile[0].header['DETNAM']
        except KeyError:
            try:
                self.numEnergyBins = rspFile[1].header['NUMEBINS']
                self.numDetChans = rspFile[1].header['DETCHANS']
                self.det = rspFile[1].header['DETNAM']
            except KeyError:
            
            #print "This is an LLE response"
                self.numDetChans = rspFile[1].header['DETCHANS']
                self.numEnergyBins = rspFile[2].header["NAXIS2"]
        
        #self.beta = angle
        try:
            self.photonE = array(zip([rspFile[2].data["ENERG_LO"], rspFile[2].data["ENERG_HI"]]))
        except KeyError:
            self.photonE = array(zip([rspFile[1].data["ENERG_LO"], rspFile[1].data["ENERG_HI"]]))
        self.photonE = self.photonE.transpose()
        self.photonE = array(map(lambda x: x[0], self.photonE))

        try:
            self.channelE = (rspFile[1].data["E_MIN"], rspFile[1].data["E_MAX"])
        except KeyError:
            self.channelE = (rspFile[2].data["E_MIN"], rspFile[2].data["E_MAX"])
        #Main component of object
        self.drm = zeros((self.numEnergyBins, self.numDetChans))

        self._ConstructDRM(rspFile)
        try:
            self.binWidths = rspFile[1].data["E_MAX"]-rspFile[1].data["E_MIN"]
        except KeyError:
            self.binWidths = rspFile[2].data["E_MAX"]-rspFile[2].data["E_MIN"]


        rspFile.close()

        


    def _ConstructDRM(self,fitFile):
        '''
        Construct the drm from the fits file
        '''
    
        try:
            mData = fitFile[2].data
            self.phtBins = array(zip(mData['ENERG_LO'],mData['ENERG_HI']))
        except KeyError:
            mData = fitFile[1].data
            self.phtBins = array(zip(mData['ENERG_LO'],mData['ENERG_HI']))
        tmp1 = mData["F_CHAN"]
        tmp2 = mData["N_CHAN"]

        for fcs,ncs,i in zip( tmp1 , tmp2  ,range(self.numEnergyBins)):
            colIndx = 0
            
            try:
                for fc,nc in zip(fcs,ncs):

                    self.drm[i,fc-1:fc+nc]=mData["MATRIX"][i][colIndx:colIndx+nc]
                    colIndx+=nc
            except TypeError: #Prolly not formatted correctly
                self.drm[i,fcs-1:fcs-1+ncs]=mData["MATRIX"][i][colIndx:colIndx+ncs]
                colIndx+=ncs

        self.drm=matrix(self.drm)
        #self.drm=self.drm.transpose()
        del mData

        

    def _FindNearestPhotonBin(self, e):
        ''' 
       	Take an energy (e) and find return the 
	corresponding column from the DRM which 
	has been normalised by the geometric area
	'''
        
        idx=(abs(self.photonE[0]-e)).argmin()
        return idx
