from astropy.table import Table
from spectralTools.step import Step
from numpy import array
import matplotlib.pyplot as plt

class DataBin:
    '''

    '''
    
    def __init__(self,databin):
        '''
        ____________________________________________________________
        arguments:
        databin: path to a databin that has been saved
        '''    
        self.data = Table.read(databin,format="fits")
        self.berr = self.data['berr']
        self.bkg = self.data["bkg"]
        self.source = self.data["source"]
        self.total = self.data["total"]
        #self.chans = self.data["chans"]
        self.meanChan = self.data["meanChan"]
        self.totBkg = self.bkg.sum()
        self.totSource =self.source.sum()
        self.totTot = self.total.sum()
        self.duration = self.data.meta["DURATION"]
        self.det = self.data.meta["DET"]
        self.rsp = self.data.meta["RSP"]
        self.instrument = self.data.meta["INST"]
        self.chanMin = self.data["emin"]
        self.chanMax = self.data["emax"]
        self.binN = self.data.meta["BINN"]
        self.fileLoc = self.data.meta["FILELOC"]
        self.tmin = self.data.meta["TMIN"]
        self.tmax = self.data.meta["TMAX"]
        self.activeLoChan = 0
        self.activeHiChan = len(self.chanMax)-1
        self.file=self.data.meta["FILE"]

        #Remove the negative counts and replace with zeros || EXPERIMENTNAL

  #      self._RemoveNegativeCounts(self.source)
  #      self._RemoveNegativeCounts(self.bkg)
        


    def _GetChannel(self,energy):
        '''
        Private function that finds the channel for a given energy

        ____________________________________________________________
        arguments:
        energy: selection energy in keV

        '''

        if energy < self.chanMin[0]:
            return 0
        elif energy > self.chanMax[-1]:
            return len(self.chanMax)-1
    

        
        ch = 0
        for lo, hi in zip(self.chanMin,self.chanMax):

            if energy >= lo and energy <= hi:
                return ch
            else:
                ch+=1

    def SetLoChan(self,lo):
        '''
        Set Lo energy channel that is used when getting rates
        and counts via the Get() functions.

        ____________________________________________________________
        arguments:
        lo: energy in keV

        '''
        self.emin = lo
        self.activeLoChan = self._GetChannel(lo)

    def SetHiChan(self,hi):
        '''
        Set HI  energy channel that is used when getting rates
        and counts via the Get() functions.

        ____________________________________________________________
        arguments:
        hi: energy in keV

        '''
        self.emax = hi
        self.activeHiChan = self._GetChannel(hi)

    def _RemoveNegativeCounts(self, spectrum):


        tt = spectrum < 0.

        spectrum[tt] = 0.


    def GetTotalCounts(self):

        return (self.total[self._energySelection])*self.duration

    def GetTotalRate(self):

        return (self.total[self._energySelection])

    def GetSourceCounts(self):

        return (self.source[self._energySelection])*self.duration

    def GetSourceRate(self):

        return (self.source[self._energySelection])


    def GetBkgCounts(self):

        return (self.bkg[self._energySelection])*self.duration

    def GetBkgRate(self):

        return (self.bkg[self._energySelection])


    def GetBkgErr(self):

        return (self.berr[self._energySelection])*self.duration

    def GetBkgErrRate(self):

        return self.berr[self._energySelection]
        #return (self.berr[self.activeLoChan:self.activeHiChan+1])

    def SelectEnergies(self,selection):

        # Make sure you've got an array
        selection = array(selection)
        self.selMins = []
        self.selMaxs = []

        if len(selection.shape) == 1:

                tmp = map(self._GetChannel,selection)

                tt = [False]*len(self.meanChan)
                
                tt[tmp[0]:tmp[1]+1] = [True]*(tmp[1]-tmp[0]+1)

                self.emin = min(selection)
                self.emax = max(selection)

                self.selMins.append(self.chanMin[tmp[0]])
                self.selMaxs.append(self.chanMax[tmp[1]])

        elif len(selection.shape) == 2:

                tmp = array(map(lambda x: [self._GetChannel(x[0]),self._GetChannel(x[1])] , selection))

                tt = [False]*len(self.meanChan)

                for x in tmp:

                    tt[x[0]:x[1]+1] = [True]*(x[1]+1-x[0])

                self.selMins = self.chanMin[tmp[:,0]]
                self.selMaxs = self.chanMax[tmp[:,1]]

                self.emin = selection.min()
                self.emax = selection.max()

                    
        tt = array(tt)


        self._energySelection = tt
    
    def ViewCountSpectrum(self):


        chans = array(zip(self.chanMin,self.chanMax))
        width = self.chanMax-self.chanMin
        fig = plt.figure(181)
        ax = fig.add_subplot(111)
        Step(ax,chans,self.total/width,"k",lw=1.)
        ax.set_xscale('log')
        ax.set_yscale('log')


        rate = self.total/width
        minRate = min(rate[rate>0.])
        
        ax.vlines(self.chanMin[self.activeLoChan],minRate,max(rate),color="r")
        ax.vlines(self.chanMax[self.activeHiChan],minRate,max(rate),color="b")

        ax.set_xlabel("Energy (keV)")
