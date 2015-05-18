from DataRead import DataRead
import astropy.io.fits as fits
from numpy import array
from numpy import mean
from astropy.table import Table


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

class PHAReader(DataRead):



    def SetBAKFile(self, bakFile):



        f= fits.open(bakFile)

        try:
            self.bkg =f[1].data["COUNTS"]
            self.berr = f[1].data["COUNTS"]
        except KeyError:
            try:
                self.bkg =f[1].data["RATE"][0]*f[1].data["EXPOSURE"][0] #Turn it back into counts
                self.berr = f[1].data["STAT_ERR"][0]*f[1].data["EXPOSURE"][0]
            except KeyError:
                try:
                    self.bkg =f[1].data["RATE"]*f[1].header["EXPOSURE"] #Turn it back into counts
                    self.berr = f[1].data["STAT_ERR"]*f[1].header["EXPOSURE"]
                except KeyError:
                    self.bkg =f[2].data["RATE"]*f[2].header["EXPOSURE"] #Turn it back into counts
                    self.berr = f[2].data["STAT_ERR"]*f[2].header["EXPOSURE"]
        f.close()


    def SetRSP(self,rsp):

        self.rsp  = rsp
        
    def CreateCounts(self,ext="bin0"):


        try:
            duration = self.data[1].header["EXPOSURE"]
        except KeyError:
            
            duration = self.data[1].data["EXPOSURE"][0]
            
        self.bkg = self.bkg/duration


        totalCounts = self.data[1].data["COUNTS"]/duration

        if len(totalCounts) == 1:
            totalCounts = totalCounts[0]

        sourceCounts = totalCounts - self.bkg
        
        self.berr = self.berr/duration 


        emin = self.data[2].data["E_MIN"]
        emax = self.data[2].data["E_MAX"]
        try:
            self.instrument = self.data[0].header["INSTRUME"]
        except KeyError:
            self.instrument = self.data[1].header["INSTRUME"]
        if self.instrument == 'LAT':
            try:
                self.det=self.data[0].header["DATATYPE"]
            except KeyError:
                try:
                    self.det=self.data[1].header["DETNAM"]
                except KeyError:
                    self.det="STD"
        else:
            try:
                self.det = detLU[self.data[0].header["DETNAM"]]
            except KeyError:
                self.det = detLU[self.data[1].header["DETNAM"]]
        
        chans = array(zip(emin,emax))


        meanChan = array(map(mean,chans))

        
        directory = ""

        for x in self.dataFile.split('/')[:-1]:
            directory+=x+"/"
        
        self.directory = directory
        
        tab = Table(array(zip(totalCounts,sourceCounts,self.bkg,self.berr,emin,emax,meanChan)),names=["total","source","bkg","berr","emin","emax","meanChan"])

        
        try:
            trigTime = self.data[0].header["TRIGTIME"]
        except KeyError:
            try:
                trigTime = self.data[1].header["TRIGTIME"]
            except KeyError:
                trigTime = self.data[2].header["TRIGTIME"]
            
        try:
            tstart = self.data[1].data["TSTART"][0]-trigTime
        except KeyError:
            try:
                tstart = self.data[3].header["TSTART"]-float(trigTime)
            except IndexError:
                print "Cannot find TSTART!\n\twill be set to 0"
                tstart = 0.
        try:    
            tstop = tstart + self.data[1].data["TELAPSE"][0]
        except KeyError:
            try:
                tstop = self.data[3].header["TSTOP"]-float(trigTime)
            except IndexError:
                tstop = tstart+duration
            

        tab.meta={"duration":duration,"INST":self.instrument,"DET":self.det,"RSP":self.rsp,"TMIN":tstart,"TMAX":tstop}

        directory = ext
        tab.meta["file"]=self.dataFile.split('/')[-1]
        
                

                
        tab.meta["binN"]=directory
        tab.meta["fileLoc"] = self.directory+directory+"/"

        
        
        self.binDict[directory] = tab
