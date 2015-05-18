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

class MultiPHAReader(DataRead):



    def SetBAKFile(self, bakFile):



        f= fits.open(bakFile)

        specItr = f[1].data["SPEC_NUM"]-1
        self.bkg = []
        self.berr = []
        for i in specItr:
        

            self.bkg.append(f[1].data["RATE"][i]*f[1].data["EXPOSURE"][i]) #Turn it back into counts
            self.berr.append(f[1].data["RATE"][i]*f[1].data["EXPOSURE"][i])
            
        f.close()


    def SetRSP(self,rsp):

        self.rsp  = rsp
        
    def CreateCounts(self):

        specItr = self.data[1].data["SPEC_NUM"]-1

        for i in specItr:
            try:
                duration = self.data[1].header["EXPOSURE"]
            except KeyError:
            
                duration = self.data[1].data["EXPOSURE"][i]
            
            self.bkg[i] = self.bkg[i]/duration


            totalCounts = self.data[1].data["COUNTS"][i]/duration

        

            sourceCounts = totalCounts - self.bkg[i]
        
            self.berr[i] = self.berr[i]/duration 


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
        
            tab = Table(array(zip(totalCounts,sourceCounts,self.bkg[i],self.berr[i],emin,emax,meanChan)),names=["total","source","bkg","berr","emin","emax","meanChan"])



            trigTime = self.data[0].header["TRIGTIME"]
            tstart = self.data[1].data["TSTART"][i]-trigTime
            tstop = tstart + self.data[1].data["TELAPSE"][i]
        


            tab.meta={"duration":duration,"INST":self.instrument,"DET":self.det,"RSP":self.rsp,"TMIN":tstart,"TMAX":tstop}

            directory = "bin%d"%i

            

            tab.meta["file"]=self.dataFile.split('/')[-1]
                

                
            tab.meta["binN"]=directory
            tab.meta["fileLoc"] = self.directory+directory+"/"

        
        
            self.binDict[directory] = tab
