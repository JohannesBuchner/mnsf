from multiFit.mnfit import mnfit
from Model import Model
from models import models
from numpy import array
from astropy.table import Table
import json

class mnSpecFit(mnfit):


    def LoadData(self,dataBins):
        '''
        This member aranges the DataBins and prepares them
        for sampling. The DataBin objects are astropy Tables
        that are created by separate objects which ultimately 
        make a uniform interface for reading here.

        ___________________________________________________
        arguments: 
        databin:  DataBin or list of DataBins

        '''


        
        
        if type(dataBins)==list:
        
            self.detList = dataBins

        else:
            self.detList = [dataBins]

        self._dataLoaded = True #Mark that data are loaded
            
    def SetLikelihood(self,lhType):
        '''
        Pass a likelihood class. This is important because the method
        instantiates a likelihood object for each detector

        ____________________________________________________________
        arguments:
        lhType: a likelihood class

        '''
        
        self.lhs = []
        for x in self.detList:

            lh = lhType()  #Instantiate the likelihood
            lh.tb = x.duration  #Set the background duration
            lh.ts = x.duration  #Set the source duration
            self.lhs.append(lh)
       


    def SetEnergyBounds(self,detector, selection):
        '''
        Set the energy bounds of a detector

        ____________________________________________________________
        arguments:
        detector: str() ex "n6"
        lo: lo channel in keV
        hi: hi channel in keV

        '''


        
        for db in self.detList:

            if db.det == detector:

                db.SelectEnergies(selection)
                #db.SetHiChan(hi)   # Call DataBin method for channel 
                #db.SetLoChan(lo)   # selection

                #print "Detector %s ignoring channels (0-%d) and (%d-%d)"%(detector,db.activeLoChan,db.activeHiChan,len(db.total))
                print "%d active PHA channels"%len(db.GetTotalCounts())
                print
                return
        print "\n Detector: %s has not been loaded!  \n"%detector
            


    def SetSaveFile(self,savefile):
        '''
        Set the name of the json file to be created
        after the fit is made

        ____________________________________________________________
        arguments:
        savefile: str() file name


        '''
        self.savefile = savefile
        self._saveFileSet = True

    def SetModel(self, model):
        '''
        Pass a model class which will be instantiated
        and the rsp loaded for each DataBin

        _____________________________________________________________
        arguments:
        model: a derived Model class

        '''

        if type(model) == list:

            listModels = True
        else:
            listModels = False
            

        self.models = []

        for db in self.detList:


            if listModels:

                tmp1 = model[0]()
                for mod in model[1:]:
                    tmp2 = mod()
                    tmp1 = tmp1+tmp2

                myModel = tmp1
            else:
                myModel = model()
            
            #mod = model()   # Instantiate a model #NEW Now the model should be an object NOT a class!
            myModel.SetRSP(db.rsp) #Pass the rsp to model
            self.models.append(myModel)
        self.n_params = self.models[0].n_params

        pass

    def ConstructLikelihood(self):
        '''
        Provides a likelihood function based off the data and
        model provided. This function is fed to MULTINEST.

        '''

        # The Likelihood function for MULTINEST
        def likelihood(cube, ndim, nparams):


            params = array([cube[i] for i in range(ndim)])

            logL = 0. # This will be -2. * log(L)
            for det,mod,lh in zip(self.detList,self.models,self.lhs):
                
                
                # Here the model's params are set and the
                # and the matrix is convolved to get
                # model counts
                
                mod.SetParams(params) #set params for the models
                
                modCnts = mod.GetModelCnts() # convolve the matrix and return counts
                
                lh.SetModelCounts(modCnts[det._energySelection]) #pass model counts to lh

                # Here the DataBin objects' source and background
                # counts are sent to the likelihood object

                lh.SetBackGround(det.GetBkgCounts(),det.GetBkgErr()) # The Get functions automatically make channel
                lh.SetCounts(det.GetTotalCounts())   # selection. Selection is performed before!

                #This is log(L) so the joint stat is an addition
    

                logL+=lh.ComputeLikelihood()
            #print "logL = %f"%logL
            
            jointLH = -0.5*(logL)
            #jointLH = logL    
            return jointLH
        
        # Becuase this is inside a class we want to create a
        # likelihood function that does not have an object ref
        # as an argument, so it is created here as a callback

        self.likelihood = likelihood  #likelihood callback
        self.prior = self.models[0].prior  #prior callback



    def _WriteFit(self):
        '''
        Private function that is called after running MULTINEST.
        It saves relevant information from the fits

        '''


        detectors = []
        rsps = []
        loChans = []
        hiChans = []
        chans = []
        
        dof = -self.n_params
        for det in self.detList:

            loChans.append(det.emin)
            hiChans.append(det.emax)
            chans.append(det._energySelection.tolist())
            detectors.append(det.instrument+"_"+det.det)
            rsps.append(det.rsp)
            dof += len(det.GetTotalCounts())


       # Construct the dictionary that will be read by
       # SpecFitView.
        out = {"outfiles":self.outfilesDir,\
               "basename":self.basename,\
               "duration":self.detList[0].duration,\
               "params":self.models[0].parameters,\
               "detectors":detectors,\
               "rsps":rsps,\
               "dataBinExt":self.detList[0].fileLoc,\
               "model":self.models[0].modName,\
               "stat":self.lhs[0].statName,\
               "dof":dof,\
               "chans":chans,\
               "loEne":loChans,\
               "hiEne":hiChans,\
               "tmin":self.detList[0].tmin,\
               "tmax":self.detList[0].tmax\
        }

        f = open(self.outfilesDir+self.savefile,'w')
        
        json.dump(out,f) # Write to a JSON file
        print
        print "Wrote "+self.outfilesDir+self.savefile
        print
        print
        
        f.close()

        self._MakeXspecTemplate()

    def _PreFitInfo(self):

        print "Starting fit of model:"
        print "\t%s"%self.models[0].modName
        print

    def _MakeXspecTemplate(self):


        s="from xspec import *\n"
        s+="Plot.device=\"/xs\"\n"
        s+="Plot.xAxis = \"keV\"\n"
        s+="Plot.xLog = True\n"
        s+="Plot.yLog = True\n\n"
        sNum = 0
        for det in self.detList:
            s+="s%d = Spectrum(\"%s\")\n"%(sNum,det.fileLoc+det.det+".pha")

            s+="s%d.ignore(\"**-%.2f "%(sNum,det.selMins[0])
            if len(det.selMins)>1:
                
                for x,y in zip(det.selMins[1:],det.selMaxs[:-1]):

                    s+=" %.2f-%.2f "%(y,x)

            
            s+="%.2f-**\")\n\n"%(det.selMaxs[-1])
            sNum+=1

        s+="m = Model(\" \")\n"
        s+="\n\nFit.statMethod = \"%s\"\n"%self.lhs[0].statName
        s+="Fit.nIterations = 10000\n\n"
        
        s+="\nFit.perform()\n"

        f = open("xspecTemplate.py","w")
        f.write(s)
        f.close()
