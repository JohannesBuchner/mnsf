from multiFit.FitView import FitView
from astropy.table import Table
from DataBin import DataBin
from models.models import models

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import Grid


from numpy import array, cumsum, linspace, sqrt, logspace, log10, mean
from scipy.stats import ks_2samp

import json


class SpecFitView(FitView):


    def _LoadData(self,data):



        f = open(data)

        fit = json.load(f)

        self.modName = fit["model"]

        self._composite = False
        #Check to see if this is a composite model

        test = self.modName.split("+")
        if len(test)>1:
            self._composite = True
            compositeModels = test
            
            for i in range(len(compositeModels)):

                tmp = compositeModels[i].split('_')[0] #Hack off the duplicate model number tag
                compositeModels[i] = tmp
            
        self.parameters = fit["params"]
        self.n_params = len(self.parameters)

        self.detectors = array(fit["detectors"])
        self.dataBinExt = fit["dataBinExt"]
        self.duration = fit["duration"]
        self.sourceCounts = []
        self.rsps = fit["rsps"]
        self.basename = fit["basename"]
        self.meanChan = []
        self.chanWidths = []
        self.tmin = fit["tmin"]
        self.tmax = fit["tmax"]


        self.meanTime = mean([self.tmin,self.tmax])

                
        if self._composite:

            model = (models[compositeModels[0]])()
            for mod in compositeModels[1:]:

                tmp = (models[mod])()
                model = model + tmp

            self._componentLU = model.componentLU
            self._componentModel = model
        else:
            model = (models[fit["model"]])()
        self.model = model.model


        self.stat = fit["stat"]
        self.dof = fit["dof"]

        self.cntMods = []
        self.activeLos = []
        self.activeHis = []
        #load counts and model counts
        for det ,lo, hi  in zip(self.detectors,fit["loEne"],fit["hiEne"]):

            db = DataBin(self.dataBinExt+det+".fits")

            if self._composite:

                mod = (models[compositeModels[0]])()
                for m in compositeModels[1:]:

                    tmp = (models[m])()
                    mod = mod + tmp
            else:
                mod = (models[fit["model"]])()

            #mod = (models[fit["model"]])()
            mod.SetRSP(db.rsp)

            chanWidth = db.chanMax - db.chanMin
            self.chanWidths.append(chanWidth)

            self.cntMods.append(mod)

            self.meanChan.append(db.meanChan)
            
            self.sourceCounts.append(db.source)

            db.SetLoChan(lo)
            db.SetHiChan(hi)
            loIndex = db.activeLoChan
            hiIndex = db.activeHiChan
            self.activeLos.append(loIndex)
            self.activeHis.append(hiIndex+1)

        #Move all of these to arrays
        self.sourceCounts = array(self.sourceCounts)
        self.meanChan = array(self.meanChan)
        self.cntMods = array(self.cntMods)


        
        
        
        
        
        self.xlabel = "Energy [keV]"

        minE = min(fit["loEne"])
        maxE = max(fit["hiEne"])
        
        self.dataRange = logspace(log10(minE),log10(maxE),700)


    def _CustomInfo(self):

        print

        print "Model:\n\t%s"%self.modName
        print "T:%.2f - %.2f"%(self.tmin,self.tmax)
        print "Detectors:"
        print self.detectors
        print "\nBest Fit Parameters (1-sigma err):"

        marg = self.anal.get_stats()["marginals"]

        for params,val,err in zip(self.parameters,self.bestFit,marg):

            err = err["1sigma"]
            
            print "\t%s:\t%.2f\t+%.2f -%.2f"%(params,val,err[1]-val,val-err[0])

        print
        print "%s per d.o.f.:\n\t %.2f/%d"%(self.stat,-2.*self.loglike,self.dof)


    def PlotvFv(self,fignum = 100):
        '''
        Plots the best fit and the surrounding likelihood space
        but in vFv space instead of the standard photon space.

        self.dataRange must be set!

        '''



        fig = plt.figure(fignum)
        ax = fig.add_subplot(111)

        bfColor="#e41a1c"
        posColor="#377eb8"
        
        yData = []


        for params in self.anal.get_equal_weighted_posterior()[::50,:-1]:

            tmp = []

            

            tmp = self.dataRange**2*self.model(self.dataRange,*params)
                
                
            yData.append(tmp)



        
       


        for y in yData:

            ax.plot(self.dataRange,y,posColor,alpha=.2,zorder=-30,linewidth=.9) ## modify later

        

        
        bfModel = self.dataRange**2*self.model(self.dataRange,*self.bestFit)
        
            
        ax.plot(self.dataRange,bfModel,bfColor,linewidth=1.6) #modify later
        ax.set_xscale('log')
        ax.set_yscale('log')

        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(r"$\nu F_{\nu}$ [keV$^2$ s$^{-s}$ cm$^{-2}$ keV$^{-1}$]")
        ax.grid(False)
        ax.set_xlim(min(self.dataRange),max(self.dataRange))
        ax.set_ylim(bottom=10)
        return ax
        
            
        
            


    def PlotvFvComponents(self,fignum=400):

        '''
        Plots the best fit and the surrounding likelihood space
        but in vFv space instead of the standard photon space.

        self.dataRange must be set!

        '''

        if not self._composite:

            print "This is not a composite model!"
            return

        fig = plt.figure(fignum)
        ax = fig.add_subplot(111)
        #model = models[self.modName]() #Remember that models must be instantiated!!
        
        bfColors = ["#FF0040","#2E2EFE","#01DF3A"]
        #contourColors = ["#AC58FA","#FF00FF","#088A29"]
        #First get the components
        
        components = self._componentLU.keys()
        
            
        colorIndex = 0
        for comp in components:
            
            thisComp= self._componentModel.SelectComponent(comp)
            
            
            #First the best fit params
            tt = self.GetParamIndex(thisComp["params"])
            bfParams = self.bestFit[tt]

            #Plot the best fit component
        
            yData = self.dataRange**2*thisComp["model"](self.dataRange, *bfParams) #Computes vFv


            ax.loglog(self.dataRange,yData,color="k")
            

            #Now Plot the contours

            yData = []


            for params in self.anal.get_equal_weighted_posterior()[::100,:-1]:

                
                params = params[tt]
                
                tmp = self.dataRange**2*thisComp["model"](self.dataRange, *params) #Computes vFv
                yData.append(tmp)
            

            for y in yData:

                ax.loglog(self.dataRange,y,color=bfColors[colorIndex],alpha=.2) ## modify later

            colorIndex+=1    
        ax.set_xlim(min(self.dataRange),max(self.dataRange))
        ax.set_ylim(bottom=10)
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(r"$\nu F_{\nu}$ [keV$^2$ s$^{-s}$ cm$^{-2}$ keV$^{-1}$]")
        return ax
    

    def PlotCounts(self,fignum=140):
        '''
        Plots the data counts along with the deconvolved models
        for each detector used in the fit. 

        '''

        
        
        fig = plt.figure(fignum,(9,5))

        ax = fig.add_subplot(111)
        
        colorLU=["#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","#999999"]
        #colorLU=["#a6cee3", "#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6"]
        #colorLU = ["#FF0000","#01DF01","#DA81F5","#0101DF","#F781D8","#58D3F7","#FFFF00"]

        for c,chan, color,cw in zip(self.sourceCounts,self.meanChan,colorLU,self.chanWidths):

            ax.errorbar(chan,c/cw,yerr = sqrt(c/cw),fmt="o", color=color, elinewidth=.6,capsize=.4,alpha=.75,markersize=4)
            #ax.errorbar(chan,c/cw,fmt="+", color=color)

        
        #ax.legend(self.detectors,loc="lower left")


        # Here the model's params are set and the
        # and the matrix is convolved to get
        # model counts

        for mod,chan,cw,lo,hi in zip(self.cntMods,self.meanChan,self.chanWidths,self.activeLos,self.activeHis):

            yData = []


            for params in self.anal.get_equal_weighted_posterior()[::50,:-1]:

                mod.SetParams(params)
                    
                yData.append(mod.GetModelCnts()[lo:hi]/cw[lo:hi])

            for y  in yData:


                ax.loglog(chan[lo:hi],y,"k",alpha=.1,zorder=-30)



        
        ax.set_xlim(left=min(self.dataRange), right=max(self.dataRange))
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(r"cnts s$^{-1}$ keV$^{-1}$")
        ax.set_yscale("log", nonposy= "clip")
        
        
        return ax

        

    def QQ(self,detector,numberOfEnergyPoints=3):
        '''
        Plot data counts against model counts for each detector
        '''

        
        fig  = plt.figure(180,(8,6))
        ax = fig.add_subplot(111)

        #grid = Grid(fig,111,nrows_ncols=(3,4),axes_pad=0.1)
        if detector not in self.detectors:
            print "Invalid Detector"
            return

        i = self.detectors == detector

        
        


        counts = self.sourceCounts[i][0]


        mod = self.cntMods[i][0]

        mod.SetParams(self.bestFit)

        modelCounts = mod.GetModelCnts()

        cumCounts = cumsum(counts)
        cumModel = cumsum(modelCounts)

        # Plot the center line
        dataMax = max(max(cumModel),max(cumCounts))
        line = linspace(0,dataMax,1000)

        ks, pvalue = ks_2samp(cumModel,cumCounts)

        ax.plot(line,line,"--",color="grey")

        # Plot QQ
        ax.plot(cumModel,cumCounts,"r-", linewidth = 1.2)

        ax.set_xlabel("Integrated Model")
        ax.set_ylabel("Integrated Counts")

        #Plot energies

        minE = min(self.meanChan[i][0])
        maxE = max(self.meanChan[i][0])

        enePoints = linspace(minE,maxE,numberOfEnergyPoints)

        for i,ep in zip(linspace(0.,1.,numberOfEnergyPoints),enePoints):

            if i<.5:

                ax.text(i,i+.05,"%.1f keV"%ep,color='k',transform=ax.transAxes,horizontalalignment="right",fontsize=7)
            else:
                 ax.text(i,i-.05,"%.1f keV"%ep,color='k',transform=ax.transAxes,fontsize=7)


        ax.text(0.55,0.5,"Model Excess",fontsize=8,color="grey",verticalalignment="center",horizontalalignment='center',transform=ax.transAxes,rotation=45)
        ax.text(0.5,0.55,"Data Excess",fontsize=8 ,color="grey",verticalalignment="center",horizontalalignment='center',transform=ax.transAxes,rotation=45)


        ax.text(.9,.1,"K-S=%.3f"%ks,transform=ax.transAxes,horizontalalignment="right",color="grey")

        ax.set_ylim(0,dataMax)
        ax.set_xlim(0,dataMax)

        ax.grid(False)

        return ax
        
        

             




        
