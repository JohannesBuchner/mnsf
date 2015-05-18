from mnSpecFit.Model import Model
from mnSpecFit.synchrotron_glue import synchrotronPy
from numpy import power
from multiFit.priorGen import *


class Synchrotron(Model):


    def __init__(self):


        def TotalSynchrotron(ene,norm,estar,index):

            norm = power(10.,norm)
            
            return synchrotronPy(ene, norm, estar,index)


        ###################################################    
        self.paramsRanges = [[1.E-15,1.,"J"],[0.,2.,"U"],[2.,5.,"U"]]#Changed the prior... should change back
                            

      
        def SynchPrior(params, ndim, nparams): 

            for i in range(ndim):
                params[i] = priorLU[self.paramsRanges[i][-1]](params[i],self.paramsRanges[i][0],self.paramsRanges[i][1])
         

       


        self.modName = "Synchrotron"
        self.model=TotalSynchrotron
        self.prior=SynchPrior
        self.n_params = 3
        self.parameters = ["logNorm",r"E$_{\star}$",r"$p$"]

        self._modelDict = {"params":self.parameters,\
                            "model":TotalSynchrotron\
                        }
        self._composite = False


