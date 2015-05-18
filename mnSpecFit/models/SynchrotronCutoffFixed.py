from mnSpecFit.Model import Model
from mnSpecFit.synchrotron_glue import synchrotron_CO_Py
from multiFit.priorGen import *
from numpy import power
class SynchrotronCutoffFixed(Model):


    def __init__(self):


        def TotalSynchrotron(ene,norm,estar,gammaMax):


            index = 3.5
            norm = power(10.,norm)

            return synchrotron_CO_Py(ene, norm, estar,index,gammaMax)



        self.paramsRanges = [[1.E-15,1.,"J"],[0.,2.,"U"],[90,10000.,"U"]]
                            

      
        def SynchPrior(params, ndim, nparams): 

            for i in range(ndim):
                params[i] = priorLU[self.paramsRanges[i][-1]](params[i],self.paramsRanges[i][0],self.paramsRanges[i][1])
        

       


        self.modName = "SynchrotronCutoffFixed"
        self.model=TotalSynchrotron
        self.prior=SynchPrior
        self.n_params = 3
        self.parameters = ["logNorm",r"E$_{\star}$",r"$\gamma_{\rm max}$"]

        self._modelDict = {"params":self.parameters,\
                            "model":TotalSynchrotron\
                        }
        self._composite = False



