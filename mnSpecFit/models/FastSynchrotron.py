from mnfit.mnSpecFit.Model import Model
from mnfit.mnSpecFit.synchrotron_glue import synchrotronFastPy
from numpy import power
from mnfit.priorGen import *


class FastSynchrotron(Model):


    def __init__(self):


        def TotalSynchrotron(ene,norm,estar,index):

            norm = power(10.,norm)

            gammaMin = 90. #Assuming that I cannot find this
            
            return synchrotronFastPy(ene, norm, estar,index, gammaMin)


        def SynchPrior(params, ndim, nparams):
         

            params[0] = jefferysPrior(params[0],1E-15,1.)
            params[1] = uniformPrior(params[1], 0., 3.)
            params[2] = uniformPrior(params[2], 2., 12.)#Must be positive!
             
            pass

       


        self.modName = "FastSynchrotron"
        self.model=TotalSynchrotron
        self.prior=SynchPrior
        self.n_params = 3
        self.parameters = ["logNorm",r"E$_{crit}$",r"$\delta$"]



