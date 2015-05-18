from mnSpecFit.Model import Model
from mnSpecFit.synchrotron_glue import synchrotronPLPy
from numpy import power
from multiFit.priorGen import *

class PLSynchrotron(Model):


    def __init__(self):


        def TotalSynchrotron(ene,norm,estar,index):
            #For now, gammaMin is an undeterminable
            # So set it to 900
            norm = power(10.,norm)

            return synchrotronPLPy(ene, norm, estar, index, 90.)


        def SynchPrior(params, ndim, nparams):
         
            params[0] = jefferysPrior(params[0],1E-15,1.)
            params[1] = uniformPrior(params[1], 0., 3.)
            params[2] = uniformPrior(params[2], 2., 12.)#Must be positive!
             
            pass

       


        self.modName = "PLSynchrotron"
        self.model=TotalSynchrotron
        self.prior=SynchPrior
        self.n_params = 3
        self.parameters = ["logNorm",r"E$_{crit}$",r"$\delta$"]



