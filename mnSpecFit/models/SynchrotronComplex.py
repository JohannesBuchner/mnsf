from mnSpecFit.Model import Model
from mnSpecFit.synchrotron_glue import synchrotronComplexPy
from numpy import power
from multiFit.priorGen import *


class SynchrotronComplex(Model):


    def __init__(self):


        def TotalSynchrotron(ene,norm,estar,gammaMin,gammaTH,index):

            norm = power(10.,norm)

            return synchrotronComplexPy(ene, norm, estar, gammaMin, gammaTH, index)


        def SynchPrior(params, ndim, nparams):

            params[0] = jefferysPrior(params[0],1E-15,1.)
            params[1] = uniformPrior(params[1], 0., 3.)
            params[2] = uniformPrior(params[2], 1., 1800.)
            params[3] = uniformPrior(params[3],1,900.)
            params[4] = uniformPrior(params[4], 2., 12.)#Must be positive!
             
            pass

       


        self.modName = "SynchrotronComplex"
        self.model=TotalSynchrotron
        self.prior=SynchPrior
        self.n_params = 5
        self.parameters = ["logNorm",r"E$_{crit}$", r"$\gamma_{min}$", r"$\gamma_{th}$", r"$\delta$"]



