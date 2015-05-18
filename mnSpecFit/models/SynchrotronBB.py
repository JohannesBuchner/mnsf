from mnfit.mnSpecFit.Model import Model
from mnfit.mnSpecFit.synchrotron_glue import synchrotronPy
from numpy import power, exp
from mnfit.priorGen import *


class SynchrotronBB(Model):


    def __init__(self):


        def TotalSynchrotronBB(ene,norm,estar,index,bblogN,kT):

            norm = power(10.,norm)




            val = synchrotronPy(ene, norm, estar,index)
            #BB
            val += power(10.,bblogN)*power(ene,2.)*power( exp(ene/float(kT)) -1.,-1.)
            return val
            
        def SynchBBPrior(params, ndim, nparams):
         
            params[0] = jefferysPrior(params[0],1E-15,1.)
            params[1] = uniformPrior(params[1], 0., 3.)
            params[2] = uniformPrior(params[2], 2., 12.)#Must be positive!
            params[3] = jefferysPrior(params[3], 1E-15,1E-0)
            params[4] = uniformPrior(params[4], 5., 500.)#keV
            pass

       


        self.modName = "SynchrotronBB"
        self.model=TotalSynchrotronBB
        self.prior=SynchBBPrior
        self.n_params = 5
        self.parameters = ["logNorm",r"E$_{crit}$",r"$\delta$","logNormBB","kT"]



