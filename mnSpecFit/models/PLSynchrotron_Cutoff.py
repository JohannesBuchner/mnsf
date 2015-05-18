from mnfit.mnSpecFit.Model import Model
from mnfit.mnSpecFit.synchrotron_glue import synchrotronPL_CO_Py
from numpy import power
from mnfit.priorGen import *
class PLSynchrotron_Cutoff(Model):


    def __init__(self):


        def TotalSynchrotron(ene,norm,estar,index,gammaMax):
            #For now, gammaMin is an undeterminable
            # So set it to 900
            norm = power(10.,norm)

            return synchrotronPL_CO_Py(ene, norm, estar, index, 900.,gammaMax)


        def SynchPrior(params, ndim, nparams):
         

            params[0] = jefferysPrior(params[0],1E-15,1) #Check on this
            params[1] = uniformPrior(params[1],0.,3.)
            params[2] = uniformPrior(params[2],2.,12.)#Must be positive!
            params[3] = uniformPrior(params[3],900.,6000.)
            pass

       


        self.modName = "PLSynchrotron_Cutoff"
        self.model=TotalSynchrotron
        self.prior=SynchPrior
        self.n_params = 4
        self.parameters = ["logNorm",r"E$_{crit}$",r"$\delta$",r"$\gamma_{max}$"]



