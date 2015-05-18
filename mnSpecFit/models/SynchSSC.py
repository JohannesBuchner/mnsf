from mnfit.mnSpecFit.Model import Model
from mnfit.mnSpecFit.synchrotron_glue import synchrotronPLPy, SSCpy
from numpy import power
from mnfit.priorGen import *

class SynchSSC(Model):


    def __init__(self):


        def SynSSCmod(ene,norm1,estar,norm2,chi,index):
            #For now, gammaMin is an undeterminable
            # So set it to 900
            norm1 =  power(10.,norm1)
            norm2 =  power(10.,norm2)
            val   =  synchrotronPLPy(ene, norm1, estar, index, 10.)
            val   += SSCpy(ene,norm2,chi,index)
            return val
        def SynchPrior(params, ndim, nparams):
         
            params[0] = jefferysPrior(params[0],1E-20,1.E1)
            params[1] = uniformPrior(params[1], 1.E-2, 1E4)
            params[2] = jefferysPrior(params[2],1E-15,1.E1)
            params[3] = uniformPrior(params[3], 1E-7,1.)
            params[4] = uniformPrior(params[4], 2., 15.)#Must be positive!
             
            pass





        #Component definitions
        def synch(ene,norm,estar,index):
            norm = power(10.,norm)
            return synchrotronPLPy(ene,norm,estar,index,10.)

        def ssc(ene,norm,chi,index):
            norm = power(10.,norm)
            return SSCpy(ene,norm,chi,index)

        synchDict={"params":\
                   ["logNorm1",r"E$_{crit}$",r"$\delta$"],\
                   "model":synch\
        }
        sscDict = {"params":\
              ["logNorm2",r"$\chi$",r"$\delta$"],\
                   "model":ssc\
        }    


        self.componentLU={"Synchrotron":synchDict,\
                          "SSC":sscDict}

        self.modName = "SynchSSC"
        self.model=SynSSCmod
        self.prior=SynchPrior
        self.n_params = 5
        self.parameters = ["logNorm1",r"E$_{crit}$","logNorm2",r"$\chi$",r"$\delta$"]

