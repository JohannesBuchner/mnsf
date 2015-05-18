from mnfit.mnSpecFit.Model import Model
from mnfit.mnSpecFit.synchrotron_glue import synchrotronPLPy, SSCpy
from numpy import power
from mnfit.priorGen import *

class SynchSSC_BB(Model):


    def __init__(self):



        def bb(x,logA,kT):
         
            val = power(10.,logA)*power(x,2.)*power( exp(x/float(kT)) -1., -1.)
            return val


        def SynSSCmod(ene,norm1,estar,norm2,chi,index,norm3,kT):
            #For now, gammaMin is an undeterminable
            # So set it to 900
            norm1 =  power(10.,norm1)
            norm2 =  power(10.,norm2)
            val   =  synchrotronPLPy(ene, norm1, estar, index, 10.)
            val   += SSCpy(ene,norm2,chi,index)
            val   += bb(ene, norm3, kT)
            return val

        

        def SynchPrior(params, ndim, nparams):
         
            params[0] = jefferysPrior(params[0],1E-20,1.E1)
            params[1] = uniformPrior(params[1], 1.E-2, 1E4)
            params[2] = jefferysPrior(params[2],1E-15,1.E1)
            params[3] = uniformPrior(params[3], 1E-7,1.)
            params[4] = uniformPrior(params[4], 2., 15.)#Must be positive!
            params[5] = jefferysPrior(params[5], 1E-15,1E-0)
            params[6] = uniformPrior(params[6], 5., 500.)#keV 
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

        bbDict = {"params":\
                [r"logN$_{\rm BB}$","kT"],\
                "model":bb\
            }

        self.componentLU={"Synchrotron":synchDict,\
                          "SSC":sscDict,\
                          "Blackbody":bb\
        }

        self.modName = "SynchSSC_BB"
        self.model=SynSSCmod
        self.prior=SynchPrior
        self.n_params = 7
        self.parameters = ["logNorm1",r"E$_{crit}$","logNorm2",r"$\chi$",r"$\delta$",r"logN$_{\rm BB}$","kT"]

