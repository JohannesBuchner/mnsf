from mnfit.mnSpecFit.Model import Model
from mnfit.mnSpecFit.synchrotron_glue import synchrotronFastPy
from numpy import power, exp
from mnfit.priorGen import *


class FastSynchrotronBB(Model):


    def __init__(self):



        def bb(x,logA,kT):
         
            val = power(10.,logA)*power(x,2.)*power( exp(x/float(kT)) -1., -1.)
            return val

        
        def TotalSynchrotron(ene,norm,estar,index):

            norm = power(10.,norm)

            gammaMin = 900. #Assuming that I cannot find this
            
            return synchrotronFastPy(ene, norm, estar,index, gammaMin)


        def fastBB(x,logA1,estar,index, logA2, kT):

            val =  TotalSynchrotron(x,logA1,estar,index)
            val += bb(x,logA2,kT)

            return val
        

        def FastBBPrior(params, ndim, nparams):
         

            params[0] = jefferysPrior(params[0],1E-15,1.)
            params[1] = uniformPrior(params[1], 0., 3.)
            params[2] = uniformPrior(params[2], 2., 12.)#Must be positive!
            params[3] = jefferysPrior(params[3], 1E-15,1E-0)
            params[4] = uniformPrior(params[4], 5., 500.)#keV
            pass

       
        
        bbDict = {"params":\
                [r"logN$_{\rm BB}$","kT"],\
                "model":bb\
            }

        fastDict = {"params":\
              ["logNorm",r"E$_{crit}$",r"$\delta$"],\
                "model":TotalSynchrotron\
               }

        self.componentLU={"Blackbody":bbDict,\
                        "FastSynchrotron":fastDict\
                      }


        self.modName = "FastSynchrotronBB"
        self.model=fastBB
        self.prior=FastBBPrior
        self.n_params = 5
        self.parameters = ["logNorm",r"E$_{crit}$",r"$\delta$",r"logN$_{\rm BB}$","kT"]



