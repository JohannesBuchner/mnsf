from mnSpecFit.Model import Model
from numpy import power
from multiFit.priorGen import *


class TsviSlow(Model):


    def __init__(self):


        def tsvinchrotron(ene,norm,vm,vc,index):

            norm = power(10.,norm)

            if ene <= vm:

                val = power(ene/vm,-.6666)
            elif vm < ene and ene <= vc:

                val = power(ene/vm,-(index+1)*.5)
            else:

                val = power(vc/vm,-(index+1)*.5)*power(ene/vc,-(index+2.)*.5)
            val=val*norm
            return val
            
            


        def TsviPrior(params, ndim, nparams):
         

            params[0] = jefferysPrior(params[0],1E-15,1.)
            params[1] = uniformPrior(params[1], 10., 20000.)
            params[2] = uniformPrior(params[2], 10.,20000)
            params[3] = uniformPrior(params[3], 2., 5.)
            
             
            pass

       


        self.modName = "TsviSlow"
        self.model=tsvinchrotron
        self.prior=TsviPrior
        self.n_params = 4
        self.parameters = ["logNorm",r"$\nu_m$", r"$\nu_c$",r"$\delta$"]



