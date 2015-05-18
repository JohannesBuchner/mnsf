from mnSpecFit.Model import Model
from numpy import power
from multiFit.priorGen import *


class TsviFast(Model):


    def __init__(self):


        def tsvinchrotron(ene,norm,vm,vc,index):

            norm = power(10.,norm)

            if ene <= vc:

                val = power(ene/vc,-.6666)
            elif vc < ene and ene <= vm:

                val = power(ene/vc,-1.5)
            else:

                val = power(vm/vc,-1.5)*power(ene/vm,-(index+2.)*.5)
            val=val*norm
            return val
            
            


        def TsviPrior(params, ndim, nparams):
         

            params[0] = jefferysPrior(params[0],1E-15,1.)
            params[1] = uniformPrior(params[1], 10., 20000.)
            params[2] = uniformPrior(params[2], 10.,20000)
            params[3] = uniformPrior(params[3], 2., 5.)
            
             
            pass

       


        self.modName = "TsviFast"
        self.model=tsvinchrotron
        self.prior=TsviPrior
        self.n_params = 4
        self.parameters = ["logNorm",r"$\nu_m$", r"$\nu_c$",r"$\delta$"]



