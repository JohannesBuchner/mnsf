from mnSpecFit.Model import Model
from numpy import log, exp, log10, power, logical_or, zeros
from multiFit.priorGen import *

import numexpr


class SBPL(Model):



    def __init__(self):

        def sbpl(ene, logN, indx1, breakE, breakScale, indx2):

            pivot =300. #keV
            B = (indx1 + indx2)/2.0
            M = (indx2 - indx1)/2.0

            arg_piv = log10(pivot/breakE)/breakScale

            if arg_piv < -6.0:
                pcosh_piv = M * breakScale * (-arg_piv-log(2.0))
            elif arg_piv > 4.0:    

                pcosh_piv = M * breakScale * (arg_piv - log(2.0))
            else:
                pcosh_piv = M * breakScale * (log( (exp(arg_piv) + exp(-arg_piv))/2.0 ))


            arg = numexpr.evaluate('log10(ene/breakE)/breakScale')
            idx1 =  arg < -6.0
            idx2 =  arg >  4.0
            idx3 =  ~logical_or(idx1,idx2)

            pcosh = zeros(ene.flatten().shape[0])

            pcosh[idx1] = M * breakScale *(-arg[idx1]-log(2.0))
            pcosh[idx2] = M * breakScale *(arg[idx2] - log(2.0))
            pcosh[idx3] = M * breakScale *(log( (exp(arg[idx3]) + exp(-arg[idx3]))/2.0 ))

            return numexpr.evaluate( '10.**logN * (ene/pivot)**B * 10.**(pcosh-pcosh_piv) ')





        self.paramsRanges = [[1.E-15,1.E3,"J"],[-5.,1.,"U"],[1E1,1E7,"U"],[.001,9.,"U"],[-5.,-1.5,"U"]]
                            

      
        def SBPLPrior(params, ndim, nparams):

            for i in range(ndim):
                params[i] = priorLU[self.paramsRanges[i][-1]](params[i],self.paramsRanges[i][0],self.paramsRanges[i][1])
         


        self.modName = "SBPL"

        self.model=sbpl
        self.prior=SBPLPrior
        self.n_params = 5
        self.parameters = [r"logN$_{\rm sbpl}$",r"indx1",r"breakE",r"breakScale","indx2"]

        self._modelDict = {"params":self.parameters,\
                           "model":sbpl\
        }
        self._composite = False

