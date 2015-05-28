from mnSpecFit.Model import Model
from numpy import exp, power, maximum, zeros
from multiFit.priorGen import *








class Band(Model):

   def __init__(self):

      

      
      
      def band(x,logA,Ep,alpha,beta):

         val = zeros(x.flatten().shape[0])

         
         
         A = power(10.,logA)
         idx  = (x < (alpha-beta)*Ep/(2+alpha))
         nidx = ~idx
         

         val[idx]  = A*( power(x[idx]/100., alpha) * exp(-x[idx]*(2+alpha)/Ep) )
         
         val[nidx] = A*power((alpha -beta)*Ep/(100.*(2+alpha)),alpha-beta)*exp(beta-alpha)*power(x[nidx]/100.,beta)
         return val




      self.paramsRanges = [[1.E-6,1.E3,"J"],[1.E2,1.E5,"U"],[-2.,1.,"U"],[-5.5,-1.9,"U"]]#Changed the prior... should change back
                            

      
      def BandPrior(params, ndim, nparams):

         for i in range(ndim):
            params[i] = priorLU[self.paramsRanges[i][-1]](params[i],self.paramsRanges[i][0],self.paramsRanges[i][1])
         

       
      
      self.modName = "Band"
      self.model=band
      self.prior=BandPrior
      self.n_params = 4
      self.parameters = ["logNorm",r"E$_{\rm p}$",r"$\alpha$",r"$\beta$"]

      self._modelDict = {"params":self.parameters,\
                         "model":band\
                      }
      self._composite = False
