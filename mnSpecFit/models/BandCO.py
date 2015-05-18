from mnSpecFit.Model import Model
from numpy import exp, power, maximum, zeros
from multiFit.priorGen import *








class BandCO(Model):

   def __init__(self):

      

      
      
      def bandco(x,logA,Ep,alpha,beta,eFolding):

         val = zeros(x.flatten().shape[0])

         
         
         A = power(10.,logA)
         idx  = (x < (alpha-beta)*Ep/(2+alpha))
         nidx = ~idx
         

         val[idx]  = A*( power(x[idx]/100., alpha) * exp(-x[idx]*(2+alpha)/Ep) )
         
         val[nidx] = A*power((alpha -beta)*Ep/(100.*(2+alpha)),alpha-beta)*exp(beta-alpha)*power(x[nidx]/100.,beta)
         val = val * exp(-x/eFolding)
         return val




      self.paramsRanges = [[1.E-6,1.E3,"J"],[1.E2,1.E5,"U"],[-2.,1.,"U"],[-10.,-2,"U"],[500.,1E7,"U"]]
                            

      
      def BandCOPrior(params, ndim, nparams):

         for i in range(ndim):
            params[i] = priorLU[self.paramsRanges[i][-1]](params[i],self.paramsRanges[i][0],self.paramsRanges[i][1])
         

       
      
      self.modName = "Cutoffband"
      self.model=bandco
      self.prior=BandCOPrior
      self.n_params = 5
      self.parameters = ["logNorm",r"E$_{\rm p}$",r"$\alpha$",r"$\beta$","eFolding"]

      self._modelDict = {"params":self.parameters,\
                         "model":bandco\
                      }
      self._composite = False
