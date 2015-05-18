from mnSpecFit.Model import Model
from numpy import exp, power, zeros
from multiFit.priorGen import *








class PL(Model):

   def __init__(self):



      def pl(x,logA,index):
         
         val = power(10.,logA)*power(x/300.,index)
         return val

      

      self.paramsRanges = [[1.E-15,1.E3,"J"],[-2.,2,"U"]]
                            

      
      def PLPrior(params, ndim, nparams):

         for i in range(ndim):
            params[i] = priorLU[self.paramsRanges[i][-1]](params[i],self.paramsRanges[i][0],self.paramsRanges[i][1])
         


      


      
      self.modName = "powerlaw"
      self.model=pl
      self.prior=PLPrior
      self.n_params = 2
      self.parameters = [r"logN$_{\rm PL}$",r"$\delta$"]

      self._modelDict = {"params":self.parameters,\
                         "model":pl\
                      }
      self._composite = False
