from mnfit.mnSpecFit.Model import Model
from numpy import exp, power, zeros
from subPhotoInterp.subPhotoInterp import subPhotoInterp
from mnfit.priorGen import *








class SubPhoto(Model):

   def __init__(self):



      sub = subPhotoInterp("/Users/jburgess/Research/mnfit/mnSpecFit/models/TableModel_grid500_res8_Nr1.fits",silent=True)
      
      

      self.paramsRanges = [[1.01,34.999,"U"],[51.,500.,"U"],[.11,299.,"U"],[.11,.499,"U"]]
                            

      
      def Prior(params, ndim, nparams):

         for i in range(ndim):
            params[i] = priorLU[self.paramsRanges[i][-1]](params[i],self.paramsRanges[i][0],self.paramsRanges[i][1])
         


      


      
      self.modName = "SubPhoto"
      self.model=sub
      self.prior=Prior
      self.n_params = 4
      self.parameters = [r"$\tau$",r"$\Gamma$",r"$L_{\rm GRB}$",r"$\epsilon_{\rm d}$"]

      self._modelDict = {"params":self.parameters,\
                         "model":sub\
                      }
      self._composite = False
