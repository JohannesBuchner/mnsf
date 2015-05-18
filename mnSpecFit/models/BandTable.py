from mnfit.mnSpecFit.Model import Model
from numpy import exp, power, zeros
from mnfit.priorGen import *
import scipy.interpolate
import numpy
import copy





class BandTable(Model):

   def __init__(self):


    npz = numpy.load("/Users/jburgess/Research/mnfit/mnSpecFit/models/band.npz")
    self._evalEnergies = npz["energy"]
    self._tableParams  = npz["params"].tolist()
    self._tableFluxes  = npz["fluxes"]
    tmp = copy.deepcopy(self._tableParams)
        
    tmp.append(self._evalEnergies.tolist())
    tmp = map(lambda x: numpy.array(x,dtype=numpy.float32),tmp)
    interpGrid = tuple(tmp)
#    self._tableFluxes.dtype = numpy.float32
    zero = numpy.float32(0.)
        
    self._interpFunc = scipy.interpolate.RegularGridInterpolator(interpGrid,self._tableFluxes,method="linear",fill_value=zero,bounds_error=False)

      
    self.paramsRanges = [[1.E-6,1.E3,"J"],[1.E2,1.E5,"U"],[-2.,1.,"U"],[-5.5,-2,"U"]]#Changed the prior... should change back

    def func(ene,A,Ep,alpha,beta):

        return self._interpFunc((10.**A,Ep,alpha,beta,ene))
      
    def Prior(params, ndim, nparams):

        for i in range(ndim):
            params[i] = priorLU[self.paramsRanges[i][-1]](params[i],self.paramsRanges[i][0],self.paramsRanges[i][1])
         


      


      
    self.modName = "BandTable"
    self.model=func
    self.prior=Prior
    self.n_params = 4
    self.parameters = ["logNorm",r"E$_{\rm p}$",r"$\alpha$",r"$\beta$"]
      

    self._modelDict = {"params":self.parameters,\
                        "model":func\
                      }
    self._composite = False
