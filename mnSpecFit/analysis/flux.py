from mnfit.mnSpecFit.SpecFitView import SpecFitView as sfv
from scipy.integrate import quad
from mnfit.mnSpecFit.models.models import models
keV2erg = 1.602E-9

class flux(object):


    def __init__(self, databin):

        self.fit = sfv(databin,silent=True)
        self.model = models[self.fit.modName]() #This is the class not the function!
        self.compFlux = {}
    def CalculateTotalFlux(self,emin=10.,emax=40000.):

        
        
        bfParams = self.fit.bestFit
        
        # Construct a callback that is only a
        # function of energy for the bestFit params
        
        def tmpModel(ene):
            ene = keV2erg*ene
            return ene*self.model.model(ene,*bfParams)

        bfFlux = self._fluxIntergral(tmpModel,emin,emax)


        # Use the propagate function to propagate the fluc
        # first ew Construct a fucntion to send to it

        def propFunc(params):

            def tmpModel(ene):
                ene = keV2erg*ene
                return ene*self.model.model(ene,*params)

            return self._fluxIntergral(tmpModel,emin,emax)

        fluxDist = self.fit.Propagate(propFunc,self.model.parameters,direct=False)
        self.flux =  {"flux":bfFlux,"distribution":fluxDist}

    def CalculateComponentFlux(self,comp,emin=10.,emax=40000.):
        compName = comp
        #Get the component
        comp = self.model.SelectComponent(comp)

        #First the best fit flux
        tt = self.fit.GetParamIndex(comp["params"])
        bfParams = self.fit.bestFit[tt]

        # Construct a callback that is only a
        # function of energy for the bestFit params
        
        def tmpModel(ene):
            
            return ene*comp["model"](ene,*bfParams)

        bfFlux = keV2erg * self._fluxIntergral(tmpModel,emin,emax)


        # Use the propagate function to propagate the fluc
        # first ew Construct a fucntion to send to it

        def propFunc(params):

            def tmpModel(ene):
                
                return ene*comp["model"](ene,*params)

            return keV2erg * self._fluxIntergral(tmpModel,emin,emax)

        fluxDist = self.fit.Propagate(propFunc,comp["params"],direct=False)
        fluxDict = {"flux":bfFlux,"distribution":fluxDist}
        self.compFlux[compName] = fluxDict
        
        
    def _fluxIntergral(self,model,emin, emax):


        
        f = quad(model,emin,emax)[0]

        return f
