cimport cython
cimport numpy as np
from scipy.integrate import quad, quadrature
from multiFit.Likelihood import Likelihood
from RSPconvolve import RSPconvolve
from numpy import array, zeros
from multiFit.priorGen import *
import numpy as np


DTYPE = np.double
ctypedef np.double_t DTYPE_t

class Model:



    def __init__(self):
        

        self.prior = 0
        self.n_params = 0
        self.likelihood = 0

    def __add__(self,other):





        #Instantiate a new Model()
        NewModel = Model()
        NewModel._composite = True #This will be a composite model so flag it

        duplicate = False #So far we don't know if the two models we are adding together are the same
        
        

        
        #Check if the current model is a composite
        if self._composite: #If it is then we can simply add the added models dictionary to this model

            NewModel.componentLU = self.componentLU

            #Need to check of this model already exists in the LU
            newKey = other.modName
            dupes = 0
            for key in NewModel.componentLU.keys():
                if newKey in key:
                    dupes+=1

            if dupes>0:
                newKey+="_%d"%(dupes+1)
                duplicate = True
                dupeParams = map(lambda x: x+"_%d"%(dupes+1),other.parameters)
                newOtherDict = other._modelDict
                newOtherDict["params"] = dupeParams
                NewModel.componentLU[newKey] = newOtherDict
            else:
                NewModel.componentLU[newKey]=other._modelDict
            

        
        else: #If not then we will have to create one

            if self.modName == other.modName:

                
                dupes=1
                newKey=self.modName+"_%d"%(dupes+1)
                dupeParams = map(lambda x: x+"_%d"%(dupes+1),other.parameters)
                
                newOtherDict = other._modelDict
                newOtherDict["params"] = dupeParams
                
                
                NewModel.componentLU={self.modName:self._modelDict,\
                                      other.modName+"_%d"%(dupes+1):other._modelDict\
                                  }
                NewModel.componentLU[newKey] = newOtherDict
                duplicate = True
            
            else:
                NewModel.componentLU={self.modName:self._modelDict,\
                                      other.modName:other._modelDict\
                                  }
                


        
        self.orig_n_params = self.n_params
        NewModel.modName=self.modName+"+"+other.modName
        NewModel.paramsRanges = self.paramsRanges +other.paramsRanges
        NewModel.n_params =  self.n_params + other.n_params

        if duplicate:
            dupeParams = map(lambda x: x+"_%d"%(dupes+1),other.parameters)
            
            NewModel.parameters = self.parameters + dupeParams
        else:
            NewModel.parameters = self.parameters + other.parameters
        def newPrior(params,ndim,nparams):
            for i in range(ndim):
                params[i] = priorLU[NewModel.paramsRanges[i][-1]](params[i],NewModel.paramsRanges[i][0],NewModel.paramsRanges[i][1])

        NewModel.prior = newPrior
        NewModel.origMod = self.model
        NewModel.otherMod = other.model
        def newModel(*args):

#            print self.orig_n_params
#            print NewModel.n_params
#            print args
            ene = args[0]
            thisParam = args[1:self.orig_n_params+1]
            otherParam = args[1+self.orig_n_params:]
#            print thisParam
#            print otherParam
            val = NewModel.origMod(ene,*thisParam)
            val +=NewModel.otherMod(ene,*otherParam)
            return val

        NewModel.model = newModel

        return NewModel

    def SetParams(self, params):
        '''
        Set the parameters of the model

        ____________________________________________
        arguments:
        params: numpy array of paramerters

        '''


        self.params = params

        ####New code for faster convolution
        self._EvalModel()

        

    def integrate(self, lims):
        '''
        Intergrate the model via scipy.intergrate.quad
        over the limits. Note model must braodcasting aware

        ____________________________________________
        arguments:
        lims: [lo/hi]

        '''

        cdef double lowE = lims[0]
        cdef double highE = lims[1]

        cdef double result

        tmp = []

        for i in range(self.n_params):

            tmp.append(self.params[i])
            
        tmp = tuple(tmp)
        
        result = quad(self.model,lowE,highE,args=tmp,full_output=0,limit=100)[0]

        return result


    def _EvalModel(self):

        tmpCounts = zeros(len(self.rsp.photonE))

        #Low res bins

        
        
        lowRes = self.model(self.rsp.lowEne,*self.params) 



        medRes = array(map(lambda x: sum(self.model(x,*self.params))/3.,self.rsp.medEne))

        hiRes =  array(map(lambda x: sum(self.model(x,*self.params))/7.,self.rsp.highEne))


        tmpCounts[self.rsp.lowEval]=lowRes
        tmpCounts[self.rsp.medEval]=medRes
        tmpCounts[self.rsp.highEval]=hiRes

        self.rsp.SetModelVec(tmpCounts)

    def _EvalModelSlow(self):

        tmpCounts = zeros(len(self.rsp.photonE))

        #Low res bins

        lowRes = array( map(lambda e: self.model(e,*self.params),self.rsp.lowEne)) 

        medRes = array(map(lambda x: sum(map(lambda e: self.model(e,*self.params), x    ))/3.,self.rsp.medEne         ))

        hiRes = array(map(lambda x: sum(map(lambda e: self.model(e,*self.params), x    ))/7.,self.rsp.highEne         ))

        tmpCounts[self.rsp.lowEval]=lowRes
        tmpCounts[self.rsp.medEval]=medRes
        tmpCounts[self.rsp.highEval]=hiRes

        tmpCounts = np.asfarray(tmpCounts, dtype="float")

        
        self.rsp.SetModelVec(tmpCounts)
    


    def __repr__(self):



        for x,y in zip(self.parameters,self.paramsRanges):


            print x,y
        return ""

                
    def SetRSP(self,rsp):
        '''
        Set the instrument response matrix for the model.

        ________________________________________________
        arguments:
        rsp: path to a fits rsp file
        
        sets self.rsp

        ##############
        29/6/2014

        Improving for a faster convolution process

        '''

        rsp = RSPconvolve(rsp)

        rsp.SetModel(self)

        self.rsp=rsp
    


        
    def GetModelCnts(self):
        '''
        Convolves the set model with the RSP and creates the 
        model counts

        ____________________________________________
        returns:
        modelCnts: numpy array of Convolved model counts.

        '''


        self.rsp.CreateModelVector()
        modelCnts = self.rsp.GetCounts()
        
        return array(modelCnts)[0]

        

    def SelectComponent(self,comp):
        '''
        Grabs the component parameters from the component
        dictionary that must be created in subclasses

        '''

        comp = self.componentLU[comp]

        return comp
        
        
