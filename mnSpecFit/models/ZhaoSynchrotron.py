from mnfit.mnSpecFit.Model import Model
from numpy import exp, power
from mnfit.priorGen import *


import subprocess
import numpy as np
from numpy import array, zeros
from decimal import Decimal
from subprocess import Popen, PIPE





class ZhaoSynchrotron(Model):

   def __init__(self):



      

      def sync_flux(E,gamma_m,p,y0,kai,alpha_B,const):
         nu=array(map(float,E*2.41789894010E17))
         E =array(map(float,E))
         gamma=300.
         gamma_max=1.e7
         q0=1.e15
         Bm=1.e4
         const = power(10.,const)
         gamma_m = power(10.,gamma_m)
         #gamma_m=P[0]
         #p=P[1]
         # q0=1.e15
         #y0=P[2]
         #kai=P[3]
         #alpha_B=P[4]
         # Bm=1.e4
         # CONST=P[5]
         #const=P[5]
         if not ('redshift_zxh' in globals()):
            redshift_zxh=1.0
         z=redshift_zxh
         
         gamma_str=str(Decimal(gamma))
         gamma_max_str=str(Decimal(gamma_max))
         gamma_m_str=str(Decimal(gamma_m))
         p_str=str(Decimal(p))
         q0_str=str(Decimal(q0))
         y0_str=str(Decimal(y0))
         z_str=str(Decimal(z))
         kai_str=str(Decimal(kai))
         alpha_B_str=str(Decimal(alpha_B))
         Bm_str=str(Decimal(Bm))
         nn_str=str(Decimal(len(nu)))
         
         nu_min_str=str(Decimal(min(nu)))
         nu_max_str=str(Decimal(max(nu)))

         nu_big_str=""
         for i in range(len(nu)):
            nu_now=nu[i]
            nu_str_tmp=str(Decimal(nu_now))
            nu_big_str=nu_big_str+' '+nu_str_tmp


         #cmd_zxh='~/CloudStation/ZXH_Fortran_Code/zxhflux_quick2'
         cmd_zxh='/Users/jburgess/Research/mnfit/mnSpecFit/models/syn_flux'
         # cmd_zxh='~/CloudStation/ZXH_Fortran_Code/zxhflux_bkpl'
         para_str=' '+gamma_str+" "+gamma_max_str+" "+gamma_m_str+' '+\
            p_str+' '+q0_str+' '+y0_str+' ' +z_str+' '+kai_str+' '+\
            alpha_B_str+' '+Bm_str+' '+nn_str+' '+nu_min_str+' '+nu_max_str+' '+nu_big_str
         whole_cmd=cmd_zxh+' '+para_str
         # print(whole_cmd)
         
         p = Popen(whole_cmd,stdout=PIPE,shell=True,stderr=PIPE)
         (out,err) = p.communicate()
         
         tmp_str_list=out.split()

         result_tmp=np.ndarray(len(out.split()))*0.
         for i in range(len(result_tmp)):
            result_tmp[i]=np.float(tmp_str_list[i])

         result_tmp=result_tmp
         ph_flux=const*np.array(result_tmp)/E
         return ph_flux


      def ZhaoPrior(params, ndim, nparams):
         
         params[0] = jefferysPrior(params[0],1E2,1E5)
         params[1] = uniformPrior(params[1], 0.5, 4.52)
         params[2] = uniformPrior(params[2], 0.4, 1.5)
         params[3] = uniformPrior(params[3], 0.5, 1.5)
         params[4] = uniformPrior(params[4], 0.02, 3.0)
         params[5] = jefferysPrior(params[5], 1E-10, 1E20)
         

         pass

       

      self.modName = "ZhaoSynchrotron"
      self.model=sync_flux
      self.prior=ZhaoPrior
      self.n_params = 6
      self.parameters = [r"log$\gamma$","p","y0","kai",r"$\alpha_B$","logConst"]




#   def _EvalModel(self):
#
#      tmpCounts = zeros(len(self.rsp.photonE))
#
      #Low res bins
#
#      lowRes = array(self.model(self.rsp.lowEne,*self.params)) 
#      
#      medRes = array(map(lambda x: sum(self.model(x,*self.params))/3.,self.rsp.medEne))

#      hiRes = array(map(lambda x: sum(self.model(x,*self.params))/7.,self.rsp.highEne))

#      tmpCounts[self.rsp.lowEval]=lowRes
#      tmpCounts[self.rsp.medEval]=medRes
#      tmpCounts[self.rsp.highEval]=hiRes
      
#      self.rsp.SetModelVec(tmpCounts)
    
