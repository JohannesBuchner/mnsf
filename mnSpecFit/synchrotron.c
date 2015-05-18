#include "synchrotron.h"



double intergrand(double gamma, void *p)
{

  double f;
  struct synch_params *params = (struct synch_params *)p;
  double energy = (params->energy);
  double norm = (params->norm);
  double estar = (params->estar);
  double index = (params->index);

  
  double gammaMin = 90.;
  double gammaTH = 30.;


  f = electronDist(gamma,norm,index, gammaMin, gammaTH)*gsl_sf_synchrotron_1(energy/(estar*gamma*gamma));
  

  return f;
}


double electronDist(double gamma, double norm, double index, double gammaMin, double gammaTH)
{

  double ed;
  double ratio = gammaMin/gammaTH;

  double epsilon = pow(ratio,2.+index)*exp(-ratio);
  
    if (gamma<=gammaMin)
      {
	ed = norm * pow(gamma/gammaTH,2.) * exp(-(gamma/gammaTH));
      }
    else
      {
	ed = norm * epsilon * pow(gamma/gammaTH,-index);
      }	
	

  return ed;


}

double synchrotron(double energy, double norm, double estar, double index)
{
  gsl_set_error_handler_off();

  double result, error;

  double epsabs = 0;
  double epsrel = 1e-5;
  double abserr;
  size_t limit = 10000;


  
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);

  struct synch_params p = {energy, norm, estar, index};

  gsl_function F;
  F.function = &intergrand;
  F.params=&p;

  gsl_integration_qagiu(&F, 1., epsabs, epsrel,limit, w, &result, &abserr);
  
  
  gsl_integration_workspace_free(w);
  
  result/=energy;
  return result;


}


/*

  More complex synchrotron model

 */

double synchrotronComplex(double energy, double norm, double estar, double gammaMin, double gammaTH, double index)
{
  gsl_set_error_handler_off();

  double result, error;

  double epsabs = 0;
  double epsrel = 1e-5;
  double abserr;
  size_t limit = 10000;


  
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);

  struct synch_params_complex p = {energy, norm, estar, gammaMin, gammaTH, index};

  gsl_function F;
  F.function = &intergrandComplex;
  F.params=&p;

  gsl_integration_qagiu(&F, 1., epsabs, epsrel,limit, w, &result, &abserr);
  
  
  gsl_integration_workspace_free(w);
  
  result/=energy;
  return result;


}

double intergrandComplex(double gamma, void *p)
{

  double f;
  struct synch_params_complex *params = (struct synch_params_complex *)p;
  double energy = (params->energy);
  double norm = (params->norm);
  double estar = (params->estar);
  double index = (params->index);

  
  double gammaMin = (params->gammaMin);
  double gammaTH = (params->gammaTH);


  f = electronDist(gamma,norm,index, gammaMin, gammaTH)*gsl_sf_synchrotron_1(energy/(estar*gamma*gamma));
  

  return f;
}


/* PL Synch*/



double electronPL(double gamma, double norm, double index, double gammaMin)
{

  double ed;


  ed = norm * (index - 1.) * pow(gammaMin,index-1.)  * pow(gamma/gammaMin,-index);
	
	

  return ed;


}

double intergrandPL(double gamma, void *p)
{

  double f;
  struct synch_params_pl *params = (struct synch_params_pl *)p;
  double energy = (params->energy);
  double norm = (params->norm);
  double estar = (params->estar);
  double index = (params->index);
  double gammaMin = (params->gammaMin);
  


  f = electronPL(gamma,norm,index, gammaMin)*gsl_sf_synchrotron_1(energy/(estar*gamma*gamma));
  

  return f;
}

double synchrotronPL(double energy, double norm, double estar, double index, double gammaMin)
{
  gsl_set_error_handler_off();

  double result, error;

  double epsabs = 0;
  double epsrel = 1e-5;
  double abserr;
  size_t limit = 10000;


  
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);

  struct synch_params_pl p = {energy, norm, estar, index, gammaMin};

  gsl_function F;
  F.function = &intergrandPL;
  F.params=&p;

  gsl_integration_qagiu(&F, gammaMin, epsabs, epsrel,limit, w, &result, &abserr);
  
  
  gsl_integration_workspace_free(w);
  
  result/=energy;
  return result;



}




/**********

	  Fast cooled synchrotron


**************/

double intergrandFast(double gamma, void *p)
{

  double f;
  struct synch_params_fast *params = (struct synch_params_fast *)p;
  double energy = (params->energy);
  double norm = (params->norm);
  double estar = (params->estar);
  double index = (params->index);
  double gammaMin = (params->gammaMin);
  

  f = electronDistFast(gamma,norm,index, gammaMin)*gsl_sf_synchrotron_1(energy/(estar*gamma*gamma));
  

  return f;
}


double electronDistFast(double gamma, double norm, double index, double gammaMin)
{

  double ed;
  

  double epsilon = norm*gammaMin/(gamma*gamma);
  
    if (gamma<=gammaMin)
      {
	ed = epsilon;
      }
    else
      {
	ed = epsilon*pow(gamma/gammaMin,1-index);
      }	
	

  return ed;


}

double synchrotronFast(double energy, double norm, double estar, double index, double gammaMin)
{
  gsl_set_error_handler_off();

  double result, error;

  double epsabs = 0;
  double epsrel = 1e-5;
  double abserr;
  size_t limit = 10000;


  
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);

  struct synch_params_fast p = {energy, norm, estar, index, gammaMin};

  gsl_function F;
  F.function = &intergrandFast;
  F.params=&p;

  gsl_integration_qagiu(&F, 1., epsabs, epsrel,limit, w, &result, &abserr);
  
  
  gsl_integration_workspace_free(w);
  
  result/=energy;
  return result;


}

/**************************************/
/**************************************/
/**************************************/
/**************************************/
/**************************************/

/* CUT OFFS */

double synchrotronPL_cutoff(double energy, double norm, double estar, double index, double gammaMin, double gammaMax)
{
  gsl_set_error_handler_off();

  double result, error;

  double epsabs = 0;
  double epsrel = 1e-5;
  double abserr;
  size_t limit = 10000;


  
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);

  struct synch_params_pl p = {energy, norm, estar, index, gammaMin};

  gsl_function F;
  F.function = &intergrandPL;
  F.params=&p;

  gsl_integration_qags(&F, gammaMin, gammaMax, epsabs, epsrel,limit, w, &result, &abserr);
  
  
  gsl_integration_workspace_free(w);
  
  result/=energy;
  return result;



}



double synchrotron_cutoff(double energy, double norm, double estar, double index, double gammaMax)
{
  gsl_set_error_handler_off();

  double result, error;

  double epsabs = 0;
  double epsrel = 1e-5;
  double abserr;
  size_t limit = 10000;


  
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);

  struct synch_params p = {energy, norm, estar, index};

  gsl_function F;
  F.function = &intergrand;
  F.params=&p;

  gsl_integration_qags(&F, 1., gammaMax ,epsabs, epsrel,limit, w, &result, &abserr);
  
  
  gsl_integration_workspace_free(w);
  
  result/=energy;
  return result;


}




/* Synchrotron Self Compton

 */



double SSC(double energy, double normalization, double chi, double delta)
{

  gsl_set_error_handler_off();

  double result, resultP, resultQ, error;

  double epsabs = 0;
  double epsrel = 1e-5;
  double abserr;
  size_t limit = 10000;


  
  gsl_integration_workspace *w1 = gsl_integration_workspace_alloc(10000);
  gsl_integration_workspace *w2 = gsl_integration_workspace_alloc(10000);

  
  chi = chi*energy;
  
  

  struct ssc_params p = {energy, chi, delta};

  /* P intergral  */

  gsl_function F1;
  F1.function = &Pintergrand;
  F1.params=&p;

  gsl_integration_qags(&F1, 0, chi,epsabs, epsrel,limit, w1, &resultP, &abserr);
  
  
  gsl_integration_workspace_free(w1);

  /* Q intergral  */

  
  gsl_function F2;
  F2.function = &Qintergrand;
  F2.params=&p;

  gsl_integration_qagiu(&F2, chi, epsabs, epsrel,limit, w2, &resultQ, &abserr);

  result = resultP + resultQ;

  result*=normalization*(delta-1.)*(delta-1.)/(delta+1.)/energy;

  gsl_integration_workspace_free(w2);


  return result;


}

double Pintergrand(double y, void *p)
{

  struct ssc_params *params = (struct ssc_params *)p;
  double energy = (params->energy);
  double delta = (params->delta);
  double chi = (params->chi);
  
  
  
  double order = 5./3.;
  double result;

  result = gsl_sf_bessel_Knu(order,y);
  result *= P(delta,y/chi);
  return result;
    

}


double Qintergrand(double y, void *p)
{

  struct ssc_params *params = (struct ssc_params *)p;
  double energy = (params->energy);
  double delta = (params->delta);
  double chi = (params->chi);
  
  
  
  double order = 5./3.;
  double result;

  result = gsl_sf_bessel_Knu(order,y);
  result *= P(delta,chi/y);
  return result;
    

}


double P(double delta, double z)
{
  double result;

  result = pow(z,(delta+1.)/2.);
  result*=((2./(delta+1.)-log(z))*G(delta,1.)+H(delta,1.));
  return result;

}


double Q(double delta, double z)
{
  double result;

  result = -2.*z*(delta+1.)/((delta+3.)*(delta+3.));
  result *= (2.*log(z)- (delta+11.)/(delta+3.));
  result -= 4.*(delta-1.)/((delta+1.)*(delta+1.));
  result -= 2.*log(z)/(delta+1.);
  result += 2.*z*z*(delta+1.)/((delta+5.)*(delta+5.));

  return result;


}


double G(double delta, double z)
{
  double result;

  result =  2.*z/(delta+3.);
  result *= (2.*log(z)+(delta-1.)/(delta+3.));
  result += 2./(delta+1.);
  result -= 4.*z*z/(delta+5.);
  
  return result;

}


double H(double delta, double z)
{
  double result;

  result =  4.*z/((delta+3.)*(delta+3.));
  result *= (2.*log(z)+(delta-5.)/(delta+3.));
  result += 4./((delta+1.)*(delta+1.));
  result -= 8.*z*z/((delta+5.)*(delta+5.));

  return result;

}

