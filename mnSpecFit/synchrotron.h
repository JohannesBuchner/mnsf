#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_errno.h>





double intergrand(double gamma, void *p);


double synchrotron(double energy, double norm, double estar, double index);


double electronDist(double gamma, double norm, double index, double gammaMin, double gammaTH);
  
  
struct synch_params {double energy; double norm; double estar; double index;};

/*         */


double intergrandComplex(double gamma, void *p);


double synchrotronComplex(double energy, double norm, double estar, double gammaMin, double gammaTH, double index);


struct synch_params_complex {double energy; double norm; double estar; double gammaMin; double gammaTH; double index;};

/*         */
double synchrotronPL(double energy, double norm, double estar, double index, double gammaMin);



double intergrandPL(double gamma, void *p);

double electronPL(double gamma, double norm, double index, double gammaMin);


struct synch_params_pl {double energy; double norm; double estar; double index; double gammaMin;};

/*         */

double synchrotronFast(double energy, double norm, double estar, double index, double gammaMin);

double electronDistFast(double gamma, double norm, double index, double gammaMin);

double intergrandFast(double gamma, void *p);

struct synch_params_fast {double energy; double norm; double estar; double index;  double gammaMin;};

/*         */

double synchrotronPL_cutoff(double energy, double norm, double estar, double index, double gammaMin, double gammaMax);


double synchrotron_cutoff(double energy, double norm, double estar, double index, double gammaMax);



/*         */

double SSC(double energy, double normalization, double chi, double delta);

struct ssc_params {double energy; double chi; double delta;};


double Qintergrand(double y, void *p);
double Pintergrand(double y, void *p);

double P(double delta, double z);
double Q(double delta, double z);
double G(double delta, double z);
double H(double delta, double z);



