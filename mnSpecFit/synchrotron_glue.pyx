cimport cython
import numpy as np
cimport numpy as np

cdef extern from "/usr/local/include/synchrotron.h":
    double synchrotron(double,double,double,double)
    double synchrotronComplex(double,double,double,double,double,double)
    double synchrotronPL(double, double, double, double, double)
    double synchrotronFast(double, double, double, double, double)
    double synchrotronPL_cutoff(double, double, double, double, double, double)
    double synchrotron_cutoff(double, double, double, double, double)
    double SSC(double, double, double , double )



DTYPE = np.double
ctypedef np.double_t DTYPE_t

#Synchrotron from shock accelerated electrons
@cython.boundscheck(False)
@cython.wraparound(False)
def synchrotronPy(energy, norm, estar, index):
    val = np.zeros(len(energy)) 
    for i in range(len(energy)):
        val[i]  = synchrotron(energy[i],norm,estar,index)

    return val



# cpdef np.ndarray[np.double_t,ndim=1] synchrotronPy(np.ndarray[DTYPE_t, ndim=1] energy, double norm, double estar, double index):
#     cdef np.ndarray[np.double_t,ndim=1] val = np.zeros(len(energy)) 
#     for i in range(len(energy)):
#         val[i]  = synchrotron(energy[i],norm,estar,index)

#     return val

#return synchrotron(energy, norm, estar, index)





#Synchrotron from shock accelerated electrons
@cython.boundscheck(False)
@cython.wraparound(False)
def synchrotron_CO_Py(energy, norm, estar, index, gammaMax):

    val = np.zeros(len(energy)) 
    
    for i in range(len(energy)):

        val[i]  = synchrotron_cutoff(energy[i], norm, estar, index, gammaMax)

    return val


# cpdef np.ndarray[np.double_t,ndim=1] synchrotron_CO_Py(np.ndarray[DTYPE_t, ndim=1] energy, double norm, double estar, double index, double gammaMax):
#     cdef np.ndarray[np.double_t,ndim=1] val = np.zeros(len(energy)) 
#     cdef double ene
#     for i in range(len(energy)):

#         ene = energy[i]
        
#         val[i]  = synchrotron_cutoff(ene, norm, estar, index, gammaMax)

#     return val



   



#Synchrotron from shock accelerated electrons but with lots of parameters
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef double synchrotronComplexPy(double energy, double norm, double estar, double gammaMin, double gammaTH, double index):

    return synchrotronComplex(energy, norm, estar, gammaMin, gammaTH, index)


#Synchrotron from a power law of electrons
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef np.ndarray[np.double_t,ndim=1] synchrotronPLPy(np.ndarray[DTYPE_t, ndim=1] energy, double norm, double estar, double index, double gammaMin):

    cdef np.ndarray[np.double_t,ndim=1] val = np.zeros(len(energy)) 
    for i in range(len(energy)):
        val[i]  = synchrotronPL(energy[i], norm, estar, index, gammaMin)

    return val


    

#Synchrotron from a power-law with a high-energy cutoff due to a max gamma
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef double synchrotronPL_CO_Py(double energy, double norm, double estar, double index, double gammaMin, double gammaMax):

    return synchrotronPL_cutoff(energy, norm, estar, index, gammaMin, gammaMax)



#Simple Fast cooled synchrotron
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef np.ndarray[np.double_t,ndim=1] synchrotronFastPy(np.ndarray[DTYPE_t, ndim=1] energy, double norm, double estar, double index, double gammaMin):

        cdef np.ndarray[np.double_t,ndim=1] val = np.zeros(len(energy)) 
        for i in range(len(energy)):
            val[i] = synchrotronFast(energy[i], norm, estar, index, gammaMin)

        return val
#Synchrotron Self Compton
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef SSCpy(double energy, double norm, double chi, double delta):

    return SSC(energy, norm, chi, delta)



