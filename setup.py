#!/usr/bin/env python
 
from distutils.core import setup
from distutils.extension import Extension
from distutils.command.install_headers import install_headers
from distutils.command.build_clib import build_clib
from Cython.Distutils import build_ext
from Cython.Build import cythonize



import os

import numpy



#This is a big kludge
os.system("gcc -shared -o libsynchrotron.so -fPIC mnSpecFit/synchrotron.c -I/usr/local/include  -L/usr/local/lib -lgsl -lgslcblas")
os.system("mv libsynchrotron.so /usr/local/lib/")
os.system("cp mnSpecFit/synchrotron.h /usr/local/include/")

#libsynchrotron = ('synchrotron',{'sources':['mnSpecFit/synchrotron.c'],'build_dir':['/usr/local/lib']})


ext_modules = [Extension("mnSpecFit/Model",["mnSpecFit/Model.pyx"],include_dirs = [numpy.get_include()]),
               Extension("mnSpecFit/synchrotron_glue",["mnSpecFit/synchrotron_glue.pyx"],
            library_dirs=['/usr/local/lib'],
            libraries=["synchrotron"],include_dirs = [numpy.get_include()])]




setup(
    
    name="mnSpecFit",

    
    packages = ['mnSpecFit','mnSpecFit.models','mnSpecFit.binning'],

    include_dirs = [numpy.get_include()],                
    
    version = 'v0.0.1',
    
    description = "Fermi Bayesian Spectral Analysis",
    
    author = 'J. Michael Burgess',
    
    author_email = 'jmichaelburgess@gmail.com',
    
    url = 'https://github.com/drJfunk/mnsf',
    
#    download_url = 'https://github.com/giacomov/3ML/archive/v0.0.5',
    
    keywords = ['Likelihood',"GBM", "Spectral"],
    
    #classifiers = [],

 #   libraries = [libsynchrotron],
    
    ext_modules=cythonize(ext_modules),
        
             
    install_requires=[
          'numpy',
          'scipy',
          'numexpr',
           #'emcee',
          'astropy',
          'matplotlib',
          'ipython>=2.0.0, <3.0.0'          
      ])

