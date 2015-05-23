import astropy.io.fits as fits

from numpy import zeros



class PHAMaker(object):


    def __init__(self,dir,det,duration,rsp,emin,emax):


        self.dir = dir
        self.det = det
        self.duration = duration
        self.rsp = rsp

        self.emin = emin
        self.emax = emax
        
        
        self.bkfile = self.dir+self.det+".bak"
        self.phafile = self.dir+self.det+".pha"



    def SetCounts(self,counts):
        self.counts = counts

    def SetBak(self,bkg,bkgErr):

        self.bkg = bkg
        self.bkgErr = bkgErr 

    def Write(self):

        self._MakePHA()
        self._MakeBAK()
        print
        print "Writing PHA files for XSPEC"
        print
        self.pha.writeto(self.phafile,clobber=True)
        self.bak.writeto(self.bkfile,clobber=True)

        

    def _MakePHA(self):

        # Initialize 
        priHDU = fits.PrimaryHDU()


        
        # Counts ext
        col1 = fits.Column(name='CHANNEL', format='I', array=range(1,len(self.counts)+1))
        col2 = fits.Column(name='COUNTS', format='J', array=self.counts)
        col3 = fits.Column(name='STAT_ERR', format='J', array=zeros(len(self.counts)))

        cols = fits.ColDefs([col1, col2, col3])

        countsHDU = fits.TableHDU.from_columns(cols)

        countsHDU.header['EXTNAME'] = 'SPECTRUM'                    
        countsHDU.header['TNULL1']  =                -9999                                     
        countsHDU.header['TNULL2']  =                -9999                                     
        countsHDU.header['HDUCLASS']= 'OGIP    '                        
        countsHDU.header['HDUCLAS1']= 'SPECTRUM'                    
        countsHDU.header['HDUVERS1']= '1.1.0   '           
        countsHDU.header['HDUVERS'] = '1.1.0   '                
        countsHDU.header['HDUCLAS2']= 'TOTAL   '     #/ Gross PHA Spectrum (source + bkgd)             
        countsHDU.header['HDUCLAS3']= 'COUNT   '     #/ PHA data stored as Counts (not count/s)        
        countsHDU.header['TLMIN1']  =              1 #/ Lowest legal channel number                    
        countsHDU.header['TLMAX1']  =              128 #/ Highest legal channel number                   
        countsHDU.header['TELESCOP']= 'GLAST   '           #/ mission/satellite name                         
        countsHDU.header['INSTRUME']= 'GBM     '           #/ instrument/detector name                       
        countsHDU.header['DETNAM']  = self.det           #/ specific detector name in use                  
        countsHDU.header['FILTER']  = 'none    '           #/ filter in use                                  
        countsHDU.header['EXPOSURE']=         self.duration #/ exposure (in seconds)                          
        countsHDU.header['AREASCAL']=         1.000000E+00 #/ area scaling factor                            
        countsHDU.header['BACKFILE']= self.bkfile #/ associated background filename       
        countsHDU.header['BACKSCAL']=         1.000000E+00 #/ background file scaling factor                 
        countsHDU.header['CORRFILE']= 'NONE    '           #/ associated correction filename                 
        countsHDU.header['CORRSCAL']=         1.000000E+00 #/ correction file scaling factor                 
        countsHDU.header['RESPFILE']= self.rsp #/ associated redistrib matrix filename 
        countsHDU.header['ANCRFILE']= 'NONE    '           #/ associated ancillary response filename         
        countsHDU.header['PHAVERSN']= '1992a   '           #/ obsolete                                       
        countsHDU.header['DETCHANS']=                  128 #/ total number possible channels                 
        countsHDU.header['CHANTYPE']= 'PHA     '           #/ channel type (PHA, PI etc)                     
        countsHDU.header['POISSERR']=                    True #/ Poissonian errors not applicable               
        countsHDU.header['SYS_ERR'] =                    0 #/ no systematic error specified                  
        countsHDU.header['GROUPING']=                    0 #/ no grouping of the data has been defined       
        countsHDU.header['QUALITY'] =                    0 #/ no data quality information specified

        # Ebounds ext

        col1 = fits.Column(name='CHANNEL', format='I', array=range(len(self.counts)))
        col2 = fits.Column(name='E_MAX', format='1E', unit="keV",array=self.emax)
        col3 = fits.Column(name='E_MIN', format='1E', unit="keV",array=self.emin)

        cols = fits.ColDefs([col1, col2, col3])

        ebHDU = fits.TableHDU.from_columns(cols)

        

        ebHDU.header['EXTNAME'] = 'EBOUNDS '           #/ Extension name                                 
        ebHDU.header['TELESCOP']= 'GLAST   '           #/ Name of mission/satellite                      
        ebHDU.header['INSTRUME']= 'GBM     '           #/ Specific instrument used for observation       
        ebHDU.header['FILTER']  = 'none    '           #/ Instrument filter in use (if any)              
        ebHDU.header['CHANTYPE']= 'PHA     '           #/ No corrections have been applied               
        ebHDU.header['DETCHANS']=                  128 #/ Total number of detector channels available.   
        ebHDU.header['HDUCLASS']= 'OGIP    '           #/ Format confirms to OGIP standard               
        ebHDU.header['HDUCLAS1']= 'RESPONSE'           #/ Extension contains response data               
        ebHDU.header['HDUCLAS2']= 'EBOUNDS '           #/ Extension contains response matrix             
        ebHDU.header['HDUVERS'] = '1.2.0   '           #/ Version number of the format                   
        ebHDU.header['DATE-OBS']= '2011-07-21T04:47:11' #/ Date of start of observation                  
        ebHDU.header['DATE-END']= '2011-07-21T04:52:46' #/ Date of end of observation                    
        ebHDU.header['TIMEUNIT']= 's       '           #/ Time since MJDREF, used in TSTART and TSTOP    
        ebHDU.header['TIMEZERO']=                   0. #/ clock correction                               
        ebHDU.header['TIMESYS'] = 'TT      '           #/ Time system used in time keywords              
        ebHDU.header['TIMEREF'] = 'LOCAL   '           #/ reference frame used for times                 
        ebHDU.header['CLOCKAPP']=                    False #/ was a clock drift correction applied           
        ebHDU.header['GPS_OUT'] =                    False #/ whether GPS time was unavailable at any time

        self.pha = fits.HDUList([priHDU,countsHDU,ebHDU])

    def _MakeBAK(self):

        # Initialize 
        priHDU = fits.PrimaryHDU()


        
        # Counts ext
        col1 = fits.Column(name='CHANNEL', format='I', array=range(1,len(self.counts)+1))
        col2 = fits.Column(name='RATE', format='E', array=self.bkg)
        col3 = fits.Column(name='STAT_ERR', format='E', array=self.bkgErr)
        col4 = fits.Column(name='SYS_ERR', format='E', array=zeros(len(self.counts)))
        
        cols = fits.ColDefs([col1, col2, col3, col4])

        countsHDU = fits.TableHDU.from_columns(cols)

        countsHDU.header['EXTNAME'] = 'SPECTRUM'                    
        countsHDU.header['TNULL1']  =                -9999                                     
        countsHDU.header['TNULL2']  =                -9999                                     
        countsHDU.header['HDUCLASS']= 'OGIP    '                        
        countsHDU.header['HDUCLAS1']= 'SPECTRUM'                    
        countsHDU.header['HDUVERS1']= '1.1.0   '           
        countsHDU.header['HDUVERS'] = '1.1.0   '                
        countsHDU.header['HDUCLAS2']= 'TOTAL   '     #/ Gross PHA Spectrum (source + bkgd)             
        countsHDU.header['HDUCLAS3']= 'COUNT   '     #/ PHA data stored as Counts (not count/s)        
        countsHDU.header['TLMIN1']  =              1 #/ Lowest legal channel number                    
        countsHDU.header['TLMAX1']  =              128 #/ Highest legal channel number                   
        countsHDU.header['TELESCOP']= 'GLAST   '           #/ mission/satellite name                         
        countsHDU.header['INSTRUME']= 'GBM     '           #/ instrument/detector name                       
        countsHDU.header['DETNAM']  = self.det           #/ specific detector name in use                  
        countsHDU.header['FILTER']  = 'none    '           #/ filter in use                                  
        countsHDU.header['EXPOSURE']= self.duration          #/ exposure (in seconds)                          
        countsHDU.header['AREASCAL']=         1.000000E+00 #/ area scaling factor                            
        countsHDU.header['BACKSCAL']=         1.000000E+00 #/ background file scaling factor                 
        countsHDU.header['CORRFILE']= 'NONE    '           #/ associated correction filename                 
        countsHDU.header['CORRSCAL']=         1.000000E+00 #/ correction file scaling factor                 
        countsHDU.header['ANCRFILE']= 'NONE    '           #/ associated ancillary response filename         
        countsHDU.header['PHAVERSN']= '1992a   '           #/ obsolete                                       
        countsHDU.header['DETCHANS']=                  128 #/ total number possible channels                 
        countsHDU.header['CHANTYPE']= 'PHA     '           #/ channel type (PHA, PI etc)                     
        countsHDU.header['POISSERR']=                    False #/ Poissonian errors not applicable               
        countsHDU.header['GROUPING']=                    0 #/ no grouping of the data has been defined       
        countsHDU.header['QUALITY'] =                    0 #/ no data quality information specified

        # Ebounds ext

        col1 = fits.Column(name='CHANNEL', format='I', array=range(len(self.counts)))
        col2 = fits.Column(name='E_MAX', format='1E', unit="keV",array=self.emax)
        col3 = fits.Column(name='E_MIN', format='1E', unit="keV",array=self.emin)

        cols = fits.ColDefs([col1, col2, col3])

        ebHDU = fits.TableHDU.from_columns(cols)

        

        ebHDU.header['EXTNAME'] = 'EBOUNDS '           #/ Extension name                                 
        ebHDU.header['TELESCOP']= 'GLAST   '           #/ Name of mission/satellite                      
        ebHDU.header['INSTRUME']= 'GBM     '           #/ Specific instrument used for observation       
        ebHDU.header['FILTER']  = 'none    '           #/ Instrument filter in use (if any)              
        ebHDU.header['CHANTYPE']= 'PHA     '           #/ No corrections have been applied               
        ebHDU.header['DETCHANS']=                  128 #/ Total number of detector channels available.   
        ebHDU.header['HDUCLASS']= 'OGIP    '           #/ Format confirms to OGIP standard               
        ebHDU.header['HDUCLAS1']= 'RESPONSE'           #/ Extension contains response data               
        ebHDU.header['HDUCLAS2']= 'EBOUNDS '           #/ Extension contains response matrix             
        ebHDU.header['HDUVERS'] = '1.2.0   '           #/ Version number of the format                   
        ebHDU.header['DATE-OBS']= '2011-07-21T04:47:11' #/ Date of start of observation                  
        ebHDU.header['DATE-END']= '2011-07-21T04:52:46' #/ Date of end of observation                    
        ebHDU.header['TIMEUNIT']= 's       '           #/ Time since MJDREF, used in TSTART and TSTOP    
        ebHDU.header['TIMEZERO']=                   0. #/ clock correction                               
        ebHDU.header['TIMESYS'] = 'TT      '           #/ Time system used in time keywords              
        ebHDU.header['TIMEREF'] = 'LOCAL   '           #/ reference frame used for times                 
        ebHDU.header['CLOCKAPP']=                    False #/ was a clock drift correction applied           
        ebHDU.header['GPS_OUT'] =                    False #/ whether GPS time was unavailable at any time

        self.bak = fits.HDUList([priHDU,countsHDU,ebHDU])

        
