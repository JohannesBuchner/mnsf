from Band import Band
from SBPL import SBPL 
from Synchrotron import Synchrotron
#from SynchrotronComplex import SynchrotronComplex
#from SynchrotronBB import SynchrotronBB
from SynchrotronCutoff import SynchrotronCutoff
from SynchrotronCutoffFixed import SynchrotronCutoffFixed
#from PLSynchrotron import PLSynchrotron
#from PLSynchrotron_Cutoff import PLSynchrotron_Cutoff
#from FastSynchrotron import FastSynchrotron
from TsviSlow import TsviSlow
from TsviFast import TsviFast
#from SynchSSC import SynchSSC
from CPL import CPL
from BandCO import BandCO
#from ZhaoSynchrotron import ZhaoSynchrotron
#from FastSynchrotronBB import FastSynchrotronBB
#from SynchSSC_BB import SynchSSC_BB
from BB import BB
from PL import PL
#from SubPhoto import SubPhoto
#from BandTable import BandTable

models = {"Band":Band,\
          "Synchrotron":Synchrotron,\
 #         "SynchrotronComplex":SynchrotronComplex,\
 #         "SynchrotronBB":SynchrotronBB ,\
          "SBPL":SBPL,\
          "SynchrotronCutoff":SynchrotronCutoff,\
          "SynchrotronCutoffFixed":SynchrotronCutoffFixed,\
 #         "PLSynchrotron":PLSynchrotron,\
 #         "PLSynchrotron_Cutoff":PLSynchrotron_Cutoff,\
          "TsviSlow":TsviSlow,\
          "TsviFast":TsviFast,\
 #         "FastSynchrotron":FastSynchrotron,\
 #         "SynchSSC":SynchSSC,\
          "CPL":CPL,\
          "BandCO":BandCO,\
#          "ZhaoSynchrotron":ZhaoSynchrotron,\
#          "FastSynchrotronBB":FastSynchrotronBB,\
 #         "SynchSSC_BB":SynchSSC_BB,\
          "BB":BB,\
          "PL":PL,\
#          "SubPhoto":SubPhoto,\
#          "BandTable":BandTable,\
}
