from mnfit.mnSpecFit.mnSpecFit import mnSpecFit

from mnfit.mnSpecFit.models.Band import Band
from mnfit.mnSpecFit.models.BandPL import BandPL
from mnfit.mnSpecFit.models.BandBB import BandBB
from mnfit.mnSpecFit.models.BandBBPL import BandBBPL


from mnfit.mnSpecFit.pgstat import pgstat
from mnfit.mnSpecFit.cstat import cstat
from mnfit.mnSpecFit.DataBin import DataBin

from mnfit.mnSpecFit.SpecFitView import SpecFitView as SFV
from mnfit.mnSpecFit.SpecCompare import SpecCompare as SC

class standardAnalysis(object):

    def __init__(self, files, stat="pgstat",ext="pdf"):
        '''
        This makes it easier to initiate a simple model comparison
        using a set of Fermi data. The files argument is intended to
        be a set of DataBin objects that already have their emin and 
        emax keywords set.


        '''

        self.figExt=ext
        if stat == "pgstat":
            self.stat = pgstat
        else:
            self.stat = cstat

        self.files = files    
        self.dmc = mnSpecFit(silent=True, live_points=300, verbose = False, sampling_efficiency="model")

        self._setupData()
        self._run()
        self._analyze()
        


    def _setupData(self):


        self.dmc.LoadData(self.files)

        print
        print "Setting energies....."
        print "\t(Ignore warnings!)"
        print
        
        self.dmc.SetEnergyBounds("LLE",20000,300000)
        self.dmc.SetEnergyBounds("n0",8.1,930)
        self.dmc.SetEnergyBounds("n1",8.1,950)
        self.dmc.SetEnergyBounds("n2",8.1,930)
        self.dmc.SetEnergyBounds("n3",8.1,950)
        self.dmc.SetEnergyBounds("n4",8.1,930)
        self.dmc.SetEnergyBounds("n5",8.1,950)
        self.dmc.SetEnergyBounds("n6",8.1,930)
        self.dmc.SetEnergyBounds("n7",8.1,950)
        self.dmc.SetEnergyBounds("n8",8.1,950)
        self.dmc.SetEnergyBounds("n9",8.1,950)
        self.dmc.SetEnergyBounds("na",8.1,950)
        self.dmc.SetEnergyBounds("nb",8.1,950)
        self.dmc.SetEnergyBounds("b0",256.,39000.)
        self.dmc.SetEnergyBounds("b1",256.,39000.)


    def _run(self):

        self.dmc.SetLikelihood(self.stat)

        #Band Fit
        self.dmc.SetBasename("standAnalysis/band_")
        self.dmc.SetSaveFile("band_params.json")


        self.dmc.SetModel(Band)

        self.dmc.ConstructLikelihood()
        self.dmc.Explore()


        #Band BB Fit
        self.dmc.SetBasename("standAnalysis/bandBB_")
        self.dmc.SetSaveFile("bandBB_params.json")


        self.dmc.SetModel(BandBB)

        self.dmc.ConstructLikelihood()
        self.dmc.Explore()

        #Band PL Fit
        self.dmc.SetBasename("standAnalysis/bandPL_")
        self.dmc.SetSaveFile("bandPL_params.json")


        self.dmc.SetModel(BandPL)

        self.dmc.ConstructLikelihood()
        self.dmc.Explore()


        #Band BB PL Fit
        self.dmc.SetBasename("standAnalysis/bandBBPL_")
        self.dmc.SetSaveFile("bandBBPL_params.json")


        self.dmc.SetModel(BandBBPL)

        self.dmc.ConstructLikelihood()
        self.dmc.Explore()

    def _analyze(self):

        print
        print "Analyzing MC chains....."
        print

        #Band
        sfv = SFV("standAnalysis/band_params.json")
        ax = sfv.ViewMarginals(9)
        ax.get_figure().savefig("band_marg."+self.figExt,bbox_inches="tight")
        ax = sfv.PlotvFv(10)
        ax.get_figure().savefig("band_vFv."+self.figExt,bbox_inches="tight")

        #Band BB
        sfv = SFV("standAnalysis/bandBB_params.json")
        ax = sfv.ViewMarginals(11)
        ax.get_figure().savefig("bandBB_marg."+self.figExt,bbox_inches="tight")
        ax = sfv.PlotvFv(12)
        ax.get_figure().savefig("bandBB_vFv."+self.figExt,bbox_inches="tight")
        ax = sfv.PlotvFvComponents(13)
        ax.get_figure().savefig("bandBB_vFv_comps."+self.figExt,bbox_inches="tight")

        #Band PL
        sfv = SFV("standAnalysis/bandPL_params.json")
        ax = sfv.ViewMarginals(14)
        ax.get_figure().savefig("bandPL_marg."+self.figExt,bbox_inches="tight")
        ax = sfv.PlotvFv(15)
        ax.get_figure().savefig("bandPL_vFv."+self.figExt,bbox_inches="tight")
        ax = sfv.PlotvFvComponents(16)
        ax.get_figure().savefig("bandPL_vFv_comps."+self.figExt,bbox_inches="tight")

        #Band BB PL
        sfv = SFV("standAnalysis/bandBBPL_params.json")
        ax = sfv.ViewMarginals(17)
        ax.get_figure().savefig("bandBBPL_marg."+self.figExt,bbox_inches="tight")
        ax = sfv.PlotvFv(18)
        ax.get_figure().savefig("bandBBPL_vFv."+self.figExt,bbox_inches="tight")
        ax = sfv.PlotvFvComponents(19)
        ax.get_figure().savefig("bandBBPL_vFv_comps."+self.figExt,bbox_inches="tight")

        sc = SC(["standAnalysis/band_params.json","standAnalysis/bandBB_params.json","standAnalysis/bandPL_params.json","standAnalysis/bandBBPL_params.json"])
        
