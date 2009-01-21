from Timba.TDL import *
from Timba.Meq import meq
from Lions.PiercePoints.PiercePoints import *
import Meow


def compile_options():
    return [TDLCompileOption("beta","Beta",[5./3.],more=float,doc="""Beta"""),
            TDLCompileOption("N","N",[20],more=int,doc="""half the gridsize"""),
            TDLCompileOption("speedx","Speed X (/s)",[100.],more=float,doc="""Speed in x direction"""),
            TDLCompileOption("speedy","Speed Y (/s)",[100.],more=float,doc="""Speed in y direction"""),
            TDLCompileOption("scale","Scale",[100.],more=float,doc="""Scale size of the grid"""),
            TDLCompileOption("amp_scale","Amplitude Scale",[1.e-5],more=float,doc="""Scale of the TEC values"""),
            TDLCompileOption("seed_nr","Seed",[None],more=int,doc="""Seeding of the random generator""")]
##    return [TDLCompileOption("scale","Scale",[100.],more=float,doc="""Scale size of the grid"""),];



class MIM(PiercePoints):
    """Create MIM_model with Kolmogorov phase screen"""
    def __init__(self,ns,name,sources,stations=None,height=300,ref_station=None,tags="iono",make_log=False):
        PiercePoints.__init__(self,ns,name,sources,stations,height,make_log);
        self.ref_station=ref_station;
        #        init_phasescreen(N,beta);
        
    def make_tec(self):
        # fit a virtual TEC value, this makes life easier (eg. include freq. and sec dependence)
        pp=self.make_pp(ref_station=self.ref_station);
        n=0;
        for src in self.src:
            for station in self.stations:
                Kol_node=self.create_Kol_node(pp,src,station);
                sec = self.ns['sec'](src,station);
                if not self.ns['tec'](src,station).initialized():
                      self.ns['tec'](src,station)<<Kol_node*sec;
                n=n+1;
 
      
        return self.ns['tec'];
       


    def create_Kol_node(self,pp,src,station):
        Kol_node=self.ns['Kol_node'](src,station);
        if Kol_node.initialized():
            return Kol_node;
        kl=self.ns['Kol_node'](src,station)  <<  Meq.PyNode(children=(pp(src,station),),class_name="KolmogorovNode",module_name="Lions.PiercePoints.modules.KolmogorovNode",grid_size=N,beta=beta,scale=scale,speedx=speedx,speedy=speedy,amp_scale=amp_scale,seed_nr=seed_nr);
        return kl;
                
