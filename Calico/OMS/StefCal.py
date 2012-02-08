# -*- coding: utf-8 -*-

from Timba import pynode
from Timba.Meq import meq
import numpy
import numpy.ma
import math
import Kittens.utils
import time
import cPickle
import os.path
import traceback
_verbosity = Kittens.utils.verbosity(name="stefcal");
dprint = _verbosity.dprint;
dprintf = _verbosity.dprintf;

def GCD (a,b):
  """Return greatest common divisor using Euclid's Algorithm."""
  while b:      
    a,b = b,a%b
  return a;

def LCM (a,b,*args):
  """Return lowest common multiple of two arguments."""
  if not args:
    return a*b//GCD(a,b);
  else:
    return reduce(LCM,[a,b]+list(args));

identity_function = lambda x:x;

def eqkey_direct_conjugate (pp,qq):
  """Helper function:
  Converts p,i,q,j tuple into equation key and conversion function for direct and conjugate values.
  If p<q, then this is just an identity relation
  If p>q, then swaps components around, and uses conjugate
  """;
  if pp<qq:
    return (pp,qq),identity_function,numpy.conj
  else:
    return (qq,pp),numpy.conj,identity_function;

class SubtiledGain (object):
  """Support class to handle a set of subtiled gains.
    For a data shape of N1,N2,..., and a subtiling of M1,M2,..., we have the following shapes in play:
      gainparm_shape = N1/M1,N2/M2,... Let's call this K1,K2,...
    Then, data is reshaped into
      subtiled_shape =  K1,M1,K2,M2,...
    to which we can apply parms with an index of [:,numpy.newaxis,:,numpy.neqaxis] etc.,
    and then collapse every second axis.
    Thus we define the following methods:
      tile_data(x):        reshapes datashape into subtiled_shape 
      untile_data(x):      reshapes subtiled_shape into datashape
      tile_gain(x):        reshapes gain into subtiled shape (by inserting a numpy.newaxis for every second axis)
      reduce_subtiles(x):  reduces subtiled_shape to gain_shape by summing every second axis
    
    If the common case of a 1,1,... subtiling, all these can be an identity
  """;
  def __init__ (self,datashape,subtiling,solve_ifrs,epsilon,conv_quota):
    self._solve_ifrs = set(solve_ifrs);
    self._epsilon = epsilon;
    self.datashape = datashape;
    self.subtiling = subtiling;
    self.gainshape = tuple([ nd/nt for nd,nt in zip(datashape,subtiling) ]);
    # if subtiling is 1,1,... then override methods with identity relations
    if max(subtiling) == 1:
      self.tile_data = self.untile_data = self.tile_gain = self.reduce_subtiles = identity_function;
    else:
      self.subtiled_shape = [];
      self.gain_expansion_slice = [];
      for ng,nt in zip(self.gainshape,subtiling):
        self.subtiled_shape += [ng,nt];
        self.gain_expansion_slice += [slice(None),numpy.newaxis];
      self.subtiled_axes = range(1,len(datashape)*2,2)[-1::-1];
    # init empty parms
    parms = set();
    for pp,qq in self._solve_ifrs:
      parms.update([pp,qq]);
    self._unity = numpy.ones(self.gainshape,dtype=complex);
    self.gain = dict([ (pp,self._unity) for pp in parms ]);
    # setup various counts and convergence targets
    self.total_parms = len(parms)*reduce(lambda a,b:a*b,self.gainshape);
    self.convergence_target = int(self.total_parms*conv_quota);
    self._residual = {};
    self._corrupt = {};
    self._gpgq = {};
    
  # define methods
  def tile_data (self,x):
    return x.reshape(self.subtiled_shape);
  def untile_data (self,x):
    return x.reshape(self.datashape);
  def tile_gain (self,x):
    return x[self.gain_expansion_slice];
  def reduce_subtiles (self,x):
    for ax in self.subtiled_axes:
      x = x.sum(ax);
    return x;
    
  def iterate (self,data,model,verbose=0):
    self._residual = {};
    self._corrupt = {};
    self._gpgq = {};
    gain1 = {};  # dict of gain parm updates from this step
    #
#    if verbose:
#      print [ (data[key].min(),model[key].min()) for key in self._solve_ifrs ];
    # converge from one side (solve for Gp)
    for pp in self.gain.keys():
      # build up sums
      sum_reim = numpy.zeros(self.gainshape,dtype=complex); 
      sum_sq = 0;
      for qq in self.gain.keys():
        eqkey,direct,conjugate = eqkey_direct_conjugate(pp,qq);
        if eqkey in self._solve_ifrs and eqkey in data:
          mh = self.tile_data(conjugate(model[eqkey]))*self.tile_gain(self.gain[qq]);
          dmh = self.tile_data(direct(data[eqkey]))*mh;
          sum_reim += self.reduce_subtiles(dmh);
          sum_sq += self.reduce_subtiles(abs(mh)**2);
      # generate update
      if verbose:
        print pp,sum_sq.min(),sum_sq.max();
      mask = sum_sq==0;
      gain1[pp] = (sum_reim/sum_sq + self.gain[pp])/2;
      gain1[pp][mask] = self.gain[pp][mask];
    if verbose:
      print "done left iter";
    gain2 = {};  # dict of gain parm updates from this step
    # then solve for Gq
    for qq in self.gain.keys():
      # build up sums
      sum_reim = numpy.zeros(self.gainshape,dtype=complex); 
      sum_sq = 0;
      for pp in self.gain.keys():
        eqkey,direct,conjugate = eqkey_direct_conjugate(pp,qq);
        if eqkey in self._solve_ifrs and eqkey in data:
          mh = self.tile_data(direct(model[eqkey]))*self.tile_gain(gain1[pp]);
          dmh = self.tile_data(conjugate(data[eqkey]))*mh;
          sum_reim += self.reduce_subtiles(dmh);
          sum_sq += self.reduce_subtiles(abs(mh)**2);
      # generate update
      if verbose:
        print qq,sum_sq.min(),sum_sq.max();
      mask = sum_sq==0;
      gain2[qq] = (sum_reim/sum_sq + gain1[qq])/2;
      gain2[qq][mask] = gain1[qq][mask];
    # check for convergence
    if verbose:
      print "done right iter";
    self.gaindiff = {};
    self.num_converged = 0;
    for pp in self.gain.iterkeys():
      delta = abs(gain2[pp] - self.gain[pp]);
      self.num_converged += (delta < self._epsilon).sum();
      self.gaindiff[pp] = abs(delta).max();
    self.maxdiff = max(self.gaindiff.itervalues());
    self.gain = gain2;
    return self.num_converged >= self.convergence_target;
    
  def gpgq (self,eqkey):
    """Returns Gp*conj(Gq), tiled into subtile shape. 
    Computes it on-demand, if not already cached""";
    g = self._gpgq.get(eqkey);
    if g is None:
      g = self._gpgq[eqkey] = self.tile_gain( self.gain.get(eqkey[0],self._unity)*
                                        numpy.conj(self.gain.get(eqkey[1],self._unity)) );
    return g;

  def reset_residuals (self):
    self._residual = {};

  def residual (self,d,m,eqkey):
    """Returns residual R = D - Gp*M*conj(Gq), tiled into subtile shape. 
    Computes it on-demand, if not already cached""";
    r = self._residual.get(eqkey);
    if r is None:
      r = self._residual[eqkey] = self.untile_data(self.tile_data(d[eqkey]) - self.tile_data(m[eqkey])*self.gpgq(eqkey));
    return r;
    
  def correct (self,d,eqkey,index=True):
    """Returns corrected data Gp^{-1}*D*Gq^{H-1}."""
    if index:
      d = d[eqkey];
    return self.untile_data(self.tile_data(d)/self.gpgq(eqkey));
    
  def corrupt (self,m,eqkey,cache=False):
    """Returns corrupted model Gp*M*Gq."""
    mc = self._corrupt.get(eqkey); 
    if mc is None or not cache:
      mc = self.untile_data(self.tile_data(m[eqkey])*self.gpgq(eqkey));
      if cache:
        self._corrupt[eqkey] = mc;
    return mc;


class StefCalNode (pynode.PyNode):
  def __init__ (self,*args):
    pynode.PyNode.__init__(self,*args);
    self._dataset_id = None;
    self.ifr_gain = {};

  def update_state (self,mystate):
    """Standard function to update our state""";
    # list of ifr pairs (as p:q strings) corresponding to first axis of tensor input 
    mystate('ifrs',[]);
    # list of ifr pairs which we use to solve for gains. If empty, all ifrs are used
    mystate('solve_ifrs',[]);
    # correlation names
    mystate('corr_names',["x","y"]);
    # labels for gain, ifr gain and differential gain parameters
    mystate('gain_parm_label',"G");
    mystate('ifr_gain_parm_label',"IG");
    mystate('diffgain_parm_label',"dE");
    # convergence criteria
    mystate('epsilon',1e-5);            # updates <epsilon are considered converged
    mystate('max_iter',50);            
    mystate('diffgain_max_iter',5);            
    mystate('max_major',10);            
    mystate('convergence_quota',0.9);   # what percentage of parms should converge
    # subtiling for gains
    mystate('gain_subtiling',[1,1]);
    mystate('diffgain_subtiling',[]);
    # return residuals (else data)
    mystate('residuals',True);
    # return corrected residuals/data (else uncorrected)
    mystate('correct',True);
    # solve for ifr gains as we go along
    mystate('solve_ifr_gains',True);
    # apply previous ifr gain solution, if available
    mystate('apply_ifr_gains',True);
    # name of ifr gain tables
    mystate('ifr_gain_table','ifrgains.cp');
    # verbosity level
    mystate('verbose',0);
    # lis of all ifrs, as p,q pairs
    self._ifrs = [ tuple(x.split(':')) for x in self.ifrs ];
    # parse set of solvable ifrs
    self._solvable_ifrs = set([ tuple(x.split(":")) for x in (self.solve_ifrs or self.ifrs) ]);
    # and update with q:p for every p:q
    self._solvable_ifrs.update([ tuple(pq[-1::-1]) for pq in self._solvable_ifrs ]);
    # other init
    _verbosity.set_verbose(self.verbose);

  def get_result (self,request,*children):
    timestamp0 = time.time();
    # get dataset ID from request
    dataset_id,domain_id = meq.split_request_id(request.request_id);
    # if new dataset ID, do setup for start of new dataset
    if dataset_id != self._dataset_id:
      self._dataset_id = dataset_id;
      dprint(1,"new dataset id",dataset_id);
      # if asked to solve for IFR gains, set up dicts for collecting stats
      if self.solve_ifr_gains:
        self.ig_sum_reim = dict([(((p,i),(q,j)),0j) for p,q in self._ifrs for i in range(2) for j in range(2) ]);
        self.ig_sum_sq = dict([(((p,i),(q,j)),0.) for p,q in self._ifrs for i in range(2) for j in range(2) ]);
        self.ifr_gain_update = {};
      # read previous IFR gains from table, if asked to apply them
      self.ifr_gain = {};
      if self.apply_ifr_gains and os.path.exists(self.ifr_gain_table):
        try:
          self.ifr_gain = cPickle.load(file(self.ifr_gain_table));
          dprint(1,"loaded %d ifr gains from %s"%(len(self.ifr_gain),self.ifr_gain_table));
        except:
          traceback.print_exc();
          dprint(1,"error loading gains from",self.ifr_gain_table);
          
    # child 0 is data
    # child 1 is direction-independent model
    # children 2 and on are models subject to dE terms
    num_diffgains = len(children)-2;
    if num_diffgains < 0:
      raise TypeError,"StefCalNode: at least 2 children (data, model) must be provided";

    # check inputs and populate mappings
    piqj_all = [];      # list of all (p,i),(q,j) pairs 
    piqj_data = [];     # subset of (p,i),(q,j) pairs for which we have non-null input
    piqj_solvable = []; # subset of (p,i),(q,j) pairs for which we solve for gains
    data  = {};         # mapping from (p,i),(q,j) to a data time-freq plane 
    model0 = {};        # mapping from (p,i),(q,j) to an M0 model time-freq plane 
    dgmodel = [ {} for i in range(num_diffgains) ]; 
                        # for each diff gain, mapping from (p,i),(q,j) to an M1,M2,... model time-freq plane 
    model = {};         # this is the full model, M0+M1+M2
    
    #
    datares = children[0]
    modelres = children[1];
    if any( [ ch.dims != datares.dims for ch in children[1:] ] ):
      raise TypeError,"tensor dimensions of data and model(s) must match";
    # expecting Nx2x2 matrices
    if len(datares.dims) == 3:
      if datares.dims[1] != 2 or datares.dims[2] != 2: 
        raise TypeError,"data and model must be of rank Nx2x2";
      nifrs = datares.dims[0];
      # setup antenna names
      if nifrs != len(self.ifrs):
        raise TypeError,"first dimension of data and model must match the number of interferometers in the ifrs field";
      # setup list of data, values and parameter names
      piqj_all = [ ((p,i),(q,j)) for p,q in self._ifrs for i in range(2) for j in range(2) ];
      for nvells,((p,i),(q,j)) in enumerate(piqj_all):
        eqkey = (p,i),(q,j);
        # get data, and apply ifr gains if we have them
        d = getattr(datares.vellsets[nvells],'value',0)*self.ifr_gain.get(eqkey,1);
        # get model
        m = getattr(modelres.vellsets[nvells],'value',0);
        if hasattr(datares.vellsets[nvells],'flags'):
          fl = datares.vellsets[nvells].flags != 0;
          dprint(2,eqkey,"has",fl.sum(),"flags");
          d[fl!=0] = 0;
          m[fl!=0] = 0;
        else:
          fl = None;
        # if model or data is null, then we're unpolarized, so skip this from the matrix entirely
        if ( numpy.isscalar(m) and not m ) or ( m.size == 1 and not m.ravel()[0] ) or \
            ( numpy.isscalar(d) and not d ):
          pass;
        else:
          # if this is the first datum, then check shape, and prepare subtilings etc.
          if not nvells:
            # this is the basic time-frequency shape
            datashape = tuple(d.shape);
            # figure out subtiling
            # if not specified, use whole tile as solution interval
            gain_subtiling = lcm_subtiling = self.gain_subtiling or datashape;
            # replace nulls in subtiling with solution interval
            gain_subtiling = [ gs or ds for gs,ds in zip(gain_subtiling,datashape) ];
            if len(gain_subtiling) != len(datashape):
              raise ValueError,"gain_subtiling vector must have the same length as the data shape";
            if min(gain_subtiling) < 1:
              raise ValueError,"invalid gain_subtiling %s"%self.gain_subtiling;
            # if diffgains are also present, then work out the least-common-multiple subtiling
            if num_diffgains:
              dg_subtiling = self.diffgain_subtiling or datashape;
              dg_subtiling = [ gs or ds for gs,ds in zip(dg_subtiling,datashape) ];
              if len(dg_subtiling) != len(datashape):
                raise ValueError,"diffgain_subtiling vector must have the same length as the data shape";
              if min(dg_subtiling) < 1:
                raise ValueError,"invalid diffgain_subtiling %s"%dg_subtiling;
              lcm_subtiling = [ LCM(a,b) for a,b in zip(gain_subtiling,dg_subtiling) ];
            # data must be expanded to match the LCM subtiling
            expanded_datashape = tuple([ (nd/np+(1 if nd%np else 0))*np for nd,np in zip(datashape,lcm_subtiling) ]);
            dprint(1,"gain parm LCM subtiling is",lcm_subtiling);
            # if tiling does not tile the data shape perfectly, we'll need to expand the input arrays
            # Define pad_array() as a function for this: it will be identity if no expansion is needed
            if datashape != expanded_datashape:
              dprint(1,"input arrays will be expanded to shape",expanded_datashape);
              expanded_dataslice = tuple([ slice(0,nd) for nd in datashape ]);
              def pad_array (x):
                x1 = numpy.zeros(expanded_datashape,dtype=x.dtype);
                x1[expanded_dataslice] = x;
                return x1;
            else:
              expanded_dataslice = None;
              pad_array = identity_function;
          # now check inputs and add them to data and model dicts
          if d.shape != datashape:
            raise TypeError,"data shape mismatch at %s:%s:%s:%s"%(p,q,self.corr_names[i],self.corr_names[j]);
          if m.shape != datashape:
            raise TypeError,"model shape mismatch at %s:%s:%s:%s"%(p,q,self.corr_names[i],self.corr_names[j]);
          piqj_data.append(eqkey);
          if (p,q) in self._solvable_ifrs:
            piqj_solvable.append(eqkey);
          # add to data/model matrices, applying the padding function defined above
          m0 = model0[eqkey] = pad_array(m);
          data[eqkey]  = pad_array(d);
          # also accumulate initial model, as M0+M1+M2
          if num_diffgains:
            m0 = model[eqkey] = m0.copy();
            for i in range(num_diffgains):
              m1 = children[2+i].vellsets[nvells].value
              if fl is not None:
                m1[fl] = 0;
              m1 = pad_array(m1);
              dgmodel[i][eqkey] = m1;
              m0 += m1;
              
    else:
      # in principle could also handle [N], but let's not bother for now
      raise TypeError,"data and model must be of rank Nx2x2";
    # init gain parms object
    gain = SubtiledGain(expanded_datashape,gain_subtiling,piqj_solvable,self.epsilon,self.convergence_quota);
    dprintf(0,"solving for %d gain parms using %d of %d inteferometers\n",
      len(gain.gain),
      len(self._solvable_ifrs)/2,len(self.ifrs));
    dprint(1,"convergence target",gain.convergence_target,"of",gain.total_parms,"parms");

    # init diffgains
    if num_diffgains:
      diffgains = [ SubtiledGain(expanded_datashape,dg_subtiling,piqj_solvable,self.epsilon,self.convergence_quota) 
        for i in range(num_diffgains) ];
      dg0 = diffgains[0];
      dprintf(0,"also solving for %dx%d differential gains\n",num_diffgains,len(dg0.gain));
      dprint(1,"convergence target for each is ",dg0.convergence_target,"of",dg0.total_parms,"parms");
    else:
      diffgains = [];
      model = model0;
    
    def compute_chisq ():
      chisq = 0;
      nterms = 0;
      for eqkey in piqj_solvable:
        r = gain.residual(data,model,eqkey);
        chisq += (r*numpy.conj(r)).sum();
        nterms += r.size;
      return chisq/nterms;
    
    # start major loop -- alternates over gains and diffgains
    for nmajor in range(self.max_major+1):
      # first iterate normal gains to convergence
      for niter in range(self.max_iter):
        # iterate over normal gains
        converged = gain.iterate(data,model);
        # check chi-square
        if ( niter and not niter%10 ) or niter >= self.max_iter-1 or converged:
          chisq = compute_chisq();
          dprint(3,"iter %d max gain update is %g converged %.2f chisq is %g"%(niter+1,
                    gain.maxdiff,gain.num_converged/float(gain.total_parms),chisq));
        # break out if converged
        if converged:
          break;
      dprint(1,"gains converge to chisq %g (last G update %g) after %d iterations"%(chisq,gain.maxdiff,niter+1));
      # break out if no diffgains to iterate over, or if we're on the last major cycle
      if not num_diffgains or nmajor >= self.max_major:
        break;
      else:
        # model is the full model, M0+corrupt(M1)+corrupt(M2)+.... 
        # subtract this from corrected data: D1 = correct(D)-M0-corrupt(M1)-corrupt(M2)-...
        data1 = dict([ (eqkey,gain.correct(data,eqkey)-model[eqkey]) for eqkey in piqj_solvable ]);
        # now loop over all diffgains and iterate each set once
        for i,dg in enumerate(diffgains):
          # add current estimate of corrupt(Mi) back into data1, and subtract from model
          for eqkey in piqj_solvable:
            corr = dg.corrupt(dgmodel[i],eqkey,cache=True);
            data1[eqkey] += corr;
            model[eqkey] -= corr;
          # iterate this diffgain solution
          for niter in range(self.diffgain_max_iter):
            if dg.iterate(data1,dgmodel[i]):
              break;
          dprint(2,"diffgain #%d converged after %d iterations, value is %s"%(i,niter,dg.gain.values()[0][0,0]));
          # add back to model, and subtract from data1 if needed
          for eqkey in piqj_solvable:
            corr = dg.corrupt(dgmodel[i],eqkey,cache=True);
            model[eqkey] += corr;
            if i<num_diffgains-1:
              data1[eqkey] -= corr;
      # done iterating over diffgains, reset cached residuals and compute chisq once again
      gain.reset_residuals();
      chisq = compute_chisq();
      dprint(1,"diffgains converge to chisq %g after %d major cycles"%(chisq,nmajor+1));
      
    # if we were solving for diffgains, then model is not completely up-to-date, since the non-solvable baselines
    # have been ignored. Fill them in here. Also, reset residuals
    if num_diffgains:
      for eqkey in set(piqj_data)-set(piqj_solvable):
        model[eqkey] = model0[eqkey];
        for i,dg in enumerate(diffgains):
          model[eqkey] += dg.corrupt(dgmodel[i],eqkey,cache=True);

    # work out result -- residual or corrected visibilities, depending on our state
    nvells = maxres = 0;
    for nvells,eqkey in enumerate(piqj_all):
      m = model.get(eqkey,None);
      if m is not None:
        r = gain.residual(data,model,eqkey);
        out = r if self.residuals else data[eqkey];
        if self.correct:
          out = gain.correct(out,eqkey,index=False);
        datares.vellsets[nvells].value[...] = out[expanded_dataslice] if expanded_dataslice else out;
        # compute stats
        maxres = max(maxres,abs(r).max());

    # figure out domain ID
    time0,time1,timestep,numtime,freq0,freq1,freqstep,numfreq = request.cells.domain.domain_id;

    # update IFR gain solutions, if asked to
    if self.solve_ifr_gains:
      for eqkey in piqj_data:
        d = data[eqkey];
        dh = numpy.conj(d);
        # work out update to ifr gains
        sri = self.ig_sum_reim[eqkey] = self.ig_sum_reim[eqkey] + (gain.corrupt(model,eqkey)*dh).sum();
        ssq = self.ig_sum_sq[eqkey]   = self.ig_sum_sq[eqkey] + (d*dh).sum();
        self.ifr_gain_update[eqkey] = sri/ssq;
    
      # if last domain, then write to file
      if time1 >= numtime:
        # apply updates
        for eqkey in piqj_data:
          self.ifr_gain[eqkey] = self.ifr_gain.get(eqkey,1)*self.ifr_gain_update[eqkey];
        dprint(2,"IFR gain solutions update: ",", ".join(
              ["%s%s:%s%s %s"%(p,self.corr_names[i],q,self.corr_names[j],
              self.ifr_gain_update[(p,i),(q,j)])
              for (p,i),(q,j) in piqj_data[0:3]]));
        dprint(2,"IFR gain solutions: ",", ".join(
              ["%s%s:%s%s %s"%(p,self.corr_names[i],q,self.corr_names[j],
              self.ifr_gain[(p,i),(q,j)])
              for (p,i),(q,j) in piqj_data[0:3]]));
        # save
        try:
          cPickle.dump(self.ifr_gain,file(self.ifr_gain_table,'w'));
          dprint(1,"saved %d ifr gains to %s"%(len(self.ifr_gain_update),self.ifr_gain_table));
        except:
          traceback.print_exc();
          dprint(0,"error saving ifr gains to",self.ifr_gain_table);
          
    dprint(0,"%s residual max %g chisq %g (last G update %g) after %d iterations and %.2f seconds"%(
              request.request_id,maxres,chisq,gain.maxdiff,niter,time.time()-timestamp0));

    return datares;

