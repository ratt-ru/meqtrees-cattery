# -*- coding: utf-8 -*-

from Timba import pynode
from Timba.Meq import meq
import numpy
import math
import Kittens.utils
import time
import cPickle
import os.path
import traceback

from MatrixOps import *
from SubtiledDiagGain import SubtiledDiagGain
from Subtiled2x2Gain import Subtiled2x2Gain

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
    # full polarization or diagonal
    mystate('full_polarization',False);
    # convergence criteria
    mystate('epsilon',1e-5);            # updates <epsilon are considered converged
    mystate('max_iter',50);
    mystate('diffgain_max_iter',5);
    mystate('max_major',10);
    mystate('convergence_quota',0.9);   # what percentage of parms should converge
    # subtiling for gains
    mystate('gain_subtiling',[1,1]);
    mystate('diffgain_subtiling',[]);
    # use stored solution (if available) as starting guess
    mystate('init_from_table',True);
    # use previous tile (timeslot) as starting guess -- if table not available
    mystate('init_from_previous',True);
    # use this value as starting guess -- if previous two not available
    mystate('init_value',1);
    # regularization factor applied to gain solutions for correction
    mystate('regularization_factor',0);
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
    # other init
    _verbosity.set_verbose(self.verbose);
    _verbosity.enable_timestamps(True,modulo=6000);
    # initial value from which to start iterating
    self._init_value_gain = self.init_value;
    self._init_value_dg = {};

  def get_result (self,request,*children):
    timestamp0 = time.time();
    # get dataset ID from request
    dataset_id,domain_id = meq.split_request_id(request.request_id);
    # get domain ID from request
    time0,time1,timestep,numtime,freq0,freq1,freqstep,numfreq = request.cells.domain.domain_id;

    # if new dataset ID, do setup for start of new dataset
    if dataset_id != self._dataset_id:
      self._dataset_id = dataset_id;
      self._init_value_gain = self.init_value;
      self._init_value_dg = {};
      dprint(1,"new dataset id",dataset_id);
      # if asked to solve for IFR gains, set up dicts for collecting stats
      if self.solve_ifr_gains:
        self.ig_sum_reim = dict([ (pq,[0j]*4) for pq in self._ifrs ]);
        self.ig_sum_sq   = dict([ (pq,[0.]*4) for pq in self._ifrs ]);
        self.ifr_gain_update = dict([ (pq,[0.]*4) for pq in self._ifrs ]);
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
    pqij_all = [];      # list of all (p,q),i,j tuples
    pqij_data = [];     # subset of (p,q),i,j tuples for which we have non-null input
    pqij_solvable = []; # subset of (p,q),i,j tuples for which we solve for gains
    data  = {};         # mapping from (p,q) to four data time-freq planes
    model0 = {};        # mapping from (p,q) to four model (M0) time-freq planes
    dgmodel = [ {} for i in range(num_diffgains) ];
                        # for each diff gain, mapping from (p,q) to M1,M2,... model time-freq planes (4 each)
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
      nvells = 0;
      for pq in self._ifrs:
        # get IFR gain for this p,q
        ifrgain = self.ifr_gain.get(pq,[1,1,1,1]);
        flag4 = False;
        # now loop over the 4 matrix elements
        for num,(i,j) in enumerate(IJ2x2):
          # get data, and apply ifr gains if we have them
          d = getattr(datares.vellsets[nvells],'value',0);
          if not is_null(d):
            d *= ifrgain[num];
          # get model
          m = getattr(modelres.vellsets[nvells],'value',0);
          if hasattr(datares.vellsets[nvells],'flags'):
            flag4 += (datares.vellsets[nvells].flags != 0);
          # if model and/or data is null, then we're unpolarized, so skip this from the matrix entirely
          skip = is_null(d) if self.full_polarization else (is_null(m) or is_null(d));
          if skip:
            pass;
          else:
            # if this is the first datum, then check shape, and prepare subtilings etc.
            if not nvells:
              # this is the basic time-frequency shape
              datashape = tuple(d.shape);
              # figure out subtiling
              # if not specified, use whole tile as solution interval
              if self.gain_subtiling:
                # replace nulls in subtiling with solution interval
                lcm_subtiling = gain_subtiling = [ min(gs,ds) or ds for gs,ds in zip(self.gain_subtiling,datashape) ];
              else:
                lcm_subtiling = gain_subtiling = datashape;
              if len(gain_subtiling) != len(datashape):
                raise ValueError,"gain_subtiling vector must have the same length as the data shape";
              if min(gain_subtiling) < 1:
                raise ValueError,"invalid gain_subtiling %s"%self.gain_subtiling;
              # if diffgains are also present, then work out the least-common-multiple subtiling
              if num_diffgains:
                dg_subtiling = self.diffgain_subtiling or datashape;
                dg_subtiling = [ min(gs,ds) or ds for gs,ds in zip(dg_subtiling,datashape) ];
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
              print d.shape,datashape,d;
              raise TypeError,"data shape mismatch at %s:%s:%s:%s"%(pq[0],pq[1],self.corr_names[i],self.corr_names[j]);
            if not is_null(m) and m.shape != datashape:
              raise TypeError,"model shape mismatch at %s:%s:%s:%s"%(pq[0],pq[1],self.corr_names[i],self.corr_names[j]);
            # add to data/model matrices, applying the padding function defined above
            m0 = model0.setdefault(pq,[0,0,0,0])[num] = pad_array(m);
            data.setdefault(pq,[0,0,0,0])[num]  = pad_array(d);
            # also accumulate initial model, as M0+M1+M2
            # if max_major==0, then we don't solve for diff
            if num_diffgains:
              m0 = model.setdefault(pq,[0,0,0,0])[num] = m0.copy();
              for k in range(num_diffgains):
                m1 = children[2+k].vellsets[nvells].value
                m1 = pad_array(m1);
                dgmodel[k].setdefault(pq,[0,0,0,0])[num] = m1;
                m0 += m1;
            else:
              model[pq] = model0[pq];
          # increment vells #
          nvells += 1;
        # do we have any flags? apply them to all 4 corrs
        if flag4 is not False:
          dprint(4,pq,"has",flag4.sum(),"flags");
          for collection in [data,model,model0] + dgmodel:
            for x in collection[pq]:
              if not is_null(x):
                x[pad_array(flag4)] = 0;
    else:
      # in principle could also handle [N], but let's not bother for now
      raise TypeError,"data and model must be of rank Nx2x2";

    GainClass = Subtiled2x2Gain if self.full_polarization else SubtiledDiagGain;
    # init gain parms object
    gain = GainClass(expanded_datashape,gain_subtiling,self._solvable_ifrs,
              self.epsilon,self.convergence_quota,
              regularization_factor=self.regularization_factor,
              init_value=self._init_value_gain);
    dprintf(0,"solving for %d gain parms using %d of %d inteferometers\n",len(gain.gain),
      len(self._solvable_ifrs)/2,len(self.ifrs));
    dprint(1,"convergence target",gain.convergence_target,"of",gain.total_parms,"parms");
    dprint(1,"initial gain value is",self._init_value_gain.values()[0].flat[0] if isinstance(self._init_value_gain,dict) else
      self._init_value_gain);

    # init diffgains
    if num_diffgains:
      diffgains = [
        GainClass(expanded_datashape,dg_subtiling,self._solvable_ifrs,self.epsilon,self.convergence_quota,
          init_value=self._init_value_dg.get(i,self.init_value) )
        for i in range(num_diffgains) ];
      dg0 = diffgains[0];
      dprintf(0,"also solving for %dx%d differential gains\n",num_diffgains,len(dg0.gain));
      dprint(1,"convergence target for each is ",dg0.convergence_target,"of",dg0.total_parms,"parms");
      for i in range(num_diffgains):
        initval = self._init_value_dg.get(i,self.init_value);
        dprint(1,"initial gain value #%d is"%i,initval.items()[0] if isinstance(initval,dict) else initval);
    else:
      diffgains = [];

    def compute_chisq ():
      chisq = 0;
      nterms = 0;
      for pq in self._solvable_ifrs:
        for r in gain.residual(data,model,pq):
          if numpy.isscalar(r):
            chisq1 = (r*numpy.conj(r));
            nterms += 1;
          else:
            chisq1 = (r*numpy.conj(r)).sum();
            nterms += r.size;
#          print pq,"chi-sqare contribution is",chisq1,"r",r if numpy.isscalar(r) else r[5,0];
          chisq += chisq1;
      return chisq/nterms;

    # start major loop -- alternates over gains and diffgains
    for nmajor in range(self.max_major+1):
      # first iterate normal gains to convergence
      for niter in range(self.max_iter):
        # iterate over normal gains
        converged = gain.iterate(data,model,first_iter=not niter);
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
        data1 = dict([ (pq,matrix_sub(gain.correct(data,pq),model[pq])) for pq in self._solvable_ifrs ]);
#        for pq in list(self._solvable_ifrs)[:1]:
#          print [ (pq,[ d1 if is_null(d1) else abs(d1).max() for d1 in data1[pq] ]) for pq in self._solvable_ifrs ];
        # now loop over all diffgains and iterate each set once
        for i,dg in enumerate(diffgains):
          # add current estimate of corrupt(Mi) back into data1, and subtract from model
          for pq in self._solvable_ifrs:
            corr = dg.corrupt(dgmodel[i],pq,cache=True);
            mm = model[pq];
            for n,(d,m,c) in enumerate(zip(data1[pq],mm,corr)):
              d += c;
              mm[n] = -c if is_null(m) else m-c;
#          print data1['0','A'][0][0,0],model0['0','A'][0][0,0],model['0','A'][0][0,0],dgmodel[i]['0','A'][0][0,0];
#          print "d1",data1['0','A'][0],"m1",dgmodel[i]['0','A'][0];
          # iterate this diffgain solution
          for niter in range(self.diffgain_max_iter):
            if dg.iterate(data1,dgmodel[i],first_iter=not niter):
              break;
          dprint(2,"diffgain #%d converged after %d iterations"%(i,niter));
#          print dg.gain.keys();
#          print "dE(0)",dg.gain['0',0];
#          print "dgcorrupt",dg.corrupt(dgmodel[i],('0','A'),cache=True)[0];
          # add back to model, and subtract from data1 if needed
          for pq in self._solvable_ifrs:
            corr = dg.corrupt(dgmodel[i],pq,cache=True);
            for d,m,c in zip(data1[pq],model[pq],corr):
              m += c;
              if i<num_diffgains-1:
                d -= c;
#      # done iterating over diffgains, reset cached residuals and compute chisq once again
#      gain.reset_residuals();
#      chisq = compute_chisq();
#      dprint(1,"diffgains converge to chisq %g after %d major cycles"%(chisq,nmajor+1));

    # if we were solving for diffgains, then model is not completely up-to-date, since the non-solvable baselines
    # have been ignored. Fill them in here. Also, reset residuals
    if num_diffgains:
      for pq in set(self._ifrs)-set(self._solvable_ifrs):
        model[pq] = model0[pq];
        for i,dg in enumerate(diffgains):
          for m,c in zip(model[pq],dg.corrupt(dgmodel[i],pq,cache=True)):
            m += c;

    # remember init value for next tile
    if self.init_from_previous:
      self._init_value_gain = gain.get_last_timeslot();
      for i,dg in enumerate(diffgains):
        self._init_value_dg[i] = dg.get_last_timeslot();

    # update IFR gain solutions, if asked to
    if self.solve_ifr_gains:
      for pq in self._ifrs:
        dd = data[pq];
        mm = gain.corrupt(model,pq,cache=True);
        for num,(d,m) in enumerate(zip(dd,mm)):
          # work out update to ifr gains
          if numpy.isscalar(m):
            m = numpy.array(m);
          if numpy.isscalar(d):
            d = numpy.array(d);
          dh = numpy.conj(d);
          sri = self.ig_sum_reim[pq][num] = self.ig_sum_reim[pq][num] + (m*dh).sum();
          ssq = self.ig_sum_sq[pq][num]   = self.ig_sum_sq[pq][num] + (d*dh).sum();
          if ssq != 0:
            self.ifr_gain_update[pq][num] = sri/ssq;
#          if num == 0 and pq[0] == '0':
#           print m[0,0],d[0,0],dh[0,0]
#            print pq,(m*dh).sum(),(d*dh).sum(),sri/ssq;

    # work out result -- residual or corrected visibilities, depending on our state
    nvells = maxres = 0;
    for pq in self._ifrs:
      m = model.get(pq);
      r = gain.residual(data,model,pq);
      out = r if self.residuals else data[pq];
      if self.correct:
        out = gain.correct(out,pq,index=False);
      for n,x in enumerate(out):
        val = getattr(datares.vellsets[nvells],'value',None);
        if val is not None:
          try:
            val[...] = x[expanded_dataslice] if expanded_dataslice \
              and not is_null(x) else x;
          except:
            print x,getattr(x,'shape',None);
        # compute stats
        maxres = max(maxres,abs(r[n]).max() if not numpy.isscalar(r[n]) else abs(r[n]));
        nvells += 1;

    # if last domain, then write ifr gains to file
    if self.solve_ifr_gains and time1 >= numtime:
      # apply updates
      for pq in self._ifrs:
        self.ifr_gain[pq] = [ g*g1 for g,g1 in zip(self.ifr_gain.get(pq,[1,1,1,1]),self.ifr_gain_update[pq]) ];
      dprint(2,"IFR gain solutions update: ",", ".join(
            ["%s%s:%s%s %s"%(p,self.corr_names[i],q,self.corr_names[j],
            self.ifr_gain_update[(p,q),i,j])
            for (p,q),i,j in pqij_data[0:3]]));
      dprint(2,"IFR gain solutions: ",", ".join(
            ["%s%s:%s%s %s"%(p,self.corr_names[i],q,self.corr_names[j],
            self.ifr_gain[(p,q),i,j])
            for (p,q),i,j in pqij_data[0:3]]));
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

