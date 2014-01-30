#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# This shows how to run a TDL script in a pipeline (aka batch-mode, aka headless mode)
#
if __name__ == '__main__':
  from Timba.Apps import meqserver
  from Timba.TDL import Compile
  from Timba.TDL import TDLOptions

  # This starts a meqserver. Note how we pass the "-mt 2" option to run two threads.
  # A proper pipeline script may want to get the value of "-mt" from its own arguments (sys.argv).
  print "Starting meqserver";
  mqs = meqserver.default_mqs(wait_init=10,extra=["-mt","2"]);

  # Once we're connected to a server, some cleanup is required before we can exit the script.
  # Since we want to perform this cleanup regardless of whether the script ran to completion
  # or was interrupted by an exception midway through, we use a try...finally block.
  try:

    TDLOptions.config.read("batch_sim_example.tdl.conf");

    script = "example-sim.py";
    print "========== Compiling batch job 1";
    mod,ns,msg = Compile.compile_file(mqs,script,config="batch job 1");
    print "========== Running batch job 1";
    mod._tdl_job_1_simulate_MS(mqs,None,wait=True);

    print "========== Compiling batch job 2";
    mod,ns,msg = Compile.compile_file(mqs,script,config="batch job 2");
    print "========== Running batch job 2";
    mod._tdl_job_1_simulate_MS(mqs,None,wait=True);

    print "========== Compiling batch job 2 with modified config";
    TDLOptions.init_options("batch job 2",save=False);
    TDLOptions.set_option("me.enable_G",False);
    mod,ns,msg = Compile.compile_file(mqs,script,config=None);
    print "========== Running batch job 2";
    mod._tdl_job_1_simulate_MS(mqs,None,wait=True);

  ### Cleanup time
  finally:
    print "Stopping meqserver";
    # this halts the meqserver
    meqserver.stop_default_mqs();
    # now we can exit
    print "Bye!";
