from Timba.TDL import *
from Timba.Meq import meq
import math

import Meow
import Meow.StdTrees
from Meow import ParmGroup, Bookmarks


Settings.forest_state.cache_policy = 20
# MS options first
mssel = Meow.MSUtils.MSSelector(has_input=True, tile_sizes=[1, 10, 50], flags=False)
# MS compile-time options
TDLCompileOptions(*mssel.compile_options())
# MS run-time options
TDLRuntimeMenu("MS/data selection options", *mssel.runtime_options())

# output mode menu
# CORRELATIONS = [ "XX","XY","YX","YY" ];
# CORR_INDICES = dict([(corr,index) for index,corr in enumerate(CORRELATIONS)]);
ALL_CORRS = "XX XY YX YY"
DIAG_CORRS = "XX YY"
CROSS_CORRS = "YX XY"
SINGLE_CORR = "XX"

# UVW (new) this line seems to be critical to include the default options for MS
TDLCompileOptions(*Meow.IfrArray.compile_options())


run_model = TDLCompileOption('run_option', "What do we want to do", ['calibrate', 'simulate'])
#(new)
# run_model.when_changed(lambda mode:mssel.enable_input_column(mode='simulate'))

calibrate_options = TDLCompileMenu("Existing UV data Options",
                                   TDLOption('do_solve', "Calibrate", True),
                                   TDLOption('cal_corr', "Use correlations",
                                             [ALL_CORRS, DIAG_CORRS, CROSS_CORRS, SINGLE_CORR]),
                                   TDLOption('do_subtract', "Subtract sky model and generate residuals", True),
                                   TDLMenu("Correct the data or residuals",
                                           TDLMenu("include sky-Jones correction",
                                                   TDLOption('src_name', "name of the source (if None, first source in list is used)", [None, "S+0+0"], more=str),
                                                   toggle='do_correct_sky', open=True),
                                           toggle='do_correct', open=True))

simulate_options = TDLCompileMenu("Simulation Options",
                                  # noise option
                                  TDLOption('do_add', "Add sky model to existing data to increase number of sources", False),
                                  TDLOption("noise_stddev", "Add noise, Jy", [None, 1e-6, 1e-3], more=float))


do_not_simulate = (run_option == 'calibrate')


def change_run_option(model):
    global do_not_simulate
    simulate_options.show(model == 'simulate')
    calibrate_options.show(model == 'calibrate')
    do_not_simulate = (model == 'calibrate')

run_model.when_changed(change_run_option)


# now load optional modules for the ME maker
from Meow import MeqMaker
meqmaker = MeqMaker.MeqMaker(solvable=do_solve and run_option == 'calibrate')

# specify available sky models
# these will show up in the menu automatically
from Siamese.OMS import fitsimage_sky
from Lions import gridded_sky
import Meow.LSM
lsm = Meow.LSM.MeowLSM(include_options=False)
meqmaker.add_sky_models([lsm, fitsimage_sky, gridded_sky])

from Lions import ZJones
meqmaker.add_sky_jones('Z', 'iono', [ZJones.ZJones()])
from Calico.OMS import solvable_sky_jones

meqmaker.add_sky_jones('GD', 'directional receiver gains/phases',
                       [solvable_sky_jones.DiagAmplPhase("GD"),
                        solvable_sky_jones.FullRealImag("GD")])

from Calico.OMS import solvable_jones
meqmaker.add_uv_jones('G', 'receiver gains/phases',
                      [solvable_jones.DiagAmplPhase("G"),
                       solvable_jones.FullRealImag("G")])

meqmaker.add_uv_jones('B', 'bandpass',
                      [solvable_jones.DiagAmplPhase("B")])

# simulate G Jones
from Siamese.OMS import oms_gain_models
meqmaker.add_uv_jones('G_sim', 'simulate receiver gains/phases',
                      [oms_gain_models])

# very important -- insert meqmaker's options properly
TDLCompileOptions(*meqmaker.compile_options())


def _define_forest(ns):
    # make pynodes, xyzcomponent for sources
    ANTENNAS = mssel.get_antenna_set(range(1, 15))
    array = Meow.IfrArray(ns, ANTENNAS, mirror_uvw=False)
    observation = Meow.Observation(ns)
    Meow.Context.set(array, observation)
    # make a predict tree using the MeqMaker
    if do_solve or do_subtract or not do_not_simulate:
        outputs = predict = meqmaker.make_tree(ns)

    # make a list of selected corrs
    selected_corrs = cal_corr.split(" ")

    # make spigot nodes
    if not do_not_simulate and do_add:
        spigots = spigots0 = outputs = array.spigots()
        sums = ns.sums
        for p, q in array.ifrs():
            sums(p, q) << spigots(p, q) + predict(p, q)
        outputs = sums

    # make spigot nodes
    if do_not_simulate:
        spigots = spigots0 = outputs = array.spigots()
        # make nodes to compute residuals

        # make nodes to compute residuals
        if do_subtract:
            residuals = ns.residuals
            for p, q in array.ifrs():
                residuals(p, q) << spigots(p, q) - predict(p, q)
            outputs = residuals

        # and now we may need to correct the outputs
        if do_correct:
            if do_correct_sky:
                if src_name:
                    sky_correct = src_name
                else:
                    srcs = meqmaker.get_source_list(ns)
                    sky_correct = srcs and srcs[0]

            else:
                sky_correct = None
            outputs = meqmaker.correct_uv_data(ns, outputs, sky_correct=sky_correct)

        # make solve trees
        if do_solve:
            # extract selected correlations
            if cal_corr != ALL_CORRS:
                index = [CORR_INDICES[c] for c in selected_corrs]
                for p, q in array.ifrs():
                    ns.sel_predict(p, q) << Meq.Selector(predict(p, q), index=index, multi=True)
                    ns.sel_spigot(p, q) << Meq.Selector(spigots(p, q), index=index, multi=True)
                spigots = ns.sel_spigot
                predict = ns.sel_predict
            model = predict
            observed = spigots

            # make a solve tree
            solve_tree = Meow.StdTrees.SolveTree(ns, model)
            # the output of the sequencer is either the residuals or the spigots,
            # according to what has been set above
            outputs = solve_tree.sequencers(inputs=observed, outputs=outputs)

    # throw in a bit of noise
    if not do_not_simulate and noise_stddev:
        # make two complex noise terms per station (x/y)
        noisedef = Meq.GaussNoise(stddev=noise_stddev)
        noise_x = ns.sta_noise('x')
        noise_y = ns.sta_noise('y')
        for p in array.stations():
            noise_x(p) << Meq.ToComplex(noisedef, noisedef)
            noise_y(p) << Meq.ToComplex(noisedef, noisedef)
        # now combine them into per-baseline noise matrices
        for p, q in array.ifrs():
            noise = ns.noise(p, q) << Meq.Matrix22(
                noise_x(p) + noise_x(q), noise_x(p) + noise_y(q),
                noise_y(p) + noise_x(q), noise_y(p) + noise_y(q)
            )
            ns.noisy_predict(p, q) << outputs(p, q) + noise
        outputs = ns.noisy_predict
    # make sinks and vdm.
    # The list of inspectors comes in handy here
    Meow.StdTrees.make_sinks(ns, outputs, spigots=None, post=meqmaker.get_inspectors())

    if not do_not_simulate:
        # add simulate job
        TDLRuntimeJob(job_simulate, "Simulate")
    if do_not_simulate and not do_solve:
        # add subtract or correct job
        TDLRuntimeJob(job_subtract, "Subtract or Correct the data")

    if do_not_simulate and do_solve:
        pg_iono = ParmGroup.ParmGroup("Z_iono",
                                      outputs.search(tags="solvable Z"),
                                      table_name="iono.mep", bookmark=4)
        ParmGroup.SolveJob("cal_iono", "Calibrate Ionosphere parameters ", pg_iono)

    # very important -- insert meqmaker's runtime options properly
    # this should come last, since runtime options may be built up during compilation.
    # TDLRuntimeOptions(*meqmaker.runtime_options(nest=False));
    # and insert all solvejobs
    TDLRuntimeOptions(*ParmGroup.get_solvejob_options())
    # finally, setup imaging options
    imsel = mssel.imaging_selector(npix=512, arcmin=meqmaker.estimate_image_size())
    TDLRuntimeMenu("Imaging options", *imsel.option_list())


def job_simulate(mqs, parent, wait=False):
    mqs.execute('VisDataMux', mssel.create_io_request(), wait=wait)


def job_subtract(mqs, parent, wait=False):
    mqs.execute('VisDataMux', mssel.create_io_request(), wait=wait)

# def _tdl_job_1_simulate_MS (mqs,parent,wait=False):
#  mqs.execute('VisDataMux',mssel.create_io_request(),wait=wait);


if __name__ == '__main__':
    # You can run the script in headless / batch mode from the
    # command line by saying something like
    #  python demo_script.py -run
    #
    # If you want to keep track of what's happening, use
    #
    #  python demo_script.py -run -dmeqserver=3
    #
    # This dumps various messages to stdout, which let you keep
    # track of how the script is progressing.
    #
    # Here's the code required to handle the '-run' flag
    if '-run' in sys.argv:
        from Timba.Apps import meqserver
        from Timba.TDL import Compile

        # you may need the following line for more complicated scripts
        # that use TDL options
        # from Timba.TDL import TDLOptions

        # this starts a kernel.
        mqs = meqserver.default_mqs(wait_init=10)

        # more complicated scripts might want to invoke TDLOptions here ...
        # e.g. the next line of (commented out) python code loads a tdl.conf file.
        # Note that it may be better to use a separate config file, rather
        # than the default .tdl.conf that the browser creates.
        TDLOptions.config.read(".tdl.conf")

        # Now compile a script as a TDL module. Any errors will be thrown as
        # an exception, so this always returns successfully. We pass in
        # __file__ so as to compile ourselves.
        (mod, ns, msg) = Compile.compile_file(mqs, __file__)

        # This next call runs the _test_forest job.
        # Note that wait should be set to True for batch processing
        # so that in _test_forest the request is executed with wait=True.
        # This makes sure that commands are executed in order.
        mod.job_simulate(mqs, None, wait=True)

    else:
        Timba.TDL._dbg.set_verbose(5)
        ns = NodeScope()
        _define_forest(ns)
        # resolves nodes
        ns.Resolve()
