#!/bin/bash
# -*- coding: utf-8 -*-

TMPDIR=.scripter.$$.tmp
GLOBAL_LOGFILE=.scripter.log
LOGFILE=.scripter.$$.log

# if -t is specified on command-line, engage test mode -- no PID in files
for arg in $*; do
  if [ "$arg" == "-t" ]; then
    TMPDIR=.scripter.tmp
    LOGFILE=$GLOBAL_LOGFILE
  fi
done


# save scripts and funcs in subdirectory, lest they be modified while we run
if [ ! -d $TMPDIR ]; then
  mkdir $TMPDIR
fi
cp scripter.*.{conf,funcs} $TMPDIR


# load configs
for conf in $TMPDIR/scripter.*.conf; do
  echo "::: Loading configuration $conf"
  source $conf
done

# parse args
ms_names=""
ms_specs=""
processing_steps=""
args="$*"

while [ "$1" != "" ]; do
  # *.{ms,MS}, with an optional var=value suffix, is an MS name
  if echo $1 | egrep '^[^=]*\.(ms|MS)(:.*)?$' >/dev/null; then
    ms_specs="$ms_specs $1"
  # -f skips initial confirmation
  elif [ "$1" == "-f" ]; then
    unset CONFIRMATION_PROMPT
  # -i enables interactive mode
  elif [ "$1" == "-i" ]; then
    INTERACTIVE=1
    echo "::: Enabling interactive mode"
  # -d sets the default processing sequence
  elif [ "${1#-d}" != "$1" ]; then
    eval ps="\$DEFAULT_PROCESSING_SEQUENCE_${1#-d}"
    processing_steps="$processing_steps $ps"
  # -a sets the MS list to ALL_MS
  elif [ "$1" == "-a" ]; then
    ms_specs="$ALL_MS"
  # -t is test mode -- handled above
  elif [ "$1" == "-t" ]; then
    /bin/true;
  # +oper specifies a global operation, same with operations that start with _
  else
    processing_steps="$processing_steps $1"
  fi
  shift;
done

# if no steps specified, use default list
if [ "$processing_steps" == "" ]; then
  echo "No processing steps specified. Use -d to run default sequence:"
  echo "$DEFAULT_PROCESSING_SEQUENCE_"
  exit 1
fi

# make list of MS names (by stripping of :var=value from MS specs)
ms_names=""
for ms in $ms_specs; do
  ms_names="$ms_names ${ms%%:*}"
done

echo "::: MSs: $ms_specs"
echo "::: processing sequence: $processing_steps";

if [ "$CONFIRMATION_PROMPT" != "" ]; then
  echo -n "$CONFIRMATION_PROMPT"
  read
fi

echo `date "+%D %T"` $BASHPID: $0 $args >>$LOGFILE
# duplicate our command-line arguments to global logfile, if needed
if [ "$LOGFILE" != "$GLOBAL_LOGFILE" ]; then
  echo `date "+%D %T"` $BASHPID: $0 $args >>$GLOBAL_LOGFILE
fi

# this displays a prompt, if interactive mode is enabled with -i
interactive_prompt ()
{
  if [ "$INTERACTIVE" != "" ]; then
    echo -n "Press Enter to continue, or enter 'go' to run non-interactively from now... "
    read answer
    if [ "$answer" == "go" ]; then
      echo "::: Disabling interactive mode"
      unset INTERACTIVE
    fi
  fi
}

# function to iterate over a series of steps
iterate_steps ()
{
  echo "::: Iterating over steps: $*";
  if [ "$*" == "" ]; then
    return
  fi
  for oper in $*; do
    # var=value argument: directly assign to local variable
    if echo $oper | egrep '^[[:alnum:]_]+=.*$' >/dev/null; then
      varname="${oper%%=*}"
      varvalue="${oper#*=}"
      echo "::: Changing variable: $varname = $varvalue"
      eval $varname='"$varvalue"'
      # if a step= is explicitly specified, do reset the step counter to this in the next per_ms call
      if [ "${oper#step=}" != "$oper" ]; then
        reset_step_counter="$step"
      fi
    # var+=value argument: append to local variable
    elif echo $oper | egrep '^[[:alnum:]_]+[+]=.*$' >/dev/null; then
      varname="${oper%%+=*}"
      varvalue="${oper#*+=}"
      echo "::: Appending to variable: $varname += $varvalue"
      eval $varname+='" $varvalue"'
      # if a step= is explicitly specified, do reset the step counter to this in the next per_ms call
      if [ "${oper#step=}" != "$oper" ]; then
        reset_step_counter="$step"
      fi
    else
      # does oper have arguments, as oper[args]?
      if [ "${oper#*[}" != "$oper" -a "${oper:${#oper}-1:1}" == "]" ]; then
        args=${oper#*[}
        args=${args%]}
        oper=${oper%%[*}
        if [ "${args#*,,}" != "$args" ]; then
          args=${args//,,/ }
        else
          args=${args//,/ }
        fi
      else
        args="";
      fi
      # reset MSNAME if not set
      if [ "$MSNAME" == "" ]; then
        MSNAME="$FULLMS"
      fi
      # run operation
      echo "::: Running $oper $args (step=$step)"
      interactive_prompt
      echo `date "+%D %T"` $BASHPID: step $oper $args >>$LOGFILE
      if ! eval $oper $args; then
        echo "::: Operation '$oper $args' (step $step) failed";
        exit 1;
      fi
    fi
  done
}

# load functions
for func in $TMPDIR/scripter.*.funcs; do
  echo "::: Loading function set $func"
  source $func
done

ddid=0
field=0
unset MSNAME
reset_step_counter=1

per_ms ()
{
  # now start outer loop over MSs
  for MSNAMESPEC in $ms_specs; do
    # reset step to 0 on each new MS, unless we happen to hit
    # a step=N assignment on the command line (see above)
    # (user can also override this with reset_step_counter=0 on the command line)
    if [ "$reset_step_counter" != "-" ]; then
      step="$reset_step_counter"
    fi
    # setup variables based on MS
    MSNAME="${MSNAMESPEC%%:*}"
    echo `date "+%D %T"` $BASHPID: MS $MSNAME >>$LOGFILE
    vars=${MSNAMESPEC#*:}
    if [ "$vars" != "$MSNAMESPEC" ]; then
      echo "::: Changing variables: $vars"
      eval ${vars//:/;}
    fi
    MS="ms_sel.msname=$MSNAME ms_sel.ddid_index=$ddid ms_sel.field_index=$field"
    CHANS="ms_sel.ms_channel_start=${CHAN0[$ddid]} ms_sel.ms_channel_end=${CHAN1[$ddid]}"
    FQSLICE=${CHAN0[$ddid]}~${CHAN1[$ddid]}:2
    FIELD="field=$field"
    msbase=`basename ${MSNAME} .MS`
    msbase="${msbase%.ms}"
    eval msbase=$FILENAME_PATTERN

    echo "::: MS $MSNAME ddid=$ddid field=$field steps:$*"
    alias
    interactive_prompt

    # load functions
    for func in $TMPDIR/scripter.*.funcs; do
      echo "::: Loading function set $func"
      source $func
    done

    iterate_steps $*
  done
  unset MSNAME
}

iterate_steps $processing_steps

echo `date "+%D %T"` $BASHPID: successful, exiting with status 0 >>$LOGFILE
if [ "$LOGFILE" != "$GLOBAL_LOGFILE" ]; then
  echo `date "+%D %T"` $BASHPID: successful, exiting with status 0 >>$GLOBAL_LOGFILE
fi
rm -fr $TMPDIR
exit 0;
