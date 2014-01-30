#!/bin/bash
# -*- coding: utf-8 -*-

#TMPDIR=.scripter.$$.tmp
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
#if [ ! -d $TMPDIR ]; then
#  mkdir $TMPDIR
#fi
#cp scripter.*.{conf,funcs} $TMPDIR

# init defaults
ddid=0
field=0

# load configs
for conf in scripter.*.conf; do
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
  if echo $1 | egrep '^[^=]*\.(ms|MS)/?(:.*)?$' >/dev/null; then
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
    _var=DEFAULT_PROCESSING_SEQUENCE_${1#-d}
    echo $_var
    eval ps="'${!_var}'"
    echo $ps
    processing_steps="$processing_steps $ps"
  # -t is test mode -- handled above
  elif [ "$1" == "-t" ]; then
    /bin/true;
  # else an operation
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
set_msnames ()
{
  MSNAMES=""
  for ms in $*; do
    MSNAMES="$MSNAMES ${ms%%:*}"
  done
  echo "::: MS names: $MSNAMES"
}

echo "::: MS specs: $ms_specs"
set_msnames $ms_specs
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

# Assigns templated variables: for each variable named VAR_Template,
# creates the variable VAR and assigns it the expansion of $VAR_Template
# A list of template names may be passed in -- otherwise it uses compgen
# to get a list of all variables ending with _Template
assign_templates ()
{
# do variable substitution
  for tt in ${*:-`compgen -v | grep _Template | sort`}; do
    varname=${tt%_Template}
    eval $varname="\"${!tt}\""
    echo "::::: Assigned $varname=${!varname}"
  done
}


# load functions
for func in scripter.*.funcs; do
  echo "::: Loading function set $func"
  source $func
done

unset MSNAME
reset_step_counter=1

# function to iterate over a series of steps
iterate_steps ()
{
  echo "::: Iterating over steps: $*";
  if [ "$*" == "" ]; then
    return
  fi
  for oper in $*; do
    # -a sets the MS list to ALL_MS
    if [ "$oper" == "-a" ]; then
      ms_specs="`eval echo $ALL_MS`"
      echo "::: Set MS list to $ms_specs"
      # make list of MS names (by stripping of :var=value from MS specs)
      set_msnames $ms_specs
    # var=value argument: directly assign to local variable
    elif echo $oper | egrep '^[[:alnum:]_]+=.*$' >/dev/null; then
      varname="${oper%%=*}"
      varvalue="${oper#*=}"
      eval $varname="\"$varvalue\""
      echo "::: Changed variable: $varname = ${!varname}"
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
      # assign templates
      assign_templates
      # reset MSNAME if not set
      if [ "$MSNAME" == "" ]; then
        MSNAME="$FULLMS"
      fi
      # remove trailing slash from MS name
      MSNAME=${MSNAME%/}
      # create DESTDIR if needed
      if [ "$DESTDIR" != "" -a ! -d "$DESTDIR" ]; then
        mkdir $DESTDIR
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
    # setup variables based on MS name spec (i.e. one of the form "msname:var1=foo:var2=bar")
    echo `date "+%D %T"` $BASHPID: MS $MSNAME >>$LOGFILE
    vars=${MSNAMESPEC#*:}
    if [ "$vars" != "$MSNAMESPEC" ]; then
      echo "::: Changing variables: $vars"
      eval ${vars//:/;}
    fi
    MSNAME="${MSNAMESPEC%%:*}"
    # assign $msbase parameter
    MSBASENAME=`basename ${MSNAME} .MS`
    MSBASENAME="${MSBASENAME%.ms}"
    # assign templates
    assign_templates
    if [ "$DESTDIR" != "" -a ! -d "$DESTDIR" ]; then
      mkdir $DESTDIR
    fi
    # go
    echo "::: MS $MSNAME ddid=$ddid field=$field steps:$*"
    alias
    interactive_prompt

#    # load functions
#    for func in $TMPDIR/scripter.*.funcs; do
#      echo "::: Loading function set $func"
#      source $func
#    done

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
