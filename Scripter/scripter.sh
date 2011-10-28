# -*- coding: utf-8 -*-


# load configs
for conf in scripter.*.conf; do
  echo "::: Loading configuration $conf"
  source $conf

# parse args
ms_names=""
ms_specs=""
per_ms_steps=""
global_pre_steps=""
global_post_steps=""

while [ "$1" != "" ]; do
  # *.{ms,MS}, with an optional var=value suffix, is an MS name
  if echo $1 | egrep '.*\\.(ms|MS)(:.*)?$' >/dev/null; then
    ms_specs="$ms_specs $1"
    ms_names="$ms_names ${1%%:*}"
  # -f skips initial confirmation
  elif [ "$1" == "-f" ]; then
    unset CONFIRMATION_PROMPT
  # -d enables debugging
  elif [ "$1" == "-d" ]; then
    DEBUG_MODE=1
  # +oper specifies a global operation
  elif [ "${1#+}" != "$1" ]; then
    oper="${1#+}"
    if [ -z "$per_ms_steps" ]; then
      global_pre_steps="$global_pre_steps $oper"
    else
      global_post_steps="$global_post_steps $oper"
    fi
  else
    per_ms_steps="$per_ms_steps $oper"
  fi
  shift;
done

if [ -z "$per_ms_steps" ]; then
  per_ms_steps="$DEFAULT_PER_MS_STEPS"
fi

echo "::: MSs are $ms_specs"
echo "::: preprocessing steps are $global_pre_steps";
echo "::: per-MS steps are $per_ms_steps";
echo "::: post-processing steps are $global_post_steps";

if [ -n "$CONFIRMATION_PROMPT" ]; then
  echo -n "$CONFIRMATION_PROMPT"
  read
fi

# function to iterate over a series of steps
iterate_steps ()
{
  echo "::: Iterating over steps: $*";
  if [ -z "$*" ]; then
    return
  fi
  # loop over global operations
  for oper in $*; do
    # var=value argument: directly assign to local variable
    if echo $1 | egrep '^[[:alnum:]]+=.*$' >/dev/null; then
      echo "::: Changing variable: $1"
      eval $1
    else
      oper_args=${oper//:/ }
      echo "::: Running $oper_args"
      if ! eval $oper_args; then
        echo "::: Operation '$oper_args' (step $step) failed";
        exit 1;
      fi
    fi
  done
}

iterate_steps $global_pre_steps


ddid=0
field=0
# now start outer loop over MSs
for MSNAMESPEC in $ms_specs; do
  # setup variables based on MS
  MSNAME="${MSNAMESPEC%%:*}"
  vars=${MSNAMESPEC#*:}
  eval ${vars//:/;}
  MS="ms_sel.msname=$MSNAME ms_sel.ddid_index=$ddid ms_sel.field_index=$field"
  CHANS="ms_sel.ms_channel_start=${CHAN0[$ddid]} ms_sel.ms_channel_end=${CHAN1[$ddid]}"
  FQSLICE=${CHAN0[$ddid]}~${CHAN1[$ddid]}:2
  FIELD="field=$field"
  msbase=`basename ${MSNAME} .MS`
  eval msbase=$FILENAME_PATTERN
  alias imager="$RUNIMAGER ms=$MSNAME field=$field spwid=$ddid"
  alias plot-ms="$PLOTMS $MSNAME -F $field -D $ddid -I '$IFRS' -L $FQSLICE"

  echo "::: MS $MSNAME ddid=$ddid field=$field"

  # load functions
  for func in scripter.*.funcs; do
    echo "::: Loading function set $func"
    source $func


  iterate_steps $per_ms_steps
done

iterate_steps $global_post_steps

exit 0;
