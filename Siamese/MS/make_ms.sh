#!/bin/bash
print_templates ()
{
  echo -n "Available templates are: "
  for t in `ls *_makems.cfg`; do
    echo -n "${t%_makems.cfg} "
  done
  echo " "
}

if [ "$1" == "" ]; then
  echo "Usage $0 <template>"
  print_templates
  exit 1
fi
configfile="$1_makems.cfg"

if [ ! -f $configfile ]; then
  echo "No such template: $1"
  print_templates
  exit 1
fi 

rm -f makems.cfg
ln -s $configfile makems.cfg
makems
# find out MS name from config file
msname="`grep MSName= $configfile | cut -d = -f 2`"
if [ "$msname" != "" ]; then
  rm -fr $msname
  mv ${msname}_p1 $msname
  chmod -R a+rX $msname
  echo "Inserting MODEL_DATA and CORRECTED_DATA columns";
  lwimager ms=$msname data=CORRECTED_DATA mode=channel weight=natural npix=1 &>/dev/null
  echo "Created $msname"
else
  echo "Couldn't find MS name in config file $configfile"
fi
