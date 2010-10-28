# -*- coding: utf-8 -*-


def substitute_pattern (filename_pattern,**substitutions):
  """Substitutes $(key) and $var instances in the given filename pattern, with values from the
  substitutions dict.
  """
  import re
  filename = filename_pattern;
  # loop over substitutions, longest to shortest
  # delimit substitutions with "\n" so that they don't interfere with subsequent word boundaries
  for key,value in sorted(substitutions.items(),lambda a,b:cmp(len(b[0]),len(a[0]))):
    filename = re.sub("\\$(%s\\b|\\(%s\\))"%(key,key),"\n"+value+"\n",filename);
  filename = re.sub("\\$\\$","$",filename);
  filename = filename.replace("\n","");
  return filename;
