# -*- coding: utf-8 -*-
print "Caught Qt3 import";
raise RuntimeError,"""It appears a TDL script is trying to import PyQt3.
PyQt3 is incompatible with the current version of MeqTrees. Please update
your scripts.""";