# -*- coding: utf-8 -*-

# This file is a guard against import of PyQt3 modules.
# Since we have switched to using PyQt4, any scripts that import PyQt3 (moudle 'qt')
# produce an unseemly segfault. To avoid the segfault, this qt.py module should
# intercept such import statements.

print("Caught Qt3 import");
raise RuntimeError("""It appears a TDL script is trying to import PyQt3.
PyQt3 is incompatible with the current version of MeqTrees. Please update
your scripts.""");