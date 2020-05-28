# The following line helps with future compatibility with Python 3
# print must now be used as a function, e.g print('Hello','World')
from __future__ import (absolute_import, division, print_function, unicode_literals)

import matplotlib
import scipy

import numpy as np
import mslice.cli as m
import matplotlib.pyplot as plt
import mslice.plotting.pyplot as mplt
from mslice.models.workspacemanager.workspace_provider import get_visible_workspace_names

# Sets paths
import os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
wd = parent_dir + '/processed/'

# Loads data in subsidiary file
import sys
sys.path.append(parent_dir + '/scripts')

import fit_bi2fe4o9_powder
if sys.version_info > (3,):
    if sys.version_info > (3, 4):
        from importlib import reload
    else:
        from imp import reload
reload(fit_bi2fe4o9_powder)

wsd = fit_bi2fe4o9_powder.load_data(wd)

cts, labs = fit_bi2fe4o9_powder.generate_cuts(wsd)
for ct in cts:
    m.SaveMatlab(ct, '{}/calculations/matlab/{}_cut.mat'.format(parent_dir, ct.name.split('.')[0]))
bkg_cuts = fit_bi2fe4o9_powder.get_backgrounds(wsd, cts)
bkg_cuts_d = {}
for ctb, ctd in zip(bkg_cuts, cts):
    bkg_cuts_d['{}_cut'.format(ctd.name.split('.')[0])] = ctb
scipy.io.savemat(parent_dir + '/calculations/matlab/bfo_cuts_bkg.mat', bkg_cuts_d)

