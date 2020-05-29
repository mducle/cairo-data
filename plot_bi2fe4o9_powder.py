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
parent_dir, _ = os.path.split(__file__)
wd = parent_dir + '/processed/'

# MERLIN data
# 40873	Bi2Fe4O9 120 meV 600Hz Gd Cooling	Mon Jul 9 18:42:14 2018	"40m 32s"	"115.4842"	Le	1820598	"23.592"	
# 40874	Bi2Fe4O9 120 meV 600Hz Gd 5K	Mon Jul 9 19:31:04 2018	"4h 40m 59s"	"800.0784"	Le	1820598	"164.077"	
# 40875	Bi2Fe4O9 180 meV 600Hz Gd 5K	Tue Jul 10 00:15:12 2018	"4h 41m 4s"	"800.05"	Le	1820598	"127.497"	
# 40876	Bi2Fe4O9 120 meV 600Hz Gd 100K	Tue Jul 10 04:59:10 2018	"1h 10m 30s"	"200.072"	Le	1820598	"40.901"	
# 40877	Bi2Fe4O9 120 meV 600Hz Gd 200K	Tue Jul 10 06:16:13 2018	"1h 14m 48s"	"200.0701"	Le	1820598	"40.2558"	
# 40878	Bi2Fe4O9 120 meV 600Hz Gd 300K	Tue Jul 10 07:54:17 2018	"16m 1s"	"45.703"	Le	1820598	"9.1442"	

# MAPS data
# 30473	Bi2Fe4O9 300meV 12S T=300K 25x45	Thu Sep 20 20:09:10 2018	"32m 29s"	"94.0521"	Le	1820598	"136.389"	
# 30474	Bi2Fe4O9 300meV 12S Cooling 25x45	Thu Sep 20 20:48:03 2018	"2h 12m 24s"	"382.9842"	Le	1820598	"556.936"	
# 30475	Bi2Fe4O9 300meV 12S T=5K 25x45	Thu Sep 20 23:00:42 2018	"2h 52m 47s"	"500.0803"	Le	1820598	"727.313"	
# 30476	Bi2Fe4O9 300meV 12S T=5K 25x45	Fri Sep 21 01:53:37 2018	"2h 52m 46s"	"500.076"	Le	1820598	"727.694"	
# 30477	Bi2Fe4O9 300meV 12S T=150K 25x45	Fri Sep 21 08:26:42 2018	"2h 32m 8s"	"429.0359"	Le	1820598	"624.146"	

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

###################################################################################
# Plot slice data and calculations
###################################################################################
sls = []
for ei, estp in zip([180, 120, 62, 38, 25], [0.9, 0.6, 0.31, 0.19, 0.125]):
    sls.append(m.Slice(wsd['mer_5K_Ei%d' % (ei)], '|Q|', 'DeltaE,{},{},{}'.format(-10, ei*0.7, estp)))

cmp = matplotlib.cm.get_cmap('cividis').reversed()
fig1, (ax1, ax2) = plt.subplots(ncols=2, subplot_kw={'projection': 'mslice'})
for sl in sls:
    ax1.pcolormesh(sl, vmin=0., vmax=200., cmap=cmp)
ax1.set_ylim(0, 120)
ax1.set_xlim(0, 5)
ax2.pcolormesh(m.Slice(wsd['map_5K'], '|Q|,0.5,5,0.075', 'DeltaE,0, 120, 1'), vmin=0, vmax=200, cmap=cmp)
fig1.show()

spinw_calc = scipy.io.loadmat(parent_dir + '/calculations/bfo_powspec.mat')
fig2, (ax1, ax2, ax3) = plt.subplots(ncols=3, subplot_kw={'projection': 'mslice'})
for sl in sls:
    ax1.pcolormesh(sl, vmin=0., vmax=200., cmap=cmp)
ax1.set_ylim(0, 100)
ax1.set_xlim(0, 5)
ax1.set_title('MERLIN data')
ax2.pcolormesh(m.Slice(wsd['map_5K'], '|Q|,0.5,5,0.075', 'DeltaE,0, 120, 1'), vmin=0, vmax=200, cmap=cmp)
ax2.set_ylim((0,100))
ax2.set_xlim(0, 5)
ax2.set_ylabel('')
ax2.set_title('MAPS data')
hkl = [op(spinw_calc['hklA']) for op in [np.min, np.max, lambda x: np.mean(np.diff(x))/2., np.shape]]
hkl = np.linspace(hkl[0]-hkl[2], hkl[1]+hkl[2], hkl[3][1]+1)
ax3.pcolormesh(hkl, np.squeeze(spinw_calc['Evect']), spinw_calc['swConv'], vmin=0, vmax=0.1, cmap=cmp)
ax3.set_xlim(0, 5)
ax3.set_xlabel('$|Q| (\mathrm{\AA}^{-1})$')
ax3.set_title('SpinW calc.')
fig2.show()

###################################################################################
# Plots cut data and calculations
###################################################################################

cts, labs = fit_bi2fe4o9_powder.generate_cuts(wsd)

fig = mplt.figure()
data_cuts = []
ax = fig.add_subplot(111, projection="mslice")
for ct, lab in zip(cts, labs):
    ax.errorbar(ct, label=lab, fmt='', marker='o', ls='-')
    data_cuts.append([ct.get_coordinates(), ct.get_signal(), ct.get_error(), lab])
ax.set_ylim(bottom=0.0, top=500.0)
ax.set_xlim(left=-20, right=130)
ax.set_waterfall(True, x_offset=0.0, y_offset=50.0)
fig.show()
m.KeepFigure()

bkg_cuts = fit_bi2fe4o9_powder.get_backgrounds(wsd, cts)
spinw_cuts = scipy.io.loadmat(parent_dir + '/calculations/bfo_powcut.mat')
offset = 100
mrk = ['o', '^', 's', 'P', 'X', 'D']
fig, ax = plt.subplots()
for (ii, dat_cut) in enumerate(data_cuts):
    ax.errorbar(dat_cut[0][list(dat_cut[0].keys())[0]], dat_cut[1] + ii*offset, dat_cut[2], label=dat_cut[3], marker=mrk[ii], ls='')
scale_fac = [1000]*5+[2500]
cc = plt.rcParams['axes.prop_cycle'].by_key()['color']
for (ii, sw_cut) in enumerate(spinw_cuts['bfo_powcut'][0,:]):
    xx = sw_cut[0,0].T
    yy = (sw_cut[0,1]*scale_fac[ii]).T
    bkg_x = bkg_cuts[ii][0]
    bkg_y = np.interp(xx, (bkg_x[:-1] + bkg_x[1:])/2., bkg_cuts[ii][1])
    ax.plot(xx, yy + bkg_y + ii*offset, ls='-', color=cc[ii])
ax.set_ylim(bottom=0.0, top=700.0)
ax.set_xlim(left=0, right=130)
ax.set_xlabel('Energy Transfer (meV)')
ax.set_ylabel('Intensity (arb. units)')
ax.legend(loc='lower right')
fig.show()

fig = mplt.figure()
ax = fig.add_subplot(111, projection="mslice")
ct5k = m.Cut(wsd['map_5K'], 'DeltaE,-10,120,1', '|Q|,0,3,3')
ct5k = ct5k / 200
ax.errorbar(ct5k, label='5K', marker='o', ls='-')
ax.errorbar(m.Cut(wsd['map_150K'], 'DeltaE,-10,120,1', '|Q|,0,3,3'), label='150K', marker='s', ls='-')
ax.errorbar(m.Cut(wsd['map_300K'], 'DeltaE,-10,120,1', '|Q|,0,3,3'), label='300K', marker='^', ls='-')
ax.set_ylim(bottom=0.0, top=2.0)
ax.set_xlim(left=-20, right=130)
fig.show()
m.KeepFigure()

