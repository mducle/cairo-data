# The following line helps with future compatibility with Python 3
# print must now be used as a function, e.g print('Hello','World')
from __future__ import (absolute_import, division, print_function, unicode_literals)

import mslice.cli as m
import matplotlib.pyplot as plt
import mslice.plotting.pyplot as mplt
from mslice.models.workspacemanager.workspace_provider import add_workspace, get_visible_workspace_names

# Sets paths
import os
parent_dir, _ = os.path.split(__file__)
wd = parent_dir + '/processed/'

# MARI data
# 25968	Bi4Fe5O13F - 66meV 350Hz Gd - T=5K	           Wed Jun 26 22:51:23 2019	"9h 15m 30s"	"1624.719"	NTC	0	"23.0677"	
# 26014	Bi4Fe5O13F - T=5K - Ei=160meV 300Hz Gd single	Fri Jul 5 20:19:05 2019	"6h 11m 13s"	"1000.0889"	 	0	"16.2536"	
# 26015	Bi4Fe5O13F - T=80K - Ei=160meV 300Hz Gd single	Sat Jul 6 02:31:20 2019	"5h 30m 14s"	"926.9114"	 	0	"14.9314"

loadedws = get_visible_workspace_names()
wsd = {}

# Loads MARI data
for sts in [[[16, 66], [25968], [5]], [[11, 160], [26014, 26015], [5, 80]]]:
    for rr, tt in zip(sts[1], sts[2]):
        for ei in sts[0]:
            wsn = 'MAR{}_Ei{:.2f}meV'.format(rr, ei)
            if wsn not in loadedws:
                wsd['mar_{}K_Ei{}'.format(tt,ei)] = m.Load(wd+wsn+'.nxspe')
            else:
                wsd['mar_{}K_Ei{}'.format(tt,ei)] = m.get_workspace_handle(wsn)

# Loads IN4 data
t2s = {'1p7':2, '40':40, '65':65, '80':80}
l2s = {'0p74':0.74, '1p11':1.11, '2p22':2.22}
l2e = {'0p74':150, '1p11':66, '2p22':16}
for sts in [['1p7', ['0p74', '1p11', '2p22']],
            ['40', ['2p22']],
            ['65', ['2p22']],
            ['80', ['1p11', '2p22']]]:
    for ll in sts[1]:
        wsn = 'T{}_{}'.format(sts[0], ll)
        if wsn not in loadedws:
            wsd['in4_{}K_Ei{}'.format(t2s[sts[0]], l2e[ll])] = m.Load(wd+wsn+'.nxs')
        else:
            wsd['in4_{}K_Ei{}'.format(t2s[sts[0]], l2e[ll])] = m.get_workspace_handle(wsn)
            
print(wsd.keys())

if True:#'mar_5K_Ei160_scaled' not in loadedws:
    wsd['mar_5K_Ei160_scaled'] = m.Scale(wsd['mar_5K_Ei160'], 0.0001, OutputWorkspace='MAR25968_Ei160.00meV_scaled')
else:
    wsd['mar_5K_Ei160_scaled'] = m.get_workspace_handle('MAR25968_Ei160.00meV_scaled')
if 'mar_5K_Ei66_scaled' not in loadedws:
    wsd['mar_5K_Ei66_scaled'] = m.Scale(wsd['mar_5K_Ei66'], 0.0002, OutputWorkspace='MAR25968_Ei66.00meV_scaled')
else:
    wsd['mar_5K_Ei66_scaled'] = m.get_workspace_handle('MAR25968_Ei66.00meV_scaled')
if 'mar_5K_Ei16_scaled' not in loadedws:
    wsd['mar_5K_Ei16_scaled'] = m.Scale(wsd['mar_5K_Ei16'], 0.00025, OutputWorkspace='MAR26014_Ei16.00meV_scaled')
else:
    wsd['mar_5K_Ei16_scaled'] = m.get_workspace_handle('MAR26014_Ei16.00meV_scaled')

###################################################################################

ct1 = []; lb1 = []; ct2 = []; lb2 = []
for wsn, inst in zip(['in4_2K_Ei16', 'mar_5K_Ei16_scaled'], ['IN4', 'MARI']):
    ct1.append(m.Cut(wsd[wsn], 'DeltaE,-5,13,0.1', '|Q|,1,1.5,0'))
    lb1.append('{} Ei=16 meV'.format(inst))
for wsn, inst in zip(['in4_2K_Ei66', 'mar_5K_Ei66_scaled'], ['IN4', 'MARI']):
    ct2.append(m.Cut(wsd[wsn], 'DeltaE,-15,60,0.5', '|Q|,1,4,0'))
    lb2.append('{} Ei=66 meV'.format(inst))
    
fig1, (ax1p1, ax1p2) = plt.subplots(nrows=2, subplot_kw={'projection': 'mslice'})
for ax, sts in [[ax1p1, [ct1, lb1]], [ax1p2, [ct2, lb2]]]:
    for ct, lb in zip(sts[0], sts[1]):
        ax.errorbar(ct, label=lb, marker='o', ls='-')
ax1p1.set_ylim(0, 0.03)
ax1p2.set_ylim(0, 0.005)
plt.show()

###################################################################################

sls = []
for ei, estp, qq in zip([160, 66], [0.2, 0.1], [[0.8, 12.6, 0.075], [0.5, 19, 0.1]]):
    sls.append(m.Slice(wsd['mar_5K_Ei{}_scaled'.format(ei)], '|Q|,{},{},{}'.format(qq[0],qq[1],qq[2]),
        'DeltaE,{},{},{}'.format(-40, ei*0.7, estp)))
sls.append(m.Slice(wsd['in4_2K_Ei16'], '|Q|', 'DeltaE,{},{},{}'.format(-10, 16*0.7, 0.025)))
fig2, (ax2r1, ax2r2) = plt.subplots(ncols=2, nrows=2, subplot_kw={'projection': 'mslice'})
ax2r1[0].pcolormesh(sls[2], vmin=0., vmax=0.025)
ax2r1[0].set_ylim(-2, 11)
ax2r1[0].set_xlim(0, 2)
for sl in sls:
    ax2r2[0].pcolormesh(sl, vmin=0., vmax=0.005)
ax2r2[0].set_ylim(-40, 80)
ax2r2[0].set_xlim(0, 5)

ct3 = m.Cut(wsd['in4_2K_Ei16'], 'DeltaE,-5,13,0.05', '|Q|,1,1.5,0')
ct3x = ct3.get_coordinates()
ct3x = ct3x[ct3x.keys()[0]]
ax2r1[1].errorbar(ct3.get_signal(), ct3x, xerr=ct3.get_error(), marker='.', ls='-')
ax2r1[1].set_xlim(0, 0.025)
ax2r1[1].set_ylim(-2, 11)

ct4 = [m.Cut(wsd['mar_5K_Ei160_scaled'], 'DeltaE,-40,120,1', '|Q|,1,4,0'),
       m.Cut(wsd['mar_5K_Ei66_scaled'], 'DeltaE,-40,60,0.5', '|Q|,1,4,0')]
lb4 = ['E$_i$=160 meV', 'E$_i$=66 meV']
for ct in ct4:
    ctx = ct.get_coordinates()
    ctx = ctx[ctx.keys()[0]]
    ax2r2[1].errorbar(ct.get_signal(), ctx, xerr=ct.get_error(), marker='.', ls='-')
ax2r2[1].set_xlim(0, 0.005)
ax2r2[1].set_ylim(-40, 80)

plt.show()

###################################################################################

sls2 = []
for ei, estp, qq in zip([160, 66], [0.2, 0.1], [[0.8, 12.6, 0.075], [0.5, 19, 0.1]]):
    sls2.append(m.Slice(wsd['mar_5K_Ei{}_scaled'.format(ei)], Axis2='|Q|,{},{},{}'.format(qq[0],qq[1],qq[2]),
        Axis1='DeltaE,{},{},{}'.format(-40, ei*0.7, estp)))
sls2.append(m.Slice(wsd['in4_2K_Ei16'], 'DeltaE,{},{},{}'.format(-10, 16*0.7, 0.025), '|Q|'))

fig3, (ax3r1, ax3r2) = plt.subplots(ncols=2, nrows=2, subplot_kw={'projection': 'mslice'})
ax3r1[0].pcolormesh(sls2[2], vmin=0., vmax=0.025)
ax3r1[0].set_xlim(-2, 11)
ax3r1[0].set_ylim(0, 2)
for sl in sls2:
    ax3r1[1].pcolormesh(sl, vmin=0., vmax=0.005)
ax3r1[1].set_xlim(-40, 100)
ax3r1[1].set_ylim(0, 5)

ax3r2[0].errorbar(ct3x, ct3.get_signal(), ct3.get_error(), marker='.', ls='-')
ax3r2[0].set_ylim(0, 0.025)
ax3r2[0].set_xlim(-2, 11)

for ct in ct4:
    ctx = ct.get_coordinates()
    ctx = ctx[ctx.keys()[0]]
    ax3r2[1].errorbar(ctx, ct.get_signal(), ct.get_error(), marker='.', ls='-')
ax3r2[1].set_ylim(0, 0.005)
ax3r2[1].set_xlim(-40, 100)

plt.show()

###################################################################################

fig4 = plt.figure()
sls3 = []
for id, tt in enumerate([2, 40, 65, 80]):
    sls3.append(m.Slice(wsd['in4_{}K_Ei16'.format(tt)], '|Q|', 'DeltaE,{},{},{}'.format(-2, 16*0.7, 0.025)))
    ax = fig4.add_subplot(4, 2, 2*id+1, projection='mslice')
    ax.pcolormesh(sls3[-1], vmin=0., vmax=0.025)
    ax.set_xlim(0, 2)
ax = fig4.add_subplot(1, 2, 2, projection='mslice')
ct3 = []
for id, tt in enumerate([2, 40, 65, 80]):
    ct3.append(m.Cut(wsd['in4_{}K_Ei16'.format(tt)], 'DeltaE,-5,13,0.1', '|Q|,1,1.5,0'))
    ax.errorbar(ct3[-1], marker='.', ls='-', label='{}K'.format(tt))
ax.set_ylim(0, 0.03)

plt.show()