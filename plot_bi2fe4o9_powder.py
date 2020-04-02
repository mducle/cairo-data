# The following line helps with future compatibility with Python 3
# print must now be used as a function, e.g print('Hello','World')
from __future__ import (absolute_import, division, print_function, unicode_literals)

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

loadedws = get_visible_workspace_names()
wsd = {}
for (rr, tt) in ([4, 5], [6, 100], [7, 200], [8, 300]):
    for ei in [25, 38, 62, 120]:
        wsn = 'MER4087%d_Ei%.2fmeV_Rings' % (rr, ei)
        if wsn not in loadedws:
            wsd['mer_%dK_Ei%d' % (tt,ei)] = m.Load(wd+wsn+'.nxspe')
        else:
            wsd['mer_%dK_Ei%d' % (tt,ei)] = m.get_workspace_handle(wsn)
for ei in [30, 47, 82, 180]:
    wsn = 'MER40875_Ei%.2fmeV_Rings' % (ei)
    if wsn not in loadedws:
        wsd['mer_5K_Ei%d' % (ei)] = m.Load(wd+wsn+'.nxspe')
    else:
        wsd['mer_5K_Ei%d' % (ei)] = m.get_workspace_handle(wsn)
for rr in [3]+range(5,8):
    wsn = 'MAP3047%d_Ei300.00meV_Rings' % (rr)
    if wsn not in loadedws:
        m.Load(wd+wsn+'.nxspe')
wsmaps1 = m.get_workspace_handle('MAP30475_Ei300.00meV_Rings')
wsmaps2 = m.get_workspace_handle('MAP30476_Ei300.00meV_Rings')
wsd['map_5K'] = (wsmaps1 + wsmaps2) * 100
wsd['map_150K'] = m.get_workspace_handle('MAP30477_Ei300.00meV_Rings')
wsd['map_300K'] = m.get_workspace_handle('MAP30473_Ei300.00meV_Rings')

sls = []
for ei, estp in zip([180, 120, 62, 38, 25], [0.9, 0.6, 0.31, 0.19, 0.125]):
    sls.append(m.Slice(wsd['mer_5K_Ei%d' % (ei)], '|Q|', 'DeltaE,{},{},{}'.format(-10, ei*0.7, estp)))

fig1, (ax1, ax2) = plt.subplots(ncols=2, subplot_kw={'projection': 'mslice'})
for sl in sls:
    ax1.pcolormesh(sl, vmin=0., vmax=200.)
ax1.set_ylim(0, 120)
ax1.set_xlim(0, 5)
ax2.pcolormesh(m.Slice(wsd['map_5K'], '|Q|,0.5,5,0.075', 'DeltaE,0, 120, 1'), vmin=0, vmax=200)
plt.show()

cts = []
labs = []
for ei, estp, maxE in zip([25, 38, 62, 120, 180], [0.9, 0.6, 0.31, 0.19, 0.125], [0.7, 0.75, 0.75, 0.6, 0.5]):
    cts.append(m.Cut(wsd['mer_5K_Ei%d' % (ei)], 'DeltaE,{},{},{}'.format(-10, ei*maxE, estp), '|Q|,0,3,3'))
    labs.append('Merlin Ei=%d meV' % (ei))
cts.append(m.Cut(wsd['map_5K'], 'DeltaE,-10,120,1', '|Q|,0,3,3'))
labs.append('MAPS Ei=300 meV')

fig = mplt.figure()
ax = fig.add_subplot(111, projection="mslice")
for ct, lab in zip(cts, labs):
    ax.errorbar(ct, label=lab, marker='o', ls='-')
ax.set_ylim(bottom=0.0, top=500.0)
ax.set_xlim(left=-20, right=130)
ax.set_waterfall(True, x_offset=0.0, y_offset=50.0)
m.KeepFigure()
m.Show()

fig = mplt.figure()
ax = fig.add_subplot(111, projection="mslice")
wsd['map_5K'] = wsd['map_5K'] / 200
ax.errorbar(m.Cut(wsd['map_5K'], 'DeltaE,-10,120,1', '|Q|,0,3,3'), label='5K', marker='o', ls='-')
ax.errorbar(m.Cut(wsd['map_150K'], 'DeltaE,-10,120,1', '|Q|,0,3,3'), label='150K', marker='s', ls='-')
ax.errorbar(m.Cut(wsd['map_300K'], 'DeltaE,-10,120,1', '|Q|,0,3,3'), label='300K', marker='^', ls='-')
ax.set_ylim(bottom=0.0, top=2.0)
ax.set_xlim(left=-20, right=130)
m.Show()
