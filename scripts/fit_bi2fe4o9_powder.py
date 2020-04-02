# The following line helps with future compatibility with Python 3
# print must now be used as a function, e.g print('Hello','World')
from __future__ import (absolute_import, division, print_function, unicode_literals)

from mantid.simpleapi import mtd, EvaluateFunction
import numpy as np
import mslice.cli as m
from mslice.models.workspacemanager.workspace_provider import get_visible_workspace_names
from mslice.models.workspacemanager.workspace_algorithms import export_workspace_to_ads
from mslice.models.workspacemanager.workspace_provider import add_workspace as import_workspace_from_ads
from CrystalField import Function

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

def load_data(wd):
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
    for rr in [3]+list(range(5,8)):
        wsn = 'MAP3047%d_Ei300.00meV_Rings' % (rr)
        if wsn not in loadedws:
            m.Load(wd+wsn+'.nxspe')
    wsmaps1 = m.get_workspace_handle('MAP30475_Ei300.00meV_Rings')
    wsmaps2 = m.get_workspace_handle('MAP30476_Ei300.00meV_Rings')
    wsd['map_5K'] = (wsmaps1 + wsmaps2) * 100
    wsd['map_150K'] = m.get_workspace_handle('MAP30477_Ei300.00meV_Rings')
    wsd['map_300K'] = m.get_workspace_handle('MAP30473_Ei300.00meV_Rings')
    return wsd

def generate_cuts(wsd):
    cts = []
    labs = []
    for ei, estp, maxE in zip([25, 38, 62, 120, 180], [0.9, 0.6, 0.31, 0.19, 0.125], [0.7, 0.75, 0.75, 0.6, 0.5]):
        cts.append(m.Cut(wsd['mer_5K_Ei%d' % (ei)], 'DeltaE,{},{},{}'.format(-10, ei*maxE, estp), '|Q|,0,3,3'))
        labs.append('Merlin Ei=%d meV' % (ei))
        if cts[-1].name not in mtd:
            export_workspace_to_ads(cts[-1])
    cts.append(m.Cut(wsd['map_5K'], 'DeltaE,-10,120,1', '|Q|,0,3,3'))
    labs.append('MAPS Ei=300 meV')
    if cts[-1].name not in mtd:
        export_workspace_to_ads(cts[-1])
    return (cts, labs)

def generate_cuts_phonons(wsd):
    cts = []
    labs = []
    for ei, estp, maxE, qmin, qmax in zip([25, 38, 62, 120, 180], [0.9, 0.6, 0.31, 0.19, 0.125], [0.7, 0.75, 0.75, 0.6, 0.5],
            [3, 5, 7, 12, 15], [6, 8, 10, 15, 18]):
        cts.append(m.Cut(wsd['mer_5K_Ei%d' % (ei)], 'DeltaE,{},{},{}'.format(-10, ei*maxE, estp), '|Q|,{},{},3'.format(qmin,qmax)))
        labs.append('Merlin Ei=%d meV' % (ei))
        if cts[-1].name not in mtd:
            export_workspace_to_ads(cts[-1])
    cts.append(m.Cut(wsd['map_5K'], 'DeltaE,-10,120,1', '|Q|,9,12,3'))
    labs.append('MAPS Ei=300 meV')
    if cts[-1].name not in mtd:
        export_workspace_to_ads(cts[-1])
    return (cts, labs)

def _pv(pars):
    if not isinstance(pars[0], list):
        pars = [pars]
    funs = [Function(str('PseudoVoigt'), Intensity=pp[0], FWHM=pp[1], PeakCentre=pp[2], Mixing=pp[3]) for pp in pars]
    return ';'.join([str(fn.function) for fn in funs])
            
def get_backgrounds(wsd, data_cuts):
    phonon_cuts, labs = generate_cuts_phonons(wsd)
    elastic_pars = [[10441.119, 1.3369, 0.25035, 0.51308],
                    [17530.788, 1.6257, 0.57291, 0.29784],
                    [[21992.176, 2.7587, 0.90112, 0.16357], [1312.669, 5.6191, 4.8713, 0.0370]],
                    [[20369.103, 5.5513, 0.51717, 0.04549], [1876.243 * (1.5/13.5), 10.7003, 16.1717, 1.0000]],
                    [20532.998, 9.7188, -0.20442, 0.16401],
                    [67295.989, 15.7065, -3.29658, 0.84432]]
    elastic_str = [_pv(pp) for pp in elastic_pars]
    phon_qs = [4.5, 6.5, 8.5, 13.5, 16.5, 10.5]  # Average |Q| of phonon cuts
    eis = [25, 38, 62, 120, 180, 1]
    emins = [0.15, 0.15, 0.2, 0.07, 0.07, 15.0]
    cts = []
    for qs, phonon_ct, data_ct, fn, ei, emin in zip(phon_qs, phonon_cuts, data_cuts, elastic_str, eis, emins):
        scale_fac = 1.5 / qs
        ws_elastic = EvaluateFunction(fn, mtd[data_ct.name], OutputWorkspace='ws_elastic', IgnoreInvalidData=True)
        xx = ws_elastic.extractX()[0]
        yy = ws_elastic.extractY()[1]
        yy_phon = mtd[phonon_ct.name].extractY()[0]
        id_min = np.where(xx < (ei * emin))[0][-1]
        yy_phon[:id_min] = yy_phon[id_min]
        cts.append([xx, yy + (yy_phon * scale_fac)])
    return cts
