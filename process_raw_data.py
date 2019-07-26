# The following line helps with future compatibility with Python 3
# print must now be used as a function, e.g print('Hello','World')
from __future__ import (absolute_import, division, print_function, unicode_literals)

# import mantid algorithms, numpy and matplotlib
from mantid import config
from mantid.simpleapi import SaveNexus, \
    DirectILLCollectData, DirectILLDiagnostics, DirectILLIntegrateVanadium
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# Imports reduction scripts
parent_dir, _ = os.path.split(__file__)
sys.path.append(parent_dir + '/scripts')
from MAPSReduction_Sample import MAPSReduction
from MERLINReduction_2018_1 import MERLINReduction
from MARIReduction_2019_1 import iliad_mari
import DirectILLReduction

# Sets paths
rawdatadir = parent_dir + '/raw'
datadir = parent_dir + '/processed'
config.appendDataSearchDir(rawdatadir)
old_save_dir = config['defaultsave.directory']
config['defaultsave.directory'] = datadir

# Reduce MERLIN data
ei = {120:[120,62,38,25], 180:[180,82, 47,30]}
runs = {120:[40874, 40876, 40877, 40878], 180:[40875]}
rd = MERLINReduction()
rd.def_advanced_properties()
rd.def_main_properties()
for efoc in [120, 180]:
    rd.reducer.prop_man.incident_energy = ei[efoc]
    rd.reducer.prop_man.sample_run = runs[efoc]
    rd.run_reduction()

# Reduce MAPS data
rd = MAPSReduction()
rd.def_advanced_properties()
rd.def_main_properties()
rd.run_reduction()

# Reduce MARI data
for runs, ei in [(25968, [66, 16]), ([26014, 26015], [160, 11])]:
    iliad_mari(runno=runs, ei=ei, wbvan=25779, monovan=None, sam_mass=0, sam_rmm=0, sum_runs=False,
        check_background=False, hard_mask_file='mari_mask_lab_spurion.msk')

# Reduce IN4 data
old_facility, old_inst = (config['default.facility'], config['default.instrument'])
config['default.facility'] = 'ILL'
config['default.instrument'] = 'IN4'
DirectILLCollectData(Run='082419.nxs', OutputWorkspace='vanadium', OutputEPPWorkspace='vanadium-epps', 
    OutputRawWorkspace='vanadium-raw')
DirectILLIntegrateVanadium(InputWorkspace='vanadium', OutputWorkspace='integrated', EPPWorkspace='vanadium-epps')
DirectILLDiagnostics(InputWorkspace='vanadium-raw', OutputWorkspace='diagnostics', EPPWorkspace='vanadium-epps')
runs = {
    'T1p7_1p22':'082349-082352',
    'T1p7_2p22':'082353-082362',
    'T1p7_1p11':'082363-082374+082379-082396',
    'T1p7_0p74':'082376-082377',
    'T80_2p22':'082398-082403',
    'T80_1p11':'082404-082408',
    'T40_2p22':'082409-082412',
    'T65_2p22':'082413-082417' }
wavelength = {'T1p7_1p22':1.22, 'T1p7_2p22':2.22, 'T1p7_1p11':1.11, 'T1p7_0p74':0.74, 
    'T80_2p22':2.22, 'T80_1p11':1.11, 'T40_2p22':2.22,  'T65_2p22':2.22}
for wsstr, numstr in runs.items():
    DirectILLCollectData(Run=rawdatadir+'/'+numstr+'.nxs', OutputWorkspace='sample', OutputIncidentEnergyWorkspace='Ei', OutputElasticChannelWorkspace='Elp')
    ei = convertUnits.input2energy(wavelength[wsstr], 'Wavelength (Angstroms)', theta=None, flightpath=None)
    DirectILLReduction.DirectILLReduction(InputWorkspace='sample', OutputWorkspace=wsstr+'_SofQW', 
        EnergyRebinning='{},{},{}'.format(-ei/2, ei/250, ei*0.9), GroupingAngleStep=0.6,
        IntegratedVanadiumWorkspace='integrated', DiagnosticsWorkspace='diagnostics',
        Cleanup='Cleanup OFF', SubalgorithmLogging='Logging ON')
    SaveNexus(wsstr+'_SofQW_grouped_detectors_', wsstr+'.nxs')

# Clean up
config['default.facility'] = old_facility
config['default.instrument'] = old_inst
config['defaultsave.directory'] = old_save_dir
