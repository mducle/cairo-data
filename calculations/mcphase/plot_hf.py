import numpy as np
import matplotlib.pyplot as plt

dat110 = np.loadtxt('bfo_hf_mag_110.dat', delimiter=',')
dat001 = np.loadtxt('bfo_hf_mag_001.dat', delimiter=',')
cal001 = np.loadtxt('bfo_001_mcphas.fum')
cal110 = np.loadtxt('bfo_110_mcphas.fum')
sus001 = np.loadtxt('bfo_susc_001_mcphas.fum')
sus110 = np.loadtxt('bfo_susc_110_mcphas.fum')

cc = plt.rcParams['axes.prop_cycle'].by_key()['color']
fig, ax = plt.subplots()
ax.plot(dat110[:,0], dat110[:,1], ls='-', color=cc[0], label='Data H||[110]')
ax.plot(dat001[:,0], dat001[:,1], ls='-', color=cc[1], label='Data H||[001]')
# Scale by the measured magnetic moment from neutron diffraction
# Calculation is in uB/Fe, data is uB/formula unit so there is a factor of 4 (Bi2Fe4O9).
scalfac = 4*((3.52/5+3.73/5)/2)
ax.plot(cal110[:,1], cal110[:,9] * scalfac, ls=':', color=cc[0], label='Calc. H||[110]')
ax.plot(cal001[:,1], cal001[:,9] * scalfac, ls=':', color=cc[1], label='Calc. H||[001]')
ax.set_xlim(0, 65)
ax.set_xlabel('Applied Field (Tesla)')
ax.set_ylim(0, 1.3)
ax.set_ylabel('Magnetisation ($\mu_B$ / f.u.)')
ax.legend(loc='upper left')
fig.tight_layout()

ax2 = fig.add_axes([0.69, 0.22, 0.25, 0.25])
# Calculation is in uB/Fe at 0.1T==1000G. To get in emu/mol-Fe: mult: NA * muB * 1000 / field
scalfac = 5.5849397 * 4*((3.52/5+3.73/5)/2)
ax2.plot(sus110[:,0], sus110[:,9] * scalfac, ls=':', color=cc[0])
ax2.plot(sus001[:,0], sus001[:,9] * scalfac, ls=':', color=cc[1])
ax2.set_xlim(0, 1000)
ax2.set_xlabel('Temperature (K)')
ax2.set_ylim(0, 0.03)
ax2.set_ylabel('Susceptibility (emu / mol)')

plt.show(block=True)
