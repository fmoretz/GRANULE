import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import NullFormatter,MultipleLocator, FormatStrFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as ticker


plt.figure(num=None, dpi=80, facecolor='w', edgecolor='k')
matplotlib.rcParams.update({'font.size': 18, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
matplotlib.rcParams['font.sans-serif'] = "Arial"
R = np.linspace(1,5,5)

Sout = [0.36, 0.22, 0.12, 0.08, 0.04]
Sout_exp = [0.343, 0.213, 0.104, 0.089, 0.043]


fig, axs = plt.subplots(1, 1, figsize=(10, 10))
plt.subplots_adjust(
		left   = 0.125,
		bottom = 0.071,
		right  = 0.9,
		top    = 0.95, 
		wspace = 0.45,
		hspace = 0.75
		)
	
axs.plot(R, Sout, 'k', linewidth=2, marker='.', markersize=10, linestyle='solid', label='model')
axs.plot(R, Sout_exp, 'r^', markersize=10, label='experimental')
axs.set_xlabel('Test n°')
axs.set_ylabel('S$_{out}$ - [kg/m$^3$]')
axs.grid(color='k', alpha=0.5, linestyle='dashed', linewidth=0.8)
axs.set_title('Effluent substrate concentration', y=1.025, fontsize=18)
axs.legend(loc='best')

fig_pathway = '../output/'
fig.savefig(fig_pathway + 'Effluent substrate_Leitao.svg')



R2 = np.linspace(1,4,4)

Sout2 = [0.35, 0.42, 0.48, 0.58]
Sout_exp2 = [0.343, 0.416, 0.436, 0.469]

fig2, axs2 = plt.subplots(1, 1, figsize=(10, 10))
plt.subplots_adjust(
		left   = 0.125,
		bottom = 0.071,
		right  = 0.9,
		top    = 0.95, 
		wspace = 0.45,
		hspace = 0.75
		)
	
axs2.plot(R2, Sout2, 'k', linewidth=2, marker='.', markersize=10, linestyle='solid', label='model')
axs2.plot(R2, Sout_exp2, 'r^', markersize=10, label='experimental')
axs2.set_xlabel('Test n°')
axs2.set_ylabel('S$_{out}$ - [kg/m$^3$]')
axs2.grid(color='k', alpha=0.5, linestyle='dashed', linewidth=0.8)
axs2.set_title('Effluent substrate concentration', y=1.025, fontsize=18)
axs2.legend(loc='best')

fig_pathway = '../output/'
fig2.savefig(fig_pathway + 'Effluent substrate_Test2_Leitao.svg')


