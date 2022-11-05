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
Sout = [0.36, ]
Sout_exp= [0.37, ]


fig, axs = plt.subplots(1, 1, figsize=(10, 10))
plt.subplots_adjust(
		left   = 0.125,
		bottom = 0.071,
		right  = 0.9,
		top    = 0.95, 
		wspace = 0.45,
		hspace = 0.75
		)
	
axs.plot(R, E, 'k', linewidth=2.5, marker='.', markersize=10, linestyle='solid', label='model')
axs.errorbar(R, E_exp, yerr=E_exp_err, marker='^', color='r',linewidth=0.01, ecolor='r', elinewidth= 1.5, markersize=10,label='experimental')
axs.set_xlabel('Test n°')
axs.set_ylabel('Efficiency - [%]')
axs.grid(color='k', alpha=0.5, linestyle='dashed', linewidth=0.8)
axs.set_title('Efficiency of substrate removal', y=1.025, fontsize=18)
axs.legend(loc='best')

fig_pathway = '../output/'
fig.savefig(fig_pathway + 'Efficiency of removal_Torkian.svg')



