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
R = np.linspace(1,6,6)

M = [10, 23.5, 30.6, 36.3, 45.2, 61.8] #L/d
M_exp =[0.68*18, 0.67*33, 0.66*43, 56.8*0.6, 83*0.54, 116*0.48] 



fig, axs = plt.subplots(1, 1, figsize=(10, 10))
plt.subplots_adjust(
		left   = 0.125,
		bottom = 0.071,
		right  = 0.9,
		top    = 0.95, 
		wspace = 0.45,
		hspace = 0.75
		)
	
axs.plot(R, M, 'k', linewidth=2, marker='.', markersize=10, linestyle='solid', label='model')
axs.plot(R, M_exp, 'r^', markersize=10, label='experimental')
axs.set_xlabel('Test n°')
axs.set_ylabel('M - [L/d]')
axs.grid(color='k', alpha=0.5, linestyle='dashed', linewidth=0.8)
axs.set_title('Methane produced', y=1.025, fontsize=18)
axs.legend(loc='best')

fig_pathway = '../output/'
fig.savefig(fig_pathway + 'Methane produced_Shanmugam.svg')





