import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import NullFormatter,MultipleLocator, FormatStrFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as ticker

plt.figure(num=None, dpi=80, facecolor='w', edgecolor='k')
matplotlib.rcParams.update({'font.size': 12, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
matplotlib.rcParams['font.sans-serif'] = "Arial"




_m  = (0.013*(To)-0.129)/24  				# h-1
