import numpy as np
import pandas as pd

dg = pd.read_excel("../input/input_data.xlsx", sheet_name='influent input', header=None)

#Sinitial = dg[5][2]
#Xinitial = dg[5][3]
Eo = 0.01
Mo = 1e-6
Ro = dg[5][4]
Tw = dg[5][5]  + 273.15
To = dg[5][6] + 273.15
day = dg[5][7]

Sinitial = 0.77 #, 0.555, 0.298, 0.195, 0.092]
Xinitial = 0.77/0.1#, 0.555/0.1, 0.298/0.1, 0.195/0.1, 0.092/0.1]


So =Sinitial/Xinitial
Xo = Xinitial/Xinitial






print('\nInfluent Data of reactor:')
print(f'Substrate: {Sinitial} kg/m3')
print(f'Biomass: {Xinitial} kg/m3')
print(f'Radius: {Ro} mm')
print(f'Temperature Jacket: {Tw} K')
print(f'Temperature Influent: {To} K')
print(f'Simulation time: {day} days')
print('=======================\n')
