import numpy as np
import pandas as pd

dg = pd.read_excel("input\input_data_NugrohoSantoso.xlsx", sheet_name='influent input', header=None)

Sinitial = dg[5][2]
Xinitial = dg[5][3]
Eo = 0.01
Mo = 1e-6
Ro = dg[5][4]
Tw = dg[5][5]  + 273.15
To = dg[5][6] + 273.15

So = Sinitial/Xinitial
Xo = Xinitial/Xinitial

print('\nInfluent Data of reactor:')
print(f'Substrate: {So}')
print(f'Biomass: {Xo}')
print(f'Radius: {Ro}')
print(f'Temperature Jacket: {Tw}')
print(f'Temperature Influent: {To}')
print('=======================\n')
