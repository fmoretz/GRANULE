import numpy as np
import pandas as pd

dg = pd.read_excel("../input/input_data.xlsx", sheet_name='kinetic coefficients', header=None)

Ks = dg[5][3]   # kg/m3 
Kd = dg[5][4]   # h
B0 = dg[5][5]   # L CH4 /g COD added
Y  = 0.27        # gCH4/gS -- from stoichiometry

print('\nKinetic coefficients of reactor:')
print(f'Half saturation constant (Monod kinetics): {Ks} kg/m3')
print(f'Decay constant: {Kd} [h]')
print(f'B0: {B0} [LCH4/g COD added]')
print(f'Yield (from stoichiometry): {Y} [-]')
print('=======================\n')