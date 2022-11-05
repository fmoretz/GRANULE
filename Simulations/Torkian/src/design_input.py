import numpy as np
import pandas as pd

df = pd.read_excel("../input/input_data.xlsx", sheet_name='design input', header=None)

V    = df[4][2]
Dc   = df[4][3]
H    = df[4][4]
H_ef = df[4][5]
Q    = df[4][6]
u    = df[4][7]

A = 3.1415 * Dc**2 / 4

Vef = A * H_ef
#
print('\nDesign Data of reactor:')
print(f'Volume: {V} m3')
print(f'Diameter: {Dc} m')
print(f'Height: {H} m')
print(f'Effective Height: {H_ef} m')
print(f'Flow Rate: {Q} m3/h')
print(f'Upflow Velocity: {u} m/h')
print('=======================\n')
