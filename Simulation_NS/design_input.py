import numpy as np
import pandas as pd

df = pd.read_excel("input\input_data_NugrohoSantoso.xlsx", sheet_name='design input', header=None)
print(df[4][2])

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
print(f'Volume: {V}')
print(f'Diameter: {Dc}')
print(f'Height: {H}')
print(f'Effective Height: {H_ef}')
print(f'Flow Rate: {Q}')
print(f'Upflow Velocity: {u}')
print('=======================\n')
