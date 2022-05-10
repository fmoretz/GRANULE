import numpy as np

with open('input/influent_input.txt') as g:
    """
    Inlet substrate concentration [kg/m3]
    Inlet biomass concentration [kg/m3]
    Initial radius of the granules (average) [mm]:
    Temperature of water jacket [°C]
    Influent temperature [°C] 
    """
    contents = g.readlines()

Substrate           = contents[2]
Biomass             = contents[3]
Radius              = contents[4]
TemperatureJacket   = contents[5]
TemperatureInfluent = contents[6]


So = np.double(Substrate[47:-1])
Xo = np.double(Biomass[47:-1])
Eo = 0.01
Mo = 1e-6
Ro = np.double(Radius[47:-1])
Tw = np.double(TemperatureJacket[47:-1]) + 273.15
To = np.double(TemperatureInfluent[47:-1]) + 273.15

print('\nInfluent Data of reactor:')
print(f'Substrate: {So}')
print(f'Biomass: {Xo}')
print(f'Radius: {Ro}')
print(f'Temperature Jacket: {Tw}')
print(f'Temperature Influent: {To}')
print('=======================\n')

g.close()
