import numpy as np

with open('influent_input.txt') as g:
    contents = g.readlines()

Substrate = contents[2]
Biomass = contents[3]
InactiveBiomass = contents[4]
Methane = contents[5]
Radius = contents[6]
TemperatureJacket = contents[7]
TemperatureInfluent = contents[8]


So = np.double(Substrate[47:-1])
Xo = np.double(Biomass[47:-1])
Eo = np.double(InactiveBiomass[47:-1])
Mo = np.double(Methane[47:-1])
Ro = np.double(Radius[47:-1])
Tw = np.double(TemperatureJacket[47:-1])
To = np.double(TemperatureInfluent[47:-1])

g.close()
