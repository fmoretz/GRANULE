import numpy as np

with open('design_input.txt') as f:
    contents = f.readlines()

Volume = contents[2]
Diameter = contents[3]
Height = contents[4]
EffectiveHeight = contents[5]
FlowRate = contents[6]
UpflowVelocity = contents[7]

V = np.double(Volume[38:-1])
Dc = np.double(Diameter[38:-1])
H = np.double(Height[38:-1])
H_ef = np.double(EffectiveHeight[38:-1])
_Q = np.double(FlowRate[38:-1])
u = np.double(UpflowVelocity[38:-1])

A = 3.1415 * Dc**2 / 4

Vef = A * H_ef
#print(V)
#print(Dc)
#print(H)
#print(H_ef)
#print(_Q)
#print(u)

f.close()
