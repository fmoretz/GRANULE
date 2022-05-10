import numpy as np

with open('input/design_input.txt') as f:
    """
    Volume of the reactor [m3]
    Diameter of the reactor [m]
    Height of the reactor [m]
    Effective height of the reactor [m]
    influent flow in the reactor [m3/h]
    Upflow velocity [m/h]
    """
    contents = f.readlines()

Volume          = contents[2]
Diameter        = contents[3]
Height          = contents[4]
EffectiveHeight = contents[5]
FlowRate        = contents[6]
UpflowVelocity  = contents[7]

V    = np.double(Volume[38:-1])
Dc   = np.double(Diameter[38:-1])
H    = np.double(Height[38:-1])
H_ef = np.double(EffectiveHeight[38:-1])
Q    = np.double(FlowRate[38:-1])
u    = np.double(UpflowVelocity[38:-1])

A = 3.1415 * Dc**2 / 4

Vef = A * H_ef

print('\nDesign Data of reactor:')
print(f'Volume: {V}')
print(f'Diameter: {Dc}')
print(f'Height: {H}')
print(f'Effective Height: {H_ef}')
print(f'Flow Rate: {Q}')
print(f'Upflow Velocity: {u}')
print('=======================\n')

f.close()
