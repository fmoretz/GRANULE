import numpy as np 

def Diameter(dp, rho_biomass, rho_H2O, g, viscosity_H2O): 
    

    vp = (4*g*(rho_biomass - rho_H2O)*dp/3/rho_H2O/(24/Re + 6/(1+(Re)**0.5)+0.4))**0.5
    
    f = 