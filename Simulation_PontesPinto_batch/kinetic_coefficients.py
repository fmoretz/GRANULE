# Kinetic coefficients

Y  	= 0.1		# gCH4/gS -- from stoichiometry
Ks	= 0.002		# kg/m3   #0.93
K   = 0.046		# g/g    
b   = K/Y
B0  = 0.516 	  # L CH4 /g COD added
bMA = 1.5 * 10**-4  #h^-1
YMA = 0.04  
#Kd  = 0.0007      # h
Kd = YMA *bMA
