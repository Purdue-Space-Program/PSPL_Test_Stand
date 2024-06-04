from CoolProp import CoolProp as cp
import numpy as np
import matplotlib.pyplot as plt

psi_to_pascal = 6894.76 # [Pa/psi] Conversion Factor
kg_to_lb = 2.20462 # [lbm/kg] Conversion Factor
m3_to_ft3 = 35.3147 # [ft^3/m^3] Conversion Factor

p_fu = 1000 * psi_to_pascal # [Pa] Fuel pressure
t_fu = 120 # [K] Fuel temperature
fuel = 'methane'

p_ox = 300 * psi_to_pascal # [Pa] Oxidizer pressure
t_ox = 116 # [K] Oxidizer temperature
oxidizer = 'oxygen'

rho_fu = cp.PropsSI('D', 'P', p_fu, 'T', t_fu, fuel) # [kg/m^3] Fuel density
rho_ox = cp.PropsSI('D', 'P', p_ox, 'T', t_ox, oxidizer) # [kg/m^3] Oxidizer density

print("Fuel density (lbm/ft^3): " + str(round(rho_fu * kg_to_lb / m3_to_ft3, 2)))
print("Oxidizer density (lbm/ft^3): " + str(round(rho_ox * kg_to_lb / m3_to_ft3, 2)))
