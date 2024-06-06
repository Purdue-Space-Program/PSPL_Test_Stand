from press_master import *
from pint import UnitRegistry
import numpy as np
import matplotlib.pyplot as plt

'''
Size for big COPV: 40L, 5500 psi
'''

S_G_GN2 = Q_(0.9669,'dimensionless')
P_bottle_min = 1401 * u('psi') # Determine this with isentropic expansion

# Choked flow regulator Cv calculations
C_v_fu = (scfm_gn2_fu.to('ft^3/min').m * 2 * np.sqrt(S_G_GN2.m)) / P_BOTTLE.to('psi').m  # Note:units in psia and SCFM
C_v_lox = (scfm_gn2_lox.to('ft^3/min').m * 2 * np.sqrt(S_G_GN2.m)) / P_BOTTLE.to('psi').m  # Note:units in psia and SCFM

# print(scfm_gn2_fu.to('ft^3/hr').m * (np.sqrt(S_g_gn2.m * T_STD.to('rankine')) / (816 * P_BOTTLE.to('psi')))) Other formula I found?

# Assume choked flow when p1 > 2 * p2


'''
C_v_list_fu = []
lowP_fu = range(1401,6001)
for i in range(1401,6001):
    P_bottle_min = i * u('psi')
    if i <= 2800:
        C_v_lowP_fu = (scfm_gn2_fu.to('ft^3/hr').m / 962) * np.sqrt((S_g_gn2.m * T_STD.to('rankine')) / ((P_bottle_min.to('psi').m)**2 - (P_MAWP.to('psi').m)**2))
    else:
        C_v_slowP_fu = scfm_gn2_fu.to('ft^3/hr').m * (np.sqrt(S_g_gn2.m * T_STD.to('rankine')) / (816 * P_BOTTLE.to('psi')))

    C_v_list_fu.append(C_v_lowP_fu.m)

# Assume choked flow when p1 > 2 * p2
C_v_list_lox = []
lowP_lox = range(1401,6001)
for i in range(1401,6001):
    P_bottle_min = i * u('psi')
    if i <= 2800:
        C_v_lowP_lox = (scfm_gn2_lox.to('ft^3/hr').m / 962) * np.sqrt((S_g_gn2.m * T_STD.to('rankine')) / ((P_bottle_min.to('psi').m)**2 - (P_MAWP.to('psi').m)**2))
    else:
        C_v_lowP_lox = scfm_gn2_lox.to('ft^3/hr').m * (np.sqrt(S_g_gn2.m * T_STD.to('rankine')) / (816 * P_BOTTLE.to('psi')))

    C_v_list_lox.append(C_v_lowP_lox.m)
'''

# print("Formula 1:")
print("Required Cv for Methane: ", f"{C_v_fu: .4f}")
print("Required Cv for LOX: ", f"{C_v_lox: .4f}")

'''
print("Formula 2 fu Cv:", C_v_lowP_fu)
print("Formula 2 lox Cv:", C_v_lowP_lox)

plt.figure(1)
plt.plot(lowP_fu, C_v_list_fu)
plt.axvline(x = 2800, color = 'r', linestyle = '--')
plt.xlabel('Bottle pressure (psia)')
plt.ylabel('Flow Coefficient (Cv)')
plt.title('Bottle Pressure vs. Required Cv to Match Propellant Flow Rate (Fuel Press)')

plt.figure(2)
plt.plot(lowP_lox, C_v_list_lox)
plt.axvline(x = 2800, color = 'r', linestyle = '--')
plt.xlabel('Bottle pressure (psia)')
plt.ylabel('Flow Coefficient (Cv)')
plt.title('Bottle Pressure vs. Required Cv to Match Propellant Flow Rate (LOX Press)')
plt.show()
'''