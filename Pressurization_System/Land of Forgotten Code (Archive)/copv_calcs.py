from CoolProp import CoolProp as cp
from pint import UnitRegistry
import numpy as np
import matplotlib.pyplot as plt
u = UnitRegistry()
Q_ = u.Quantity

# INPUTS
# Constants
T_STD = Q_(60, 'degF').to('kelvin')  # Standard temperature
P_STD = 14.696 * u('psi')  # Standard pressure
R_N2 = 0.2968 * u('kJ/(kg*K)')  # R constant for N2
GAMMA_N2 = 1.47 * u('dimensionless')  # Gamma constant for N2
S_G_GN2 = Q_(0.9669,'dimensionless')

# Defined Constant Quantities
P_FILL = 40 * u('psi')  # Fill pressure; maybe should assume atmospheric

P_BOTTLE = 6000 * u('psi')  # Bottle max pressure
P_BOTTLE_MIN = 1400 * u('psi')  # Bottle min pressure
V_BOTTLE = 43.3 * u('L')  # Bottle volume (standard GN2 cylinder)

P_COPV = 6000 * u('psi')  # COPV max pressure
P_COPV_MIN = 1400 * u('psi')  # COPV min pressure
V_COPV = 40.0 * u('L')  # COPV volume (big CMS)

P_mawp = 1400 * u('psi')  # Propellant tank MAWP (maximum MAWP)
V_TANK = np.pi*72*(5.761/2)**2 * u('in^3')  # Propellant tank volume | ID = 5.761"
V_TANK = 1738.368 * u('in^3')
print("Tank volume: ", V_TANK)
T_ISENT_TANK = T_STD * (P_mawp / P_BOTTLE)**(1 - 1 / GAMMA_N2)  # Temperature of pressurant gas after isentropically expanding

P_choke = P_mawp * u('psi') * 1.895 # For N2 only

copv_pressures = []
tank_pressures = range(1,1400)
for i in range(1,1400):
    # COPV Calcs
    P_mawp = i * u('psi')
    P_choke = P_mawp * 1.895 # For N2 only

    v_1_copv_isentropic = V_TANK * (P_mawp / P_COPV)**(1 / GAMMA_N2)  # Volume at P_COPV that expands to fill tanks
    m_1_copv_isentropic = (P_COPV * v_1_copv_isentropic) / (R_N2 * (T_STD))  # Mass at P_COPV that fills tanks

    m_copv_postburn = (P_choke * V_COPV) / (R_N2 * T_STD)  # Mass required to keep flow choked

    m_total = (m_copv_postburn + m_1_copv_isentropic * 2).to('kg')
    P_copv_init = (m_total * R_N2 * T_STD) / V_COPV

    if P_copv_init.to('psi').m <= 4350:
        max_tank_press = P_mawp

    copv_pressures.append(P_copv_init.to('psi').m)

print("Max achievable tank pressure at 4350 copv press: ", max_tank_press)
plt.figure(1)
plt.plot(tank_pressures, copv_pressures)
plt.axhline(y = 4350, color = 'r', linestyle = '--')
plt.xlabel('Tank Pressure (psia)')
plt.ylabel('COPV Starting Pressure (psia)')
plt.title('Initial COPV pressure as a function of tank pressure')

plt.show()