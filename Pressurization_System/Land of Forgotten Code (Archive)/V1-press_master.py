from CoolProp import CoolProp as cp
from pint import UnitRegistry
import numpy as np
u = UnitRegistry()
Q_ = u.Quantity

# ------- CONSTANTS -------

# Standard Temperature and Pressure
T_STD = Q_(60, 'degF').to('kelvin')  # Standard temperature
P_STD = 14.696 * u('psi')  # Standard pressure

# Nitrogen Constants
R_N2 = 0.2968 * u('kJ/(kg*K)')  # R constant for N2
GAMMA_N2 = 1.47 * u('dimensionless')  # Gamma constant for N2
S_G_GN2 = Q_(0.9669,'dimensionless')

# Helium Constants
R_HE = 0.0 * u('kJ/(kg*K)')
GAMMA_N2 = 0.0 * u('dimensionless')  # Gamma constant for N2
S_G_GN2 = Q_(0.0,'dimensionless')

# ------- INPUTS -------

# Operations
P_FILL = 40 * u('psi')  # Fill pressure; maybe should assume atmospheric

# Bottle
P_BOTTLE = 6000 * u('psi')  # Bottle max pressure
P_BOTTLE_MIN = 1400 * u('psi')  # Bottle min pressure
V_BOTTLE = 43.3 * u('L')  # Bottle volume (standard GN2 cylinder)

# COPV
P_COPV = 6000 * u('psi')  # COPV max pressure
P_COPV_MIN = 1400 * u('psi')  # COPV min pressure
V_COPV = 40.0 * u('L')  # COPV volume (big CMS)

P_MAWP = 980 * u('psi')  # Propellant tank MAWP (maximum MAWP)
V_TANK = np.pi*72*(5.761/2)**2 * u('in^3')  # Propellant tank volume | ID = 5.761"
V_TANK = 1738.368 * u('in^3')
print("Tank volume: ", V_TANK)
T_ISENT_TANK = T_STD * (P_MAWP / P_BOTTLE)**(1 - 1 / GAMMA_N2)  # Temperature of pressurant gas after isentropically expanding
P_CHOKE = P_MAWP * 1.895 # For N2 only

# Liquid Oxygen
MDOT_LOX = 8 * u('kg/s')  # LOX mass flow out of tank

# Methane
MDOT_FU = 3 * u('kg/s')  # Methane mass flow out of tank 

# CALCS
# Liquid Oxygen
rho_lox = cp.PropsSI('D', 'P', P_FILL.to('Pa').m, 'Q', 0, 'Oxygen') * u('kg/m^3')
Tsat_lox = cp.PropsSI('T', 'P', P_FILL.to('Pa').m, 'Q', 0, 'Oxygen') * u('kelvin')
vdot_lox = MDOT_LOX / rho_lox  # LOX volumetric flow out of tank 

# Methane
rho_fu = cp.PropsSI('D', 'P', P_FILL.to('Pa').m, 'Q', 0, 'Methane') * u('kg/m^3')
Tsat_fu = cp.PropsSI('T', 'P', P_FILL.to('Pa').m, 'Q', 0, 'Methane') * u('kelvin')
vdot_fu = MDOT_FU / rho_fu  # Methane volumetric flow out of tank

# GN2
scfm_gn2_fu = vdot_fu * (P_MAWP / P_STD) * (T_STD / T_ISENT_TANK)
scfm_gn2_fu.ito('ft^3/min')
mdot_gn2_fu = (P_MAWP * vdot_fu) / (R_N2 * T_ISENT_TANK)

scfm_gn2_lox = vdot_lox * (P_MAWP / P_STD) * (T_STD / T_ISENT_TANK)
scfm_gn2_lox.ito('ft^3/min')
mdot_gn2_lox = (P_MAWP * vdot_lox) / (R_N2 * T_ISENT_TANK)   # CHANGE TO ISENTROPIC TEMP

throat_area_lox = ((mdot_gn2_lox * np.sqrt(T_STD)) / (P_BOTTLE * np.sqrt(GAMMA_N2 / R_N2) * np.power(((GAMMA_N2 + 1) / 2), (-(GAMMA_N2 + 1)/(2 * (GAMMA_N2 - 1))))))
throat_dia_lox = 2 * np.sqrt(throat_area_lox / np.pi)
throat_area_lox.ito('in^2')
throat_dia_lox.ito('in')

throat_area_fu = ((mdot_gn2_fu * np.sqrt(T_STD)) / (P_BOTTLE * np.sqrt(GAMMA_N2 / R_N2) * np.power(((GAMMA_N2 + 1) / 2), (-(GAMMA_N2 + 1)/(2 * (GAMMA_N2 - 1))))))
throat_dia_fu = 2 * np.sqrt(throat_area_fu / np.pi)
throat_area_fu.ito('in^2')
throat_dia_fu.ito('in')

print("Isentropic temp:", f"{T_ISENT_TANK.to('kelvin'): ~P.4f}")
print("T_sat (lox, fu):", f"{Tsat_lox.to('fahrenheit'): ~P.4f},  " f"{Tsat_fu.to('fahrenheit'): ~P.4f}")
print("rho (lox, fu):", f"{rho_lox.to('kg/m^3'): ~P.4f},  " f"{rho_fu.to('kg/m^3'): ~P.4f}")
print("m_dot (lox, fu):", f"{mdot_gn2_lox.to('kg/s'): ~P.4f},  " f"{mdot_gn2_fu.to('kg/s'): ~P.4f}")
print("SCFM (lox, fu):", f"{scfm_gn2_lox: ~P.4f},  " f"{scfm_gn2_fu: ~P.4f}")
print("Area (lox, fu):", f"{throat_area_lox: ~P.4f},        ", f"{throat_area_fu: ~P.4f}")
print("Dia (lox, fu): ", f"{throat_dia_lox: ~P.4f},         ", f"{throat_dia_fu: ~P.4f}")

'''
print("Area, one 3/8 line:", f"{np.pi*((0.378-0.049*2)/2)**2 * u('in^2'): ~P.4f}")
print("Area, one 1/4 line:", f"{np.pi*((0.25-0.049*2)/2)**2 * u('in^2'): ~P.4f}")
'''
# Cv Calcs


# Choked flow regulator Cv calculations
C_v_fu = (scfm_gn2_fu.to('ft^3/min').m * 2 * np.sqrt(S_G_GN2.m)) / P_BOTTLE.to('psi').m  # Note:units in psia and SCFM
C_v_lox = (scfm_gn2_lox.to('ft^3/min').m * 2 * np.sqrt(S_G_GN2.m)) / P_BOTTLE.to('psi').m  # Note:units in psia and SCFM

# 5828.5310116569635
# COPV Calcs
v_1_copv_isentropic = V_TANK * (P_MAWP / P_COPV)**(1 / GAMMA_N2)  # Volume at P_COPV that expands to fill tanks
m_1_copv_isentropic = (P_COPV * v_1_copv_isentropic) / (R_N2 * T_STD)  # Mass at P_COPV that fills tanks

m_copv_postburn = (P_CHOKE * V_COPV) / (R_N2 * T_STD)  # Mass required to keep flow choked

m_total = (m_copv_postburn + m_1_copv_isentropic * 2).to('kg')
P_copv_init = (m_total * R_N2 * T_STD) / V_COPV
print(P_copv_init.to('psi'))
