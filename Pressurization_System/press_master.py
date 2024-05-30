from CoolProp import CoolProp as cp
from pint import UnitRegistry
import numpy as np
u = UnitRegistry()
Q_ = u.Quantity

# INPUTS
# Common properties
P_FILL = 40 * u('psi')  # Fill pressure; maybe should assume atmospheric
P_BOTTLE = 6000 * u('psi')  # Upstream (press tank) pressure
P_MAWP = 1400 * u('psi')  # Downstream (propellant tank) pressure
T_STD = Q_(60, 'degF').to('kelvin')  # Standard temperature
P_STD = 14.696 * u('psi')  # Standard pressure
R_N2 = 0.2968 * u('kJ/(kg*K)')  # R constant for N2
GAMMA_N2 = 1.47 * u('dimensionless')  # Gamma constant for N2

# Liquid Oxygen
mdot_lox = 8 * u('kg/s')  # LOX mass flow out of tank

# Methane
mdot_fu = 3 * u('kg/s')  # Methane mass flow out of tank 

# CALCS
# Liquid Oxygen
rho_lox = cp.PropsSI('D', 'P', P_FILL.to('Pa').m, 'Q', 0, 'Oxygen') * u('kg/m^3')
Tsat_lox = cp.PropsSI('T', 'P', P_FILL.to('Pa').m, 'Q', 0, 'Oxygen') * u('kelvin')
vdot_lox = mdot_lox / rho_lox  # LOX volumetric flow out of tank 

# Methane
rho_fu = cp.PropsSI('D', 'P', P_FILL.to('Pa').m, 'Q', 0, 'Methane') * u('kg/m^3')
Tsat_fu = cp.PropsSI('T', 'P', P_FILL.to('Pa').m, 'Q', 0, 'Methane') * u('kelvin')
vdot_fu = mdot_fu / rho_fu  # Methane volumetric flow out of tank

# GN2
scfm_gn2_fu = vdot_fu * (P_MAWP / P_STD) * (T_STD / Tsat_fu)
scfm_gn2_fu.ito('ft^3/min')
mdot_gn2_fu = (P_MAWP * vdot_fu) / (R_N2 * Tsat_fu)

scfm_gn2_lox = vdot_lox * (P_MAWP / P_STD) * (T_STD / Tsat_lox)
scfm_gn2_lox.ito('ft^3/min')
mdot_gn2_lox = (P_MAWP * vdot_lox) / (R_N2 * Tsat_lox)

throat_area_lox = ((mdot_gn2_lox * np.sqrt(T_STD)) / (P_BOTTLE * np.sqrt(GAMMA_N2 / R_N2) * np.power(((GAMMA_N2 + 1) / 2), (-(GAMMA_N2 + 1)/(2 * (GAMMA_N2 - 1))))))
throat_dia_lox = 2 * np.sqrt(throat_area_lox / np.pi)
throat_area_lox.ito('in^2')
throat_dia_lox.ito('in')

throat_area_fu = ((mdot_gn2_fu * np.sqrt(T_STD)) / (P_BOTTLE * np.sqrt(GAMMA_N2 / R_N2) * np.power(((GAMMA_N2 + 1) / 2), (-(GAMMA_N2 + 1)/(2 * (GAMMA_N2 - 1))))))
throat_dia_fu = 2 * np.sqrt(throat_area_fu / np.pi)
throat_area_fu.ito('in^2')
throat_dia_fu.ito('in')

print("Area:",throat_area_lox, throat_area_fu,3)
print("Dia:", throat_dia_lox, throat_dia_fu)

print("One 3/8 line:", np.pi*((0.378-0.049*2)/2)**2)
print("One 1/4 line:", np.pi*((0.25-0.049*2)/2)**2)