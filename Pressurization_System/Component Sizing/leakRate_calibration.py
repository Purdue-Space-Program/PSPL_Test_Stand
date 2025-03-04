# IMPORT FUNCTIONS
import pint
import math as m
import numpy as np
from CoolProp.CoolProp import PropsSI as p
from press_funcs import *

# SETUP PINT
#u = pint.UnitRegistry()

# FUNCTIONS
# Calculate ratio of specific heats; return float
def calcGamma(P, T, working_fluid):
    P.to_base_units() # [Pa]
    T.to_base_units() # [K]

    c_p = p('Cpmass', 'T', T.magnitude, 'P', P.magnitude, working_fluid) * (u.joule / (u.kilogram * u.kelvin)) # [J/(kg*K)]
    c_v = p('Cvmass', 'T', T.magnitude, 'P', P.magnitude, working_fluid) * (u.joule / (u.kilogram * u.kelvin)) # [J/(kg*K)] 
    
    gamma = (c_p / c_v) # [unitless], ratio of specific heats

    return gamma

# Calculate temperature due to isentropic expansion; return float
def calcIsentropicExpansion(P1, P2, T1, gamma):
    P1.to_base_units() # [Pa]
    P2.to_base_units() # [Pa]
    T1.to_base_units() # [K]

    T2 = T1 * ((P2 / P1) ** (1 - (1 / gamma))) # [K], temperature after isentropic expansion

    return T2

# Calculate density at given conditions
def calcDensity(P, T, R):
    P.to_base_units() # [Pa]
    T.to_base_units() # [K]

    rho = P / (R * T) # [kg/m^3]

    return rho

# Calculate leak rate; return float
def calculateLeakRate(V, MM, R, T, rho, P1, P2, t, P_amb, T_amb):
    LR = (((V * MM)  / (R * T * rho)) * ((P1 - P2) / t) * (P1 / P_amb) * (T_amb / T)) # [standard m^3 / s]
    LR = LR.to((u.cm ** 3) / u.minute) # [standard cm^3 / min] or [scc/min]

    return LR


# GIVEN INPUT PARAMETERS
fluid = "nitrogen" # testing with N2
time = 1 * (u.seconds) # [sec]; time intervals at each pressure increment
V = ((3) * ((u.inch ** 3))).to(u.meter ** 3) # [m^3], estimation as of now; need to get real volume

MM = ((14.007 * 2) / 1000) * (u.kg / u.mol)  # [kg/mol]; molar mass of N2
R_mol = (8.314) * (((u.kg * (u.meter ** 2)) / (u.second ** 2)) / (u.mol * u.kelvin))# [J/(K-mol)]; ideal gas constant in terms of mols
R_kg = (R_mol / MM) # [J/(K-kg)]; ideal gas constant in terms of kg

P1 = (2000 * (u.psi)).to(u.kg / (u.meter * u.second ** 2)) # [Pa]; Pressure at time 1
P2 = (1999 * (u.psi)).to(u.kg / (u.meter * u.second ** 2)) # [Pa]; Pressure at time 2
P_amb = (14.7 * (u.psi)).to(u.kg / (u.meter * u.second ** 2)) # [Pa]; ambient pressure

T1 =  (240) * (u.K) # [K]; temperature inside tube volume
T_amb = (240) * (u.K) # [K], ambient temperature


# CALCULATE INPUT PARAMETERS
gamma = (calcGamma(P1, T1, fluid)) # [unitless]; ratio of specific heats

T2 = (calcIsentropicExpansion(P1, P_amb, T1, gamma)) # [K]; temperature after expanding out of the leak

density = min((calcDensity(P2, T2, R_kg)), calcDensity(P1, T1, R_kg)) # [kg/s]; smaller density will result in a more conservative leak rate value


# CALCULATE OUTPUTS
leak_rate = calculateLeakRate(V, MM, R_mol, T2, density, P1, P2, time, P_amb, T_amb) # [scc/sec]; leak rate through component


m_dot = scfm_to_mdot(leak_rate.to('ft^3/min'), "Nitrogen")
# OUTPUT RESULTS
print("Leak rate: " + str(leak_rate.magnitude) + " sccm")
print(f"Leak rate: {leak_rate.to('cm^3/min'): ~P.4f}")
print(f"Leak rate: {m_dot.to('g/s'): ~P.4f}")
