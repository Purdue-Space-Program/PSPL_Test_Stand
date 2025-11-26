## Zoha Peterson
## PSP-Liquids, Test Stand
## 11/8/24
## Calculate minimum GOx and H2 pressures to deliver a pressure of 1000psi and 22g/s to the igniter
## Assumes a pressure of 1000psi is delivered to the igniter chamber; assumes flow is choked over the orifice; 
## assumes theres no choking through regulator; assumes no pressure drop across check valve or pneumatic valve; assumed ideal gas

import pint
import math as m
from CoolProp.CoolProp import PropsSI as p
from scipy.optimize import fsolve

## ------------------------ FUNCTIONS ------------------------
# Calculate standard density
def calcStandardDensity(P, T, R):
    rho = P / (R * T) # [kg/m^3]
    return rho

# Calculate Gamma; return float
def calcGamma(P, T, working_fluid):
    P.to_base_units() # [Pa]
    T.to_base_units() # [K]
    c_p = p('Cpmass', 'T', T.magnitude, 'P', P.to(u.pascal).magnitude, working_fluid) * (u.joule / (u.kilogram * u.kelvin)) # [J/(kg*K)]
    c_v = p('Cvmass', 'T', T.magnitude, 'P', P.to(u.pascal).magnitude, working_fluid) * (u.joule / (u.kilogram * u.kelvin)) # [J/(kg*K)] 
    
    return (c_p / c_v) # [unitless], ratio of specific heats

# Calcualte Specifc Volume; return float
def calcSpecifcGravity(P, T, working_fluid, ref_fluid):
    P.to_base_units() # [Pa]
    T.to_base_units() # [K]
    workingFluid_molarMass = p("M", "T", T.magnitude, "P", (P.to(u.pascal)).magnitude, working_fluid) # [kg/mol]
    referenceFluid_molarMass = p("M", "T", T.magnitude, "P", (P.to(u.pascal)).magnitude, ref_fluid) # [kg/mol]

    return (workingFluid_molarMass / referenceFluid_molarMass) # [unitless], density of gas over density of reference fluid, in this case (for gases) air is used

# Calculate Critical Pressure, return float
def calcChokedUpstreamPressure(downstream_P, gamma):
    upstream_pressure = downstream_P / ((2 / (gamma + 1)) ** (gamma / (gamma - 1))) # [Pa]
    return upstream_pressure # [Pa]

# Equation for pressure drop through regulator, not choking through the regulator
def pressureDrop_equation(P1, q, Cv, G, T, P2, N2):
    return (-1 * q) + (N2 * Cv * P1) * (1 - (2 * (P1 - P2)) / (3 * P1)) * (((P1 - P2) / (P1 * G * T)) ** (0.5)) # per swage handbook - pressure drop for gases

# Calcualte regulator Upstream Pressure; return float
def calcComponentUpstreamPressure(q, Cv, G, T, P2, N2):
    initial_guess = 2000 # [psi], upstream pressure initial guess for fsolve
    upstream_presure = fsolve(lambda P1: pressureDrop_equation(P1, q.magnitude, Cv.magnitude, G, T.to(u.rankine).magnitude, P2.to(u.lbf / (u.inch ** 2)).magnitude, N2), initial_guess) # find the root of the equation
    return upstream_presure[0] # [psi]

## ------------------------ DEFINE UNITS ------------------------
u = pint.UnitRegistry()

## ------------------------ INPUTS ------------------------
# General
standardPressure = 101325 * (u.pascal) # [Pa], standard atmospheric pressure
standardTemperature = 300 * (u.kelvin) # [K], standard atmospheric temperature
reference_fluid = "Air"
N2 = 22.67 # [-], swage constant

# GOx
oxidizer = "Oxygen"
GOx_massFlow = 12.5 * (u.grams / u.second) # [g/s], requirement: GOx mass flow
GOx_regulatorCv = 0.06 * (u.dimensionless) # [unitless], GOx regulator Cv
GOx_targetPressure = 750 * (u.lbf / (u.inch ** 2)) # [lbf/in^2], requirement
GOx_standardR = 259.84 * ((u.joule) / (u.kilogram * u.kelvin)) # [J/kgK]; assumes ideal gas and standard conditions
GOx_temperature = 280 * (u.kelvin) # [K], estimate

# H2
fuel = "Hydrogen"
H2_massFlow = 5 * (u.grams / u.second) # [g/s], requirement: GOx mass flow
H2_regulatorCv = 0.06 * (u.dimensionless) # [unitless], H2 regulator Cv
H2_targetPressure = 750 * (u.lbf / (u.inch ** 2)) # [lbf/in^2], requirement
H2_standardR = 4124.2 * ((u.joule) / (u.kilogram * u.kelvin)) # [J/kgK]; assumes ideal gas and standard conditions
H2_temperature = 280 * (u.kelvin) # [K], estimate

# Unit Conversion
# GOx
GOx_targetPressure = GOx_targetPressure.to_base_units().to(u.pascal) # [Pa]
GOx_massFlow = GOx_massFlow.to_base_units() # [kg/s]

# H2
H2_targetPressure = H2_targetPressure.to_base_units().to(u.pascal) # [Pa]
H2_massFlow = H2_massFlow.to_base_units() # [kg/s]

## ------------------------ CALCULATIONS ------------------------
# GOx
GOx_standardDensity = (calcStandardDensity(standardPressure, standardTemperature, GOx_standardR)).to(u.kilogram / (u.meter ** 3)) # [kg/m^3], oxygen density at atmospheric conditions
GOx_scfm = (GOx_massFlow / GOx_standardDensity).to((u.feet ** 3) / (u.minute)) # [SCFM], GOx
GOx_gamma = calcGamma(GOx_targetPressure, GOx_temperature, oxidizer) # [unitless], ratio of specific heats
GOx_specificGravity = calcSpecifcGravity(GOx_targetPressure, GOx_temperature, oxidizer, reference_fluid) # [unitless], GOx specific gravity
GOx_upstreamOrficePressure = calcChokedUpstreamPressure(GOx_targetPressure, GOx_gamma) # [Pa], minimum pressure needed to check at orifice
GOx_upstreamRegulatorPressure = calcComponentUpstreamPressure(GOx_scfm, GOx_regulatorCv, GOx_specificGravity, GOx_temperature, GOx_upstreamOrficePressure, N2) * (u.lbf / (u.inch ** 2)) # [psi], upstream pressure

# H2
H2_standardDensity = calcStandardDensity(standardPressure, standardTemperature, H2_standardR).to(u.kilogram / (u.meter ** 3)) # [kg/m^3], hydrogen density at atmospheric conditions
H2_scfm = (H2_massFlow / H2_standardDensity).to((u.feet ** 3) / (u.minute)) # [SCFM], H2
H2_gamma = calcGamma(H2_targetPressure, H2_temperature, fuel) # [unitless], ratio of specific heats
H2_specificGravity = calcSpecifcGravity(H2_targetPressure, H2_temperature, fuel, reference_fluid) # [unitless], H2 specific gravity
H2_upstreamOrficePressure = calcChokedUpstreamPressure(H2_targetPressure, H2_gamma) # [Pa], minimum pressure needed to check at orifice
H2_upstreamRegulatorPressure = calcComponentUpstreamPressure(H2_scfm, H2_regulatorCv, H2_specificGravity, H2_temperature, H2_upstreamOrficePressure, N2) * (u.lbf / (u.inch ** 2)) # [psi], upstream pressure

## ------------------------ OUTPUT ------------------------
print("Minimum Source Pressures to reach 1000psi chamber pressure:")
print("GOx: " + str(GOx_upstreamRegulatorPressure))
print("H2: " + str(H2_upstreamRegulatorPressure))