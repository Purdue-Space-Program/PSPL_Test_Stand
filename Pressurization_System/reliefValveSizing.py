## Zoha Peterson
## PSP-Liquids, Test Stand
## 10/25/24
## Calcualte Cv needed for pressure relief valves 
## BELOW MR-N2-301, -302, -321 on Master P&ID
## Analysis assumes regulator is fully open and relief valve is open; 6ksi tank pressure, 1.5ksi relief pressure, 0.08 Regulator Cv, 14.7psi ambient pressure, 298.15K ambient temperature

import pint
from CoolProp.CoolProp import PropsSI
import math as m

## ------------------------ FUNCTIONS ------------------------
# Calculate isentropic temperature due to isentropic expansion
def calcIsentropicExpansion(pressure1, pressure2, tankTemperature):
    gamma = calcGamma(tankTemperature, pressure1, working_fluid) # [unitless], ratio of specific heats

    return tankTemperature * ((pressure2 / pressure1) ** (1 - (1 / gamma))) # [K], temperature after isentropic expansion

# Calculate gamma at a given pressure and temp, return float
def calcGamma(upstream_T, upstream_P, working_fluid):
    c_p = PropsSI('Cpmass', 'T', upstream_T.magnitude, 'P', upstream_P.to(units.pascal).magnitude, working_fluid) * (units.joule / (units.kilogram * units.kelvin)) # [J/(kg*K)]
    c_v = PropsSI('Cvmass', 'T', upstream_T.magnitude, 'P', upstream_P.to(units.pascal).magnitude, working_fluid) * (units.joule / (units.kilogram * units.kelvin)) # [J/(kg*K)] 
    
    return (c_p / c_v) # [unitless], ratio of specific heats

# Get specific volume, return float
def calcSpecificGravity(working_T, working_P, working_fluid, ref_T, ref_P, ref_fluid):
    density_N2 = PropsSI("D", "T", working_T.magnitude, "P", (working_P.to(units.pascal)).magnitude, working_fluid) # [kg/m^3]
    density_water = PropsSI("D", "T", ref_T.magnitude, "P", (ref_P.to(units.pascal)).magnitude, ref_fluid) # [kg/m^3]

    return (density_N2 / density_water) # [unitless], density of gas over density of reference fluid, in this case (for gases) air is used

# Check choked flow, return boolean
def checkChokedFlow(upstream_P, downstream_P, gamma):
    critical_pressure = upstream_P * ((2 / (gamma + 1)) ** (gamma / (gamma - 1)))

    if (downstream_P <= critical_pressure):
        return True
    else:
        return False

# Calculate the Volumetric flow rate through the Regulator, return float
def calcVolumetricFlow(constant_N2, cv, upstream_P, downstream_P, upstream_T, working_fluid, ambient_pressure, ambient_temp, reference_fluid):
    gamma = calcGamma(upstream_T, upstream_P, working_fluid) # [unitless], ratio of specific heats
    choked = checkChokedFlow(upstream_P, downstream_P, gamma) # [unitless]; boolean, stores whether flow is choked or not
    SG = calcSpecificGravity(upstream_T, upstream_P, working_fluid, ambient_temp, ambient_pressure, reference_fluid) # [unitless]


    delta_P = upstream_P - downstream_P # [psi], difference in pressure
    upstream_T = upstream_T.to(units.rankine) # [Ranking], formula calls for rankine temperature (Swagelok)

    if (not choked):
        Q = constant_N2 * cv * upstream_P * (1 - ((2 * delta_P) / (3 * upstream_P))) * ((delta_P / (upstream_P * SG * upstream_T)) ** (1/2)) # [ft^3/min], volumetric flow rate
    else:
        Q = 0.471 * constant_N2 * cv * upstream_P * ((1 / (SG * upstream_T)) ** (1/2)) # [ft^3/min], volumetric flow rate

    return Q # [ft^3/min], volumetric flow rate

# Calculate the needed Cv for a valve
def calcValveCv(constant_N2, Q, upstream_P, downstream_P, upstream_T, working_fluid, ambient_pressure, ambient_temp, reference_fluid):
    gamma = calcGamma(upstream_T, upstream_P, working_fluid) # [unitless], ratio of specific heats
    choked = checkChokedFlow(upstream_P, downstream_P, gamma) # [unitless]; boolean, stores whether flow is choked or not
    SG = calcSpecificGravity(upstream_T, upstream_P, working_fluid, ambient_temp, ambient_pressure, reference_fluid) # [unitless]

    delta_P = upstream_P - downstream_P # [psi], difference in pressure
    upstream_T = upstream_T.to(units.rankine) # [Ranking], formula calls for rankine temperature (Swagelok)

    if (not choked):
        cv = Q / (constant_N2 * upstream_P * (1 - ((2 * delta_P) / (3 * upstream_P))) * ((delta_P / (upstream_P * SG * upstream_T)) ** (1/2))) # [unitless], flow coefficient
    else:
        cv = Q / (0.471 * constant_N2 * upstream_P * ((1 / (SG * upstream_T)) ** (1/2))) # [unitless], flow coefficient

    return cv # [unitless], flow coefficient

## ------------------------ DEFINE UNITS ------------------------
units = pint.UnitRegistry()

## ------------------------ INPUTS ------------------------
# Hardware
regulator_cv = 0.08 * units.dimensionless # [unitless], regulator upstream of relief valves

# Fluids
working_fluid = "nitrogen" # working fluid
reference_fluid = "air" # for calculating specific volume

# Pressures
relief_pressure = 1500 * (units.lbf / (units.inch ** 2)) # [psi]
tank_pressure = 6000 * (units.lbf / (units.inch ** 2)) # [psi]
ambient_pressure = 14.7 * (units.lbf / (units.inch ** 2)) # [psi]

# Temperatures
ambient_temp = 298.15 * (units.kelvin) # [K]
tankTemp_N2 = 250 * (units.kelvin) # [K], estimate - need to revise; note that the Zucrow tank volume is 250 cubic feet

# Constants
constant_N2 = 22.67 * units.dimensionless # [unitless], constant for calculating the proper units; per Swagelok Valve Sizing

## ------------------------ CALCULATIONS ------------------------
volumetricFlow = (calcVolumetricFlow(constant_N2, regulator_cv, tank_pressure, ambient_pressure, tankTemp_N2, working_fluid, ambient_pressure, ambient_temp, reference_fluid)).magnitude * ((units.foot ** 3) / (units.minute)) # [ft^3/min]
relief_temp = calcIsentropicExpansion(tank_pressure, relief_pressure, tankTemp_N2)
valve_cv = calcValveCv(constant_N2, volumetricFlow, relief_pressure, ambient_pressure, relief_temp, working_fluid, ambient_pressure, ambient_temp, reference_fluid) # [unitless]

## ------------------------ OUTPUT ------------------------
print("Volumetric flow: " + str(volumetricFlow))
print("Cv, Relief Valve: " + str(valve_cv.magnitude))