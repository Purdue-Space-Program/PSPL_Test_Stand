from CoolProp import CoolProp as cp
import numpy as np
import matplotlib.pyplot as plt

###########################################################################
# Ullage Volume Calcs (Test Stand Tanks)
#
# Author: Max Van de Graaff
# Contact me on Slack @ Max Van de Graaff
#
# Description: This code calculates required ullage volume based on 
# anticipated bang-bang pressure spike
#
# Notes: You will need the Python CoolProp, NumPy, and MatPlotLib libraries 
# to run this script.
#
# Key Assumptions:
# - Assume tank ullage volume contains only gaseous pressurant at start, not 
#   evaporated propellants
# - Assume isentropic conditions (same specific entropy) within the tank
#   and COPV relative to starting enthalphy values for each
# - Assume choked flow (not a valid assumption for lower COPV pressures)
#
###########################################################################

###########################################################################
# Function Definitions
def getChoked(p_tank, p_He):
    press_ratio = p_He / p_tank # [unitless ratio]
    if press_ratio > 1.9:
        choked = 1 # [binary]
    else:
        choked = 0 # [binary]
    return(choked)

###########################################################################
# Define inputs
in_to_m = (1/39.37); # [in/m] The conversion factor between in and m.
psi_to_pascal = 6894.76 # [Pa/psi]
sec_to_millisec = 10 ** -3
tank_diameter = 6.625; # [in] The outer diameter of the tanks.
p_He = 5500 # [psi] Helium pressure in COPV initially
p_tank = 1650 # [psi] Gas pressure in tank ullage volume initially
t_He = 300 # [K] Helium temperature
t_tank = 117 # [K] Tank temperature
pressurant = 'helium'
cycle_time = 125 # [millisec] Bang-bang valve cycle time
sol_Cd = 0.6 # [unitless] Solenoid discharge coefficient
sol_ESEOD = 0.058 # [in] Solenoid equivalent sharp edge orifice diameter
vol_He = 0.009 # [m^3] Helium COPV volume
vol_tank = 0.0045 # [m^3] Tank ullage volume
num_valves = 2 # [unitless] Number of valves flowing

###########################################################################
# Convert to metric units
tank_diameter *= in_to_m; # [m]
p_He *= psi_to_pascal # [Pa]
p_tank *= psi_to_pascal # [Pa]
sol_ESEOD *= in_to_m

###########################################################################
# Setup calcs
p_tank_init = p_tank # [Pa]
times = range(0, cycle_time)
tankPressures = []
HePressures = []
mDots = []
sol_area = np.pi * (sol_ESEOD / 2) ** 2 # [m^2] Equivalent valve orifice area
s_He = cp.PropsSI('S', 'P', p_He, 'T', t_He, pressurant) # [kJ/kgK] Helium specific entropy
s_tank = cp.PropsSI('S', 'P', p_tank, 'T', t_tank, pressurant) # [kJ/kgK] Tank specific entropy

for time in range(0, cycle_time):
    # choked = getChoked(p_tank, p_He)
    rho_He = cp.PropsSI('D', 'S', s_He, 'P', p_He, pressurant) # [kg/m^3] Helium density
    rho_tank = cp.PropsSI('D', 'S', s_tank, 'P', p_tank, pressurant) # [kg/m^3] Helium density
    Cp_He = cp.PropsSI('C', 'S', s_He, 'P', p_He, pressurant) # [kJ/kgK] Specific heat capacity for helium at contant pressure
    Cv_He = cp.PropsSI('O', 'S', s_He, 'P', p_He, pressurant) # [kJ/kgK] Specific heat capacity for helium at contant volume
    k_ratio_He = Cp_He / Cv_He # [unitless] Heat capacity ratio 
    m_He = rho_He * vol_He # [kg] Mass of pressurant in tank ullage
    m_tank = rho_tank * vol_tank # [kg] Mass of pressurant in tank ullage

    mDot = num_valves * sol_Cd * sol_area * np.sqrt(rho_He * k_ratio_He * p_He * (2 / (k_ratio_He + 1)) ** ((k_ratio_He + 1)/(k_ratio_He - 1)))
    # [kg/sec] Choked mass flow of helium from COPV to tank ullage volume

    # Update Masses
    mDot *= sec_to_millisec # [kg/ms]
    m_He -= mDot # [kg]
    m_tank += mDot # [kg]

    # Update Density, Temp, Press
    rho_He = m_He / vol_He # [kg/m^3]
    p_He = cp.PropsSI('P', 'D', rho_He, 'S', s_He, pressurant) # [Pa] Helium pressure
    t_He = cp.PropsSI('T', 'D', rho_He, 'S', s_He, pressurant) # [K] Helium temperature
    rho_tank = m_tank / vol_tank # [kg/m^3]
    p_tank = cp.PropsSI('P', 'D', rho_tank, 'S', s_tank, pressurant) # [Pa] Tank pressure
    t_tank = cp.PropsSI('T', 'D', rho_tank, 'S', s_tank, pressurant) # [K] Tank temperature

    tankPressures.append(p_tank / psi_to_pascal)
    HePressures.append(p_He / psi_to_pascal)
    mDots.append(mDot / sec_to_millisec)

    time += 1

###########################################################################
# Return outputs
print("##########################################################################")
print("Final tank pressure after one cycle (psi): ", round((p_tank / psi_to_pascal), 2))
print("Net tank pressure increase (psi): ", round(((p_tank - p_tank_init) / psi_to_pascal), 2))
print("##########################################################################")

plt.figure() # Tank Pressure Plot
plt.plot(times, tankPressures, 'r-')
plt.title("Tank Pressure over Time")
plt.xlabel("Time (Milliseconds)")
plt.ylabel("Tank Pressure (psi)")

plt.figure() # Mass Flow Plot
plt.plot(times, mDots, 'b-')
plt.title("Mass Flow Rate over Time")
plt.xlabel("Time (Milliseconds)")
plt.ylabel("Mass Flow Rate (kg/s)")
plt.show()