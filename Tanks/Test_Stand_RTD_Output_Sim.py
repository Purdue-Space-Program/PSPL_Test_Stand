import numpy as np
import matplotlib.pyplot as plt

###########################################################################
# RTD Current Script
#
# Author: Max Van de Graaff (Based on work by Hugo Filmer)
# Contact me on Slack @ Max Van de Graaff
#
# Description: This script varies the convection coefficient to determine
# the current across an RTD assuming a constant voltage.
#
# Last updated 05/31/2023
#
# Notes:
# - You will need the NumPy and MatPlotLib libraries to run this script.
# - RTD1 is the upper RTD, and RTD 2 is the lower RTD
# - RTD Properties given by 
#   https://evosensors.com/products/cryogenic-2-wire-pt100-rtd-probe-2-inches-long-x-1-8-diameter-with-40-inches-of-pfa-insulated-lead-wire-and-miniature-male-flat-pin-connector
#
###########################################################################

###########################################################################
# Function Definitions

def RTD_Resistance(temp):
    # https://www.engineeringtoolbox.com/pt100-electrical-resistance-d_1651.html

    temp = temp - 273 # Convert from K to C for the fit
    R_0 = 100 # [Ohm] The resistance at 0 C

    # Linearization constants
    A = 3.9083 * 10 ** -3
    B = -5.775 * 10 ** -7

    if temp > 273:
        C = 0
    else:
        C = -4.183 * 10 ** -12

    resistance = R_0 * (1 + A * temp + (B*temp ** 2) + C * (temp - 100) * (temp ** 3)) # [Ohm]
    return(resistance)

###########################################################################
# Define inputs
in_to_m = (1/39.37) # [in/m] The conversion factor between in and m.

voltage_in = -5 # [V] The voltage into the RTD circuit

RTD_length = 1 * in_to_m # [m] The length of the RTD
RTD_diam = 0.125 * in_to_m # [m] The diameter of the RTD
RTD_thick = 0.05 * in_to_m # [m] The thickness of the sheath (estimated)
RTD_area = np.pi * RTD_diam * RTD_length # [m^2] The exposed surface area, assuming RTD is a cylinder
RTD_volume = (((RTD_diam/2) ** 2) - (((RTD_diam - 2 * RTD_thick) / 2) ** 2)) * np.pi * RTD_length # [m^3]

RTD_density = 7930 # [kg/m^3] The density of the RTD, typical value for stainless steel
RTD_mass = RTD_volume * RTD_density; # [kg] The mass of the RTD

RTD_specific_heat = 500 # [J/Kg K] The specific heat of the material, typical value for stainless steel
# RTD_specific_heat = 125 # [J/Kg K] The specific heat of the material, typical value for platinum
RTD_heat_cap = RTD_specific_heat * RTD_mass # [J/K] The heat capacity of the RTD


###########################################################################
# Define fluid properties
temp = 100 # [K] The temp of the fluids, arbitrary value
gas_coeff = 10 # [W/m^2 K] The convection coefficient of the gas, typical value for free convection
liquid_coeff = 100 # [W/m^2 K] The convection coefficient of the liquid, typical value for free convection

###########################################################################
# RTD initial conditions
RTD1_initial_temp = 300 # [K] The initial temperature of the upper RTD
RTD2_initial_temp = 300 # [K] The initial temperature of the lower RTD
RTD3_initial_temp = 300 # [K] The initial temperature of the lower RTD

###########################################################################
# Simulation Parameters
timestep = 0.01 # [s] The timestep that the simulation uses
stop_time = 500 # [s] The time the simulation runs to

# Note: Gas time must be greater than liquid time!
RTD1_liquid_time = 100 # [s] The time that RTD1 (upper) switches from gas to liquid
RTD1_gas_time = 500 # [s] The time that RTD1 (upper) switches from liquid to gas

RTD2_liquid_time = 200 # [s] The time that RTD2 (lower) switches from gas to liquid
RTD2_gas_time = 500 # [s] The time that RTD2 (lower) switches from liquid to gas

RTD3_liquid_time = 375 # [s] The time that RTD2 (lower) switches from gas to liquid
RTD3_gas_time = 500 # [s] The time that RTD2 (lower) switches from liquid to gas

###########################################################################
# Create Data Lists
times = []
RTD1_temps = []
RTD2_temps = []
RTD3_temps = []
RTD1_resistances = []
RTD2_resistances = []
RTD3_resistances = []
voltages_measured = []
currents_measured = []

###########################################################################
# Create Data Lists

RTD1_temp = RTD1_initial_temp
RTD2_temp = RTD2_initial_temp
RTD3_temp = RTD3_initial_temp
index = 1

for t in range (0, stop_time):
    # Determine RTD heat transfer coefficients
    if (t >= RTD1_liquid_time) and (t <= RTD1_gas_time):
        RTD1_coeff = liquid_coeff
    else:
        RTD1_coeff = gas_coeff

    if (t >= RTD2_liquid_time) and (t <= RTD2_gas_time):
        RTD2_coeff = liquid_coeff
    else:
        RTD2_coeff = gas_coeff

    if (t >= RTD3_liquid_time) and (t <= RTD3_gas_time):
        RTD3_coeff = liquid_coeff
    else:
        RTD3_coeff = gas_coeff

    ###########################################################################
    # Calculate Resistance
    RTD1_resistance = RTD_Resistance(RTD1_temp) # [Ohm]
    RTD2_resistance = RTD_Resistance(RTD2_temp) # [Ohm]
    RTD3_resistance = RTD_Resistance(RTD3_temp) # [Ohm]
    
    ###########################################################################
    # Calculate electrical parameters
    voltage_measured = voltage_in*(RTD2_resistance/(RTD1_resistance + RTD2_resistance)) - voltage_in/2 # [V]
    RTD1_voltage = voltage_in*(RTD1_resistance/(RTD1_resistance + RTD2_resistance)) # [V]
    RTD2_voltage = voltage_in*(RTD2_resistance/(RTD1_resistance + RTD2_resistance)) # [V]
    # RTD3_voltage = voltage_in*(RTD3_resistance/(RTD1_resistance + RTD3_resistance)) # [V]
    current_measured = voltage_in/(RTD1_resistance + RTD2_resistance) # [Amp]

    ###########################################################################
    # Calculate heat transfer
    Q1_in = ((RTD1_voltage **2) / RTD1_resistance) * timestep # [J]
    Q1_out = (RTD1_temp - temp) * RTD_area * RTD1_coeff * timestep # [J]

    Q2_in = ((RTD2_voltage ** 2) / RTD2_resistance) * timestep # [J]
    Q2_out = (RTD2_temp - temp) * RTD_area * RTD2_coeff * timestep # [J]

    # Q3_in = ((RTD3_voltage ** 2)/RTD3_resistance) * timestep # [J]
    # Q3_out = (RTD3_temp - temp) * RTD_area * RTD3_coeff * timestep # [J]

    ###########################################################################
    # Find new values of tracked variables
    RTD1_temps.append(RTD1_temp)
    RTD2_temps.append(RTD2_temp)
    RTD3_temps.append(RTD3_temp)
    voltages_measured.append(voltage_measured)
    currents_measured.append(current_measured)
    RTD1_resistances.append(RTD1_resistance)
    RTD2_resistances.append(RTD2_resistance)
    RTD3_resistances.append(RTD3_resistance)

    ###########################################################################
    # Update temperatures
    RTD1_temp = RTD1_temp + (Q1_in - Q1_out) / RTD_heat_cap # [K]
    RTD2_temp = RTD2_temp + (Q2_in - Q2_out) / RTD_heat_cap # [K]
    # RTD3_temp = RTD3_temp + (Q3_in - Q3_out) / RTD3_heat_cap # [K]

###########################################################################
# Find voltage and current derivatives

voltage_derivatives = []

current_derivatives = []

for index in range(0, stop_time - 2):
    voltage_derivatives.append((voltages_measured[index + 1] - voltages_measured[index])/timestep)
    current_derivatives.append((currents_measured[index + 1] - currents_measured[index])/timestep)

times = range(0, 500)

###########################################################################
# Output Plots

#fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
#example_plot(ax1)
#example_plot(ax2)
#example_plot(ax3)
#example_plot(ax4)
#plt.tight_layout()

plt.plot(times, RTD1_temps)
plt.plot(times, RTD2_temps)
plt.plot(times, RTD3_temps)
plt.vlines(RTD1_liquid_time, 250, 300)
#plt.vlines(RTD1_gas_time, "m")
#plt.vlines(RTD2_liquid_time, "g")
#plt.vlines(RTD2_gas_time, "m")
#plt.vlines(RTD3_liquid_time, "b")
#plt.vlines(RTD3_gas_time, "m")
plt.title("RTD Temperature Over Time")
plt.xlabel("Time [s]")
plt.ylabel("Temperature [K]")
plt.legend(["Upper RTD", "Middle RTD", "Lower RTD"], loc="lower right")

#nexttile
plt.figure()
plt.plot(times, RTD1_resistances, 'r')
plt.plot(times, RTD2_resistances, 'g')
plt.plot(times, RTD3_resistances, 'b')
#xline(RTD1_liquid_time, "r")
#xline(RTD1_gas_time, "m")
#xline(RTD2_liquid_time, "g")
#xline(RTD2_gas_time, "m")
#xline(RTD3_liquid_time, "b")
#xline(RTD3_gas_time, "m")
#title("RTD Resistance Over Time")
#xlabel("Time [s]")
#ylabel("Resistance [Ohm]")
#legend("Upper RTD", "Middle RTD", "Lower RTD")
#grid on

#nexttile
#hold on
#plot(times, -1 .* voltages_measured, "LineWidth", 2)
#xline(RTD1_liquid_time, "r")
#xline(RTD1_gas_time, "m")
#xline(RTD2_liquid_time, "g")
#xline(RTD2_gas_time, "m")
#xline(RTD3_liquid_time, "b")
#xline(RTD3_gas_time, "m")
#title("Voltage Measured Over Time")
#xlabel("Time [s]")
#ylabel("Voltage Measured [V]")
#grid on


#nexttile
#hold on
#plot(times, -1 .* voltage_derivatives, "LineWidth", 2)
#xline(RTD1_liquid_time, "r")
#xline(RTD1_gas_time, "m")
#xline(RTD2_liquid_time, "g")
#xline(RTD2_gas_time, "m")
#xline(RTD3_liquid_time, "b")
#xline(RTD3_gas_time, "m")
#title("Voltage Measured Derivative Over Time")
#xlabel("Time [s]")
#ylabel("Voltage Measured Derivative [V/s]")

#nexttile
#hold on
#plot(times, currents_measured, "LineWidth", 2)
#xline(RTD1_liquid_time, "r")
#xline(RTD1_gas_time, "m")
#xline(RTD2_liquid_time, "g")
#xline(RTD2_gas_time, "m")
#xline(RTD3_liquid_time, "b")
#xline(RTD3_gas_time, "m")
#title("Current Measured Over Time")
#xlabel("Time [s]")
#ylabel("Current Measured [Amp]")

#hold on
#plot(times, current_derivatives, "LineWidth", 2)
#xline(RTD1_liquid_time, "r")
#xline(RTD1_gas_time, "m")
#xline(RTD2_liquid_time, "g")
#xline(RTD2_gas_time, "m")
#xline(RTD3_liquid_time, "b")
#xline(RTD3_gas_time, "m")
#title("Current Measured Derivative Over Time")
#xlabel("Time [s]")
#ylabel("Current Measured Derivative [Amp/s]")
#grid on
plt.show()