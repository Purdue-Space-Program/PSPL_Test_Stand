from math import *
print()

V_tank = pi*72*(7.625/2)**2  # in^3 | ID = 7.625"
P_tank = 1500  # psi
P_bottle = 6000 # psi
V_bottle = 43.3 * 61.024 * (1 - 1500/6000) # L -> in^3

gamma_N2 = 1.40

T_lox = 90.15 # K
T_amb = 293 # K

# Constant Temp
V_1_isothermal = (P_tank*V_tank / P_bottle)

# Isentropic
V_1_isentropic = V_tank*(P_tank / P_bottle)**(1/gamma_N2)

# Lox Temp
V_1_loxtemp = ((P_tank*V_tank*T_amb)/(T_lox)) / P_bottle

num_isothermal = V_1_isothermal / V_bottle
num_isentropic = V_1_isentropic / V_bottle
num_loxtemp = V_1_loxtemp / V_bottle

print("Initial V required: ")
print("Isothermal: ", V_1_isothermal / 61.024, " L")
print("Isentropic: ", V_1_isentropic / 61.024, " L")
print("At LOX temp: ", V_1_loxtemp / 61.024, " L")
print("# bottles (Isothermal, isentropic, LOX temp): ", num_isothermal, " ", num_isentropic, " ", num_loxtemp)

print()