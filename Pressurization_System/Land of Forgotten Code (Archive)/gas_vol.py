from math import *
from press_funcs import *
from reg_sizing import *

# Bottles
# Constant Temp
V_1_isothermal = (P_MAWP * V_TANK / P_BOTTLE)

# Isentropic
V_1_isentropic = V_TANK * (P_MAWP / P_BOTTLE)**(1 / GAMMA_N2)

# Lox Temp
V_1_loxtemp = ((P_MAWP * V_TANK * T_STD) / (Tsat_lox)) / P_BOTTLE

num_isothermal = 2 * V_1_isothermal / V_BOTTLE
num_isentropic = 2 * V_1_isentropic / V_BOTTLE
num_loxtemp = 2 * V_1_loxtemp / V_BOTTLE

print("Initial V required: ")
print("Isothermal: ", V_1_isothermal)
print("Isentropic: ", V_1_isentropic)
print("At LOX temp: ", V_1_loxtemp)
print("# bottles for two tanks (Isothermal, isentropic, LOX temp): ", num_isothermal, " ", num_isentropic, " ", num_loxtemp)

'''
# ACFM = SCFM x (PStandard/PActual) x (TActual/TStandard)
# For 1400 psi minimum flow condition
max_bottle_acfm = scfm_gn2_lox * (P_STD / (P_GN2_MIN)) * (T_GN2_BOTTLE / T_STD)
max_bottle_acfm.ito('ft^3/min')
print("Min bottle acfm: ", f"{max_bottle_acfm: ~P.4f}")

# For 6000 psi maximum flow condition
min_bottle_acfm = scfm_gn2_lox * (P_STD / (P_BOTTLE)) * (T_GN2_BOTTLE / T_STD)
min_bottle_acfm.ito('ft^3/min')
print("Max bottle acfm: ", f"{min_bottle_acfm: ~P.4f}")
'''

# COPV
# Isentropic
v_1_copv_isentropic = V_TANK * (P_MAWP / P_COPV)**(1 / GAMMA_N2)  # Volume at P_COPV that expands to fill tanks
m_1_copv_isentropic = (P_COPV * v_1_copv_isentropic) / (R_N2 * T_STD)  # Mass at P_COPV that fills tanks

m_copv_postburn = (P_CHOKE * V_COPV) / (R_N2 * T_STD)  # Mass required to keep flow choked

m_total = (m_copv_postburn + m_1_copv_isentropic * 2).to('kg')
print("m_total: ", m_total)

P_copv_init = (m_total * R_N2 * T_STD) / V_COPV
print(P_copv_init.to('psi'))

'''
m_copv_postburn = (P_CHOKE * V_COPV) / (R_N2 * T_STD)

m_1_copv_isentropic = (P_COPV * V_1_copv_isentropic) / (R_N2 * T_STD)
v_choke_press_6k = (m_postburn * R_N2 * T_STD) / P_BOTTLE

V_init_total = 2 * (v_choke_press_6k + V_1_isentropic)
print(V_init_total.to('L'))

m_choke_press_copv = (P_CHOKE * V_COPV) / (R_N2 * T_STD)
v_choke_press_copv_maxP = (m_choke_press_copv * R_N2 * T_STD) / P_COPV

V_init_total_copv = 2 * (v_choke_press_copv_maxP + V_1_isentropic)
print(V_init_total_copv.to('L'))
'''