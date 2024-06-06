from CoolProp import CoolProp as cp
from pint import UnitRegistry
import numpy as np
import matplotlib.pyplot as plt
u = UnitRegistry()
Q_ = u.Quantity

T_STD = Q_(60, 'degF').to('kelvin')  # Standard temperature
P_STD = 14.696 * u('psi')  # Standard pressure
R_N2 = 0.2968 * u('kJ/(kg*K)')  # R constant for N2
GAMMA_N2 = 1.47 * u('dimensionless')  # Gamma constant for N2
S_G_GN2 = Q_(0.9669,'dimensionless')
V_TANK = 1738.368 * u('in^3')
V_COPV = 40.0 * u('L')  # COPV volume (big CMS)

P_copv = 15000 * u('psi')

for i in range(30):
    P_mawp = 980 * u('psi')
    P_choke = P_mawp * 1.895 # For N2 only

    v_1_copv_isentropic = V_TANK * (P_mawp / P_copv)**(1 / GAMMA_N2)  # Volume at P_COPV that expands to fill tanks
    m_1_copv_isentropic = (P_copv * v_1_copv_isentropic) / (R_N2 * (T_STD))  # Mass at P_COPV that fills tanks

    m_copv_postburn = (P_choke * V_COPV) / (R_N2 * T_STD)  # Mass required to keep flow choked

    m_total = (m_copv_postburn + m_1_copv_isentropic * 2).to('kg')
    P_copv_init = (m_total * R_N2 * T_STD) / V_COPV

    P_copv = P_copv_init.to('psi')
    print(P_copv_init.to('psi'))

print(m_1_copv_isentropic.to('kg'))