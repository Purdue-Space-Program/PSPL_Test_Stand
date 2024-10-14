from press_funcs import *

'''
print(choked_mdot(0.37 * u('in'), 2392 * u('psi'), 288.7060 * u('kelvin'), 'Nitrogen'))
print(c_v_flow(1.3, 4100 * u('psi'), 'Choked', 288.7060 * u('kelvin'), 'Nitrogen'))

#P2 = 4100 * u('psi') - valve_dP(1.3 * u('dimensionless'), 3099.6308 * u('ft^3/min'), 4100 * u('psi'), 288.7060 * u('kelvin'), 'Nitrogen')
#T2 = T_isent_exp(288.7060 * u('kelvin'), )
'''

vdot_fu = prop_vol_flow('Methane', 3 * u('kg/s'), 14.7 * u('psi'))
T2 = T_isent_exp(288.7060 * u('kelvin'), 4000 * u('psi'), 1500 * u('psi'), 'Nitrogen')
scfm_fu = acfm_to_scfm(vdot_fu, 1500 * u('psi'), T2)  
mdot_gas_fu = acfm_to_mdot(vdot_fu, 1500 * u('psi'), T2, 'Nitrogen')

print("------------------------ BRYCE'S NAPKIN MATH - FUEL SIDE 4000 PSI ----------------------------")
print("m_dot to maintain tank pressure:", f"{mdot_gas_fu.to('kg/s'): ~P.4f}")
print("SCFM:", f"{scfm_fu: ~P.4f}")

mdot_flow = scfm_to_mdot((c_v_flow(1.3 * u('dimensionless'), 4000 * u('psi'), 1500 * u('psi'), 288.7060 * u('kelvin'), 'Nitrogen')), 'Nitrogen')
throat_area, throat_dia = choked_flow_dim(mdot_flow, 4000 * u('psi'), 288.7060 * u('kelvin'), 'Nitrogen')

print("mdot_max: ", f"{mdot_flow: ~P.4f}")
print("Area, Diameter: ", f"{throat_area: ~P.4f}", f"{throat_dia: ~P.4f}")

print(get_constants('Nitrogen', 288.7060 * u('kelvin'), 4000 * u('psi')))

T1, P1, T01, P01, M1, a1, v1 = isentropic_flow(288.7060 * u('kelvin'), 4000 * u('psi'), mdot_gas_lox + mdot_gas_fu, tube_dia, GAS_TYPE_LOX)
T2, P2, T02, P02, M2, a2, v2 = fanno_flow(T01, P01, T1, P1, M1, L_12, tube_dia, GAS_TYPE_LOX)