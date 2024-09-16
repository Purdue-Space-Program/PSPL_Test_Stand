from press_funcs import *

'''
mdot_fu, mdot_lox = prop_mdot_from_press_mdot(2.14793 * u('lb/s'), 4500 * u('psi'), 1350 * u('psi'), T_STD, 'Nitrogen', 'Kerosene', 2.2857)

scfm_n2 = mdot_to_scfm(2.14793 * u('lb/s'), 'Nitrogen')
print(f'SCFM N2: {scfm_n2: ~P.4}')
flow = c_v_flow(1.3 * u('dimensionless'), 2820 * u('psi'), 1350 * u('psi'), T_STD, 'Nitrogen')
print(f'Flow N2: {flow: ~P.4}')

area, dia = choked_flow_dim(2.14793 * u('lb/s'), 4000 * u('psi'), T_STD, 'Nitrogen')
print(f'{area: ~P.4}, {dia: ~P.4}')
'''
print('TEST 1:')
total_pressure(4800 * u('psi'), 277.594 * u('kelvin'), 0.041 * u('in'), (0.5 - 0.065 * 2) * u('in'), 'Nitrogen')
print()

print('TEST 2:')
total_pressure(3600 * u('psi'), 277.594 * u('kelvin'), 0.115 * u('in'), (0.5 - 0.065 * 2) * u('in'), 'Nitrogen')
print()

print('TEST 3:')
total_pressure(2460 * u('psi'), 277.594 * u('kelvin'), 0.160 * u('in'), (0.5 - 0.065 * 2) * u('in'), 'Nitrogen')
print()

print('TEST 4:')
total_pressure(1900 * u('psi'), 277.594 * u('kelvin'), 0.190 * u('in'), (0.5 - 0.065 * 2) * u('in'), 'Nitrogen')
print()

mdot_fu, mdot_lox = prop_mdot_from_press_mdot( * u('lb/s'), 4500 * u('psi'), 1350 * u('psi'), T_STD, 'Nitrogen', 'Kerosene', 2.2857)
