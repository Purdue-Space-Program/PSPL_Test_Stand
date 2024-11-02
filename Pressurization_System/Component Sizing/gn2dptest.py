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
V1, M1, P_s1, P_d1, P_t1, mdot1 = total_pressure(4800 * u('psi'), 277.594 * u('kelvin'), 0.041 * u('in'), (0.5 - 0.065 * 2) * u('in'), 'Nitrogen')
print()

print('TEST 2:')
V2, M2, P_s2, P_d2, P_t2, mdot2 = total_pressure(3600 * u('psi'), 277.594 * u('kelvin'), 0.115 * u('in'), (0.5 - 0.065 * 2) * u('in'), 'Nitrogen')
print()

print('TEST 3:')
V3, M3, P_s3, P_d3, P_t3, mdot3 = total_pressure(2460 * u('psi'), 277.594 * u('kelvin'), 0.160 * u('in'), (0.5 - 0.065 * 2) * u('in'), 'Nitrogen')
print()

print('TEST 4:')
V4, M4, P_s4, P_d4, P_t4, mdot4 = total_pressure(1880 * u('psi'), 277.594 * u('kelvin'), 0.190 * u('in'), (0.5 - 0.065 * 2) * u('in'), 'Nitrogen')
print()

mdot_fu1, mdot_lox1 = prop_mdot_from_press_mdot(mdot1, P_t1, 1470 * u('psi'), 277.594 * u('kelvin'), 'Nitrogen', 'Kerosene', 2.5)
mdot_fu2, mdot_lox2 = prop_mdot_from_press_mdot(mdot2, P_t2, 1470 * u('psi'), 277.594 * u('kelvin'), 'Nitrogen', 'Kerosene', 2.5)
mdot_fu3, mdot_lox3 = prop_mdot_from_press_mdot(mdot3, P_t3, 1470 * u('psi'), 277.594 * u('kelvin'), 'Nitrogen', 'Kerosene', 2.5)
mdot_fu4, mdot_lox4 = prop_mdot_from_press_mdot(mdot4, P_t4, 1470 * u('psi'), 277.594 * u('kelvin'), 'Nitrogen', 'Kerosene', 2.5)
print(f"MDOT: FU, LOX: {mdot_fu1.to('kg/s'): ~P.4}, {mdot_lox1.to('kg/s'): ~P.4f}")
print(f"MDOT: FU, LOX: {mdot_fu2.to('kg/s'): ~P.4}, {mdot_lox2.to('kg/s'): ~P.4f}")
print(f"MDOT: FU, LOX: {mdot_fu3.to('kg/s'): ~P.4}, {mdot_lox3.to('kg/s'): ~P.4f}")
print(f"MDOT: FU, LOX: {mdot_fu4.to('kg/s'): ~P.4}, {mdot_lox4.to('kg/s'): ~P.4f}")

# Regression equation: -1775.58x^2 + 319.202x + 4791.41 where units are x = mdot in lb/s and y = P_t in psi for Zucrow bulk at 4850 psi
#P1 = 3000 * u('psi') # Total pressure upstream
#P1 = (-1775.58 * mdot.to('lb/s').m ** 2 + 319.202 * mdot.to('lb/s').m + 4791.41) * u('psi')

    
mdot = scfm_to_mdot(c_v_flow(1.3, 1935.9830 * u('psi'), 1350 * u('psi'), 277.594 * u('kelvin'), 'Nitrogen'), 'Nitrogen')
print(f"Mdot: {mdot.to('lb/s'): ~P.4f}")

T1, P1, T01, P01, M1, a1, v1 = isentropic_flow(288.7060 * u('kelvin'), 4850 * u('psi'), 1.3452 * u('lb/s'), (0.5 - 0.083 * 2) * u('in'), 'Nitrogen')
T2, P2, T02, P02, M2, a2, v2 = fanno_flow(T01, P01, T1, P1, M1, 225 * u('ft'), (0.5 - 0.083 * 2) * u('in'), 'Nitrogen')
print("T02_fu, P02_fu, T2_fu, P2_fu, M2_fu, speed of sound 2 fu, velocity 2 fu:", f"{T02: ~P.4f}, {P02: ~P.4f}, {T2: ~P.4f}, {P2: ~P.4f}, {M2: ~P.4f}, {a2: ~P.4f}, {v2: ~P.4f},")

'''
T1, P1, T01, P01, M1, a1, v1 = isentropic_flow(288.7060 * u('kelvin'), 4850 * u('psi'), 1.5850 * u('kg/s'), (0.5 - 0.083 * 2) * u('in'), 'Nitrogen')
T2, P2, T02, P02, M2, a2, v2 = fanno_flow(T01, P01, T1, P1, M1, 300 * u('ft'), (0.5 - 0.065 * 2) * u('in'), 'Nitrogen')
print("T02_fu, P02_fu, T2_fu, P2_fu, M2_fu, speed of sound 2 fu, velocity 2 fu:", f"{T02: ~P.4f}, {P02: ~P.4f}, {T2: ~P.4f}, {P2: ~P.4f}, {M2: ~P.4f}, {a2: ~P.4f}, {v2: ~P.4f},")
'''
mdot = scfm_to_mdot(c_v_flow(1.3, 1500 * u('psi'), 1470 * u('psi'), 277.594 * u('kelvin'), 'Nitrogen'), 'Nitrogen')


x = [scfm_to_mdot(c_v_flow(1.3, x * u('psi'), 1470 * u('psi'), 277.594 * u('kelvin'), 'Nitrogen'), 'Nitrogen').to('lb/s').m for x in np.arange(4850, 1055, -10)]
y1 = [(5042*np.exp(-0.1505*x) + -92.31*np.exp(2.358*x)) for x in np.arange(0, 1.6, 0.01)]
x1 = np.arange(0, 1.6, 0.01)
y = np.arange(4850, 1055, -10)

plt.figure(1)
plt.plot(x, y, label = 'Cv = 1.3 Regulator Mass Flow', color="#252526")
plt.plot(x1, y1, label="Stand Inlet Total Pressure",  color="#daaa00")
plt.grid()
plt.xlim(0, 1.8)
plt.ylim(0, 5000)
plt.xlabel('Mass Flow Rate (lb/s)')
plt.ylabel('Stand Inlet Pressure (psi)')
plt.title('Regulator Mass Flow Rate vs. Inlet Pressure')
plt.legend()
plt.show()
plt.gca().set_aspect(1.0)
print(mdot)

#(5042*np.exp(-0.1505*1.3) + -92.31*np.exp(2.358*1.3))
'''
mdot_fu2, mdot_lox2 = prop_mdot_from_press_mdot(1.3452 * u('lb/s'), 1915.5755 * u('psi'), 1400 * u('psi'), T_STD, 'Nitrogen', 'Kerosene', 2.28)
mdot = scfm_to_mdot(c_v_flow(1.3, 1915.5755 * u('psi'), 1400 * u('psi'), T_STD, 'Nitrogen'), 'Nitrogen').to('lb/s')
print(mdot)
print(f"MDOT: FU, LOX: {mdot_fu2.to('kg/s'): ~P.4}, {mdot_lox2.to('kg/s'): ~P.4f}")
'''

inlet = 2500
for i in range(15):
    if inlet < 1470:
        inlet = 1470
    mdot = scfm_to_mdot(c_v_flow(1.3, inlet * u('psi'), 1470 * u('psi'), 277.594 * u('kelvin'), 'Nitrogen'), 'Nitrogen')
    mdot.ito('lb/s')
    inlet = 5042*np.exp(-0.1505*mdot.m) + -92.31*np.exp(2.358*mdot.m)
    print(f'Mdot: {mdot}')
    print(f'Inlet Pressure: {inlet}')

