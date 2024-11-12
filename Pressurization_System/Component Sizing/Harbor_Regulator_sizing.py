from press_funcs import *

P_tank = 100 * u('psi')
P_source = 2400 * u('psi')
mdot_prop = 4 * u('lb/s')
T_ambient = 273.15 * u('kelvin')
#tube_dia = (0.5 - 0.083*2) * u('in')

T_tank = T_isent_exp(T_ambient, P_source, P_tank, "Nitrogen")
vdot = prop_vol_flow("Ethanol", mdot_prop, 14.7 * u('psi'))
mdot_n2 = acfm_to_mdot(vdot, P_tank, T_tank, "Nitrogen")
scfm_n2 = acfm_to_scfm(vdot, P_tank, T_tank)
#scfm_n2 = mdot_to_scfm(mdot_n2, "Nitrogen")
C_v = c_v(scfm_n2, P_source, P_tank, T_ambient, "Nitrogen") * u('dimensionless ')

#tube_length = 10 * u('ft')
#T1_gox, P1_gox, T01_gox, P01_gox, M1_gox, a1_gox, v1_gox = isentropic_flow(T_tank, P_tank, mdot_n2, tube_dia, "Nitrogen")
#T2_gox, P2_gox, T02_gox, P02_gox, M2_gox, a2_gox, v2_gox = fanno_flow(T01_gox, P01_gox, T1_gox, P1_gox, M1_gox, tube_length, tube_dia, "Nitrogen")

#tube_length = 20 * u('ft')
#T1_gox, P1_gox, T01_gox, P01_gox, M1_gox, a1_gox, v1_gox = isentropic_flow(T_ambient, P_source, mdot_n2, tube_dia, "Nitrogen")
#T2_gox, P2_gox, T02_gox, P02_gox, M2_gox, a2_gox, v2_gox = fanno_flow(T01_gox, P01_gox, T1_gox, P1_gox, M1_gox, tube_length, tube_dia, "Nitrogen")

print(f"ACFM: {vdot: ~P.4f}")
print(f"mdot_n2: {mdot_n2.to('lb/s'): ~P.4f}")
print(f"T: {T_tank: ~P.4f}")
print(f"SCFM: {scfm_n2: ~P.4f}")
print(f"C_v: {C_v: ~P.4f}")
#print("dT:  ", f"      {(T2_gox - T1_gox).to('kelvin'): ~P.4f}")
#print("dP:  ", f"      {(P02_gox - P01_gox).to('psi'): ~P.4f}")
#print("Velocity:  ", f"      {v2_gox.to('ft/s'): ~P.4f}")
#print("Velocity (gox, fu):  ", f"      {vel_gox.to('ft/s'): ~P.4f},       ", f"{vel_h2.to('ft/s'): ~P.4f}")
#print("Mach:  ", f"          {M2_gox.to('dimensionless'): ~P.4f}")
#plot_cv(scfm_n2, P_source, T_tank, 2300, "Nitrogen")