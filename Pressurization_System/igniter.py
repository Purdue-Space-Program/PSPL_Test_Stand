from press_funcs import *

h2_source_press = 5000 * u('psi')
gox_source_press = 2400 * u('psi')
gox_regulated_press = 200 * u('psi')
h2_regulated_press = 200 * u('psi')
tube_dia = (0.375 - 0.049*2) * u('in')


mdot_h2 = 0.94 * u('g/s')
mdot_gox = 3 * u('g/s')

scfm_h2 = mdot_to_scfm(mdot_h2, "Hydrogen")
scfm_gox = mdot_to_scfm(mdot_gox, "Oxygen")

c_v_h2 = c_v(scfm_h2, h2_source_press, 14.7 * u('psi'), T_STD, "Hydrogen")
c_v_gox = c_v(scfm_gox, gox_source_press, 14.7 * u('psi'), T_STD, "Oxygen")

area_h2, dia_h2 = choked_flow_dim(mdot_h2, h2_source_press, T_STD, "Hydrogen")
area_gox, dia_gox = choked_flow_dim(mdot_gox, gox_source_press, T_STD, "Oxygen")

ox_tube_length = 10 * u('ft')
T1_gox, P1_gox, T01_gox, P01_gox, M1_gox, a1_gox, v1_gox = isentropic_flow(T_STD, gox_regulated_press, mdot_gox, tube_dia, "Oxygen")
T2_gox, P2_gox, T02_gox, P02_gox, M2_gox, a2_gox, v2_gox = fanno_flow(T01_gox, P01_gox, T1_gox, P1_gox, M1_gox, ox_tube_length, tube_dia, "Oxygen")

h2_tube_length = 10 * u('ft')
T1_h2, P1_h2, T01_h2, P01_h2, M1_h2, a1_h2, v1_h2 = isentropic_flow(T_STD, h2_regulated_press, mdot_h2, tube_dia, "Hydrogen")
T2_h2, P2_h2, T02_h2, P02_h2, M2_h2, a2_h2, v2_h2 = fanno_flow(T01_h2, P01_h2, T1_h2, P1_h2, M1_h2, h2_tube_length, tube_dia, "Hydrogen")

vel_h2, mach_h2 = vel_ideal_gas(T_STD, h2_regulated_press, tube_dia, mdot_h2, "Hydrogen")
vel_gox, mach_gox = vel_ideal_gas(T_STD, gox_regulated_press, tube_dia, mdot_gox, "Oxygen")

#a, b, c, d, e, f, vel_h2 = isentropic_flow(T_STD, h2_source_press, mdot_h2, tube_dia, "Hydrogen")
print("mdot (gox, fu): ", f"           {mdot_gox: ~P.4f},         " f"{mdot_h2: ~P.4f}")
print("SCFM (gox, fu): ", f"           {scfm_gox: ~P.4f},    " f"{scfm_h2: ~P.4f}")
print("Area (gox, fu): ", f"           {area_gox: ~P.4f},        ", f"{area_h2: ~P.4f}")
print("Dia (gox, fu):  ", f"           {dia_gox: ~P.4f},         ", f"{dia_h2: ~P.4f}")
print("C_v (gox, fu):  ", f"           {c_v_gox: <.4f},            ", f"{c_v_h2: <.4f}")
print("Velocity (gox, fu):  ", f"      {v1_gox.to('ft/s'): ~P.4f},       ", f"{v1_h2.to('ft/s'): ~P.4f}")
print("Velocity (gox, fu):  ", f"      {v2_gox.to('ft/s'): ~P.4f},       ", f"{v2_h2.to('ft/s'): ~P.4f}")
print("Velocity (gox, fu):  ", f"      {vel_gox.to('ft/s'): ~P.4f},       ", f"{vel_h2.to('ft/s'): ~P.4f}")