from press_funcs import *
from CoolProp.CoolProp import PropsSI

P1 = 850 * u('psi')
P2 = P_STD
purgerate = (18/5) * u('m^3/s')
tube_dia = (0.5 - 2*0.065) * u("in")
length_tubing = 200 * u('ft')

T1 = T_STD
h1 = PropsSI("H", "P", P1.to("Pa").magnitude, "T", T1.to("K").magnitude, "CO2")
T2 = PropsSI("T", "P", P2.to("Pa").magnitude, "H", h1, "CO2") * u.kelvin

mdot = acfm_to_mdot(purgerate, P2, T2, "CO2")
print("Mdot:", f"{mdot: ~P.4f}")

T1, P1, T01, P01, M1, a1, v1 = isentropic_flow(288.7060 * u('kelvin'), P1, mdot, tube_dia, 'CO2')
T1 = 250 * u('kelvin')
T2, P2, T02, P02, M2, a2, v2 = fanno_flow(T01, P01, T1, P1, M1, length_tubing, tube_dia, 'CO2')

T2, P2, T02, P02, M2, a2, v2 = run_fanno(T_STD, P1, mdot, tube_dia, length_tubing, "CO2")
#print(P02)
print("T02_fu, P02_fu, T2_fu, P2_fu, M2_fu, speed of sound 2 fu, velocity 2 fu:", f"{T02: ~P.4f}, {P02: ~P.4f}, {T2: ~P.4f}, {P2: ~P.4f}, {M2: ~P.4f}, {a2: ~P.4f}, {v2: ~P.4f},")
#print("dP:", f"{P1: ~P.4f}")