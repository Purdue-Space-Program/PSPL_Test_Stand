from press_funcs import *

P_tank = 600 * u('psi')
P_source = 6000 * u('psi')
mdot_prop = 4 * u('lb/s')


T_tank = T_isent_exp(T_STD, P_source, P_tank, "Nitrogen")
vdot = prop_vol_flow("Ethanol", mdot_prop, 14.7 * u('psi'))
mdot_n2 = acfm_to_mdot(vdot, P_tank, T_tank, "Nitrogen")
scfm_n2 = mdot_to_scfm(mdot_n2, "Nitrogen")
C_v = c_v(scfm_n2, P_source, P_tank, T_STD, "Nitrogen") * u('dimensionless ')

print(f"SCFM: {scfm_n2: ~P.4f}")
print(f"C_v: {C_v: ~P.4f}")