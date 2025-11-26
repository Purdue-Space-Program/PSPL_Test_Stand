from press_funcs import *

print('------------------------- CONSTANT PRESSURE --------------------------')
gas = 'Nitrogen'
reg_cv = 1.3 * u('dimensionless')
source_P = 4800 * u('psi')
relief_P = 650 * u('psi')
amb_P = P_STD
source_T = T_STD

flow_in = c_v_flow(reg_cv, source_P, relief_P, source_T, gas)

t_relief = source_T
mdot_out = scfm_to_mdot(flow_in, gas)
area, dia = choked_flow_dim(mdot_out, relief_P, t_relief, gas)

# Hydrogen flow exits through a valve from tank MAWP (1500 psi) to atmospheric at room temperature
c_v_out = c_v(flow_in, relief_P, amb_P, t_relief, gas) * u('dimensionless')
print(f"SCFM in: {flow_in: ~P.4f}, Required Cv out: {c_v_out: ~P.4f}, Required mdot out: {mdot_out: ~P.4f}, Required diameter out: {dia: ~P.4f}")
print(f"Relief temp: {t_relief: ~P.4f}, {amb_P: ~P.4f}")