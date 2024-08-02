from press_funcs import *

# Relief worst case scenario (sizing case):
# Regulator breaks; Hydrogen gas flows from 1.3 Cv regulator at a very high pressure drop (past choke) at very low temperatures.
# Gas flow out of relief valves is at very low pressure and room temperature. 
# Find the required Cv

# Hydrogen flow enters a 1.3 Cv regulator at 5000 psi upstream, goes to temp of methane in tank
print('------------------------- HYDROGEN --------------------------')
gas = 'Hydrogen'
flow_in = c_v_flow(1.3 * u('dimensionless'), 5000 * u('psi'), 14.7 * u('psi'), 288 * u('kelvin'), gas)

t_tank = T_isent_exp(288 * u('kelvin'), 5000 * u('psi'), 1500 * u('psi'), gas)
mdot_out = scfm_to_mdot(flow_in, gas)
area, dia = choked_flow_dim(mdot_out, 1500 * u('psi'), t_tank, gas)

# Hydrogen flow exits through a valve from tank MAWP (1500 psi) to atmospheric at room temperature
c_v_out = c_v(flow_in, 1500 * u('psi'), 14.7 * u('psi'), t_tank, gas) * u('dimensionless')
print(f"SCFM in: {flow_in: ~P.4f}, Required Cv out: {c_v_out: ~P.4f}, Required mdot out: {mdot_out: ~P.4f}, Required diameter out: {dia: ~P.4f}")
print(f"Gas temp in tank: {t_tank: ~P.4f}")
print("mdot for 0.25\" orifice:", f"{choked_mdot(0.25 * u('in'), 1500 * u('psi'), t_tank, gas): ~P.4f}")


print('------------------------- NITROGEN --------------------------')
gas = 'Nitrogen'
flow_in = c_v_flow(1.3 * u('dimensionless'), 5000 * u('psi'), 14.7 * u('psi'), 288 * u('kelvin'), gas)

t_tank = T_isent_exp(288 * u('kelvin'), 5000 * u('psi'), 1500 * u('psi'), gas)
mdot_out = scfm_to_mdot(flow_in, gas)
area, dia = choked_flow_dim(mdot_out, 1500 * u('psi'), t_tank, gas)

# Hydrogen flow exits through a valve from tank MAWP (1500 psi) to atmospheric at room temperature
c_v_out = c_v(flow_in, 1500 * u('psi'), 14.7 * u('psi'), t_tank, gas) * u('dimensionless')
print(f"SCFM in: {flow_in: ~P.4f}, Required Cv out: {c_v_out: ~P.4f}, Required mdot out: {mdot_out: ~P.4f}, Required diameter out: {dia: ~P.4f}")
print(f"Gas temp in tank: {t_tank: ~P.4f}")
print("mdot for 0.25\" orifice:", f"{choked_mdot(0.25 * u('in'), 1500 * u('psi'), t_tank, gas): ~P.4f}")