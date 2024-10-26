from press_funcs import *

# Relief worst case scenario (sizing case):
# Regulator breaks; Hydrogen gas flows from 1.3 Cv regulator at a very high pressure drop (past choke) at very low temperatures.
# Gas flow out of relief valves is at very low pressure and room temperature. 
# Find the required Cv

# Hydrogen flow enters a 1.3 Cv regulator at 5000 psi upstream, goes to temp of methane in tank
'''
print('------------------------- HYDROGEN --------------------------')
gas = 'Hydrogen'
flow_in = c_v_flow(1.3 * u('dimensionless'), 5000 * u('psi'), 1500 * u('psi'), 288 * u('kelvin'), gas)

t_tank = T_isent_exp(288 * u('kelvin'), 5000 * u('psi'), 1500 * u('psi'), gas)
mdot_out = scfm_to_mdot(flow_in, gas)
area, dia = choked_flow_dim(mdot_out, 1500 * u('psi'), t_tank, gas)

# Hydrogen flow exits through a valve from tank MAWP (1500 psi) to atmospheric at room temperature
c_v_out = c_v(flow_in, 1500 * u('psi'), 14.7 * u('psi'), t_tank, gas) * u('dimensionless')
print(f"SCFM in: {flow_in: ~P.4f}, Required Cv out: {c_v_out: ~P.4f}, Required mdot out: {mdot_out: ~P.4f}, Required diameter out: {dia: ~P.4f}")
print(f"Gas temp in tank: {t_tank: ~P.4f}")
print("mdot for 0.25\" orifice:", f"{choked_mdot(0.25 * u('in'), 1500 * u('psi'), t_tank, gas): ~P.4f}")
'''
tank_P = 1500 * u('psi')

print('------------------------- NITROGEN --------------------------')
gas = 'Nitrogen'
total_mdot = 0.1 * u('lb/s')
regflow = 1 * u('lb/s')
flag = False
while abs(total_mdot.to('lb/s').m - regflow.to('lb/s').m) > 0.00001:

    #print(".")
    if flag:
        total_mdot += (regflow - total_mdot) / 8
    flag = True

    #print(total_mdot)
    source = (5042*np.exp(-0.1505*total_mdot.to('lb/s').m) + -92.31*np.exp(2.358*total_mdot.to('lb/s').m)) * u('psi')
    flow_in = c_v_flow(1.3 * u('dimensionless'), source, tank_P, 288 * u('kelvin'), gas)
    regflow = scfm_to_mdot(flow_in, "Nitrogen")
    print(source)

t_tank = 288 * u('kelvin')
mdot_out = scfm_to_mdot(flow_in, gas)
area, dia = choked_flow_dim(mdot_out, tank_P, t_tank, gas)

# Hydrogen flow exits through a valve from tank MAWP (1500 psi) to atmospheric at room temperature
c_v_out = c_v(flow_in, tank_P, 14.7 * u('psi'), t_tank, gas) * u('dimensionless')
print(f"SCFM in: {flow_in: ~P.4f}, Required Cv out: {c_v_out: ~P.4f}, Required mdot out: {mdot_out.to('lb/s'): ~P.4f}, Required diameter out: {dia: ~P.4f}, area: {area: ~P.4f}")
print(f"Gas temp in tank: {t_tank: ~P.4f}")
print("mdot for 0.25\" orifice:", f"{choked_mdot(0.25 * u('in'), tank_P, t_tank, gas): ~P.4f}")
print(cp.PropsSI('T', 'P', (20*u('psi')).to('Pa').m, 'Q', 0, "Nitrogen"))