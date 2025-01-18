import numpy as np
import CEA_Wrap as CEA
import CoolProp.CoolProp as CP

# Unit conversions
IN2M = 0.0254
PSI2PA = 6895

# Inlet orifices
Cd_ox = 0.7
Cd_fu = 0.7

# Operating parameters

temp_ox = 290
temp_fu = 290

ox = 'Oxygen'
fu = 'Hydrogen'
oxidizer = CEA.Oxidizer('O2', temp = temp_ox)
fuel = CEA.Fuel('H2', temp = temp_fu)

chamber_pressure = 120 * PSI2PA

inlet_pressure = 1.2 * chamber_pressure

gamma_ox = CP.PropsSI('C', 'T', temp_ox, 'P', inlet_pressure, ox) / CP.PropsSI('O', 'T', temp_ox, 'P', inlet_pressure, ox)
gamma_fu = CP.PropsSI('C', 'T', temp_fu, 'P', inlet_pressure, fu) / CP.PropsSI('O', 'T', temp_fu, 'P', inlet_pressure, fu)
R_ox = CP.PropsSI('C', 'T', temp_ox, 'P', inlet_pressure, ox) - CP.PropsSI('O', 'T', temp_ox, 'P', inlet_pressure, ox)
R_fu = CP.PropsSI('C', 'T', temp_fu, 'P', inlet_pressure, fu) - CP.PropsSI('O', 'T', temp_fu, 'P', inlet_pressure, fu)

ox_set_pressure = inlet_pressure / (2 / (gamma_ox + 1))**(gamma_ox / (gamma_ox - 1))
fu_set_pressure = inlet_pressure / (2 / (gamma_fu + 1))**(gamma_fu / (gamma_fu - 1))

core_OF = 40
exit_flame_OF = 7

problem = CEA.RocketProblem(pressure = chamber_pressure, sup = 1, materials = [fuel, oxidizer], o_f = core_OF, analysis_type='frozen')
data = problem.run()
cstar = data.cstar 

throat_diameter = 0.18 * IN2M

throat_area = (np.pi/4) * throat_diameter**2

mdot_core = (throat_area * chamber_pressure) / cstar
mdot_fu = mdot_core / (core_OF + 1)
mdot_ox = mdot_fu * core_OF
mdot_fu_total = mdot_ox / exit_flame_OF
mdot_fu_cooling = mdot_fu_total - mdot_fu

ox_inlet_area = mdot_ox / (Cd_ox * inlet_pressure * gamma_ox * (np.sqrt((2 / (gamma_ox + 1))**((gamma_ox + 1) / (gamma_ox - 1))) / np.sqrt(gamma_ox * R_ox * temp_ox)))
fu_inlet_area = mdot_fu / (Cd_fu * inlet_pressure * gamma_fu * (np.sqrt((2 / (gamma_fu + 1))**((gamma_fu + 1) / (gamma_fu - 1))) / np.sqrt(gamma_fu * R_fu * temp_fu)))
fu_cooling_inlet_area = mdot_fu_cooling / (Cd_fu * inlet_pressure * gamma_fu * (np.sqrt((2 / (gamma_fu + 1))**((gamma_fu + 1) / (gamma_fu - 1))) / np.sqrt(gamma_fu * R_fu * temp_fu)))

ox_inlet_diam = 2 * np.sqrt(ox_inlet_area / np.pi)
fu_inlet_diam = 2 * np.sqrt(fu_inlet_area / np.pi)
fu_cooling_inlet_diam = 2 * np.sqrt(fu_cooling_inlet_area / np.pi)

print(f'Oxygen inlet diameter {ox_inlet_diam / IN2M:.4f} in')
print(f'Fuel inlet diameter {fu_inlet_diam / IN2M:.4f} in')
print(f'Fuel cooling inlet diameter {fu_cooling_inlet_diam / IN2M:.4f} in')
print(f'Oxygen set pressure {ox_set_pressure / PSI2PA:.2f} psi')
print(f'Fuel set pressure {fu_set_pressure / PSI2PA:.2f} psi')
print(f'Oxygen mass flow rate {mdot_ox * 1000:.2f} g/s')
print(f'Fuel core mass flow rate {mdot_fu * 1000:.2f} g/s')
print(f'Fuel cooling mass flow rate {mdot_fu_cooling * 1000:.2f} g/s')
