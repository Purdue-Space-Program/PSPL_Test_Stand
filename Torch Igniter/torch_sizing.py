import numpy as np
import CEA_Wrap as CEA
import CoolProp.CoolProp as CP

# Unit conversions
IN2M = 0.0254
PSI2PA = 6895

# Inlet orifices
ox_inlet_diam = 0.081 * IN2M # No. 46 drill
fu_inlet_diam = 0.025 * IN2M # No. 72 drill
fu_cooling_inlet_diam = 0.052 * IN2M # No. 55 drill

ox_inlet_area = (np.pi / 4) * ox_inlet_diam**2
fu_inlet_area = (np.pi / 4) * fu_inlet_diam**2
fu_cooling_inlet_area = (np.pi / 4) * fu_cooling_inlet_diam**2

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

ox_mdot = Cd_ox * ox_inlet_area * inlet_pressure * gamma_ox * (np.sqrt((2 / (gamma_ox + 1))**((gamma_ox + 1) / (gamma_ox - 1))) / np.sqrt(gamma_ox * R_ox * temp_ox))
fu_mdot = Cd_fu * fu_inlet_area * inlet_pressure * gamma_fu * (np.sqrt((2 / (gamma_fu + 1))**((gamma_fu + 1) / (gamma_fu - 1))) / np.sqrt(gamma_fu * R_fu * temp_fu))
fu_cooling_mdot = Cd_fu * fu_cooling_inlet_area * inlet_pressure * gamma_fu * (np.sqrt((2 / (gamma_fu + 1))**((gamma_fu + 1) / (gamma_fu - 1))) / np.sqrt(gamma_fu * R_fu * temp_fu))

core_OF = ox_mdot / fu_mdot
exit_flame_OF = ox_mdot / (fu_mdot + fu_cooling_mdot)

print(f'Core OF: {core_OF:.2f}')
print(f'Exit flame OF: {exit_flame_OF:.2f}')
print(f'Oxygen mass flow rate {ox_mdot * 1000:.2f} g/s')
print(f'Fuel core mass flow rate {fu_mdot * 1000:.2f} g/s')
print(f'Fuel cooling mass flow rate {fu_cooling_mdot * 1000:.2f} g/s')

problem = CEA.RocketProblem(pressure = chamber_pressure, sup = 1, materials = [fuel, oxidizer], o_f = core_OF)
data = problem.run()

cstar = data.cstar
total_mdot = ox_mdot + fu_mdot 

throat_diameter = 0.18 * IN2M

throat_area = (np.pi/4) * throat_diameter**2

chamber_pressure = (cstar * total_mdot) / throat_area
print(chamber_pressure / PSI2PA)
