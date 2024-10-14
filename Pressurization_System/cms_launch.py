from press_funcs import *
P_tank_fu = 314.78 * u('psi')
P_tank_lox = 318.07 * u('psi')
V_tank_fu = 1323.69 * u('in^3')
V_tank_lox = 1186.92 * u('in^3')
V_copv = 9.01289 * u('liter')
gas = 'Helium'

print(copv_init_press(P_tank_fu, P_tank_lox, V_tank_fu, V_tank_lox, V_copv, gas))