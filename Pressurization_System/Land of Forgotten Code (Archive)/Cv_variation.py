from press_master import run_test
from press_master import u
from press_master import copv_init_press
from press_master import V_TANK
from press_master import V_COPV
from CoolProp import CoolProp as cp
from pint import UnitRegistry
import numpy as np
import matplotlib.pyplot as plt

V_bottle = 43.3 * u('L')
gas_type = 'N2'

plot_lox = []
plot_fu = []
x = []
y = []
P_tank = 1400 * u('psi')
for j in range(int(copv_init_press(P_tank, V_TANK, V_COPV, gas_type).to('psi').m), 6000):
    P_bottle = j * u('psi')
    output = run_test(P_tank, P_bottle, V_bottle, gas_type)
    x.append(P_bottle.to('psi').m)
    plot_lox.append(output['OX']['Cv'])
    plot_fu.append([P_tank, P_bottle, output['OX']['Cv']])
    print(j)

plt.figure(1)
plt.plot(x, plot_lox)
plt.xlabel('Bottle Pressure (psia)')
plt.ylabel('Cv Value')
plt.title('Cv as a function of bottle pressure; constant tank pressure of 1400 psi')

plt.show()