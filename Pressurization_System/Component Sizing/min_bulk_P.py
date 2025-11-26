from press_funcs import *
from scipy.interpolate import griddata

c_v = 2.6 * u('dimensionless')
T = T_STD
tube_dia = (0.5 - 2*0.083) * u('in')
tube_length = 225 * u('ft')
gas = "Nitrogen"
P2 = 1055 * u('psi')
P0_list = []
mdot_list = np.linspace(0.1, 1.5, 20)
for m in mdot_list:
    mdot = m * u('lb/s')
    scfm = mdot_to_scfm(mdot, "Nitrogen")
    P1 = upstream_P(c_v, scfm, T, P2, gas)
    P0 = P1 + 2000 * u('psi')
    P1a = 9999 * u('psi')
    while round(abs(P1a - P1), 4) > 0.001 * u('psi'):
        P1a = run_fanno(T, P0, mdot, tube_dia, tube_length, gas)[3]
        P0 -= (P1a - P1) * 0.1

    P0_list.append(P0.to('psi').m)
    print(P0)

plt.figure()
plt.plot(mdot_list, P0_list, color="#daaa00")
plt.xlabel("Total Mass Flow Rate [lb/s]")
plt.ylim(0, 6000)
plt.xlim(0, 1.5)
plt.ylabel("Required Bulk Supply Pressure [psi]")
plt.title("Required Bulk Supply Pressure for 1055 psi Tank Pressure")
plt.grid()
plt.show()

'''
mdot_lox = 0.834 * u('lb/s')
mdot_fu = 0.3 * u('lb/s')
scfm_lox = mdot_to_scfm(mdot_lox, "Nitrogen")
scfm_fu = mdot_to_scfm(mdot_fu, "Nitrogen")
P1 = max(upstream_P(c_v, scfm_lox, T, P2, gas), upstream_P(c_v, scfm_fu, T, P2, gas))
print(f"P1: {P1: ~P.4f}")

P1 = upstream_P(c_v, scfm, T, P2, gas)
flow = c_v_flow(c_v, P1, P2, T, gas)
scfm = mdot_to_scfm(mdot, "Nitrogen")
print(f"Desired flow: {scfm: ~P.4f}, Actual flow: {flow: ~P.4f}")

P0 = P1 + 2000 * u('psi')
P1a = 9999 * u('psi')
while round(abs(P1a - P1), 4) > 0.001 * u('psi'):
    P1a = run_fanno(T, P0, mdot, tube_dia, tube_length, gas)[3]
    P0 -= (P1a - P1) * 0.1
print(f"Regulator inlet pressure: {P1a: ~P.4f} vs target {P1: ~P.4f}")
print(f"Bulk pressure required: {P0: ~P.4f}")

'''
'''
mdot_list = np.linspace(0.3, 1.43, 10)
tank_P_list = np.linspace(30, 1470, 10)
P0_list = []

for m in mdot_list:
    mdot = m * u('lb/s')
    for P in tank_P_list:
        P2 = P * u('psi')
        scfm = mdot_to_scfm(mdot, "Nitrogen")
        P1 = upstream_P(c_v, scfm, T, P2, gas)
        P0 = P1 + 2000 * u('psi')
        P1a = 9999 * u('psi')
        while round(abs(P1a - P1), 4) > 0.001 * u('psi'):
            P1a = run_fanno(T, P0, mdot, tube_dia, tube_length, gas)[3]
            P0 -= (P1a - P1) * 0.1

        P0_list.append([m, P, P0.to('psi').m])
        print(P1)

coordinates = P0_list
# Convert the list to a NumPy array
points = np.array(coordinates)

# Separate x, y, z values
x = points[:, 0]
y = points[:, 1]
z = points[:, 2]

# Create a grid for x and y
x_grid = np.linspace(min(x), max(x), 50)
y_grid = np.linspace(min(y), max(y), 50)
X, Y = np.meshgrid(x_grid, y_grid)

# Interpolate z values over the grid
Z = griddata((x, y), z, (X, Y), method='linear')

# Plot the surface
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='k')

# Add labels and color bar
ax.set_xlabel('Total Mass Flow [lb/s]')
ax.set_ylabel('Propellant Tank Pressure [psi]')
ax.set_zlabel('Required Bulk Supply Pressure')
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.title('Required Bulk Supply Pressure, Varying Mass Flows and Tank Pressures')
plt.show()

plt.figure(2)
nominalP_list = [x for x in P0_list if round(x[1]) == 1055]
plt.plot(nominalP_list[0], nominalP_list[2])
plt.xlabel("Total Mass Flow [lb/s]")
plt.ylabel("Required Bulk Supply Pressure")
plt.title("Bulk Supply Pressure Required at 1055 psi Tank Pressure")
plt.show()

plt.figure(3)
nominalmdot_list = [x for x in P0_list if round(x[0], 1) == 1.4]
plt.plot(nominalP_list[0], nominalP_list[2])
plt.xlabel("Propellant Tank Pressure [psi]")
plt.ylabel("Required Bulk Supply Pressure [psi]")
plt.title("Bulk Supply Pressure Required at 1.4 lb/s flow")
plt.show()
'''