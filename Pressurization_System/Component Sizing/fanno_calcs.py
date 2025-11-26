from press_funcs import *
T = T_STD
P = 4850 * u('psi')
mdot = 0.834 * u('lb/s')
tube_dia = (0.5 - 2*0.083) * u('in')
tube_length = 225 * u('ft')

mdot = np.linspace(0.1, 1.43, 500)
P2s = []
for m in mdot:
    P2s.append(run_fanno(T, P, m * u('lb/s'), tube_dia, tube_length, "Nitrogen")[3].to('psi').m)

model = 5042*np.exp(-0.1505*mdot) - 92.31*np.exp(2.358*mdot)

plt.figure(1)
plt.plot(mdot, P2s, label = "Fanno Flow Model, 225 ft 0.5 x 0.083", color="#daaa00")
plt.plot(mdot, model, label = "Actual Test Data Regression", color="#252526")
plt.axhline(1470, label = "Tank MEOP", color="black", linestyle="dashed")
plt.title("Pressure Drop to Regulator Inlet: Model vs. Test Data")
plt.xlabel("Mass Flow Rate [lb/s]")
plt.ylabel("Regulator Inlet Pressure [psi]")
plt.legend()
plt.grid()
plt.ylim(0, 6000)
plt.show()

