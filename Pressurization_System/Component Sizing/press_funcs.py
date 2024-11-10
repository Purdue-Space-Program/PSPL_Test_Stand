from CoolProp import CoolProp as cp
from pint import UnitRegistry
import numpy as np
import matplotlib.pyplot as plt
u = UnitRegistry()
Q_ = u.Quantity
from scipy.optimize import fsolve

# Standard Temperature and Pressure
T_STD = Q_(60, 'degF').to('kelvin')  # Standard temperature
P_STD = 14.696 * u('psi')  # Standard pressure

# ------- INPUTS -------
# Big CMS COPV Parameters
PMAX_COPV = 4350 * u('psi')
V_COPV = 40.0 * u('l')

# Bottle Parameters 
PMAX_BOTTLE = 6000 * u('psi')
V_BOTTLE = 43.3 * u('L')

# Prop Tank Parameters
PMAX_TANK = 1400 * u('psi')
V_TANK = 3.40 * u('ft^3')

# ------- FUNCTIONS -------

def prop_vol_flow(prop, m_dot, P):
    '''
    Returns propellant volumetric flow rate out of tanks 

    Returns: Volumetric flow (ft^3/min)
    '''

    if prop == 'Kerosene':
        rho = 820.1 * u('kg/m^3')
    elif prop == 'RP1':
        rho = 810 * u('kg/m^3')
    elif prop == 'Isopropanol':
        rho = 785.027 * u('kg/m^3')
    else:
        rho = cp.PropsSI('D', 'P', P.to('Pa').m, 'Q', 0, prop) * u('kg/m^3')
        #print(rho)
    # T_sat = cp.PropsSI('T', 'P', P.to('Pa').m, 'Q', 0, prop) * u('kelvin')
    v_dot = m_dot / rho  # Volumetric flow out of tank 

    return v_dot.to('ft^3/min')

def prop_mdot_from_press_mdot(mdot_gas, P1, P2, T1, gas, prop, o_f):
    T2 = T_isent_exp(T1, P1, P2, gas)

    rho = cp.PropsSI('D', 'P', P2.to('Pa').m, 'T', T2.to('kelvin').m, gas) * u('kg/m^3')
    vdot_total = mdot_gas / rho
    # print(f"{vdot_total.to('ft^3/s'): ~P.4}")

    if prop == 'Kerosene':
        rho_fu = 820.1 * u('kg/m^3')
    elif prop == 'RP1':
        rho_fu = 810 * u('kg/m^3')
    elif prop == 'Isopropanol':
        rho_fu = 785.027 * u('kg/m^3')
    else:
        rho_fu = cp.PropsSI('D', 'P', P_STD.to('Pa').m, 'Q', 0, prop) * u('kg/m^3')

    rho_lox = cp.PropsSI('D', 'P', P_STD.to('Pa').m, 'Q', 0, 'Oxygen') * u('kg/m^3')

    vdot_fu = 0.000000001 * u('ft^3/s')
    mdot_lox = 0.0000000001
    mdot_fu = 0.0000000001

    while abs(o_f - mdot_lox / mdot_fu) > 0.01:
        mdot_fu = vdot_fu * rho_fu

        vdot_lox = vdot_total - vdot_fu
        mdot_lox = vdot_lox * rho_lox

        of_calc = mdot_lox / mdot_fu

        if of_calc >= o_f:
            vdot_fu = vdot_fu + 0.00001 * u('ft^3/s')
        else:
            vdot_fu = vdot_fu - 0.00001 * u('ft^3/s')

        #print(abs(o_f - mdot_lox / mdot_fu))
    
    #print(mdot_fu.to('kg/s'), mdot_lox.to('kg/s'))
    #print((mdot_lox / mdot_fu).to('dimensionless'))

    return mdot_fu, mdot_lox

# region Compressible Gas Tools
def get_constants(gas, T0=T_STD, P0=P_STD):
    '''
    Returns constant gas properties with units for various gases. 

    Returns: R, gamma, S_g, choke ratio (P1 / P2)
    '''
    if gas == 'N2':
        gas = 'Nitrogen'
    elif gas == 'He':
        gas = 'Helium'
    elif gas == 'H2':
        gas = 'Hydrogen'

    T0.ito_base_units()
    P0.ito_base_units()
    CP = cp.PropsSI('CPMASS', 'P', P0.m, 'T', T0.m, gas) * u('J/kg/K')
    CV = cp.PropsSI('CVMASS', 'P', P0.m, 'T', T0.m, gas) * u('J/kg/K')
    mm_gas = cp.PropsSI('M', 'P', P0.m, 'T', T0.m, gas) * u('kg/mol')
    mm_air = cp.PropsSI('M', 'P', P0.m, 'T', T0.m, 'Air') * u('kg/mol')
    S_g = mm_gas / mm_air
    R = CP - CV
    gamma = CP / CV
    choke = 1 / (2 / (gamma + 1)) ** (gamma / (gamma - 1))
    
    R = 8.314 * u('J/mol/K') / mm_gas

    Z = cp.PropsSI('Z', 'P', P0.m, 'T', T0.m, gas) * u('dimensionless')
    # print("compressibility ratio for ", gas, ": ", Z)
    R *= Z
    T0.ito('kelvin')
    P0.ito('psi')

    return R, gamma, S_g, choke

def acfm_to_scfm(acfm, P_act, T_act):
    '''
    Converts acfm at the given pressure and temperature to scfm

    Returns: SCFM (ft^3/min)
    '''

    return (acfm * (P_act / P_STD) * (T_STD / T_act)).to('ft^3/min')

def mdot_to_scfm(mdot, gas):
    R, gamma, S_g, choke = get_constants(gas, T_STD, P_STD)

    return ((mdot * R * T_STD) / (P_STD)).to('ft^3/min')

def scfm_to_mdot(scfm, gas):
    R, gamma, S_g, choke = get_constants(gas, T_STD, P_STD)

    return ((P_STD * scfm) / (R * T_STD)).to('kg/s')

def acfm_to_mdot(acfm, P, T, gas):
    '''
    Converts a volumetric flow rate at a given pressure and temperature into mass flow rate

    Returns: Mass flow rate (kg/s)
    '''

    R, gamma, S_g, choke = get_constants(gas, T, P)
    #print(R)
    return ((P * acfm) / (R * T)).to('kg/s')

#def vdot_out_scfm_in(acfm_out, P1, P2, T1)

def choked_flow_dim(m_dot, P0, T0, gas):
    '''
    Returns throat area and diameter for a choked flow in terms of m_dot at the given pressure and temperature. 

    Returns: Area (in^3), Diameter (in)
    '''

    R, gamma, S_g, choke = get_constants(gas, T0, P0)

    throat_area = ((m_dot * np.sqrt(T0)) / (P0 * np.sqrt(gamma / R) * np.power(((gamma + 1) / 2), (-(gamma + 1)/(2 * (gamma - 1))))))
    throat_dia = 2 * np.sqrt(throat_area / np.pi)
    return throat_area.to('in^2'), throat_dia.to('in')

def choked_mdot(dia, P0, T0, gas):
    R, gamma, S_g, choke = get_constants(gas, T0, P0)
    A = np.pi / 4 * dia ** 2
    
    return ((A * P0) / np.sqrt(T0) * np.sqrt(gamma / R) * ((gamma + 1) / 2) ** -((gamma + 1) / (2 * (gamma - 1)))).to('kg/s')

    # endregion

def total_pressure(P, T, dia_vent, dia_tube, gas):
    mdot0 = choked_mdot(dia_vent, P, T, gas)
    V0 = flow_to_vel(mdot0, dia_tube, P, T, gas)
    rho = cp.PropsSI('D', 'P', P.to_base_units().m, 'T', T.to_base_units().m, gas) * u('kg/m^3')
    P_t = P + 0.5 * rho * V0 ** 2

    for i in range(5):
        mdot = choked_mdot(dia_vent, P_t, T, gas)
        V = flow_to_vel(mdot, dia_tube, P, T, gas)
        rho = cp.PropsSI('D', 'P', P.to_base_units().m, 'T', T.to_base_units().m, gas) * u('kg/m^3')
        P_t = P + 0.5 * rho * V ** 2
        #print(f"Iteration {i}, P_t = {P_t: ~P.4f}")
    
    R, gamma, S_g, choke = get_constants(gas, T, P_t)
    a = np.sqrt(gamma * R * T)
    print(f"Velocity: {V.to('ft/s'): ~P.4f}")
    print(f"Mach: {(V / a).to('dimensionless'): ~P.4f}")
    print(f"Static Pressure: {(P).to('psi'): ~P.4f}")
    print(f"Dynamic Pressure: {(0.5 * rho * V ** 2).to('psi'): ~P.4f}")
    print(f"Total Pressure: {(P_t).to('psi'): ~P.4f}")
    print(f"Calculated Mass Flow: {mdot.to('lb/s'): ~P.4f}")
    
    return V.to('ft/s'), (V / a).to('dimensionless'), (P).to('psi'), (0.5 * rho * V ** 2).to('psi'), (P_t).to('psi'), mdot.to('lb/s')

def flow_to_vel(mdot, dia, P, T, gas):
    rho = cp.PropsSI('D', 'P', P.to_base_units().m, 'T', T.to_base_units().m, gas) * u('kg/m^3')

    return mdot / (rho * np.pi * (dia / 2) ** 2)

# region Valves
def c_v(scfm, P1, P2, T, gas):
    '''
    Returns the Cv for a choked compressible flow at the given pressure.
 

    Returns: Cv
    '''

    R, gamma, S_g, choke = get_constants(gas, T, P1)
    #print("--------------------------------------------------------------------")
    if P2 == 'choked':
        P1 = choke * P2 + 5

    if P1 < (choke * P2):
        # Sub-critical
        C_v = ((scfm.to('ft^3/min').m) / (22.67 * P1.to('psi').m * (1 - (2 * (P1.to('psi').m - P2.to('psi').m) / (3 * P1.to('psi').m))) * (np.sqrt((P1.to('psi').m - P2.to('psi').m) / (P1.to('psi').m * S_g.to('dimensionless').m * T.to('rankine').m)))))

    else:
        C_v = (scfm.to('ft^3/min').m) / (0.471 * 22.67 * P1.to('psi').m * np.sqrt(1 / (S_g.to('dimensionless').m * T.to('rankine').m)))
    return C_v

def c_v_flow(Cv, P1, P2, T1, gas):
    R, gamma, S_g, choke = get_constants(gas, T1, P1)
    if P2 == 'Choked':
        P2 = (1/choke) * P1

    if P1 < (choke * P2):
        # Sub-critical
        dP = P1 - P2
        N2 = 22.67
        q = N2 * Cv * P1.to('psi').m * (1 - (2 * dP.to('psi').m) / (3 * P1.to('psi').m)) * np.sqrt((dP.to('psi').m) / (P1.to('psi').m * S_g.to('dimensionless').m * T1.to('rankine').m))
        q *= u('ft^3/min')

    else:
        N2 = 22.67

        q = 0.471 * N2 * Cv * P1.to('psi').m * np.sqrt(1 / (S_g.to('dimensionless').m * T1.to('rankine').m))
        q *= u('ft^3/min')

    return q

def valve_dP(c_v, scfm, P0, T, gas):
    R, gamma, S_g, choke = get_constants(gas, T, P0)
    choked = False

    c_v = np.float64(c_v.m)
    scfm = np.float64(scfm.to('ft^3/min').m)
    P0 = np.float64(P0.to('psi').m)
    T = np.float64(T.to('rankine').m)
    N2 = 22.67
    eqn = lambda dP: -scfm + N2 * c_v * P0 * (1 - (2 * dP)/(3 * P0)) * np.sqrt((dP)/(P0 * S_g * T))
    dP = float(np.abs(fsolve(eqn, 0.001)[0])) * u('psi')
    P0 *= u('psi')
    c_v *= u('dimensionless')

    if P0 - dP < P0 * 1 / choke:
        choked = True
        # print("Choked!!")
        dP = P0 - (P0 * 1 / choke)

    #print(f'Pressure drop over valve of Cv {c_v.to("dimensionless"): ~P.2f}: {-dP.to("psi"): ~P.2f}')
    return -dP, choked
# endregion

# region COPV calcs
def copv_init_press(P_tank_fu, P_tank_lox, V_tank_fu, V_tank_lox, V_copv, gas):
    '''
    Calculates the initial pressure required to pressurize tanks to a given pressure throughout the burn cycle while keeping flow choked.

    Returns: Initial COPV pressure (psi)
    '''

    R, gamma, S_g, choke = get_constants(gas)
    P_choke = max(P_tank_fu * choke, P_tank_lox * choke)
    P_copv = 6000 * u('psi')

    for i in range(30):
        v_1_copv_fu = V_tank_fu * (P_tank_fu / P_copv)**(1 / gamma)
        m_press_fu = (P_copv * v_1_copv_fu) / (R * T_STD)

        v_1_copv_lox = V_tank_lox * (P_tank_lox / P_copv)**(1 / gamma)
        m_press_lox = (P_copv * v_1_copv_lox) / (R * T_STD)

        m_copv_postburn = (P_choke * V_copv) / (R * T_STD)

        m_total = (m_copv_postburn + m_press_fu + m_press_lox)
        P_copv_init = (m_total * R * T_STD) / V_copv

        P_copv = P_copv_init.to('psi')
        print(P_copv)

    return P_copv.to('psi')

def copv_max_tank_press(P_max_copv, V_tank, V_copv, gas):
    '''
    Calculates the maximum tank pressure achievable with a given COPV throughout the burn cycle while keeping flow choked. 

    Returns: Initial tank pressure (psi)
    '''

    R, gamma, S_g, choke = get_constants(gas)
    max_press = 0

    for i in range(1, 1600):
        P_tank = i * u('psi')
        if copv_init_press(P_tank, V_tank, V_copv, gas).to('psi').m <= P_max_copv.to('psi').m:
            max_press = P_tank

        else:
            return max_press.to('psi')
# endregion

def reqd_volume(T1, P1, P2, V2, gas):
    T2 = T_isent_exp(T1, P1, P2, gas)

    return P2 * V2 / T2 * T1 / P1

def T_isent_exp(T1, P1, P2, gas):
    '''
    Calculates gas temperature following isentropic expansion (T2) from P1, T1 to P2, T2 

    Returns: Gas temperature (Kelvin)
    '''
    
    R, gamma, S_g, choke = get_constants(gas, T1, P1)
    print(gamma)
    return (T1 * (P2 / P1)**(1 - 1 / gamma)).to('kelvin')

# region Fanno Flow
def isentropic_flow(T0, P0, m_dot, dia, gas):
    '''
    Calculates isentropic flow given required mass flow and stagnation conditions from T0 (stagnation, bottle), P0 to T1 and P1.

    Returns: T1, P1, T0, P0, M, a, v (Temperature post-expansion, pressure post-expansion, mach, speed of sound, velocity)
    '''
    area = np.pi * (dia / 2) ** 2  # Flow area (tube)
    R, gamma, S_g, choke = get_constants(gas, T0, P0)  # Gas constants

    if gas == 'N2':
        gas = 'Nitrogen'
    elif gas == 'He':
        gas = 'Helium'
    elif gas == 'H2':
        gas = 'Hydrogen'

    p0 = cp.PropsSI('D', 'P', P0.to_base_units().m, 'T', T0.to_base_units().m, gas) * u('kg/m^3')

    m_dot.ito_base_units()
    area.ito_base_units()
    gamma.ito_base_units()
    R.ito_base_units()
    p0.ito_base_units()
    T0.ito_base_units()

    eqn = lambda M : -np.float64(m_dot.m) + np.float64(p0.m) * (1 + (np.float64(gamma.m) - 1) / 2 * M ** 2) ** (-1 / (np.float64(gamma.m) - 1)) * np.float64(area.m) * np.sqrt(np.float64(gamma.m) * np.float64(R.m) * np.float64(T0.m) * (1 + (np.float64(gamma.m) - 1) / 2 * M ** 2) ** -1) * M
    M = fsolve(eqn, 0.01)
    M = float(np.abs(M[0])) * u('dimensionless')

    # Isentropic flow equations: https://www.grc.nasa.gov/www/k-12/airplane/Images/isentrop.gif
    T1 = T0 * (1 + (gamma - 1) / 2 * M ** 2) ** -1  # Flow temperature
    P1 = P0 * (1 + (gamma - 1) / 2 * M ** 2) ** (-gamma / (gamma - 1))  # Flow pressure
    a = np.sqrt(gamma * R * T1)  # Speed of sound
    vel = a * M  # Speed of gas

    m_dot.ito('kg/s')
    area.ito('in^2')
    p0.ito('kg/m^3')

    return T1.to('kelvin'), P1.to('psi'), T0.to('kelvin'), P0.to('psi'), M.to('dimensionless'), a.to('ft/s'), vel.to('ft/s')

def fanno_parameter(gamma, M):
    return ((1 - M**2) / (gamma * M**2) + (gamma + 1) / (2 * gamma) * np.log(((gamma + 1) * M**2) / (2 * (1 + (gamma - 1) / 2 * M**2)))).to('dimensionless')

def mach_from_fanno(gamma, fanno):
    fanno.ito_base_units()
    gamma.ito_base_units()

    gamma = np.float64(gamma.m)
    fanno = np.float64(fanno.m)

    eqn = lambda M : -fanno + (1 - M**2) / (gamma * M**2) + (gamma + 1) / (2 * gamma) * np.log(((gamma + 1) * M**2) / (2 * (1 + (gamma - 1) / 2 * M**2)))
    M = fsolve(eqn, 0.01)
    M = float(np.abs(M[0])) * u('dimensionless')

    return M

def fanno_flow(T0, P0, T1, P1, M1, L, dia, gas_type):
    R, gamma, S_g, choke = get_constants(gas_type, T0, P0)  # Gas constants

    if gas_type == "He":
        gas_type = "Helium"
    elif gas_type == "N2":
        gas_type = "Nitrogen"
    elif gas_type == "H2":
        gas_type = "Hydrogen"

    dyn_visc = cp.PropsSI('V', "P", P1.to('Pa').m, "T", T1.to('kelvin').m, gas_type) * u('Pa') * u('second')
    rho = cp.PropsSI('D', 'P', P1.to_base_units().m, 'T', T1.to_base_units().m, gas_type) * u('kg/m^3')
    a = np.sqrt(gamma * R * T1)
    vel = M1 * a
    Re = (rho * vel * dia / dyn_visc).to('dimensionless')

    # Fanning friction factor using equations from https://en.wikipedia.org/wiki/Fanning_friction_factor#Hydraulically_smooth_piping
    # Assumes hydraulically smooth piping
    '''
    if 2100 < Re < 10**5:
        fanning = 0.0791 / Re**0.25
    elif 10**4 < Re < 10**7:
        fanning = 0.0014 + 0.125 / Re**0.32
    else:
        print("Reynolds number for Fanning friction factor out of range!!")
    '''

    fanning = 0.0014 + 0.125 / Re**0.32
    
    M1_fanno = fanno_parameter(gamma, M1)
    L_fanno = (4 * fanning * L / dia).to('dimensionless')
    #print("choke length for m1", (M1_fanno * dia / (4 * fanning)))

    M2_fanno = M1_fanno - L_fanno
    M2 = mach_from_fanno(gamma, M2_fanno)
    #if M2 > 0.97: print('FLOW MAY BE CHOKED IN TUBE!!!')

    P01_P0star = 1 / M1 * ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M1 ** 2)) ** ((gamma + 1)/(2 * (gamma - 1)))
    P02_P0star = 1 / M2 * ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M2 ** 2)) ** ((gamma + 1)/(2 * (gamma - 1)))
    P02 = P0 * P02_P0star / P01_P0star 

    P1_Pstar = 1 / (M1 * np.sqrt((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M1**2)))
    P2_Pstar = 1 / (M2 * np.sqrt((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M2**2)))
    P2 = P1 * P2_Pstar / P1_Pstar
    
    T1_Tstar = 1 / ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M1 ** 2))
    T2_Tstar = 1 / ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M2 ** 2))
    T2 = T1 * T2_Tstar / T1_Tstar

    a = np.sqrt(gamma * R * T2)  # Speed of sound
    a.ito('m/s')
    v = M2 * a
    v.ito('m/s')

    #print(f"Pressure drop over {L.to('ft'): ~P.2f} long section: {(P02-P0).to('psi'): ~P.2f}")
    return T2.to('kelvin'), P2.to('psi'), T0.to('kelvin'), P02.to('psi'), M2.to('dimensionless'), a.to('ft/s'), v.to('ft/s')  

def vel_ideal_gas(T, P, dia, m_dot, gas_type):
    '''
    Finds the velocity of an ideal gas using conservaiton of mass. Depreciated in favor of idsentropic model.
    '''
    area = np.pi * (dia / 2) ** 2
    R, gamma, S_g, choke = get_constants(gas_type, T, P)
    '''
    print()
    print(f"getting ideal gas velocity for R={R}, T={T}, mdot={m_dot.to('kg/s')}, P={P.to('psi')}, area={area.to('in^2')}, gamma={gamma}")
    print()
    '''

    vel = (R * T * m_dot) / (P * area)
    c = np.sqrt(gamma * R * T)
    mach = vel / c

    return vel.to('m/s'), mach.to('dimensionless')
# endregion

# region Graphical Functions
def plot_copv():
    copv_pressures = []
    tank_pressures = range(1, 1400)

    for i in range(1, 1400):
        P_tank = i * u('psi')
        press = copv_init_press(P_tank, V_TANK, V_COPV, 'N2').to('psi').m
        copv_pressures.append(press)
        if press <= PMAX_COPV.to('psi').m:
            max_press = P_tank

    plt.figure(1)
    plt.plot(tank_pressures, copv_pressures, label = 'Initial COPV Pressure')
    plt.axhline(y = PMAX_COPV.to('psi').m, color = 'r', linestyle = '--', label = 'COPV MEOP')
    plt.xlabel('Tank Pressure (psia)')
    plt.ylabel('COPV Starting Pressure (psia)')
    plt.title('Initial COPV pressure as a function of tank pressure')
    plt.legend()

    plt.plot(max_press, PMAX_COPV, color = 'yellow')
    print(max_press)
    #plt.annotate(f"({max_press.to('psi'): ~P.2f}, {PMAX_COPV.to('psi'): ~P.2f})", (max_press, PMAX_COPV))
    plt.show()

def plot_cv(scfm, P_1, T_tank, max_P, gas):
    cv = [c_v(scfm, P_1, i * u('psi'), T_tank, gas) for i in range(max_P)]
    P = [x for x in range(max_P)]
    plt.figure(1)
    plt.plot(P, cv, label = 'C_v Value')
    plt.xlabel('Downstream Pressure')
    plt.ylabel('C_v Value')
    plt.title('C_v Value as a function of downstream pressure')
    plt.legend()
    plt.show()


def plot_cv_flow():
    flow = [c_v_flow(1.3, 6000 * u('psi'), i * u('psi'), 298 * u('K'), 'N2').to('ft^3/min').m for i in range(6000)]
    P = [x for x in range(6000)]
    plt.figure(1)
    plt.plot(flow, P, label = '6000 psi Inlet')
    ax = plt.gca()
    ax.set_xlim([0, 4000])
    ax.set_ylim([0, 6000])
    plt.xticks(np.arange(0, 4001, 400))
    plt.yticks(np.arange(0, 6001, 750))
    plt.grid()
    plt.ylabel('OUTLET PRESSURE (psi)')
    plt.xlabel('FLOW RATE (SCFM) Nitrogen')
    plt.title('Flow Curve - Cv = 1.3 Tescom 26-1100 3/8" Orifice')
    plt.legend()
    plt.show()
# endregion

def check(station, mdot, dia, gas, Cv = False, valve_choked = False):
    T = station[0]
    P = station[1]
    T0 = station[2]
    P0 = station[3]
    M = station[4]
    a = station[5]
    v = station[6]

    message = ''

    if choked_mdot(dia, P0, T0, gas) < mdot:
        message += f"\nERROR: Flow mdot too low at choked condition for given diameter. Choked: {P0: ~P.4f}, {T0: ~P.4f}, {gas}, {choked_mdot(dia, P0, T0, gas): ~P.4f}, Actual: {mdot: ~P.4f}"
    
    if Cv != False:
        scfm = mdot_to_scfm(mdot, gas)
        if c_v_flow(Cv, P0, "Choked", T0, gas) < scfm:
            message += f"\nERROR: Valve/Regulator Cv too low at given condition; choked and will not allow further. "

    if M >= 0.995:
        message += f"\nERROR: Flow is choked at given condition. "

    if valve_choked:
        message += f"\nERROR: Valve/Regulator choked at given conditions. "

    return message