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

    rho = cp.PropsSI('D', 'P', P.to('Pa').m, 'Q', 0, prop) * u('kg/m^3')
    T_sat = cp.PropsSI('T', 'P', P.to('Pa').m, 'Q', 0, prop) * u('kelvin')
    v_dot = m_dot / rho  # Volumetric flow out of tank 

    return v_dot.to('ft^3/min')

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
    # R *= Z
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
    return ((P * acfm) / (R * T)).to('kg/s')

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

# region Valves
def c_v(scfm, P1, P2, T, gas):
    '''
    Returns the Cv for a choked compressible flow at the given pressure.
 

    Returns: Cv
    '''

    R, gamma, S_g, choke = get_constants(gas, T, P1)

    if P2 != 'Choked' and P1 < (choke * P2):
        # Sub-critical
        C_v = ((scfm.to('ft^3/min').m) / (22.67 * P1.to('psi').m * (1 - (2 * (P1.to('psi').m - P2.to('psi').m) / (3 * P1.to('psi').m))) * (np.sqrt((P1.to('psi').m - P2.to('psi').m) / (P1.to('psi').m * S_g.to('dimensionless').m * T.to('rankine').m)))))

    else:
        C_v = (scfm.to('ft^3/min').m) / (0.471 * 22.67 * P1.to('psi').m * np.sqrt(1 / (S_g.to('dimensionless').m * T.to('rankine').m)))
    
    return C_v

def c_v_flow(Cv, P1, P2, T1, gas):
    R, gamma, S_g, choke = get_constants(gas, T1, P1)
    if P2 != 'Choked' and P1 < (choke * P2):
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
def copv_init_press(P_tank, V_tank, V_copv, gas):
    '''
    Calculates the initial pressure required to pressurize tanks to a given pressure throughout the burn cycle while keeping flow choked.

    Returns: Initial COPV pressure (psi)
    '''

    R, gamma, S_g, choke = get_constants(gas)

    P_choke = P_tank * choke
    P_copv = 6000 * u('psi')

    for i in range(30):
        v_1_copv = V_tank * (P_tank / P_copv)**(1 / gamma)
        m_press = (P_copv * v_1_copv) / (R * T_STD)

        m_copv_postburn = (P_choke * V_copv) / (R * T_STD)

        m_total = (m_copv_postburn + m_press * 2)
        P_copv_init = (m_total * R * T_STD) / V_copv

        P_copv = P_copv_init.to('psi')

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
   
def T_isent_exp(T1, P1, P2, gas):
    '''
    Calculates gas temperature following isentropic expansion (T2) from P1, T1 to P2, T2 

    Returns: Gas temperature (Kelvin)
    '''
    
    R, gamma, S_g, choke = get_constants(gas, T1, P1)

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

def plot_cv():
    cv = [c_v(2370 * u('ft^3/min'), 6000 * u('psi'), i * u('psi'), 193.35213377632095 * u('K'), 'N2') for i in range(5500)]
    P = [x for x in range(5500)]
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
        message += f"\nERROR: Flow mdot too low at choked condition for given diameter. "
    
    if Cv != False:
        scfm = mdot_to_scfm(mdot, gas)
        if c_v_flow(Cv, P0, "Choked", T0, gas) < scfm:
            message += f"\nERROR: Valve/Regulator Cv too low at given condition; choked and will not allow further. "

    if M >= 1:
        message += f"\nERROR: Flow is choked at given condition. "

    if valve_choked:
        message += f"\nERROR: Valve/Regulator choked at given conditions. "

    return message

'''

# region Main Test
def run_test(P_tank, P_bottle, V_bottle, tube_dia, gas_type):
    ''''''
    Inputs: P_tank, P_bottle, V_bottle, gas_type

    Outputs: m_dot, SCFM, Area, Diameter, C_v
    ''''''
    
    # Oxidizer
    OX = 'Oxygen'
    MDOT_LOX = 8 * u('kg/s')  # LOX mass flow out of tank

    # Fuel
    FU = 'Methane'
    MDOT_FU = 3 * u('kg/s')  # Methane mass flow out of tank 

    # Press Gas
    GAS_TYPE = gas_type

    # Bottle Test Case
    # Change These
    P_TANK = P_tank  # Prop tank pressure
    P_bottle = P_bottle
    V_bottle = V_bottle

    # Leave these
    T_COPV = T_STD
    T_TANK = T_isent_exp(T_STD, P_bottle, P_TANK, GAS_TYPE)
    P_FILL = 40 * u('psi')  # Propellant fill pressure, used for density
    # ------- MAIN RUN - DON'T CHANGE VARS -------

    # Liquid Oxygen
    vdot_lox = prop_vol_flow(OX, MDOT_LOX, P_FILL)
    scfm_lox = scfm(vdot_lox, P_TANK, T_TANK)
    mdot_gas_lox = m_dot(vdot_lox, P_TANK, T_TANK, GAS_TYPE)
    throat_area_lox, throat_dia_lox = throat_calcs(mdot_gas_lox, P_bottle, T_COPV, GAS_TYPE)
    Cv_lox = c_v(scfm_lox, P_bottle, P_TANK, T_COPV, GAS_TYPE)
    vel_lox, mach_lox = vel_ideal_gas(T_COPV, P_bottle, tube_dia, mdot_gas_lox, GAS_TYPE)

    # Fuel
    vdot_fu = prop_vol_flow(FU, MDOT_FU, P_FILL)
    scfm_fu = scfm(vdot_fu, P_TANK, T_TANK)
    mdot_gas_fu = m_dot(vdot_fu, P_TANK, T_TANK, GAS_TYPE)
    throat_area_fu, throat_dia_fu = throat_calcs(mdot_gas_fu, P_bottle, T_COPV, GAS_TYPE)
    Cv_fu = c_v(scfm_fu, P_bottle, P_TANK, T_COPV, GAS_TYPE)
    
    print(f"Simulation of pressurization of {FU} and {OX} test with two {V_TANK.to('L'): ~P.2f} tanks from {V_bottle.to('L'): ~P.2f} COPV.")
    print(f"Maintaining tank pressure of {P_TANK.to('psi'): ~P.2f} from {P_bottle.to('psi'): ~P.2f} initial COPV pressurization.")
    print("Isentropic temp:", f"{T_TANK.to('kelvin'): ~P.2f}")
    print("Mass Flow Rates (lox, fu): ", f"{MDOT_LOX.to('kg/s'): ~P.2f}, {MDOT_FU.to('kg/s'): ~P.2f}")
    print("Line size (lox, fu): ", f"{tube_dia.to('in'): ~P.4f}, {tube_dia.to('in'): ~P.4f}")
    print("Tank pressure and press gas temp: ", f"{P_TANK.to('psi'): ~P.2f}, {T_TANK.to('kelvin'): ~P.2f}")
    print("COPV pressure and press gas temp: ", f"{P_bottle.to('psi'): ~P.2f}, {T_COPV.to('kelvin'): ~P.2f}")

    print()
    print("m_dot (lox, fu):", f"{mdot_gas_lox.to('kg/s'): ~P.4f},       " f"{mdot_gas_fu.to('kg/s'): ~P.4f}")
    print("SCFM (lox, fu): ", f"{scfm_lox: ~P.4f},  " f"{scfm_fu: ~P.4f}")
    print("Area (lox, fu): ", f"{throat_area_lox: ~P.4f},        ", f"{throat_area_fu: ~P.4f}")
    print("Dia (lox, fu):  ", f"{throat_dia_lox: ~P.4f},         ", f"{throat_dia_fu: ~P.4f}")
    print("C_v (lox, fu):  ", f"{Cv_lox: <.4f},            ", f"{Cv_fu: <.4f}")
    #print("Velocity (lox, fu):  ", f"{vel_lox: ~P.4f},         ", f"{vel_fu: ~P.4f}")
    #print("Mach (lox, fu):  ", f"{mach_lox: ~P.4f},         ", f"{mach_fu: ~P.4f}")

    T1, P1, T01, P01, M, a1, v1 = isentropic_flow(T_COPV, P_bottle, mdot_gas_lox + mdot_gas_fu, tube_dia, "N2")
    T02, P02, T2, P2, M2, a2, v2 = fanno_flow(T01, P01, T1, P1, M, 2 * u('ft'), tube_dia, "N2")

    T3_lox, P3_lox, T03_lox, P03_lox, M3_lox, a3_lox, v3_lox = isentropic_flow(T02, P02, mdot_gas_lox, tube_dia, "N2")
    T04_lox, P04_lox, T4_lox, P4_lox, M4_lox, a4_lox, v4_lox = fanno_flow(T03_lox, P03_lox, T3_lox, P3_lox, M3_lox, 7 * u('ft'), tube_dia, "N2")

    T3_fu, P3_fu, T03_fu, P03_fu, M3_fu, a3_fu, v3_fu = isentropic_flow(T02, P02, mdot_gas_fu, tube_dia, "N2")
    T04_fu, P04_fu, T4_fu, P4_fu, M4_fu, a4_fu, v4_fu = fanno_flow(T03_fu, P03_fu, T3_fu, P3_fu, M3_fu, 7 * u('ft'), tube_dia, "N2")
    
    a2.ito('m/s')
    v2.ito('m/s')
    P1.ito('psi')
    P2.ito('psi')
    P3_fu.ito('psi')
    P3_lox.ito('psi')
    P4_fu.ito('psi')
    P4_lox.ito('psi')
    P01.ito('psi')
    P02.ito('psi')
    P03_fu.ito('psi')
    P03_lox.ito('psi')
    P04_fu.ito('psi')
    P04_lox.ito('psi')

    print("T01, P01, T1, P1, M1, speed of sound, velocity:", f"{T01: ~P.4f}, {P01: ~P.4f}, {T1: ~P.4f}, {P1: ~P.4f}, {M: ~P.4f}, {a1: ~P.4f}, {v1: ~P.4f},")
    print("T02, P02, T2, P2, M2, speed of sound, velocity:", f"{T02: ~P.4f}, {P02: ~P.4f}, {T2: ~P.4f}, {P2: ~P.4f}, {M2: ~P.4f}, {a2: ~P.4f}, {v2: ~P.4f},")

    print("T03_lox, P03_lox, T3_lox, P3_lox, M3_lox, speed of sound 3 lox, velocity 3 lox:", f"{T03_lox: ~P.4f}, {P03_lox: ~P.4f}, {T3_lox: ~P.4f}, {P3_lox: ~P.4f}, {M3_lox: ~P.4f}, {a3_lox: ~P.4f}, {v3_lox: ~P.4f},")
    print("T04_lox, P04_lox, T4_lox, P4_lox, M4_lox, speed of sound 4 lox, velocity 4 lox:", f"{T04_lox: ~P.4f}, {P04_lox: ~P.4f}, {T4_lox: ~P.4f}, {P4_lox: ~P.4f}, {M4_lox: ~P.4f}, {a4_lox: ~P.4f}, {v4_lox: ~P.4f},")

    print("T03_fu, P03_fu, T3_fu, P3_fu, M3_fu, speed of sound 3 fu, velocity 3 fu:", f"{T03_fu: ~P.4f}, {P03_fu: ~P.4f}, {T3_fu: ~P.4f}, {P3_fu: ~P.4f}, {M3_fu: ~P.4f}, {a3_fu: ~P.4f}, {v3_fu: ~P.4f},")
    print("T04_fu, P04_fu, T4_fu, P4_fu, M4_fu, speed of sound 4 fu, velocity 4 fu:", f"{T04_fu: ~P.4f}, {P04_fu: ~P.4f}, {T4_fu: ~P.4f}, {P4_fu: ~P.4f}, {M4_fu: ~P.4f}, {a4_fu: ~P.4f}, {v4_fu: ~P.4f},")

    outputs = {
        "OX": {
            "m_dot": MDOT_LOX,
            "V_dot": vdot_lox,
            "SCFM": scfm_lox,
            "Gas m_dot": mdot_gas_lox,
            "Throat Area": throat_area_lox,
            "Throat Diameter": throat_dia_lox,
            "Cv": Cv_lox
        },

        "FU": {
            "m_dot": MDOT_FU,
            "V_dot": vdot_fu,
            "SCFM": scfm_fu,
            "Gas m_dot": mdot_gas_fu,
            "Throat Area": throat_area_fu,
            "Throat Diameter": throat_dia_fu,
            "Cv": Cv_fu
        },

        "Config": {
            "Oxidizer": OX,
            "Fuel": FU,
            "Tank Pressure": P_TANK,
            "Tank Volume": V_TANK,
            "Bottle Pressure": P_bottle,
            "Bottle Volume": V_bottle,
            "Bottle Temp": T_COPV,
            "Press Gas Temp": T_TANK
        }
    }

    return outputs
# endregion

test4 = run_test(1500*u('psi'), 6000*u('psi'), V_BOTTLE*6, (0.5 - 2*0.065) * u("in"), 'N2')
'''