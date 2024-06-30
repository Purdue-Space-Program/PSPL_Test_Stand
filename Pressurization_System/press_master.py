from CoolProp import CoolProp as cp
from pint import UnitRegistry
import numpy as np
import matplotlib.pyplot as plt
u = UnitRegistry()
Q_ = u.Quantity
from sympy.solvers import nsolve
from scipy.optimize import fsolve
from sympy import Symbol
from sympy import log
from sympy import sqrt

# ------- CONSTANTS -------
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
# V_TANK = 1.009 * u('ft^3')
# V_TANK = np.pi * (7.625 / 2 * u('in')) ** 2 * 9 * u('ft')
V_TANK = 3 * u('ft^3')


# region Functions
# ------- FUNCTIONS -------
# Functions
def get_constants(gas):
    '''
    Returns constant gas properties with units for various gases. 

    Returns: R, gamma, S_g, choke ratio
    '''

    if gas == 'N2':
        R_N2 = 0.2968 * u('kJ/(kg*K)')  # R constant for N2
        GAMMA_N2 = 1.47 * u('dimensionless')  # Gamma constant for N2
        S_G_N2 = Q_(0.9669,'dimensionless')
        CHOKE_N2 = 1.895

        return R_N2, GAMMA_N2, S_G_N2, CHOKE_N2
    elif gas == 'He':
        R_HE = 0.0 * u('kJ/(kg*K)')
        GAMMA_HE = 0.0 * u('dimensionless')  # Gamma constant for N2
        S_G_HE = Q_(0.0,'dimensionless')
        CHOKE_HE = 0

        return R_HE, GAMMA_HE, S_G_HE, CHOKE_HE

def prop_vol_flow(prop, m_dot, P):
    '''
    Returns propellant volumetric flow rate out of tanks 

    Returns: Volumetric flow (ft^3/min)
    '''

    rho = cp.PropsSI('D', 'P', P.to('Pa').m, 'Q', 0, prop) * u('kg/m^3')
    T_sat = cp.PropsSI('T', 'P', P.to('Pa').m, 'Q', 0, prop) * u('kelvin')
    v_dot = m_dot / rho  # Volumetric flow out of tank 

    return v_dot.to('ft^3/min')

def scfm(acfm, P_act, T_act):
    '''
    Converts acfm at the given pressure and temperature to scfm

    Returns: SCFM (ft^3/min)
    '''

    return (acfm * (P_act / P_STD) * (T_STD / T_act)).to('ft^3/min')

def m_dot(acfm, P, T, gas):
    '''
    Converts a volumetric flow rate at a given pressure and temperature into mass flow rate

    Returns: Mass flow rate (kg/s)
    '''

    R, gamma, S_g, choke = get_constants(gas)

    return ((P * acfm) / (R * T)).to('kg/s')

def throat_calcs(m_dot, P, T, gas):
    '''
    Returns throat area and diameter for a choked flow at the given pressure and temperature. 

    Returns: Area (in^3), Diameter (in)
    '''

    R, gamma, S_g, choke = get_constants(gas)

    throat_area = ((m_dot * np.sqrt(T)) / (P * np.sqrt(gamma / R) * np.power(((gamma + 1) / 2), (-(gamma + 1)/(2 * (gamma - 1))))))
    throat_dia = 2 * np.sqrt(throat_area / np.pi)
    return throat_area.to('in^2'), throat_dia.to('in')

def c_v(scfm, P, gas):
    '''
    Returns the Cv for a choked compressible flow at the given pressure.
    UPDATES NEEDED: Need non-choked functionality

    Returns: Cv
    '''

    R, gamma, S_g, choke = get_constants(gas)
    C_v = (scfm.to('ft^3/min').m * 2 * np.sqrt(S_g.m)) / P.to('psi').m  # Note:units in psia and SCFM
    return C_v

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
        
def T_isent_exp(T1, P1, P2, gas):
    '''
    Calculates gas temperature following isentropic expansion (T2) from P1, T1 to P2, T2 

    Returns: Gas temperature (Kelvin)
    '''
    
    R, gamma, S_g, choke = get_constants(gas)

    return (T1 * (P2 / P1)**(1 - 1 / gamma)).to('kelvin')

def isentropic_flow(T0, P0, m_dot, dia, gas_type):
    '''
    Calculates orifice isentropic flow given required mass flow and stagnation conditions from T0 (stagnation, bottle), P0 to T1 and P1.

    Returns: T1, P1, T0, P0, M, a, v (Temperature post-expansion, pressure post-expansion, mach, speed of sound, velocity)
    '''
    area = np.pi * (dia / 2) ** 2  # Flow area (tube)
    R, gamma, S_g, choke = get_constants(gas_type)  # Gas constants

    p0 = P0 / (R * T0)  # Stagnation density

    m_dot.ito_base_units()
    area.ito_base_units()
    gamma.ito_base_units()
    R.ito_base_units()
    p0.ito_base_units()
    T0.ito_base_units()

    M = Symbol('M')
    #eqn = -m_dot + p0 * (1 + (gamma - 1) / 2 * M ** 2) ** (-1 / (gamma - 1)) * area * np.sqrt(gamma * R * T0 * (1 + (gamma - 1) / 2 * M ** 2) ** -1)

    #eqn = -m_dot + (area * P0) / np.sqrt(T0) * np.sqrt(gamma / R) * M * (1 + ((gamma - 1) / 2) * M ** 2) ** (-(gamma + 1) / (2 * (gamma - 1))) #  Compressible orifice flow equation: https://www.grc.nasa.gov/www/k-12/airplane/Images/mflchk.gif
    eqn = -m_dot.m + p0.m * (1 + (gamma.m - 1) / 2 * M ** 2) ** (-1 / (gamma.m - 1)) * area.m * sqrt(gamma.m * R.m * T0.m * (1 + (gamma.m - 1) / 2 * M ** 2) ** -1) * M
    M = float(nsolve(eqn, M, 0.01)) * u('dimensionless') # Solve for mach number following flow 

    # Isentropic flow equations: https://www.grc.nasa.gov/www/k-12/airplane/Images/isentrop.gif
    T1 = T0 * (1 + (gamma - 1) / 2 * M ** 2) ** -1  # Flow temperature
    P1 = P0 * (1 + (gamma - 1) / 2 * M ** 2) ** (-gamma / (gamma - 1))  # Flow pressure
    a = np.sqrt(gamma * R * T1)  # Speed of sound
    vel = a * M  # Speed of gas

    return T1, P1, T0, P0, M, a, vel

def fanno_parameter(gamma, M):
    return (1 - M**2) / (gamma * M**2) + (gamma + 1) / (2 * gamma) * log(((gamma + 1) * M**2) / (2 * (1 + (gamma - 1) / 2 * M**2)))

def mach_from_fanno(gamma, fanno):
    fanno.ito_base_units()
    gamma.ito_base_units()

    gamma = np.float64(gamma.m)
    fanno = np.float64(fanno.m)

    eqn = lambda M : -fanno + (1 - M**2) / (gamma * M**2) + (gamma + 1) / (2 * gamma) * np.log(((gamma + 1) * M**2) / (2 * (1 + (gamma - 1) / 2 * M**2)))
    M = fsolve(eqn, 0.01)
    M = float(M[0]) * u('dimensionless')

    return M

def fanno_flow(T0, P0, T1, P1, M1, L, m_dot, dia, gas_type):
    R, gamma, S_g, choke = get_constants(gas_type)  # Gas constants

    if gas_type == "He":
        gas_type = "Helium"
    elif gas_type == "N2":
        gas_type = "Nitrogen"

    dyn_visc = cp.PropsSI('V', "P", P1.to('Pa').m, "T", T1.to('kelvin').m, gas_type) * u('Pa') * u('second')
    Re = ((m_dot * dia) / (dyn_visc * (np.pi * dia ** 2 / 4))).to('dimensionless')

    # Fanning friction factor using equations from https://en.wikipedia.org/wiki/Fanning_friction_factor#Hydraulically_smooth_piping
    # Assumes hydraulically smooth piping
    if 2100 < Re < 10**5:
        fanning = 0.0791 / Re**0.25
    elif 10**4 < Re < 10**7:
        fanning = 0.0014 + 0.125 / Re**0.32
    else:
        print("Fanning error!!")

    M1_fanno = fanno_parameter(gamma, M1)
    L_fanno = 4 * fanning * L / dia

    M2_fanno = M1_fanno - L_fanno
    M2 = mach_from_fanno(gamma, M2_fanno)

    print(M2)

    P1_Pstar = 1 / M1 * np.sqrt((gamma + 1) / 2 * (1 + (gamma - 1) / 2 * M1 ** 2))
    P2_Pstar = 1 / M2 * np.sqrt((gamma + 1) / 2 * (1 + (gamma - 1) / 2 * M2 ** 2))
    P2 = P1 * 1 / P1_Pstar * P2_Pstar

    T1_Tstar = (gamma + 1) / (2 * (1 + (gamma - 1) / 2) * M1 ** 2)
    T2_Tstar = (gamma + 1) / (2 * (1 + (gamma - 1) / 2) * M2 ** 2)
    T2 = T1 * 1 / T1_Tstar * T2_Tstar

    a = np.sqrt(gamma * R * T2)  # Speed of sound
    v = M2 * a

    print(f"Pressure drop over {L.to('ft'): ~P.2f} long section: {(P2-P1).to('psi'): ~P.2f}")
    return T0, P0, T2, P2, M2, a, v
    

def vel_ideal_gas(T, P, dia, m_dot, gas_type):
    area = np.pi * (dia / 2) ** 2
    R, gamma, S_g, choke = get_constants(gas_type)
    '''
    print()
    print(f"getting ideal gas velocity for R={R}, T={T}, mdot={m_dot.to('kg/s')}, P={P.to('psi')}, area={area.to('in^2')}, gamma={gamma}")
    print()
    '''

    vel = (R * T * m_dot) / (P * area)
    c = np.sqrt(gamma * R * T)
    mach = vel / c

    return vel.to('m/s'), mach.to('dimensionless')

# Graphical Functions
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


def run_test(P_tank, P_bottle, V_bottle, tube_dia, gas_type):
    '''
    Inputs: P_tank, P_bottle, V_bottle, gas_type

    Outputs: m_dot, SCFM, Area, Diameter, C_v
    '''
    
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
    Cv_lox = c_v(scfm_lox, P_bottle, GAS_TYPE)
    vel_lox, mach_lox = vel_ideal_gas(T_COPV, P_bottle, tube_dia, mdot_gas_lox, GAS_TYPE)

    # Fuel
    vdot_fu = prop_vol_flow(FU, MDOT_FU, P_FILL)
    scfm_fu = scfm(vdot_fu, P_TANK, T_TANK)
    mdot_gas_fu = m_dot(vdot_fu, P_TANK, T_TANK, GAS_TYPE)
    throat_area_fu, throat_dia_fu = throat_calcs(mdot_gas_fu, P_bottle, T_COPV, GAS_TYPE)
    Cv_fu = c_v(scfm_fu, P_bottle, GAS_TYPE)
    # vel_fu, mach_fu = vel_ideal_gas(T_COPV, P_bottle, tube_dia, mdot_gas_fu, GAS_TYPE)

    # vel_combined, mach_combined = vel_ideal_gas(T_COPV, P_bottle, tube_dia, mdot_gas_fu + mdot_gas_lox, GAS_TYPE)
    
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

    T1, P1, T0, P0, M, a1, v1 = isentropic_flow(T_COPV, P_bottle, mdot_gas_lox, tube_dia, "N2")
    T0, P0, T2, P2, M2, a2, v2 = fanno_flow(T0, P0, T1, P1, M, 5000 * u('ft'), mdot_gas_lox, tube_dia, "N2")
    a2.ito('m/s')
    v2.ito('m/s')
    print("T1, P1, T0, P0, M, speed of sound, velocity:", f"{T1: ~P.4f}, {P1: ~P.4f}, {T0: ~P.4f}, {P0: ~P.4f}, {M: ~P.4f}, {a1: ~P.4f}, {v1: ~P.4f},")
    print("T0, P0, T2, P2, M2, speed of sound, velocity:", f"{T0: ~P.4f}, {P0: ~P.4f}, {T2: ~P.4f}, {P2: ~P.4f}, {M2: ~P.4f}, {a2: ~P.4f}, {v2: ~P.4f},")

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

# test1 = run_test(980*u('psi'), copv_init_press(980*u('psi'), V_TANK, V_COPV, 'N2'), V_COPV, 'N2')
# test2 = run_test(1400*u('psi'), copv_init_press(1400*u('psi'), V_TANK, V_BOTTLE, 'N2'), V_BOTTLE, 'N2')
# test3 = run_test(1400*u('psi'), 6000*u('psi'), V_BOTTLE, 'N2')
# test4 = run_test(1049*u('psi'), copv_init_press(980*u('psi'), V_TANK, V_COPV, 'N2'), V_COPV, 'N2')
test4 = run_test(1500*u('psi'), 6000*u('psi'), V_BOTTLE*6, (0.5 - 2*0.065) * u("in"), 'N2')
#test5 = run_test(1500*u('psi'), 6000 * u('psi'), V_TANK, V_BOTTLE*6, 'N2', V_BOTTLE*6, (3/8 - 2*0.049) * u("in"), 'N2'

'''
T0, P0, T2, P2, M2, a, v = fanno_flow(288 * u('kelvin'), 6000 * u('psi'), 270 * u('kelvin'), 5500 * u('psi'), 0.288 * u('dimensionless'), 15 * u('ft'), 1.3 * u('kg/s'), 0.375 * u('in'), "N2")
T0.ito_base_units()
P0.ito('psi')
T2.ito_base_units()
P2.ito('psi')
M2.ito_base_units()
a.ito_base_units()
v.ito_base_units()
'''
