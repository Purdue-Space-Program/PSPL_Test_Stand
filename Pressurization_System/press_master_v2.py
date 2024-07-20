from CoolProp import CoolProp as cp
from pint import UnitRegistry
import numpy as np
import matplotlib.pyplot as plt
u = UnitRegistry()
Q_ = u.Quantity
from scipy.optimize import fsolve
from press_master import *

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
V_TANK = 3.40 * u('ft^3')



def run_test(combined=False):

    # Change These
    # Oxidizer
    OX = 'Oxygen'
    #MDOT_LOX = 8 * u('kg/s')  # LOX mass flow out of tank
    MDOT_LOX = 2.916 * u('lb/s')

    # Fuel
    FU = 'Methane'
    #MDOT_FU = 3 * u('kg/s')  # Methane mass flow out of tank
    MDOT_FU = 1.188 * u('lb/s')

    # Press Gas
    GAS_TYPE_LOX = 'N2'
    GAS_TYPE_FU = 'H2'

    # Bottle Test Case
    # Change These
    P_TANK = 315 * u('psi')  # Prop tank pressure
    P_source_lox = 6000 * u('psi')
    P_source_fu = 6000 * u('psi')
    tube_dia = (0.5 - 2*0.065) * u("in")

    # Leave these
    T_SOURCE = T_STD
    T_TANK_LOX = T_isent_exp(T_STD, P_source_lox, P_TANK, GAS_TYPE_LOX)
    T_TANK_FU = T_isent_exp(T_STD, P_source_fu, P_TANK, GAS_TYPE_FU)
    P_FILL = 40 * u('psi')  # Propellant fill pressure, used for density
    # ------- MAIN RUN - DON'T CHANGE VARS -------

    '''
    # Liquid Oxygen
    vdot_lox = prop_vol_flow(OX, MDOT_LOX, P_FILL)
    scfm_lox = scfm(vdot_lox, P_TANK, T_TANK)
    mdot_gas_lox = m_dot(vdot_lox, P_TANK, T_TANK, GAS_TYPE)
    throat_area_lox, throat_dia_lox = throat_calcs(mdot_gas_lox, P_source, T_SOURCE, GAS_TYPE)
    Cv_lox = c_v(scfm_lox, P_source, P_TANK, T_SOURCE, GAS_TYPE)
    vel_lox, mach_lox = vel_ideal_gas(T_SOURCE, P_source, tube_dia, mdot_gas_lox, GAS_TYPE)

    
    
    print(f"Simulation of pressurization of {FU} and {OX} test with two {V_TANK.to('L'): ~P.2f} tanks.")
    print(f"Maintaining tank pressure of {P_TANK.to('psi'): ~P.2f} from {P_source.to('psi'): ~P.2f} initial COPV pressurization.")
    print("Isentropic temp:", f"{T_TANK.to('kelvin'): ~P.2f}")
    print("Mass Flow Rates (lox, fu): ", f"{MDOT_LOX.to('kg/s'): ~P.2f}, {MDOT_FU.to('kg/s'): ~P.2f}")
    print("Line size (lox, fu): ", f"{tube_dia.to('in'): ~P.4f}, {tube_dia.to('in'): ~P.4f}")
    print("Tank pressure and press gas temp: ", f"{P_TANK.to('psi'): ~P.2f}, {T_TANK.to('kelvin'): ~P.2f}")
    print("COPV pressure and press gas temp: ", f"{P_source.to('psi'): ~P.2f}, {T_SOURCE.to('kelvin'): ~P.2f}")

    print()
    print("m_dot (lox, fu):", f"{mdot_gas_lox.to('kg/s'): ~P.4f},       " f"{mdot_gas_fu.to('kg/s'): ~P.4f}")
    print("SCFM (lox, fu): ", f"{scfm_lox: ~P.4f},  " f"{scfm_fu: ~P.4f}")
    print("Area (lox, fu): ", f"{throat_area_lox: ~P.4f},        ", f"{throat_area_fu: ~P.4f}")
    print("Dia (lox, fu):  ", f"{throat_dia_lox: ~P.4f},         ", f"{throat_dia_fu: ~P.4f}")
    print("C_v (lox, fu):  ", f"{Cv_lox: <.4f},            ", f"{Cv_fu: <.4f}")
    #print("Velocity (lox, fu):  ", f"{vel_lox: ~P.4f},         ", f"{vel_fu: ~P.4f}")
    #print("Mach (lox, fu):  ", f"{mach_lox: ~P.4f},         ", f"{mach_fu: ~P.4f}")
    '''

    if combined:
        #                    BULK GN2
        #                       |
        #                       |
        #         -------------------------------
        #        |                              |
        #       FU                             LOX
        print('I\'m combined circuit')

    else:
        #             BULK H2                  BULK GN2
        #               |                         |
        #               |                         |
        #             CH4                        LOX

        # METHANE CIRCUIT / HYDROGEN
        # Fuel
        vdot_fu = prop_vol_flow(FU, MDOT_FU, P_FILL)
        scfm_fu = scfm(vdot_fu, P_TANK, T_TANK_FU)
        mdot_gas_fu = m_dot(vdot_fu, P_TANK, T_TANK_FU, GAS_TYPE_FU)
        throat_area_fu, throat_dia_fu = throat_calcs(mdot_gas_fu, P_source_fu, T_SOURCE, GAS_TYPE_FU)
        Cv_fu = c_v(scfm_fu, P_source_fu, P_TANK, T_SOURCE, GAS_TYPE_FU)

        L_12 = 3 * u('ft')
        L_34 = 6 * u('ft')
        L_56 = 10 * u('ft')
        L_78 = 2 * u('ft')
        Cv_1 = 8 * u('dimensionless')
        Cv_2 = 8 * u('dimensionless')
        Cv_3 = 1.3 * u('dimensionless')

        T1, P1, T01, P01, M1, a1, v1 = isentropic_flow(T_SOURCE, P_source_fu, mdot_gas_fu, tube_dia, GAS_TYPE_FU)
        T02, P02, T2, P2, M2, a2, v2 = fanno_flow(T01, P01, T1, P1, M1, L_12, tube_dia, GAS_TYPE_FU)
        dP_1 = valve_dP(Cv_1, scfm_fu, P02, T2, GAS_TYPE_FU)
        P3 = P2 + dP_1
        P03 = P02 + dP_1

        T3, P3, T03, P03, M3, a3, v3 = isentropic_flow(T02, P03, mdot_gas_fu, tube_dia, GAS_TYPE_FU)
        T04, P04, T4, P4, M4, a4, v4 = fanno_flow(T03, P03, T3, P3, M3, L_34, tube_dia, GAS_TYPE_FU)
        dP_2 = valve_dP(Cv_2, scfm_fu, P04, T4, GAS_TYPE_FU)
        P5 = P4 + dP_2
        P05 = P04 + dP_2

        T5, P5, T05, P05, M5, a5, v5 = isentropic_flow(T04, P05, mdot_gas_fu, tube_dia, GAS_TYPE_FU)
        T06, P06, T6, P6, M6, a6, v6 = fanno_flow(T05, P05, T5, P5, M5, L_56, tube_dia, GAS_TYPE_FU)
        dP_3 = valve_dP(Cv_3, scfm_fu, P06, T6, GAS_TYPE_FU)
        P7 = P6 + dP_3
        P07 = P06 + dP_3

        T7, P7, T07, P07, M7, a7, v7 = isentropic_flow(T06, P07, mdot_gas_fu, tube_dia, GAS_TYPE_FU)
        T08, P08, T8, P8, M8, a8, v8 = fanno_flow(T07, P07, T7, P7, M7, L_78, tube_dia, GAS_TYPE_FU)

        print("T08_fu, P08_fu, T8_fu, P8_fu, M8_fu, speed of sound 8 fu, velocity 8 fu:", f"{T08: ~P.4f}, {P08: ~P.4f}, {T8: ~P.4f}, {P8: ~P.4f}, {M8: ~P.4f}, {a8: ~P.4f}, {v8: ~P.4f},")

        # LOX CIRCUIT / GN2
        vdot_lox = prop_vol_flow(OX, MDOT_LOX, P_FILL)
        scfm_lox = scfm(vdot_lox, P_TANK, T_TANK_LOX)
        mdot_gas_lox = m_dot(vdot_lox, P_TANK, T_TANK_LOX, GAS_TYPE_LOX)
        throat_area_lox, throat_dia_lox = throat_calcs(mdot_gas_lox, P_source_lox, T_SOURCE, GAS_TYPE_LOX)
        Cv_lox = c_v(scfm_lox, P_source_lox, P_TANK, T_SOURCE, GAS_TYPE_LOX)

        L_12 = 3 * u('ft')
        L_34 = 6 * u('ft')
        L_56 = 10 * u('ft')
        L_78 = 2 * u('ft')
        Cv_1 = 8 * u('dimensionless')
        Cv_2 = 8 * u('dimensionless')
        Cv_3 = 1.3 * u('dimensionless')

        T1, P1, T01, P01, M1, a1, v1 = isentropic_flow(T_SOURCE, P_source_lox, mdot_gas_lox, tube_dia, GAS_TYPE_LOX)
        T02, P02, T2, P2, M2, a2, v2 = fanno_flow(T01, P01, T1, P1, M1, L_12, tube_dia, GAS_TYPE_LOX)
        dP_1 = valve_dP(Cv_1, scfm_fu, P02, T2, GAS_TYPE_LOX)
        P3 = P2 + dP_1
        P03 = P02 + dP_1

        T3, P3, T03, P03, M3, a3, v3 = isentropic_flow(T02, P03, mdot_gas_lox, tube_dia, GAS_TYPE_LOX)
        T04, P04, T4, P4, M4, a4, v4 = fanno_flow(T03, P03, T3, P3, M3, L_34, tube_dia, GAS_TYPE_LOX)
        dP_2 = valve_dP(Cv_2, scfm_fu, P04, T4, GAS_TYPE_LOX)
        P5 = P4 + dP_2
        P05 = P04 + dP_2

        T5, P5, T05, P05, M5, a5, v5 = isentropic_flow(T04, P05, mdot_gas_lox, tube_dia, GAS_TYPE_LOX)
        T06, P06, T6, P6, M6, a6, v6 = fanno_flow(T05, P05, T5, P5, M5, L_56, tube_dia, GAS_TYPE_LOX)
        dP_3 = valve_dP(Cv_3, scfm_fu, P06, T6, GAS_TYPE_LOX)
        P7 = P6 + dP_3
        P07 = P06 + dP_3

        T7, P7, T07, P07, M7, a7, v7 = isentropic_flow(T06, P07, mdot_gas_lox, tube_dia, GAS_TYPE_LOX)
        T08, P08, T8, P8, M8, a8, v8 = fanno_flow(T07, P07, T7, P7, M7, L_78, tube_dia, GAS_TYPE_LOX)
        print("T08_fu, P08_fu, T8_fu, P8_fu, M8_fu, speed of sound 8 fu, velocity 8 fu:", f"{T08: ~P.4f}, {P08: ~P.4f}, {T8: ~P.4f}, {P8: ~P.4f}, {M8: ~P.4f}, {a8: ~P.4f}, {v8: ~P.4f},")
        print("Total Pressure: ", f"{P01.to('psi'): ~P.4f}, {P02.to('psi'): ~P.4f}, {P03.to('psi'): ~P.4f}, {P04.to('psi'): ~P.4f}, {P05.to('psi'): ~P.4f}, {P06.to('psi'): ~P.4f}, {P07.to('psi'): ~P.4f}, {P08.to('psi'): ~P.4f},")
        print("dP: ", f"{dP_1: ~P.4f}, {dP_2: ~P.4f}, {dP_3: ~P.4f}")

        print("m_dot (lox, fu):", f"{mdot_gas_lox.to('kg/s'): ~P.4f},       " f"{mdot_gas_fu.to('kg/s'): ~P.4f}")
        print("SCFM (lox, fu): ", f"{scfm_lox: ~P.4f},  " f"{scfm_fu: ~P.4f}")
        print("Area (lox, fu): ", f"{throat_area_lox: ~P.4f},        ", f"{throat_area_fu: ~P.4f}")
        print("Dia (lox, fu):  ", f"{throat_dia_lox: ~P.4f},         ", f"{throat_dia_fu: ~P.4f}")
        print("C_v (lox, fu):  ", f"{Cv_lox: <.4f},            ", f"{Cv_fu: <.4f}")

test4 = run_test()