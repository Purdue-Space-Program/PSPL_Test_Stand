from pint import UnitRegistry
u = UnitRegistry()
Q_ = u.Quantity
from press_funcs import *

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
    MDOT_LOX = 7.908 * u('lb/s')  # LOX mass flow out of tank

    # Fuel
    FU = 'Kerosene'
    MDOT_FU = 3.468 * u('lb/s')  # Methane mass flow out of tank

    # Press Gas
    GAS_TYPE_LOX = 'N2'
    GAS_TYPE_FU = 'N2'

    # Change These
    P_TANK_LOX = 450 * u('psi')  # Prop tank pressure 1055 450
    P_TANK_FU = 450 * u('psi')

    total_mdot = (1) * u('lb/s') # USER INPUT
    mdot_gas_fu = 0.1 * u('lb/s')
    mdot_gas_lox = 0.1 * u('lb/s')
    flag = False

    while abs(total_mdot.to('lb/s').m - (mdot_gas_lox.to('lb/s').m + mdot_gas_fu.to('lb/s').m)) > 0.00001:

        #print(".")
        if flag:
            total_mdot += ((mdot_gas_fu + mdot_gas_lox) - total_mdot) / 4
        flag = True
        #print(total_mdot)
        source = (5042*np.exp(-0.1505*total_mdot.to('lb/s').m) + -92.31*np.exp(2.358*total_mdot.to('lb/s').m)) * u('psi')
        #print(source)
        
        if source < 0:
            raise Exception("Source dP too low!!!")
        if source < P_TANK_FU:
            raise Exception("Press line below tank pressure!!")

        P_source_lox = source
        P_source_fu = source
        tube_dia = (0.75 - 2*0.083) * u("in")

        # Leave these
        T_SOURCE = T_STD
        T_TANK_LOX = T_isent_exp(T_STD, P_source_lox, P_TANK_LOX, GAS_TYPE_LOX)
        T_TANK_FU = T_isent_exp(T_STD, P_source_fu, P_TANK_FU, GAS_TYPE_FU)
        P_FILL = 14.7 * u('psi')  # Propellant fill pressure, used for density
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
            scfm_fu = acfm_to_scfm(vdot_fu, P_TANK_FU, T_TANK_FU)
            mdot_gas_fu = acfm_to_mdot(vdot_fu, P_TANK_FU, T_TANK_FU, GAS_TYPE_FU)
            throat_area_fu, throat_dia_fu = choked_flow_dim(mdot_gas_fu, P_source_fu, T_SOURCE, GAS_TYPE_FU)
            Cv_fu = c_v(scfm_fu, P_source_fu, P_TANK_FU, T_SOURCE, GAS_TYPE_FU)

            '''
            # region 10-station
            L_12 = 20 * u('ft')
            L_34 = 6 * u('ft')
            L_56 = 10 * u('ft')
            L_78 = 5 * u('ft')
            Cv_1 = 8 * u('dimensionless')
            Cv_2 = 8 * u('dimensionless')
            Cv_3 = 1.3 * u('dimensionless')

            T1, P1, T01, P01, M1, a1, v1 = isentropic_flow(T_SOURCE, P_source_fu, mdot_gas_fu, tube_dia, GAS_TYPE_FU)
            T2, P2, T02, P02, M2, a2, v2 = fanno_flow(T01, P01, T1, P1, M1, L_12, tube_dia, GAS_TYPE_FU)
            dP_1 = valve_dP(Cv_1, scfm_fu, P02, T2, GAS_TYPE_FU)[0]
            P3 = P2 + dP_1
            P03 = P02 + dP_1

            T3, P3, T03, P03, M3, a3, v3 = isentropic_flow(T02, P03, mdot_gas_fu, tube_dia, GAS_TYPE_FU)
            T4, P4, T04, P04, M4, a4, v4 = fanno_flow(T03, P03, T3, P3, M3, L_34, tube_dia, GAS_TYPE_FU)
            dP_2 = valve_dP(Cv_2, scfm_fu, P04, T4, GAS_TYPE_FU)[0]
            P5 = P4 + dP_2
            P05 = P04 + dP_2

            T5, P5, T05, P05, M5, a5, v5 = isentropic_flow(T04, P05, mdot_gas_fu, tube_dia, GAS_TYPE_FU)
            T6, P6, T06, P06, M6, a6, v6 = fanno_flow(T05, P05, T5, P5, M5, L_56, tube_dia, GAS_TYPE_FU)
            dP_3 = valve_dP(Cv_3, scfm_fu, P06, T6, GAS_TYPE_FU)[0]
            P7 = P6 + dP_3
            P07 = P06 + dP_3

            T7, P7, T07, P07, M7, a7, v7 = isentropic_flow(T06, P07, mdot_gas_fu, tube_dia, GAS_TYPE_FU)
            T8, P8, T08, P08, M8, a8, v8 = fanno_flow(T07, P07, T7, P7, M7, L_78, tube_dia, GAS_TYPE_FU)

            print("T08_fu, P08_fu, T8_fu, P8_fu, M8_fu, speed of sound 8 fu, velocity 8 fu:", f"{T08: ~P.4f}, {P08: ~P.4f}, {T8: ~P.4f}, {P8: ~P.4f}, {M8: ~P.4f}, {a8: ~P.4f}, {v8: ~P.4f},")
            # endregion
            '''

            # LOX CIRCUIT / GN2
            vdot_lox = prop_vol_flow(OX, MDOT_LOX, P_FILL)
            scfm_lox = acfm_to_scfm(vdot_lox, P_TANK_LOX, T_TANK_LOX)
            mdot_gas_lox = acfm_to_mdot(vdot_lox, P_TANK_LOX, T_TANK_LOX, GAS_TYPE_LOX)
            throat_area_lox, throat_dia_lox = choked_flow_dim(mdot_gas_lox, P_source_lox, T_SOURCE, GAS_TYPE_LOX)
            Cv_lox = c_v(scfm_lox, P_source_lox, P_TANK_LOX, T_SOURCE, GAS_TYPE_LOX)

            throat_area_combined, throat_dia_combined = choked_flow_dim(mdot_gas_lox + mdot_gas_fu, P_source_lox, T_SOURCE, GAS_TYPE_LOX)
            '''
            # region 10-station
            L_12 = 20 * u('ft')
            L_34 = 4 * u('ft')
            L_56 = 8 * u('ft')
            L_78 = 10 * u('ft')
            Cv_1 = 8 * u('dimensionless')
            Cv_2 = 8 * u('dimensionless')
            Cv_3 = 1.3 * u('dimensionless')

            T1, P1, T01, P01, M1, a1, v1 = isentropic_flow(T_SOURCE, P_source_lox, mdot_gas_lox + mdot_gas_fu, tube_dia, GAS_TYPE_LOX)
            T2, P2, T02, P02, M2, a2, v2 = fanno_flow(T01, P01, T1, P1, M1, L_12, tube_dia, GAS_TYPE_LOX)
            
            dP_1 = valve_dP(Cv_1, scfm_lox, P02, T2, GAS_TYPE_LOX)[0]
            P3 = P2 + dP_1
            P03 = P02 + dP_1

            T3, P3, T03, P03, M3, a3, v3 = isentropic_flow(T02, P03, mdot_gas_lox, tube_dia, GAS_TYPE_LOX)
            T4, P4, T04, P04, M4, a4, v4 = fanno_flow(T03, P03, T3, P3, M3, L_34, tube_dia, GAS_TYPE_LOX)
            dP_2 = valve_dP(Cv_2, scfm_lox, P04, T4, GAS_TYPE_LOX)[0]
            P5 = P4 + dP_2
            P05 = P04 + dP_2

            T5, P5, T05, P05, M5, a5, v5 = isentropic_flow(T04, P05, mdot_gas_lox, tube_dia, GAS_TYPE_LOX)
            T6, P6, T06, P06, M6, a6, v6 = fanno_flow(T05, P05, T5, P5, M5, L_56, tube_dia, GAS_TYPE_LOX)
            dP_3 = valve_dP(Cv_3, scfm_lox, P06, T6, GAS_TYPE_LOX)[0]
            P7 = P6 + dP_3
            P07 = P06 + dP_3

            T7, P7, T07, P07, M7, a7, v7 = isentropic_flow(T06, P07, mdot_gas_lox, tube_dia, GAS_TYPE_LOX)
            T8, P8, T08, P08, M8, a8, v8 = fanno_flow(T07, P07, T7, P7, M7, L_78, tube_dia, GAS_TYPE_LOX)
            print("T08_fu, P08_fu, T8_fu, P8_fu, M8_fu, speed of sound 8 fu, velocity 8 fu:", f"{T08: ~P.4f}, {P08: ~P.4f}, {T8: ~P.4f}, {P8: ~P.4f}, {M8: ~P.4f}, {a8: ~P.4f}, {v8: ~P.4f},")
            print("Total Pressure: ", f"{P01.to('psi'): ~P.4f}, {P02.to('psi'): ~P.4f}, {P03.to('psi'): ~P.4f}, {P04.to('psi'): ~P.4f}, {P05.to('psi'): ~P.4f}, {P06.to('psi'): ~P.4f}, {P07.to('psi'): ~P.4f}, {P08.to('psi'): ~P.4f},")
            print("dP: ", f"{dP_1: ~P.4f}, {dP_2: ~P.4f}, {dP_3: ~P.4f}")
            # endregion
            '''
            '''
            print("T08_fu, P08_fu, T8_fu, P8_fu, M8_fu, speed of sound 8 fu, velocity 8 fu:", f"{T01: ~P.4f}, {P01: ~P.4f}, {T1: ~P.4f}, {P1: ~P.4f}, {M1: ~P.4f}, {a1: ~P.4f}, {v1: ~P.4f},")
            print("T08_fu, P08_fu, T8_fu, P8_fu, M8_fu, speed of sound 8 fu, velocity 8 fu:", f"{T02: ~P.4f}, {P02: ~P.4f}, {T2: ~P.4f}, {P2: ~P.4f}, {M2: ~P.4f}, {a2: ~P.4f}, {v2: ~P.4f},")
            '''
            regflow = c_v_flow(1.3 * u('dimensionless'), P_source_fu, P_TANK_FU, T_STD, "Nitrogen")
            if regflow < scfm_fu or regflow < scfm_lox:
                raise Exception("REGULATORS TOO SMALL")

            if round(P_TANK_FU.to('psi').m) == 1055.0:
                chamber = 500
            elif round(P_TANK_FU.to('psi').m) == 450.0:
                chamber = 200
            print(f"---------------{FU} and {OX}, {chamber} psi Chamber Pressure-----------------")
            print("Propellant m_dot (lox, fu):", f"{MDOT_LOX.to('lb/s'): ~P.4f},       " f"{MDOT_FU.to('lb/s'): ~P.4f}")
            print(f"Source P (lox, fu):", f"        {P_source_lox.to('psi'): ~P.4f},     " f"{P_source_fu.to('psi'): ~P.4f}")
            print("m_dot (lox, fu):", f"           {mdot_gas_lox.to('lb/s'): ~P.4f},       " f"{mdot_gas_fu.to('lb/s'): ~P.4f}"f"    Total:{(mdot_gas_fu + mdot_gas_lox).to('lb/s'): ~P.4f}")
            print("v_dot (lox, fu):", f"           {vdot_lox.to('ft^3/min'): ~P.4f},       " f"{vdot_fu.to('ft^3/min'): ~P.4f}")
            print("T tank gas (lox, fu): ", f"     {T_TANK_LOX: ~P.4f},         " f"{T_TANK_FU: ~P.4f}")
            print("SCFM (lox, fu): ", f"           {scfm_lox: ~P.4f},   " f"{scfm_fu: ~P.4f}")
            print("Area (lox, fu): ", f"           {throat_area_lox: ~P.4f},        ", f"{throat_area_fu: ~P.4f}")
            print("Dia (lox, fu):  ", f"           {throat_dia_lox: ~P.4f},         ", f"{throat_dia_fu: ~P.4f}")
            print("C_v (lox, fu):  ", f"           {Cv_lox: <.4f},            ", f"{Cv_fu: <.4f}")
            print("Comb. flow (area, dia)", f"     {throat_area_combined: ~P.4f}      ", f"   {throat_dia_combined: ~P.4f}")
            print("Regulator Flow (SCFM)", f"      {regflow: ~P.4f}")
            print(c_v(scfm_lox, 6000 * u('psi'), 1547 * u('psi'), T_STD, "Nitrogen"))

run_test()
