import pandas as pd
from openpyxl import load_workbook
from pint import UnitRegistry
u = UnitRegistry()
Q_ = u.Quantity
from press_funcs import *

file = pd.read_excel("C:/Users/brigh/OneDrive - purdue.edu/PSP/Test Stand/PSPL_Test_Stand/Pressurization_System/Press System Test Cases.xlsx", 0)
with pd.ExcelWriter('Press System Test Case Outputs.xlsx', engine='openpyxl', mode='w') as writer:
    pd.DataFrame([]).to_excel(writer)

for test_num in range(len(file)):
    # Test Title
    title = file.loc[test_num, "Test Case Description"]
    
    # Mass Flows
    LOX = file.loc[test_num, "Oxidizer"]
    MDOT_LOX = file.loc[test_num, "OX mdot (kg/s)"] * u('kg/s')

    FU = file.loc[test_num, "Fuel"]
    MDOT_FU = file.loc[test_num, "FU mdot (kg/s)"] * u('kg/s')

    # Gas Types
    GAS_TYPE_LOX = file.loc[test_num, "OX Press Gas"]
    GAS_TYPE_FU = file.loc[test_num, "FU Press Gas"]

    # Pressures
    source = file.loc[test_num, "Single/Dual Source"]
    P_TANK_LOX = file.loc[test_num, "OX Tank Pressure (psi)"] * u('psi')
    P_TANK_FU = file.loc[test_num, "FU Tank Pressure (psi)"] * u('psi')
    P_SOURCE_LOX = file.loc[test_num, "OX Source Pressure (psi)"] * u('psi')
    P_SOURCE_FU = file.loc[test_num, "FU Source Pressure (psi)"] * u('psi')
    P_FILL = file.loc[test_num, "Fill Pressure (used for prop density) (psi)"] * u('psi')

    # Temperatures
    T_SOURCE_LOX = file.loc[test_num, "OX Source Temperature (K)"] * u('kelvin')
    T_SOURCE_FU = file.loc[test_num, "FU Source Temperature (K)"] * u('kelvin')
    T_TANK_LOX = T_isent_exp(T_STD, P_SOURCE_LOX, P_TANK_LOX, GAS_TYPE_LOX)
    T_TANK_FU = T_isent_exp(T_STD, P_SOURCE_FU, P_TANK_FU, GAS_TYPE_FU)

    config_num = file.loc[test_num, "Configuration No."]
    configuration = pd.read_excel("C:/Users/brigh/OneDrive - purdue.edu/PSP/Test Stand/PSPL_Test_Stand/Pressurization_System/Press System Test Cases.xlsx", config_num)

    # FUEL CIRCUIT
    # Mass flow calculations
    vdot_fu = prop_vol_flow(FU, MDOT_FU, P_FILL) # Computes required volumetric flow into the tank to match props
    scfm_fu = acfm_to_scfm(vdot_fu, P_TANK_FU, T_TANK_FU) # Computes the SCFM required at the given tank temp and pressure to match vol flow
    mdot_gas_fu = acfm_to_mdot(vdot_fu, P_TANK_FU, T_TANK_FU, GAS_TYPE_FU) # Mass flow @ given T & P to match vdot

    # Single orifice calculations
    throat_area_fu, throat_dia_fu = choked_flow_dim(mdot_gas_fu, P_SOURCE_FU, T_SOURCE_FU, GAS_TYPE_FU)
    Cv_fu = c_v(scfm_fu, P_SOURCE_FU, P_TANK_FU, T_SOURCE_FU, GAS_TYPE_FU)

    # LOX CIRCUIT / GN2
    # Mass flow calculations
    vdot_lox = prop_vol_flow(LOX, MDOT_LOX, P_FILL)
    scfm_lox = acfm_to_scfm(vdot_lox, P_TANK_LOX, T_TANK_LOX)
    mdot_gas_lox = acfm_to_mdot(vdot_lox, P_TANK_LOX, T_TANK_LOX, GAS_TYPE_LOX)

    # Single orifice calculations
    throat_area_lox, throat_dia_lox = choked_flow_dim(mdot_gas_lox, P_SOURCE_LOX, T_SOURCE_LOX, GAS_TYPE_LOX)
    Cv_lox = c_v(scfm_lox, P_SOURCE_LOX, P_TANK_LOX, T_SOURCE_LOX, GAS_TYPE_LOX)

    stations = []
    dP_list = []
    errors = []
    station = 0
    combined = False
    
    for component in range(len(configuration)):
        if configuration.loc[component, "Component"] == "Combined":
            flowpath = "Combined"
            combined = True
            mdot = mdot_gas_fu + mdot_gas_lox
            gas_type = GAS_TYPE_LOX
            stations.append([0 * u('dimensionless')])
            station += 1
            errors.append('')
            dia = configuration.loc[component + 1, "Diameter (in)"] * u('in')
            [T, P, T0, P0, M, a, v] = isentropic_flow(T_SOURCE_LOX, P_SOURCE_LOX, mdot, dia, gas_type)

        elif configuration.loc[component, "Component"] == "Fu":
            mdot = mdot_gas_fu
            stations.append([1 * u('dimensionless')])
            station += 1
            errors.append('')
            dia = configuration.loc[component + 1, "Diameter (in)"] * u('in')
            if combined:
                gas_type = GAS_TYPE_LOX
                if flowpath == "Combined":  # Exiting combined section
                    combined_exit = s
                    [T, P, T0, P0, M, a, v] = isentropic_flow(T0, P0, mdot, dia, gas_type)
                else:  # Exiting LOX section
                    [T, P, T0, P0, M, a, v] = isentropic_flow(combined_exit[2], combined_exit[3], mdot, dia, gas_type)
            else:  # Exiting tank; non-combined setup
                gas_type = GAS_TYPE_FU
                [T, P, T0, P0, M, a, v] = isentropic_flow(T_SOURCE_FU, P_SOURCE_FU, mdot, dia, gas_type)
                

            flowpath = "Fu"

        elif configuration.loc[component, "Component"] == "LOX":
            mdot = mdot_gas_lox
            stations.append([2 * u('dimensionless')])
            station += 1
            errors.append('')
            dia = configuration.loc[component + 1, "Diameter (in)"] * u('in')
            if combined:
                if flowpath == "Combined":  # Exiting combined section
                    combined_exit = s
                    gas_type = GAS_TYPE_LOX
                    [T, P, T0, P0, M, a, v] = isentropic_flow(T0, P0, mdot, dia, gas_type)
                else:  # Exiting LOX section
                    [T, P, T0, P0, M, a, v] = isentropic_flow(combined_exit[2], combined_exit[3], mdot, dia, gas_type)
            else:  # Exiting tank; non-combined setup
                gas_type = GAS_TYPE_LOX
                [T, P, T0, P0, M, a, v] = isentropic_flow(T_SOURCE_LOX, P_SOURCE_LOX, mdot, dia, gas_type)

            flowpath = "LOX"

        elif configuration.loc[component, "Component"] == "Tube":
            L = configuration.loc[component, "Length (ft)"] * u('ft')
            dia = configuration.loc[component, "Diameter (in)"] * u('in')

            # Tube entrance
            s = [T, P, T0, P0, M, a, v] = isentropic_flow(T0, P0, mdot, dia, gas_type)
            P01 = P0
            stations.append(s)
            errors.append(check(s, mdot, dia, gas_type))
            station += 1

            # Tube exit
            s = [T, P, T0, P0, M, a, v] = fanno_flow(T0, P0, T, P, M, L, dia, gas_type)
            
            dP = P0 - P01
            stations.append(s)
            errors.append(check(s, mdot, dia, gas_type))
            station += 1
            dP_list.append([f'{configuration.loc[component, "Component"]}, stations {station-1}-{station}:', dP])

        elif configuration.loc[component, "Component"] == "Regulator" or configuration.loc[component, "Component"] == "Valve":
            # Valve Entrance
            dia = configuration.loc[component, "Diameter (in)"] * u('in')
            s = [T, P, T0, P0, M, a, v] = isentropic_flow(T0, P0, mdot, dia, gas_type)
            P01 = P0
            stations.append(s)
            errors.append(check(s, mdot, dia, gas_type))
            station += 1

            # Valve Body
            Cv = configuration.loc[component, "Cv (dimensionless)"] * u('dimensionless')
            dP = valve_dP(Cv, scfm_lox, P0, T, gas_type)
            P += dP[0]
            P0 += dP[0]

            # Valve exit
            s = [T, P, T0, P0, M, a, v] = isentropic_flow(T0, P0, mdot, dia, gas_type)

            P01 = P0
            stations.append(s)
            errors.append(check(s, mdot, dia, gas_type, Cv, dP[1]))
            station += 1
            dP_list.append([f'{configuration.loc[component, "Component"]}, stations {station-1}-{station}:', dP[0]])

        else:
            invalid = configuration.loc[component, "Component"]
            print('INVALID COMPONENT' + f'{invalid}')
            #stations.append("INVALID COMPONENT")

    info = [[mdot_gas_lox, mdot_gas_fu],
            [scfm_lox, scfm_fu],
            [Cv_lox, Cv_fu]]

    print("m_dot (lox, fu):", f"{mdot_gas_lox.to('kg/s'): ~P.4f},       " f"{mdot_gas_fu.to('kg/s'): ~P.4f}")
    print("SCFM (lox, fu): ", f"{scfm_lox: ~P.4f},  " f"{scfm_fu: ~P.4f}")
    print("Area (lox, fu): ", f"{throat_area_lox: ~P.4f},        ", f"{throat_area_fu: ~P.4f}")
    print("Dia (lox, fu):  ", f"{throat_dia_lox: ~P.4f},         ", f"{throat_dia_fu: ~P.4f}")
    print("C_v (lox, fu):  ", f"{Cv_lox: <.4f},            ", f"{Cv_fu: <.4f}")

    print(throat_dia_fu, throat_dia_lox)
    stations = [[f"{elem: ~P.4f}" for elem in row] for row in stations]
    dP_list = [[row[0], f"{row[1]: ~P.4f}"] for row in dP_list]
    dP_df = pd.DataFrame(dP_list, columns=['Component', 'dP (psi)'])
    dP_df.insert(0, '', ['' for x in range(len(dP_df))])
    output = pd.DataFrame(stations, columns=['T (K)',  'P (psi)', 'T0 (K)', 'P0 (psi)', 'Mach', 'Speed of Sound', 'Velocity'])
    errors_df = pd.DataFrame(errors, columns=['Error'])
    output.index += 1
    dP_df.index += 1
    errors_df.index += 1
    output = pd.concat([output, errors_df, dP_df], axis=1)

    with pd.ExcelWriter('Press System Test Case Outputs.xlsx', engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
        output.to_excel(writer, sheet_name=f"Test {test_num}; Config {config_num}")