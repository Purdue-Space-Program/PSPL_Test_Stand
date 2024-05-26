%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RTD Current Script
%
% Author: Hugo Filmer (Updated Max Van de Graaff
%
% Description: This script varies the convection coefficient to determine
% the current across an RTD assuming a constant voltage.
%
% Last updated 05/24/2023
%
% Contact me on Slack if you have any questions!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

%% Define system properties

% RTD1 is the upper RTD, and RTD 2 is the lower RTD

% Conversions
in2m = 1/39.37; % [m/in]

% Circuit properties
voltage_in = -5; % [V] The voltage into the RTD circuit

% RTD properties
%https://evosensors.com/products/cryogenic-2-wire-pt100-rtd-probe-2-inches-long-x-1-8-diameter-with-40-inches-of-pfa-insulated-lead-wire-and-miniature-male-flat-pin-connector
RTD1_length = 1 * in2m; % [m] The length of the RTD
RTD2_length = 1 * in2m; % [m] The length of the RTD
RTD3_length = 1 * in2m; % [m] The length of the RTD
RTD_diam = 0.125 * in2m; % [m] The diameter of the RTD
RTD_thick = 0.05 * in2m; % [m] The thickness of the sheath (estimated)
RTD1_area = pi*RTD_diam*RTD1_length; % [m^2] The exposed surface area, assuming RTD is a cylinder
RTD2_area = pi*RTD_diam*RTD2_length; % [m^2] The exposed surface area, assuming RTD is a cylinder
RTD3_area = pi*RTD_diam*RTD3_length; % [m^2] The exposed surface area, assuming RTD is a cylinder

RTD1_volume = (((RTD_diam/2)^2) - (((RTD_diam - 2*RTD_thick)/2)^2))*pi*RTD1_length; % [m^3]
RTD2_volume = (((RTD_diam/2)^2) - (((RTD_diam - 2*RTD_thick)/2)^2))*pi*RTD2_length; % [m^3]
RTD3_volume = (((RTD_diam/2)^2) - (((RTD_diam - 2*RTD_thick)/2)^2))*pi*RTD3_length; % [m^3]
RTD_density = 7930; % [kg/m^3] The density of the RTD, typical value for stainless steel
%RTD_density = 21500; % [kg/m^3] The density of the RTD, typical value for stainless steel
RTD1_mass = RTD1_volume * RTD_density; % [kg] The mass of the RTD
RTD2_mass = RTD2_volume * RTD_density; % [kg] The mass of the RTD
RTD3_mass = RTD2_volume * RTD_density; % [kg] The mass of the RTD

RTD1_specific_heat = 500; % [J/Kg K] The specific heat of the material, typical value for stainless steel
%RTD1_specific_heat = 125; % [J/Kg K] The specific heat of the material, typical value for platinum
RTD1_heat_cap = RTD1_specific_heat * RTD1_mass; % [J/K] The heat capacity of the RTD

RTD2_specific_heat = 500; % [J/Kg K] The specific heat of the material, typical value for stainless steel
%RTD2_specific_heat = 125; % [J/Kg K] The specific heat of the material, typical value for platinum
RTD2_heat_cap = RTD2_specific_heat * RTD2_mass; % [J/K] The heat capacity of the RTD

RTD3_specific_heat = 500; % [J/Kg K] The specific heat of the material, typical value for stainless steel
%RTD2_specific_heat = 125; % [J/Kg K] The specific heat of the material, typical value for platinum
RTD3_heat_cap = RTD3_specific_heat * RTD3_mass; % [J/K] The heat capacity of the RTD

% Fluid properties
temp = 100; % [K] The temp of the fluids, arbitrary value
gas_coeff = 10; % [W/m^2 K] The convection coefficient of the gas, typical value for free convection
liquid_coeff = 100; % [W/m^2 K] The convection coefficient of the liquid, typical value for free convection

%% Setup simulation

% RTD initial conditions
RTD1_initial_temp = 300; % [K] The initial temperature of the upper RTD
RTD2_initial_temp = 300; % [K] The initial temperature of the lower RTD
RTD3_initial_temp = 300; % [K] The initial temperature of the lower RTD

% Simulation parameters
timestep = 0.01; % [s] The timestep that the simulation uses
stop_time = 500; % [s] The time the simulation runs to

% Gas time must be greater than liquid time!
RTD1_liquid_time = 100; % [s] The time that RTD1 (upper) switches from gas to liquid
RTD1_gas_time = 500; % [s] The time that RTD1 (upper) switches from liquid to gas

RTD2_liquid_time = 200; % [s] The time that RTD2 (lower) switches from gas to liquid
RTD2_gas_time = 500; % [s] The time that RTD2 (lower) switches from liquid to gas

RTD3_liquid_time = 375; % [s] The time that RTD2 (lower) switches from gas to liquid
RTD3_gas_time = 500; % [s] The time that RTD2 (lower) switches from liquid to gas

% Create vectors
times = 0:timestep:stop_time;
RTD1_temps = zeros(1, length(times));
RTD2_temps = zeros(1, length(times));
RTD3_temps = zeros(1, length(times));
RTD1_resistances = zeros(1, length(times));
RTD2_resistances = zeros(1, length(times));
RTD3_resistances = zeros(1, length(times));
voltages_measured = zeros(1, length(times));
currents_measured = zeros(1, length(times));

%% Simulation

% Final setup
RTD1_temp = RTD1_initial_temp;
RTD2_temp = RTD2_initial_temp;
RTD3_temp = RTD3_initial_temp;
index = 1;

for t = times

    % Determine RTD heat transfer coefficients
    if (t >= RTD1_liquid_time) && (t <= RTD1_gas_time)
        RTD1_coeff = liquid_coeff;
    else
        RTD1_coeff = gas_coeff;
    end

    if (t >= RTD2_liquid_time) && (t <= RTD2_gas_time)
        RTD2_coeff = liquid_coeff;
    else
        RTD2_coeff = gas_coeff;
    end

    if (t >= RTD3_liquid_time) && (t <= RTD3_gas_time)
        RTD3_coeff = liquid_coeff;
    else
        RTD3_coeff = gas_coeff;
    end

    % Calculate new resistances
    RTD1_resistance = RTD_Resistance(RTD1_temp); % [Ohm]
    RTD2_resistance = RTD_Resistance(RTD2_temp); % [Ohm]
    RTD3_resistance = RTD_Resistance(RTD3_temp); % [Ohm]

    % Calculate electrical parameters !!
    voltage_measured = voltage_in*(RTD2_resistance/(RTD1_resistance + RTD2_resistance)) - voltage_in/2; % [V]
    RTD1_voltage = voltage_in*(RTD1_resistance/(RTD1_resistance + RTD2_resistance)); % [V]
    RTD2_voltage = voltage_in*(RTD2_resistance/(RTD1_resistance + RTD2_resistance)); % [V]
    RTD3_voltage = voltage_in*(RTD3_resistance/(RTD1_resistance + RTD3_resistance)); % [V]
    current_measured = voltage_in/(RTD1_resistance + RTD2_resistance); % [Amp]

    % Calculate amount of heat transferred
    Q1_in = ((RTD1_voltage^2)/RTD1_resistance) * timestep; % [J]
    Q1_out = (RTD1_temp - temp)*RTD1_area*RTD1_coeff * timestep; % [J]

    Q2_in = ((RTD2_voltage^2)/RTD2_resistance) * timestep; % [J]
    Q2_out = (RTD2_temp - temp)*RTD2_area*RTD2_coeff * timestep; % [J]

    Q3_in = ((RTD3_voltage^2)/RTD3_resistance) * timestep; % [J]
    Q3_out = (RTD3_temp - temp)*RTD3_area*RTD3_coeff * timestep; % [J]

    % Find new values of tracked variables
    RTD1_temps(index) = RTD1_temp;
    RTD2_temps(index) = RTD2_temp;
    RTD3_temps(index) = RTD3_temp;
    voltages_measured(index) = voltage_measured;
    currents_measured(index) = current_measured;
    RTD1_resistances(index) = RTD1_resistance;
    RTD2_resistances(index) = RTD2_resistance;
    RTD3_resistances(index) = RTD3_resistance;

    % Update temperatures
    RTD1_temp = RTD1_temp + (Q1_in - Q1_out)/RTD1_heat_cap; % [K]
    RTD2_temp = RTD2_temp + (Q2_in - Q2_out)/RTD2_heat_cap; % [K]
    RTD3_temp = RTD3_temp + (Q3_in - Q3_out)/RTD3_heat_cap; % [K]

    index = index + 1;

end

%% Find voltage and current derivatives

voltage_derivatives = zeros(1, length(times));

current_derivatives = zeros(1, length(times));

for index = 1:(length(voltages_measured) - 2)
    voltage_derivatives(index) = (voltages_measured(index + 1) - voltages_measured(index))/timestep;
end

for index = 1:(length(currents_measured) - 2)
    current_derivatives(index) = (currents_measured(index + 1) - currents_measured(index))/timestep;
end

%% Output

tiledlayout(3,2)

% Setup plots
nexttile
hold on
plot(times, RTD1_temps, "LineWidth", 2)
plot(times, RTD2_temps, "LineWidth", 2)
plot(times, RTD3_temps, "LineWidth", 3)
xline(RTD1_liquid_time, "r")
xline(RTD1_gas_time, "m")
xline(RTD2_liquid_time, "g")
xline(RTD2_gas_time, "m")
xline(RTD3_liquid_time, "b")
xline(RTD3_gas_time, "m")
title("RTD Temperature Over Time")
xlabel("Time [s]")
ylabel("Temperature [K]")
legend("Upper RTD", "Middle RTD", "Lower RTD")
grid on

nexttile
hold on
plot(times, RTD1_resistances, "LineWidth", 2)
plot(times, RTD2_resistances, "LineWidth", 2)
plot(times, RTD3_resistances, "LineWidth", 2)
xline(RTD1_liquid_time, "r")
xline(RTD1_gas_time, "m")
xline(RTD2_liquid_time, "g")
xline(RTD2_gas_time, "m")
xline(RTD3_liquid_time, "b")
xline(RTD3_gas_time, "m")
title("RTD Resistance Over Time")
xlabel("Time [s]")
ylabel("Resistance [Ohm]")
legend("Upper RTD", "Middle RTD", "Lower RTD")
grid on

nexttile
hold on
plot(times, -1 .* voltages_measured, "LineWidth", 2)
xline(RTD1_liquid_time, "r")
xline(RTD1_gas_time, "m")
xline(RTD2_liquid_time, "g")
xline(RTD2_gas_time, "m")
xline(RTD3_liquid_time, "b")
xline(RTD3_gas_time, "m")
title("Voltage Measured Over Time")
xlabel("Time [s]")
ylabel("Voltage Measured [V]")
grid on


nexttile
hold on
plot(times, -1 .* voltage_derivatives, "LineWidth", 2)
xline(RTD1_liquid_time, "r")
xline(RTD1_gas_time, "m")
xline(RTD2_liquid_time, "g")
xline(RTD2_gas_time, "m")
xline(RTD3_liquid_time, "b")
xline(RTD3_gas_time, "m")
title("Voltage Measured Derivative Over Time")
xlabel("Time [s]")
ylabel("Voltage Measured Derivative [V/s]")
grid on

nexttile
hold on
plot(times, currents_measured, "LineWidth", 2)
xline(RTD1_liquid_time, "r")
xline(RTD1_gas_time, "m")
xline(RTD2_liquid_time, "g")
xline(RTD2_gas_time, "m")
xline(RTD3_liquid_time, "b")
xline(RTD3_gas_time, "m")
title("Current Measured Over Time")
xlabel("Time [s]")
ylabel("Current Measured [Amp]")
grid on

nexttile
hold on
plot(times, current_derivatives, "LineWidth", 2)
xline(RTD1_liquid_time, "r")
xline(RTD1_gas_time, "m")
xline(RTD2_liquid_time, "g")
xline(RTD2_gas_time, "m")
xline(RTD3_liquid_time, "b")
xline(RTD3_gas_time, "m")
title("Current Measured Derivative Over Time")
xlabel("Time [s]")
ylabel("Current Measured Derivative [Amp/s]")
grid on

% Setup legend
qw{1} = plot(nan, 'r', "LineWidth", 2);
qw{2} = plot(nan, 'g', "LineWidth", 2);
qw{3} = plot(nan, 'b', "LineWidth", 2);
qw{4} = plot(nan, 'm', "LineWidth", 2);
l = legend([qw{:}], {'Liquid Upper','Gas Upper', 'Liquid Lower', 'Gas Lower'}, "Orientation", "Horizontal");
l.Layout.Tile = "north";

%% RTD resistance as a function of temperature
function resistance = RTD_Resistance(temp)
% https://www.engineeringtoolbox.com/pt100-electrical-resistance-d_1651.html

temp = temp - 273; % Convert from K to C for the fit

R_0 = 100; % [Ohm] The resistance at 0 C

% Linearization constants
A = 3.9083 * 10^-3;
B = -5.775 * 10^-7;

if temp > 273
    C = 0;
else
    C = -4.183 * 10^-12;
end

% Fit equation
resistance = R_0 * (1 + A*temp + (B*temp^2) + C*(temp - 100)*(temp^3)); % [Ohm]

end