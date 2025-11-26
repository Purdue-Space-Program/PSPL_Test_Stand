% Weld sizing script
% Daniel DeConti and Tyler Trostle
% 2023-04-12

% Inputs:

weld_inputs = readtable("welds.xlsx");

weld_numbers = table2array(weld_inputs(:,1));
weld_labels = table2array(weld_inputs(:,2));

reaction_force = table2array(weld_inputs(:,3));
weld_thickness = table2array(weld_inputs(:,4));
weld_diameter = table2array(weld_inputs(:,5));
yield_stress = table2array(weld_inputs(:,6));
weld_type = table2array(weld_inputs(:,7));

num_inputs = size(reaction_force, 1);

shear_knockdown_factor = sqrt(3); % 1 / (1 / 0.577) 
weld_knockdown_factor = 1 / 0.8;

% Calculations

unity = ones(num_inputs, 1);
fillet = repelem(["Fillet"], num_inputs);
area =  weld_thickness .* pi .* weld_diameter .* (unity ./ (sqrt(2) - (sqrt(2) - unity) .* (unity - (weld_type == fillet))));
stress = reaction_force ./ area;
margin = yield_stress ./((shear_knockdown_factor - (shear_knockdown_factor - 1) .* (unity - (weld_type == fillet))) .* weld_knockdown_factor .* stress) - unity;

for(i = 1:1:num_inputs)
    fprintf("#%02d %s: %f\n", weld_numbers(i), string(weld_labels(i)), margin(i))
end

% Sources:
% Al 5356 weld yield strength: dfmachinespecialties.com/downloads/One-Sheets/alloy_5356.pdf