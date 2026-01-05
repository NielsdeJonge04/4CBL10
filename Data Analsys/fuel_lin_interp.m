% ============================================================
% Load emission data from CSV files (new format)
% ============================================================
diesel_results = readtable('Results\Diesel_g_per_MJ.csv');
hvo_results    = readtable('Results\HVO_g_per_MJ.csv');
gtl_results    = readtable('Results\GTL_g_per_MJ.csv');
dieselhvo_results = readtable('Results\HVO_Diesel_g_per_MJ.csv');
dieselgtl_results = readtable('Results\GTL_Diesel_g_per_MJ.csv');

% ============================================================
% Define IMEP values, crank angles, fuel names
% ============================================================
IMEP_values   = [1.5, 2.5, 3.5];    % bar
crank_angles  = [12, 15, 18];       % deg
fuel_names    = {'Diesel', 'HVO', 'GTL', 'Diesel\HVO blend', 'Diesel\GTL blend'};

% ============================================================
% Build matrices with correct row order
% ============================================================
fuelTables = {diesel_results, hvo_results, gtl_results, dieselhvo_results, dieselgtl_results};
GHG_100 = buildFuelMatrix(fuelTables, ...
    'BS_GHG100_gCO2eq_per_kWh', IMEP_values, crank_angles);
GHG_20 = buildFuelMatrix(fuelTables, ...
    'BS_GHG20_gCO2eq_per_kWh', IMEP_values, crank_angles);
soot_bsg = buildFuelMatrix(fuelTables, ...
    'BS_soot_mg_per_kWh', IMEP_values, crank_angles);
NOx_bsg = buildFuelMatrix(fuelTables, ...
    'BS_NNOx_g_per_kWh', IMEP_values, crank_angles);

% ============================================================
% USER SELECTION: Choose IMEP and CA values to plot
% ============================================================
selected_IMEP = [1.5, 2.5, 3.5];  % Select which IMEP values to plot
selected_CA = [12, 15, 18];        % Select which CA values to plot

% Select blend type: 'HVO' or 'GTL'
selected_blend = 'HVO';  % Options: 'HVO' or 'GTL'

% ============================================================
% Plotting function for fuel blend interpolation
% ============================================================
plotFuelBlendInterpolation(GHG_100, GHG_20, soot_bsg, NOx_bsg, ...
    IMEP_values, crank_angles, selected_IMEP, selected_CA, selected_blend);

% ============================================================
% Function: Plot linear interpolation for HVO/Diesel blend
% ============================================================
function plotFuelBlendInterpolation(GHG_100, GHG_20, soot_bsg, NOx_bsg, ...
    IMEP_values, crank_angles, selected_IMEP, selected_CA, selected_blend)
    
    % Determine which fuel columns to use based on selection
    if strcmpi(selected_blend, 'HVO')
        blend_fuel_col = 2;      % Column 2: HVO
        measured_blend_col = 4;  % Column 4: Diesel\HVO blend
        blend_name = 'HVO';
        fuel2_color = 'g';       % Green for HVO
    elseif strcmpi(selected_blend, 'GTL')
        blend_fuel_col = 3;      % Column 3: GTL
        measured_blend_col = 5;  % Column 5: Diesel\GTL blend
        blend_name = 'GTL';
        fuel2_color = 'c';       % Cyan for GTL
    else
        error('Invalid blend selection. Choose ''HVO'' or ''GTL''');
    end
    
    % Emission matrices (fuel types are columns)
    emission_data = {GHG_100, GHG_20, soot_bsg, NOx_bsg};
    emission_names = {'GHG-100 (gCO_2eq/kWh)', 'GHG-20 (gCO_2eq/kWh)', ...
                      'Soot (mg/kWh)', 'NO_x (g/kWh)'};
    
    % Calculate number of subplots needed
    n_subplots = length(selected_IMEP) * length(selected_CA);
    n_cols = min(3, n_subplots);  % Max 3 columns
    n_rows = ceil(n_subplots / n_cols);
    
    % Loop through each emission type (each gets its own figure)
    for em_idx = 1:length(emission_data)
        emission_matrix = emission_data{em_idx};
        
        % Create new figure for this emission type
        figure;
        
        subplot_idx = 1;
        
        % Loop through selected IMEP values
        for imep_val = selected_IMEP
            % Find IMEP index
            imep_idx = find(IMEP_values == imep_val);
            if isempty(imep_idx)
                warning('IMEP value %.1f not found', imep_val);
                continue;
            end
            
            % Loop through selected CA values
            for ca_val = selected_CA
                % Find CA index
                ca_idx = find(crank_angles == ca_val);
                if isempty(ca_idx)
                    warning('CA value %d not found', ca_val);
                    continue;
                end
                
                % Calculate row index in matrix
                % Row structure: [CA12,CA15,CA18] at IMEP1.5, then IMEP2.5, then IMEP3.5
                row_idx = (imep_idx - 1) * length(crank_angles) + ca_idx;
                
                % Extract data for this operating point
                diesel_val = emission_matrix(row_idx, 1);              % Column 1: Diesel
                fuel2_val = emission_matrix(row_idx, blend_fuel_col);  % HVO or GTL
                measured_blend = emission_matrix(row_idx, measured_blend_col);  % Measured blend
                
                % Create interpolation line (0% fuel2 to 100% fuel2)
                blend_ratio = 0:0.01:1;  % 0 = pure diesel, 1 = pure HVO/GTL
                interp_line = diesel_val * (1 - blend_ratio) + fuel2_val * blend_ratio;
                
                % 50/50 blend estimation from interpolation
                blend_50_estimated = diesel_val * 0.5 + fuel2_val * 0.5;
                
                % Create subplot
                subplot(n_rows, n_cols, subplot_idx);
                
                % Plot interpolation line
                plot(blend_ratio * 100, interp_line, 'b-', 'LineWidth', 2);
                hold on;
                
                % Plot actual measurement points
                plot(0, diesel_val, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
                plot(100, fuel2_val, [fuel2_color 'o'], 'MarkerSize', 10, 'MarkerFaceColor', fuel2_color);
                plot(50, measured_blend, 'ms', 'MarkerSize', 10, ...
                    'MarkerFaceColor', 'm', 'LineWidth', 2);
                
                % Plot 50/50 estimated point
                plot(50, blend_50_estimated, 'bd', 'MarkerSize', 10, ...
                    'MarkerFaceColor', 'b', 'LineWidth', 2);
                
                % Formatting
                xlabel(sprintf('%s Percentage (%%)', blend_name), 'FontSize', 11);
                ylabel(emission_names{em_idx}, 'FontSize', 11);
                title(sprintf('IMEP=%.1f bar, CA=%d°', imep_val, ca_val), 'FontSize', 12);
                grid on;
                legend('Linear Interpolation', 'Diesel', blend_name, ...
                    'Measured Blend', 'Estimated 50/50', 'Location', 'best', 'FontSize', 9);
                
                subplot_idx = subplot_idx + 1;
            end
        end
        
        % Overall title for this emission type
        sgtitle(sprintf('Fuel Blend Behavior (Diesel/%s): %s', blend_name, emission_names{em_idx}), ...
            'FontSize', 14, 'FontWeight', 'bold');
    end
end

% ============================================================
% Helper function to build fuel matrix
% ============================================================
function matrix = buildFuelMatrix(fuelTables, columnName, IMEP_values, crank_angles)
    % Initialize matrix: rows = operating points, columns = fuels
    n_points = length(IMEP_values) * length(crank_angles);
    n_fuels = length(fuelTables);
    matrix = zeros(n_points, n_fuels);
    
    % Fill matrix for each fuel
    for f = 1:n_fuels
        tbl = fuelTables{f};
        
        % Check for actual column names (display for first fuel only)
        if f == 1
            fprintf('Available columns in table: %s\n', strjoin(tbl.Properties.VariableNames, ', '));
            fprintf('Sample OperatingPoint values:\n');
            disp(tbl.OperatingPoint(1:min(3, height(tbl))));
        end
        
        row_idx = 1;
        
        % Loop through IMEP values
        for imep = IMEP_values
            % Loop through crank angles
            for ca = crank_angles
                % Create operating point string in format "CA12_IMEP1_5"
                % Replace decimal point with underscore
                imep_str = strrep(sprintf('%.1f', imep), '.', '_');
                op_point_str = sprintf('CA%d_IMEP%s', ca, imep_str);
                
                % Find matching row in table
                match = strcmp(tbl.OperatingPoint, op_point_str);
                
                if any(match)
                    matrix(row_idx, f) = tbl.(columnName)(match);
                else
                    matrix(row_idx, f) = NaN;
                    warning('No match found for IMEP=%.1f, CA=%d (tried "%s") in fuel %d', ...
                        imep, ca, op_point_str, f);
                end
                
                row_idx = row_idx + 1;
            end
        end
    end
end

