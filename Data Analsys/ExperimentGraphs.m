clear; clc;

% ============================================================
% Load emission data from CSV files (new format)
% ============================================================
diesel_results = readtable('Diesel_g_per_MJ.csv');
hvo_results    = readtable('HVO_g_per_MJ.csv');
gtl_results    = readtable('GTL_g_per_MJ.csv');

% ============================================================
% Define IMEP values, crank angles, fuel names
% ============================================================
IMEP_values   = [1.5, 2.5, 3.5];    % bar
crank_angles  = [12, 15, 18];       % deg
fuel_names    = {'Diesel', 'HVO', 'GTL'};

% ============================================================
% Build matrices with correct row order:
% [CA12,CA15,CA18] at IMEP=1.5,
% [CA12,CA15,CA18] at IMEP=2.5,
% [CA12,CA15,CA18] at IMEP=3.5
% ============================================================
fuelTables = {diesel_results, hvo_results, gtl_results};

GHG_100 = buildFuelMatrix(fuelTables, ...
    'BS_GHG100_gCO2eq_per_kWh', IMEP_values, crank_angles);

GHG_20 = buildFuelMatrix(fuelTables, ...
    'BS_GHG20_gCO2eq_per_kWh', IMEP_values, crank_angles);

soot_bsg = buildFuelMatrix(fuelTables, ...
    'BS_soot_mg_per_kWh', IMEP_values, crank_angles);

NOx_bsg = buildFuelMatrix(fuelTables, ...
    'BS_NNOx_g_per_kWh', IMEP_values, crank_angles);

% ============================================================
% Plot appearance: colors per fuel, line style per IMEP
% ============================================================
colors = [
    0.74 0.84 0.66;   % Diesel   - sage green
    0.78 0.63 0.86;   % HVO      - lilac
    0.80 0.33 0.00    % GTL      - burnt orange
];

line_styles = {'-', '--', ':'};  % IMEP 1.5, 2.5, 3.5

% ============================================================
% Plots (using your existing helper)
% ============================================================
plotFuelImeps(crank_angles, GHG_100, fuel_names, IMEP_values, ...
              line_styles, colors, ...
              'Brake Specific GHG-100 (gCO2eq/kWh)', ...
              'Brake Specific GHG-100 Emissions vs Crank Angle');

plotFuelImeps(crank_angles, GHG_20, fuel_names, IMEP_values, ...
              line_styles, colors, ...
              'Brake Specific GHG-20 (gCO2eq/kWh)', ...
              'Brake Specific GHG-20 Emissions vs Crank Angle');

plotFuelImeps(crank_angles, soot_bsg, fuel_names, IMEP_values, ...
              line_styles, colors, ...
              'Brake Specific Soot (mg/kWh)', ...
              'Brake Specific Soot vs Crank Angle');

plotFuelImeps(crank_angles, NOx_bsg, fuel_names, IMEP_values, ...
              line_styles, colors, ...
              'Brake Specific NOx (g/kWh)', ...
              'Brake Specific NOx vs Crank Angle');

%% =====================================================================
% Helper: buildFuelMatrix
% ======================================================================
function Y = buildFuelMatrix(fuelTables, fieldName, imepVals, crankAngles)
% buildFuelMatrix
%   Reorders data for multiple fuels so that rows are grouped by IMEP,
%   then by crank angle, matching the expectation of plotFuelImeps.
%
%   For each fuel table, we expect rows with OperatingPoint of the form:
%       CA12_IMEP1_5
%       CA15_IMEP2_5
%       CA18_IMEP3_5
%   etc.
%
% INPUTS
%   fuelTables : 1 x nFuels cell array of tables (Diesel, HVO, GTL, ...)
%   fieldName  : name of the column in the table to extract
%   imepVals   : 1 x nImeps (e.g. [1.5 2.5 3.5])
%   crankAngles: 1 x nAngles (e.g. [12 15 18])
%
% OUTPUT
%   Y : [nAngles*nImeps x nFuels] matrix, with row order:
%       IMEP1: CA12, CA15, CA18
%       IMEP2: CA12, CA15, CA18
%       IMEP3: CA12, CA15, CA18

    nFuels   = numel(fuelTables);
    nImeps   = numel(imepVals);
    nAngles  = numel(crankAngles);

    Y = nan(nAngles * nImeps, nFuels);

    for f = 1:nFuels
        T  = fuelTables{f};
        op = string(T.OperatingPoint);   % OperatingPoint as string array

        for i = 1:nImeps
            imep = imepVals(i);

            for k = 1:nAngles
                angle = crankAngles(k);

                % Construct OperatingPoint name, e.g. 'CA12_IMEP1_5'
                pat = sprintf('CA%d_IMEP%.1f', angle, imep);
                pat = strrep(pat, '.', '_');  % turn '1.5' into '1_5'

                mask   = (op == pat);
                rowIdx = find(mask);

                if numel(rowIdx) ~= 1
                    error('Expected exactly one row for %s in fuel index %d, found %d.', ...
                          pat, f, numel(rowIdx));
                end

                % Global index in the [nAngles*nImeps] row order
                globalIdx = (i-1)*nAngles + k;

                Y(globalIdx, f) = T.(fieldName)(rowIdx);
            end
        end
    end
end

%% =====================================================================
% Your existing plotFuelImeps function (unchanged)
% ======================================================================
function plotFuelImeps(crank_angles, Y, fuel_names, imep_vals, line_styles, colors, yLabel, plotTitle)
% plotFuelImeps
%   Plots a matrix of emission data (or any quantity) versus crank angle
%   for multiple fuels (colors) and multiple IMEPs (line styles).
%
% INPUTS:
%   crank_angles : [nAngles x 1] vector of crank angles (e.g. [12 15 18])
%   Y            : [nAngles*nImeps x nFuels] matrix.
%                  For each fuel (column), rows are grouped as:
%                  IMEP1: all crank angles
%                  IMEP2: all crank angles
%                  IMEP3: all crank angles
%
%   fuel_names   : 1 x nFuels cell array of fuel names (e.g. {'Diesel','HVO','GTL'})
%   imep_vals    : 1 x nImeps vector of IMEP values (e.g. [1.5 2.5 3.5])
%   line_styles  : 1 x nImeps cell array of line styles (e.g. {'-','--',':'})
%   colors       : nFuels x 3 matrix of RGB colors for each fuel
%   yLabel       : (optional) string for y-axis label
%   plotTitle    : (optional) string for title

    % ------------------- Basic size checks -------------------
    if nargin < 3 || isempty(fuel_names)
        fuel_names = {'Fuel 1','Fuel 2','Fuel 3'};
    end

    nFuels  = size(Y, 2);
    nAngles = numel(crank_angles);

    if nargin < 4 || isempty(imep_vals)
        % Infer IMEP count from size of Y
        nImeps = size(Y, 1) / nAngles;
        imep_vals = 1:nImeps;  % generic labels if not given
    else
        nImeps = numel(imep_vals);
    end

    if nargin < 5 || isempty(line_styles)
        line_styles = repmat({'-'}, 1, nImeps);
    end

    if nargin < 6 || isempty(colors)
        colors = lines(nFuels); % default MATLAB colormap
    end

    % Sanity check on dimensions
    assert(size(Y,1) == nAngles * nImeps, ...
        'Y must have nAngles * nImeps rows.');
    assert(numel(fuel_names) == nFuels, ...
        'fuel_names length must match number of columns of Y.');
    assert(size(colors,1) == nFuels && size(colors,2) == 3, ...
        'colors must be nFuels x 3 RGB matrix.');
    assert(numel(line_styles) == nImeps, ...
        'line_styles length must match number of IMEP values.');

    % Make IMEP labels
    imep_labels = arrayfun(@(v) sprintf('IMEP = %.1f bar', v), ...
                           imep_vals, 'UniformOutput', false);

    % ------------------- Actual plotting -------------------
    figure;
    hold on; grid on;

    % Plot all real data (hidden from legend)
    for j = 1:nFuels          % loop over fuels (columns)
        for i = 1:nImeps      % loop over IMEPs (blocks in rows)
            idx_start = (i-1)*nAngles + 1;
            idx_end   =  i   *nAngles;

            plot(crank_angles, ...
                 Y(idx_start:idx_end, j), ...
                 'LineStyle',       line_styles{i}, ...
                 'Color',           colors(j,:), ...
                 'Marker',          '.', ...
                 'MarkerSize',      15, ...
                 'HandleVisibility','off'); % hide from legend
        end
    end

    % --- Fuel legend entries (color) ---
    for j = 1:nFuels
        plot(NaN, NaN, '-', ...
             'Color',      colors(j,:), ...
             'LineWidth',  1.5, ...
             'DisplayName', fuel_names{j});
    end

    % --- IMEP legend entries (line styles) ---
    for i = 1:nImeps
        plot(NaN, NaN, line_styles{i}, ...
             'Color',      'k', ...        % only style matters here
             'LineWidth',  1.5, ...
             'DisplayName', imep_labels{i});
    end

    legend('show', 'Location', 'bestoutside');

    % Labels and title
    xlabel('Crank Angle (°)');
    if nargin >= 7 && ~isempty(yLabel)
        ylabel(yLabel);
    else
        ylabel('Quantity (units)');
    end

    if nargin >= 8 && ~isempty(plotTitle)
        title(plotTitle);
    end

    hold off;
end
