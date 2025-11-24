%% ========================================================================
%  DIESEL ENGINE THERMODYNAMIC CYCLE ANALYSIS
%  ========================================================================
%  This script performs a theoretical thermodynamic analysis of a diesel
%  engine cycle and compares it with measured pressure data. The analysis
%  includes:
%  - Intake, compression, combustion, expansion, and exhaust processes
%  - Calculation of key performance indicators (efficiency, BSFC, BSCO2)
%  - Visualization of p-V and p-Ca diagrams
%  ========================================================================

%% Initialization
% Clear workspace, command window, and close all figures
clear; clc; close;

% Add paths to required functions and NASA thermodynamic database
addpath("Functions", "Nasa");

%% ========================================================================
%  UNIT DEFINITIONS
%  ========================================================================
%  Define conversion factors for various units used throughout the script
%  ========================================================================

mm      = 1e-3;         % Millimeters to meters
dm      = 0.1;          % Decimeters to meters
bara    = 1e5;          % Bar to Pascal (absolute pressure)
MJ      = 1e6;          % Megajoules to Joules
kWhr    = 1000*3600;    % Kilowatt-hour to Joules
volperc = 0.01;         % Volume percentage (for emissions)
ppm     = 1e-6;         % Parts per million (volume fraction)
g       = 1e-3;         % Grams to kilograms
s       = 1;            % Seconds (base unit)
kelvin  = 273.15;       % Celsius to Kelvin conversion offset

%% ========================================================================
%  NASA THERMODYNAMIC DATABASE SETUP
%  ========================================================================
%  Load NASA polynomial coefficients for calculating temperature-dependent
%  thermodynamic properties (Cp, Cv, enthalpy, entropy)
%  ========================================================================

% Universal gas constant [J/(mol·K)]
global Runiv 
Runiv = 8.314;

% Load NASA thermodynamic database
TdataBase = fullfile('Nasa', 'NasaThermalDatabase');
load(TdataBase);

% Select chemical species relevant to diesel combustion
iSp = myfind({Sp.Name}, {'Diesel', 'O2', 'N2', 'CO2', 'H2O'});
SpS = Sp(iSp);              % Subselection of species from database
NSp = length(SpS);          % Number of species
Mi = [SpS.Mass];            % Molar masses [g/mol]

%% ========================================================================
%  ENGINE GEOMETRY DATA
%  ========================================================================
%  Physical dimensions and configuration of the diesel engine
%  ========================================================================

Cyl.Bore                = 104*mm;       % Cylinder bore diameter [m]
Cyl.Stroke              = 85*mm;        % Piston stroke length [m]
Cyl.CompressionRatio    = 21.5;         % Compression ratio [-]
Cyl.ConRod              = 136.5*mm;     % Connecting rod length [m]
Cyl.TDCangle            = 0;            % Top Dead Center angle [°CA]

% Valve timing events [°CA] (obsolute but kept for future reference)
CaIVO = -355;           % Intake Valve Opening angle
CaIVC = -135;           % Intake Valve Closing angle
CaEVO = 149;            % Exhaust Valve Opening angle
CaEVC = -344;           % Exhaust Valve Closing angle
CaSOI = -3.2;           % Start of Injection angle

% Calculate minimum and maximum cylinder volumes
Vmin = CylinderVolume(0, Cyl);          % Volume at TDC [m³]
Vmax = CylinderVolume(180, Cyl);        % Volume at BDC [m³]

%% ========================================================================
%  BASE OPERATING CONDITIONS
%  ========================================================================
%  Define ambient conditions, fuel properties, and engine operating point
%  ========================================================================

% Ambient conditions
pamb = 1.31325*bara;            % Ambient pressure [Pa]
Tamb = 21 + kelvin;             % Ambient temperature [K]

% Fuel properties
Q_LHV = 29.6*MJ;                % Lower Heating Value of diesel [J/kg]
mfuelinj = 1.8*1.1*10^-5;       % Fuel mass injected per cycle [kg]

% Exhaust conditions
p_exhaust = pamb + 0.2*bara;    % Exhaust manifold pressure [Pa]

% Fuel elemental composition (weight percentages)
FuelWeightPercent = struct();
FuelWeightPercent.C = 0.86;     % % Carbon by weight
FuelWeightPercent.H = 0.14;     % % Hydrogen by weight
FuelWeightPercent.O = 0.00;     % % Oxygen by weight (pure hydrocarbon)

% Atomic masses [g/mol]
M_atom = struct('C', 12.01, 'H', 1.008, 'O', 16.00);

% Air composition (molar fractions)
Xair = [0 0.21 0.79 0 0];       % [Diesel, O2, N2, CO2, H2O]
MAir = Xair*Mi';                % Mean molar mass of air [g/mol]
Yair = Xair.*Mi/MAir;           % Mass fractions in air

%% ========================================================================
%  CRANK ANGLE AND VOLUME VECTORS
%  ========================================================================
%  Create arrays for crank angle sweep and corresponding cylinder volumes
%  ========================================================================

% Define crank angle range and resolution
a_start = -360;                 % Starting crank angle [°CA]
a_fin = 360;                    % Ending crank angle [°CA]
a_res = 0.2;                    % Angular resolution [°CA]
a_steps = (abs(a_start) + abs(a_fin))/a_res;    % Number of steps

% Generate crank angle vector
Ca = linspace(a_start, a_fin, a_steps);         % [°CA]

% Calculate theoretical volume as function of crank angle
V_theory = CylinderVolume(Ca, Cyl);             % [m³]

% Initialize theoretical pressure vector
p_theory = zeros(size(Ca));                     % [Pa]

%% ========================================================================
%  INTAKE STROKE ANALYSIS (-360° to -180° CA)
%  ========================================================================
%  During intake, the cylinder fills with fresh air at ambient conditions
%  ========================================================================

% State point 6: End of exhaust / Start of intake
v_6 = Vmax;                     % Volume at BDC [m³]
p_6 = pamb;                     % Pressure equals ambient [Pa]

% Set pressure for entire intake stroke (constant at ambient pressure)
intake_indices = (Ca > -360) & (Ca <= 0);
p_theory(intake_indices) = pamb;

%% ========================================================================
%  COMPRESSION STROKE ANALYSIS (-180° to 0° CA)
%  ========================================================================
%  Isentropic compression of air from BDC to TDC
%  ========================================================================

% Initial state (State 1: BDC, intake valve closes)
p1 = pamb;                      % Initial pressure [Pa]
V1 = Vmax;                      % Initial volume [m³]
T1 = Tamb;                      % Initial temperature [K]
R1 = Runiv/MAir;                % Specific gas constant for air [J/(kg·K)]

% Calculate mass of air trapped in cylinder
m_air = (p1*V1)/(R1*T1);        % [kg]

% Calculate specific heat capacities using NASA polynomials at T1
for i = 1:NSp
    Cpi1(i) = CpNasa(T1, SpS(i));       % Cp for each species [J/(mol·K)]
    Cvi1(i) = CvNasa(T1, SpS(i));       % Cv for each species [J/(mol·K)]
end
cp1 = Yair * Cpi1';             % Mass-weighted Cp for air [J/(kg·K)]
cv1 = Yair * Cvi1';             % Mass-weighted Cv for air [J/(kg·K)]
gamma1 = 1.1;                   % Specific heat ratio (approximate calculation doesnt hold due to heat losses)

% Apply isentropic compression: p*V^gamma = constant
compression_indices = (Ca > -180) & (Ca <= 0);
V_compression = V_theory(compression_indices);

% Calculate pressure during compression
p_compression = @(Vcom) ((p1 * V1.^gamma1) ./ Vcom.^gamma1);
p_theory(compression_indices) = p_compression(V_compression);

% End state (State 2: TDC, end of compression)
V2 = Vmin;                      % Volume at TDC [m³]
p2 = p_compression(V2);         % Pressure at end of compression [Pa]
T2 = (p2*V2)/(m_air*R1);        % Temperature at end of compression [K]

%% ========================================================================
%  COMBUSTION ANALYSIS (NEAR TDC)
%  ========================================================================
%  Calculate fuel-air mixture composition before and after combustion
%  ========================================================================

% --- PRE-COMBUSTION MIXTURE ---
% Calculate absolute masses of each species in the cylinder
mabsair = Yair .* m_air;        % Mass of each species in air [kg]
mtot = m_air + mfuelinj;        % Total mass in cylinder [kg]

% Add fuel to mixture (Species order: [Diesel, O2, N2, CO2, H2O])
mabsair(1) = mfuelinj;          % Mass of diesel [kg]

% Calculate mass and molar fractions of pre-combustion mixture
Yprecomb = mabsair ./ mtot;                     % Mass fractions
Mmixpre = mabsair ./ Mi;                        % Moles of each species [mol]
Xprecomb = Mmixpre ./ sum(Mmixpre);             % Molar fractions

% --- FUEL ELEMENTAL ANALYSIS ---
% Calculate moles of C, H, O in the fuel based on composition
m_C_fuel = mfuelinj * FuelWeightPercent.C;      % Mass of carbon [kg]
m_H_fuel = mfuelinj * FuelWeightPercent.H;      % Mass of hydrogen [kg]
m_O_fuel = mfuelinj * FuelWeightPercent.O;      % Mass of oxygen [kg]

n_C_fuel = m_C_fuel / M_atom.C;                 % Moles of C [mol]
n_H_fuel = m_H_fuel / M_atom.H;                 % Moles of H [mol]
n_O_fuel = m_O_fuel / M_atom.O;                 % Moles of O [mol]

% --- STOICHIOMETRIC COMBUSTION ---
% Elemental combustion reactions:
%   C + O2 → CO2        (1 mol O2 per mol C, produces 1 mol CO2)
%   H + 0.25·O2 → 0.5·H2O   (0.25 mol O2 per mol H, produces 0.5 mol H2O)
%   O in fuel reduces O2 requirement by 0.5 mol per mol O

% Calculate oxygen consumption and product formation
O2_for_C = n_C_fuel;                    % O2 needed to burn carbon [mol]
O2_for_H = n_H_fuel / 4;                % O2 needed to burn hydrogen [mol]
O2_from_fuel = n_O_fuel / 2;            % O2 credit from fuel oxygen [mol]

O2_consumed = O2_for_C + O2_for_H - O2_from_fuel;   % Net O2 consumed [mol]
CO2_produced = n_C_fuel;                            % CO2 produced [mol]
H2O_produced = n_H_fuel / 2;                        % H2O produced [mol]

% --- POST-COMBUSTION MIXTURE ---
% Update species moles after combustion
Mmixpost = Mmixpre;                     % Start with pre-combustion values
Mmixpost(1) = 0;                        % All diesel consumed
Mmixpost(2) = Mmixpost(2) - O2_consumed;        % O2 consumed
Mmixpost(4) = Mmixpost(4) + CO2_produced;       % CO2 produced
Mmixpost(5) = Mmixpost(5) + H2O_produced;       % H2O produced

% Calculate post-combustion mass and molar fractions
Xpostcomb = Mmixpost ./ sum(Mmixpost);          % Molar fractions
Mpostcomb = Xpostcomb*Mi';                      % Mean molar mass [g/mol]
Ypostcomb = Xpostcomb.*Mi/Mpostcomb;            % Mass fractions

% --- TEMPERATURE RISE FROM COMBUSTION ---
% Energy released by fuel combustion
Q_in = mfuelinj*Q_LHV;          % Heat input from fuel [J]

% Calculate temperature after combustion (constant pressure assumption)
T3 = T2 + (Q_in/(cp1*mtot));    % Temperature after combustion [K]

% Update thermodynamic properties at new temperature
for i = 1:NSp
    Cpi2(i) = CpNasa(T3, SpS(i));       % Cp at T3 [J/(mol·K)]
    Cvi2(i) = CvNasa(T3, SpS(i));       % Cv at T3 [J/(mol·K)]
end
cp2 = Yair * Cpi2';             % Updated Cp [J/(kg·K)]
cv2 = Yair * Cvi2';             % Updated Cv [J/(kg·K)]
gamma2 = gamma1;                % Assume same gamma for simplicity
R3 = Runiv/Mpostcomb;           % Gas constant for products [J/(kg·K)]

% --- STATE AFTER COMBUSTION (State 3) ---
p3 = p2;                        % Pressure remains constant during combustion [Pa]
V3 = (mtot*R3*T3)/p3;           % Volume after combustion [m³]

% --- DETERMINE END OF COMBUSTION PERIOD ---
% Find crank angle where cylinder volume reaches V3
expansion_region = (Ca > 0) & (Ca <= 180);
combustion_end_index = find((V_theory >= V3) & expansion_region, 1, 'first');

% Handle case where V3 exceeds available volume
if isempty(combustion_end_index)
    warning('V3 exceeds maximum volume in expansion stroke. Setting combustion end to 180°');
    combustion_end_angle = 180;
else
    combustion_end_angle = Ca(combustion_end_index);
end

disp(['Combustion ends at: ', num2str(combustion_end_angle), '° CA'])

% Set constant pressure during combustion period
combustion_indices = (Ca > 0) & (Ca <= combustion_end_angle);
p_theory(combustion_indices) = p3;

%% ========================================================================
%  EXPANSION STROKE ANALYSIS (0° to 180° CA)
%  ========================================================================
%  Isentropic expansion of combustion products from end of combustion to BDC
%  ========================================================================

% Final state (State 4: BDC, end of expansion)
V4 = Vmax;                      % Volume at BDC [m³]

% Apply isentropic expansion: p*V^gamma = constant
expansion_indices = (Ca > combustion_end_angle) & (Ca <= 180);
V_expansion = V_theory(expansion_indices);

% Calculate pressure during expansion
p_expansion = @(Vexp) ((p3 * V3^gamma2) ./ Vexp.^gamma2);
p_theory(expansion_indices) = p_expansion(V_expansion);

% End state properties
p4 = p_expansion(V4);           % Pressure at end of expansion [Pa]
T4 = (p4*V4)/(mtot*R3);         % Temperature at end of expansion [K]

% Heat rejected during expansion (for energy balance)
Q_out = cv2*mtot*(T4 - T1);     % Heat lost [J]

%% ========================================================================
%  EXHAUST STROKE ANALYSIS (180° to 360° CA)
%  ========================================================================
%  Expulsion of combustion products at exhaust manifold pressure
%  ========================================================================

% State point 5: End of expansion / Start of exhaust
v_5 = Vmin;                     % Volume at TDC [m³]
p_5 = p_exhaust;                % Pressure drops to exhaust pressure [Pa]

% Set pressure for entire exhaust stroke (constant at exhaust pressure)
exhaust_indices = (Ca > 180) & (Ca <= 360);
p_theory(exhaust_indices) = p_exhaust;

% Set initial point to exhaust pressure for plot continuity
p_theory(1) = p_exhaust;

%% ========================================================================
%  WORK OUTPUT CALCULATION
%  ========================================================================
%  Calculate net work done by the engine per cycle
%  ========================================================================

% Net work = Heat in - Heat out - Pumping work
W_theoretical = Q_in - Q_out - (p_exhaust - pamb)*(Vmax - Vmin);    % [J]

%% ========================================================================
%  KEY PERFORMANCE INDICATORS (KPIs)
%  ========================================================================
%  Calculate thermal efficiency, fuel consumption, and emissions
%  ========================================================================

% --- THERMAL EFFICIENCY ---
% Ratio of work output to heat input
eta_th = W_theoretical / Q_in;                  % [-]

% --- BRAKE SPECIFIC FUEL CONSUMPTION (BSFC) ---
% Fuel consumption per unit work output
BSFC_kg_per_J = mfuelinj / W_theoretical;       % [kg/J]
BSFC_g_per_kWh = BSFC_kg_per_J * 1e3 * 3.6e6;   % [g/kWh]

% --- BRAKE SPECIFIC CO2 EMISSIONS (BSCO2) ---
% CO2 emissions per unit work output
% Formula: BSEM_E = y_E_fuel * BSFC
% where y_E_fuel = mass of pollutant E per kg of fuel

% Calculate CO2 produced per kg of fuel
kg_CO2_per_kg_fuel = Mmixpost(4) * Mi(4) / mfuelinj;   % [kg_CO2/kg_fuel]

% Apply BSEM formula for CO2
BSEM_CO2_kg_per_J = kg_CO2_per_kg_fuel * BSFC_kg_per_J;        % [kg_CO2/J]
BSEM_CO2_g_per_kWh = BSEM_CO2_kg_per_J * 1e3 * 3.6e6;          % [g_CO2/kWh]

% Assign to BSCO2 variable for clarity
BSCO2_g_per_kWh = BSEM_CO2_g_per_kWh;           % [g_CO2/kWh]

% --- DISPLAY RESULTS ---
disp(['Thermal Efficiency (η) = ', num2str(eta_th)]);
disp(['BSFC (kg/J) = ', num2str(BSFC_kg_per_J)]);
disp(['BSFC (g/kWh) = ', num2str(BSFC_g_per_kWh)]);
disp(['BSCO2 (g/kWh) = ', num2str(BSCO2_g_per_kWh)]);

%% ========================================================================
%  LOAD MEASURED EXPERIMENTAL DATA
%  ========================================================================
%  Read pressure trace data from engine test measurements
%  ========================================================================

% Load data file
FullName = fullfile('Data', 'ExampleDataSet.txt');
dataIn = table2array(readtable(FullName));
[Nrows, Ncols] = size(dataIn);                  % Get data dimensions

% Parse data structure
NdatapointsperCycle = 720/0.2;                  % Points per cycle (720°/0.2° resolution)
Ncycles = Nrows/NdatapointsperCycle;            % Number of engine cycles recorded

% Reshape data into matrix format [angle × cycle]
Ca_measured = reshape(dataIn(:,1), [], Ncycles);        % Crank angles [°CA]
p_measured = reshape(dataIn(:,2), [], Ncycles)*bara;    % Pressures [Pa]

%% ========================================================================
%  PLOTTING: PRESSURE vs CRANK ANGLE
%  ========================================================================
%  Visualize measured pressure traces and theoretical cycle
%  ========================================================================

f1 = figure(1);

% Plot all measured cycles
pp = plot(Ca_measured, p_measured/bara, 'LineWidth', 1);

% Plot theoretical cycle
pt = plot(Ca, p_theory/bara, 'LineWidth', 1);

% Axis labels and limits
xlabel('Crank Angle [°CA]'); 
ylabel('Pressure [bar]');
xlim([-360 360]); 
ylim([0 50]);

% Highlight a specific cycle (cycle 10)
iselect = 10;
line(Ca_measured(:,iselect), p_measured(:,iselect)/bara, ...
    'LineWidth', 2, 'Color', 'r');

% Add reference lines for valve timing events
YLIM = ylim;
line([CaIVC CaIVC], YLIM, 'LineWidth', 1, 'Color', 'b');    % Intake valve closes
line([CaEVO CaEVO], YLIM, 'LineWidth', 1, 'Color', 'r');    % Exhaust valve opens

% Configure plot appearance
set(gca, 'XTick', -360:60:360, 'XGrid', 'on', 'YGrid', 'on');
title('All cycles in one plot');

%% ========================================================================
%  PLOTTING: p-V DIAGRAM (Linear and Log-Log)
%  ========================================================================
%  Show thermodynamic cycle in pressure-volume coordinates
%  ========================================================================

% Calculate measured volume from crank angle
V_measured = CylinderVolume(Ca_measured(:,iselect), Cyl);

f2 = figure(2);

% --- Subplot 1: Linear scale ---
subplot(2,1,1)
plot(V_measured/dm^3, p_measured(:,iselect)/bara);
xlabel('Volume [dm³]'); 
ylabel('Pressure [bar]');
xlim([0 0.8]); 
ylim([0.5 50]);
set(gca, 'XTick', 0:0.1:0.8, 'XGrid', 'on', 'YGrid', 'on');
title({'p-V Diagram', '(with wrong Volume function btw)'})

% --- Subplot 2: Log-log scale ---
subplot(2,1,2)
loglog(V_measured/dm^3, p_measured(:,iselect)/bara);
loglog(V_theory/dm^3, p_theory/bara);
xlabel('Volume [dm³]'); 
ylabel('Pressure [bar]');
xlim([0.02 0.8]); 
ylim([1 55]);
set(gca, 'XTick', [0.02 0.05 0.1 0.2 0.5 0.8], ...
    'YTick', [1 2 5 10 20 50], 'XGrid', 'on', 'YGrid', 'on');
title('p-V Diagram (Log-Log Scale)');
legend('Measured', 'Theory');

%% ========================================================================
%  PLOTTING: p-V DIAGRAM WITH COLORED STROKES
%  ========================================================================
%  Distinguish different engine strokes by color
%  ========================================================================

figure;

% --- Intake stroke: -360° to -180° (Blue) ---
intake_indices = (Ca > -360) & (Ca <= -180);
intake_indices = intake_indices(1:length(V_theory));
plot(V_theory(intake_indices), p_theory(intake_indices), ...
    'b-', 'LineWidth', 2, 'DisplayName', 'Intake');

% --- Compression stroke: -180° to 0° (Red) ---
compression_indices = (Ca > -180) & (Ca <= 0);
compression_indices = compression_indices(1:length(V_theory));
plot(V_theory(compression_indices), p_theory(compression_indices), ...
    'r-', 'LineWidth', 2, 'DisplayName', 'Compression');

% --- Expansion/Power stroke: 0° to 180° (Green) ---
expansion_indices = (Ca > 0) & (Ca <= 180);
expansion_indices = expansion_indices(1:length(V_theory));
plot(V_theory(expansion_indices), p_theory(expansion_indices), ...
    'g-', 'LineWidth', 2, 'DisplayName', 'Expansion');

% --- Exhaust stroke: 180° to 360° (Magenta) ---
exhaust_indices = (Ca > 180) & (Ca <= 360) & 1;
exhaust_indices = exhaust_indices(1:length(V_theory));
plot(V_theory(exhaust_indices), p_theory(exhaust_indices), ...
    'm-', 'LineWidth', 2, 'DisplayName', 'Exhaust');

xlabel('Volume [m³]');
ylabel('Pressure [Pa]');
title('p-V Diagram - Engine Cycle');
legend('Location', 'best');
grid on;

%% ========================================================================
%  PLOTTING: COMBINED THEORETICAL AND MEASURED p-V DIAGRAM
%  ========================================================================
%  Overlay theoretical stroke-by-stroke analysis with measured data
%  ========================================================================

figure;


% --- Theoretical cycle with color-coded strokes ---

% Intake stroke (Blue)
intake_indices = (Ca > -360) & (Ca <= -180);
intake_indices = intake_indices(1:length(V_theory));
plot(V_theory(intake_indices)/dm^3, p_theory(intake_indices)/bara, ...
    'b-', 'LineWidth', 2, 'DisplayName', 'Intake (Theory)');

% Compression stroke (Red)
compression_indices = (Ca > -180) & (Ca <= 0);
compression_indices = compression_indices(1:length(V_theory));
plot(V_theory(compression_indices)/dm^3, p_theory(compression_indices)/bara, ...
    'r-', 'LineWidth', 2, 'DisplayName', 'Compression (Theory)');

% Expansion/Power stroke (Green)
expansion_indices = (Ca > 0) & (Ca <= 180);
expansion_indices = expansion_indices(1:length(V_theory));
plot(V_theory(expansion_indices)/dm^3, p_theory(expansion_indices)/bara, ...
    'g-', 'LineWidth', 2, 'DisplayName', 'Expansion (Theory)');

% Exhaust stroke (Magenta)
exhaust_indices = (Ca > 180) & (Ca <= 360);
exhaust_indices = exhaust_indices(1:length(V_theory));
plot(V_theory(exhaust_indices)/dm^3, p_theory(exhaust_indices)/bara, ...
    'm-', 'LineWidth', 2, 'DisplayName', 'Exhaust (Theory)');

% --- Measured data overlay (Black dashed) ---
V_measured = CylinderVolume(Ca_measured(:,iselect), Cyl);
plot(V_measured/dm^3, p_measured(:,iselect)/bara, ...
    'k--', 'LineWidth', 1.5, 'DisplayName', 'Measured');

% Configure plot
xlabel('Volume [dm³]');
ylabel('Pressure [bar]');
title('p-V Diagram - Engine Cycle (Linear Scale)');
legend('Location', 'best');
xlim([0 0.8]);
ylim([0.5 50]);
set(gca, 'XTick', 0:0.1:0.8, 'XGrid', 'on', 'YGrid', 'on');
grid on;

%% ========================================================================
%  END OF SCRIPT
%  ========================================================================