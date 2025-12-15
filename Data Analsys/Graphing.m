%% ===============================================================
%  ENGINE COMBUSTION + DATA ANALYSIS (FINAL VERSION + KPIs)
%% ===============================================================

clear; clc; close all;
addpath("Functions","Nasa");

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
Mi = [SpS.Mass];            % Molar masses [kg/mol]

%% ===============================================================
%  Starting Parameters
% enter all parameters here that came from the matching test result
%% ===============================================================

% Import emission data
fdaq_file = fullfile('Data\Diesel','20251125_0000006_CA15_IMEP2.5_fdaq.txt'); % fast data file
sdaq_file = fullfile('Data\Diesel','20251125_0000006_CA15_IMEP2.5_sdaq.txt'); % slow data file
csv_file = fullfile('Data', 'Diesel_emissions.csv');              % emissions file
emission_data = readtable(csv_file);

% User inputs for selection
CA_value = 12;      % [deg] Combustion advance (12, 15, or 18)
IMEP_value = 1.5;   % [bar] IMEP (1.5, 2.5, or 3.5)
use_variable_gamma = false; % can switch of variable vs non variable gamma

% Static parameters (not in CSV)
AFR_stoich = 14.5;       % [-] stoich AFR for diesel
LHV    = 42.7e6;         % [J/kg]
nonrenfactor = 1;        % Factor from 0 to 1 of renewable carbon percentage
n_test = 1500;           % [rpm] measurement speed

% Construct the name to search for
imep_str = strrep(num2str(IMEP_value), '.', '_');
search_name = sprintf('CA%d_IMEP%s', CA_value, imep_str);

% % Debug: display what we're searching for and what's available
% fprintf('Searching for: %s\n', search_name);
% fprintf('Available names in CSV:\n');
% disp(emission_data.Name);

% Find matching row
row_idx = find(strcmp(emission_data.Name, search_name));

if isempty(row_idx)
    error('No matching data found for CA=%d and IMEP=%g', CA_value, IMEP_value);
end

% Extract data from the matched row
lambda   = emission_data.lambda(row_idx);      % [-] excess-air ratio
CO_vol   = emission_data.CO_volpct(row_idx);   % [vol-%] CO dry gas
CO2_vol  = emission_data.CO2_volpct(row_idx);  % [vol-%] CO2 dry gas
HC_ppm   = emission_data.HC_ppm(row_idx);      % [ppm] HC
NOx_ppm  = emission_data.NOx_ppm(row_idx);     % [ppm] NOx
O2_vol   = emission_data.O2_volpct(row_idx);   % [vol-%] O2 dry gas
FSN      = emission_data.FSN(row_idx);         % [-] Filter Smoke Number

% Convert HC from ppm to vol-%
HC_vol = HC_ppm / 10000;  % [vol-%]

% Display loaded values
fprintf('Loaded data for %s:\n', search_name);
fprintf('  Lambda:  %.2f\n', lambda);
fprintf('  CO:      %.2f vol-%%\n', CO_vol);
fprintf('  CO2:     %.2f vol-%%\n', CO2_vol);
fprintf('  HC:      %.2f vol-%% (%d ppm)\n', HC_vol, HC_ppm);
fprintf('  NOx:     %d ppm\n', NOx_ppm);
fprintf('  O2:      %.2f vol-%%\n', O2_vol);
fprintf('  FSN:     %.2f\n', FSN);

% Manually entered data
fprintf('Manually entered data for %s:\n', search_name);
fprintf('  AFR_stoich:  %.2f\n', AFR_stoich);
fprintf('  LHV:         %.2f vol-%%\n', LHV);
fprintf('  nonrenfactor:%.2f vol-%%\n', nonrenfactor);

%% ========================================================================
%  BASE OPERATING CONDITIONS
%  ========================================================================
%  Define ambient conditions, fuel properties, and engine operating point
%  ========================================================================

Tamb = 21 + 273.15;             % Ambient temperature [K]

% Fuel elemental composition (weight percentages)
FuelWeightPercent = struct();
FuelWeightPercent.C = 0.86;     % % Carbon by weight
FuelWeightPercent.H = 0.14;     % % Hydrogen by weight
FuelWeightPercent.O = 0.00;     % % Oxygen by weight (pure hydrocarbon)

rho_air = 1.200012; % [kg/m^3]
eta_v = 0.95;


% Atomic masses [kg/mol]
M_atom = struct('C', 12.01e-3, 'H', 1.008e-3, 'O', 16.00e-3);


% Air composition (molar fractions)
Xair = [0 0.21 0.79 0 0];       % [Diesel, O2, N2, CO2, H2O]
MAir = Xair*Mi';                % Mean molar mass of air [kg/mol]
Yair = Xair.*Mi/MAir;           % Mass fractions in air

%% ======================
% UNITS
mm   = 1e-3;
dm   = 0.1;
bara = 1e5;

%% ======================
% ENGINE GEOMETRY
Cyl.Bore             = 104*mm;
Cyl.Stroke           = 85*mm;
Cyl.CompressionRatio = 21.5;
Cyl.ConRod           = 136.5*mm;
Cyl.TDCangle         = 0;          % TDC reference
n_engine             = 1500; % RPM

%% ======================
% LOAD FAST DATA (fdaq)
fdaq = table2array(readtable(fdaq_file));

Ca_raw   = fdaq(:,1);          % crank angle [deg]
p_raw    = fdaq(:,2) * bara;   % in-cylinder pressure [Pa] (NOT pegged)
Iinj_raw = fdaq(:,3);          % injector current [A]

%% Reshape into cycles
dCa = 0.2;                     % deg/step
pts_per_cycle = 720/dCa;       % 4-stroke -> 720 deg/cycle
Ncycles = floor(length(Ca_raw)/pts_per_cycle);

Ca   = reshape(Ca_raw ,[],Ncycles);
p    = reshape(p_raw  ,[],Ncycles);
Iinj = reshape(Iinj_raw,[],Ncycles);

iselect = 10;                  % cycle index to highlight

%% ======================
% LOAD SLOW DATA (sdaq)
sdaq = table2array(readtable(sdaq_file));

mfuel = sdaq(:,1) * 1e-3;      % [kg/s]   (g/s -> kg/s)
Tint  = sdaq(:,2);             % [°C]
Texh  = sdaq(:,3);             % [°C]
Pint  = sdaq(:,4) * bara;      % [Pa]     (bar(a) -> Pa)

mfuel_mean   = mean(mfuel);    % average fuel mass flow [kg/s]
Texh_mean_C  = mean(Texh);     % mean exhaust temperature [°C]

%% =============================================================
% DRIFT CORRECTION (PEGGING at BDC) - Applied to all cycles
[~, idxBDC] = min(abs(Ca(:,1) + 180));  % BDC ~ -180 deg CA

p_ref = mean(Pint);                     % intake pressure reference [Pa]

p_corr_all = zeros(size(p));
for i = 1:Ncycles
    p_meas = p(idxBDC, i);             % measured p at BDC (cycle i)
    p_corr_all(:, i) = p(:, i) + (p_ref - p_meas);
end

%% =============================================================
% PRESSURE SMOOTHING (noise reduction)
% First average all pegged cycles:
p_avg = mean(p_corr_all, 2);           % [Pa]

% Then apply Savitzky-Golay filter:
p_smooth = sgolayfilt(p_avg, 3, 81);   % (order 3, window 81 samples)

% Extract selected cycle (for plotting / volume):
Ca_sel  = Ca(:, iselect);
p_sel   = p(:, iselect);               % raw, uncorrected
p_corr  = p_corr_all(:, iselect);      % corrected (single cycle)

%% ======================
% PLOT: p–CA (raw)
figure;
plot(Ca, p/bara,'LineWidth',1);
hold on;
plot(Ca_sel, p_sel/bara,'r','LineWidth',2);
xlabel('Crank angle [deg]');
ylabel('Pressure [bar]');
title('p–CA Diagram (all cycles, raw)');
grid on;

%% ======================
% CYLINDER VOLUME FOR SELECTED CYCLE
V_sel = CylinderVolume(Ca_sel, Cyl);   % [m^3]

%% ======================
% PLOT: p–V (raw)
figure;
subplot(2,1,1)
plot(V_sel/dm^3, p_sel/bara,'LineWidth',1.5);
xlabel('V [dm^3]');
ylabel('p [bar]');
title('p–V (Linear, raw)');
grid on;

subplot(2,1,2)
mask_raw = (V_sel>0) & (p_sel>0);
loglog(V_sel(mask_raw)/dm^3, p_sel(mask_raw)/bara,'LineWidth',1.5);
xlabel('V [dm^3]');
ylabel('p [bar]');
title('p–V (Log-Log, raw)');
grid on;

%% =============================================================
% IMEP (using smoothed average pressure)
W_ind = trapz(V_sel, p_smooth);              % [J/cycle]
IMEP  = W_ind / (max(V_sel)-min(V_sel));     % [Pa] = [J/m^3]

fprintf('Indicated work: %.2f J/cycle\n', W_ind);
fprintf('IMEP (from p-V): %.2f bar\n', IMEP/bara);



%% =============================================================
% ROHR (using smoothed pressure)
gamma  = 1.35;                               % constant gamma
dtheta = deg2rad(mean(diff(Ca_sel)));        % [rad/step]
    
dp_dtheta = gradient(p_smooth, dtheta);      % [Pa/rad]
dV_dtheta = gradient(V_sel    , dtheta);     % [m^3/rad]


ROHRnonvar = (gamma/(gamma-1))*p_smooth.*dV_dtheta + ...
           (1/(gamma-1))*V_sel.*dp_dtheta;       % [J/rad]
    
% Convert to J/deg:
ROHR_degnonvar = ROHRnonvar * (pi/180);                  % [J/deg]


%% =============================================================
% VARIABLE GAMMA ROHR CALCULATION
%% =============================================================

% 1. Initial conditions
T_intake = mean(Tint) + 273.15;  % [K]

% --- PRE-COMBUSTION MIXTURE (fresh charge) ---
% Calculate mass of fresh charge per cycle
cycles_per_sec = n_test / 120;
mfuelinj = mfuel_mean / cycles_per_sec;  % [kg/cycle]


M_mean = MAir;
R_mixpre = (Runiv / MAir);   % [J/(kg·K)]

P_ivc = p_ref;          % intake pressure
T_ivc = T_intake;     % intake temperature
V_ivc = CylinderVolume(-180, Cyl);
m_trapped = (P_ivc * V_ivc) / (R_mixpre * T_ivc);
Ttest = (P_ivc * V_ivc) / (m_trapped * R_mixpre );
% Calculate total moles of air
n_air_total = m_trapped / MAir;  % [mol] (MAir in g/mol -> kg/mol)

% --- FUEL ELEMENTAL ANALYSIS ---
m_C_fuel = mfuelinj * FuelWeightPercent.C;
m_H_fuel = mfuelinj * FuelWeightPercent.H;
m_O_fuel = mfuelinj * FuelWeightPercent.O;
n_C_fuel = m_C_fuel / M_atom.C;
n_H_fuel = m_H_fuel / M_atom.H;
n_O_fuel = m_O_fuel / M_atom.O;

% --- STOICHIOMETRIC COMBUSTION ---
O2_for_C = n_C_fuel;
O2_for_H = n_H_fuel / 4;
O2_from_fuel = n_O_fuel / 2;
O2_consumed = O2_for_C + O2_for_H - O2_from_fuel;
CO2_produced = n_C_fuel;
H2O_produced = n_H_fuel / 2;

% --- POST-COMBUSTION MIXTURE ---
m_species = Yair * m_trapped;
n_species = m_species ./ Mi;
n_post = n_species;      % start with air composition
n_post(1) = 0;           % no diesel vapor — already included in mfuelinj
n_post(2) = n_post(2) - O2_consumed;
n_post(4) = n_post(4) + CO2_produced;
n_post(5) = n_post(5) + H2O_produced;

% Mass fractions of combustion products
m_post = n_post .* Mi;
Ypost = m_post / sum(m_post);
M_mean_post = sum(n_post.*Mi) / sum(n_post);   % molar mass in kg/mol
R_mixpost = Runiv / M_mean_post;

% Calculate total mass in cylinder per cycle
m_total_pre = m_trapped;                     % Before combustion
m_total_post = m_trapped + mfuelinj;         % After combustion (includes fuel)

% Calculate variable gamma and temperature at each crank angle
gamma = zeros(size(Ca_sel));
T = zeros(size(Ca_sel));
Cpi = zeros(1, NSp);
Cvi = zeros(1, NSp);

for idx = 1:length(Ca_sel)
    % ---- DETERMINE COMPOSITION AT THIS CRANK ANGLE ----
    if Ca_sel(idx) < 0
        % Before combustion: fresh charge
        Y_current = Yair;
        R_current = R_mixpre;
        m_current = m_total_pre;

    else
        % After combustion: fully burned products
        Y_current = Ypost;
        R_current = R_mixpost;
        m_current = m_total_post;
    end

    if (Ca_sel(idx)) > -180 &&  (Ca_sel(idx)) < 180
        % ---- CALCULATE TEMPERATURE FROM IDEAL GAS LAW ----
        T(idx) = (p_smooth(idx) * V_sel(idx)) / (m_current * R_current);  % [K]
    else
        T(idx) = T_intake;
    end
    % Convert species Cp, Cv (molar) to mass-based values:
    for i = 1:NSp
        Cp_mass(i) = CpNasa(T(idx), SpS(i)) / Mi(i);  % J/kg/K
        Cv_mass(i) = CvNasa(T(idx), SpS(i)) / Mi(i);  % J/kg/K
    end

    % mass-fraction weighted mixture properties
    cp = Y_current * Cp_mass';
    cv = Y_current * Cv_mass';

    gamma(idx) = cp / cv;

end

% 6. Calculate ROHR with variable gamma
ROHRvar = gamma./(gamma-1) .* p_smooth .* dV_dtheta + ...
    1./(gamma-1) .* V_sel .* dp_dtheta;

% 7. Convert to J/deg and continue with existing analysis
ROHR_degvar = ROHRvar * (pi/180);

%gamma diagnosis

gamma_min = min(gamma);
gamma_max = max(gamma);

gamma_delta = gamma_max - gamma_min;



% Plot variable gamma for diagnostics
figure;
subplot(2,1,1)
plot(Ca_sel, gamma, 'LineWidth', 1.5);
xlabel('Crank angle [deg]');
ylabel('Gamma [-]');
title('Variable Gamma vs Crank Angle');
grid on;
yline(1.33, '--r', 'Constant \gamma = 1.33');
legend('Variable \gamma', 'Constant \gamma');

subplot(2,1,2)
plot(Ca_sel, T - 273.15, 'LineWidth', 1.5);
xlabel('Crank angle [deg]');
ylabel('Temperature [°C]');
title('In-Cylinder Temperature vs Crank Angle');
grid on;

% Plot both variable and non-variable ROHR for comparison
comb_mask   = (Ca_sel >= -40) & (Ca_sel <= 120);
Ca_c = Ca_sel(comb_mask);

% Get both ROHR versions for the combustion window
ROHR_var_c = ROHR_degvar(comb_mask);
ROHR_nonvar_c = ROHR_degnonvar(comb_mask);

% Smooth both versions
ROHR_var_smooth = sgolayfilt(ROHR_var_c, 3, 31);
ROHR_nonvar_smooth = sgolayfilt(ROHR_nonvar_c, 3, 31);

% Create figure with both plots
figure;
plot(Ca_c, ROHR_var_smooth, 'LineWidth', 1.7, 'DisplayName', 'Variable γ');
hold on;
plot(Ca_c, ROHR_nonvar_smooth, 'LineWidth', 1.7, 'DisplayName', 'Constant γ');
hold off;

xlabel('Crank angle [deg]');
ylabel('ROHR [J/deg]');
title('Corrected Smoothed Rate of Heat Release');
legend('Location', 'best');
grid on;

% Select which one to use for further analysis based on flag
if use_variable_gamma == true
    ROHR_smooth = ROHR_var_smooth; % Use variable gamma for ROHR calculation
else
    ROHR_smooth = ROHR_nonvar_smooth;  % Use constant gamma for ROHR calculation
end

%% =============================================================
% CUMULATIVE HEAT RELEASE + CA10/50/90
cumHR = cumtrapz(Ca_c, ROHR_smooth);         % [J]

% Shift so minimum = 0 (removes integration offset)
[min_val, min_idx] = min(cumHR);
cumHR = cumHR - min_val;

cumHR_from_min = cumHR(min_idx:end);
Ca_c_from_min  = Ca_c(min_idx:end);

Q_total = max(cumHR_from_min);

CA10 = interp1(cumHR_from_min, Ca_c_from_min, 0.10 * Q_total);
CA50 = interp1(cumHR_from_min, Ca_c_from_min, 0.50 * Q_total);
CA90 = interp1(cumHR_from_min, Ca_c_from_min, 0.90 * Q_total);

figure;
plot(Ca_c, cumHR,'LineWidth',1.8); hold on;

yl = ylim;
y_label_pos = yl(1) + 0.05*(yl(2)-yl(1));

xline(CA10,'--r');
xline(CA50,'--g');
xline(CA90,'--b');
text(CA10, y_label_pos, 'CA10', 'Color', 'r', 'FontWeight','bold', ...
     'HorizontalAlignment','center');
text(CA50, y_label_pos, 'CA50', 'Color', 'g', 'FontWeight','bold', ...
     'HorizontalAlignment','center');
text(CA90, y_label_pos, 'CA90', 'Color', 'b', 'FontWeight','bold', ...
     'HorizontalAlignment','center');

xlabel('Crank angle [deg]');
ylabel('Cumulative Heat Release [J]');
title('Corrected Cumulative Heat Release');
grid on;

%% =============================================================
% --------- COMBUSTION KPIs ----------
[p_max, idx_pmax]    = max(p_smooth);
CA_pmax              = Ca_sel(idx_pmax);

[ROHR_max, idx_rmax] = max(ROHR_smooth);
CA_ROHR_max          = Ca_c(idx_rmax);

comb_duration        = CA90 - CA10;    % [deg]

fprintf('\n========== COMBUSTION KPIs ==========\n');
fprintf('p_max               = %.2f bar at CA = %.1f deg\n', p_max/bara, CA_pmax);
fprintf('ROHR_max            = %.1f J/deg at CA = %.1f deg\n', ROHR_max, CA_ROHR_max);
fprintf('CA10 / CA50 / CA90  = %.1f / %.1f / %.1f deg\n', CA10, CA50, CA90);
fprintf('Combustion duration = %.1f deg (CA10–CA90)\n', comb_duration);

%% =============================================================
% INJECTOR CURRENT (cycle-resolved)
figure;
plot(Ca_sel, Iinj(:,iselect),'LineWidth',1.5);
xlabel('Crank angle [deg]');
ylabel('Injector Current [A]');
title('Injector Current Profile (selected cycle)');
grid on;

%% =============================================================
% SMOOTHED p–V AND p–CA (using p_smooth)
figure;
subplot(2,1,1)
plot(V_sel/dm^3, p_smooth/bara,'LineWidth',1.8);
xlabel('V [dm^3]');
ylabel('p [bar]');
title('Smoothed p–V (Linear)');
grid on;

subplot(2,1,2)
mask_s = (V_sel>0) & (p_smooth>0);
loglog(V_sel(mask_s)/dm^3, p_smooth(mask_s)/bara,'LineWidth',1.8);
xlabel('V [dm^3]');
ylabel('p [bar]');
title('Smoothed p–V (Log-Log)');
grid on;

figure;
plot(Ca_sel, p_sel/bara,'LineWidth',1,'DisplayName', 'Raw pressure');
hold on
plot(Ca_sel, p_smooth/bara,'LineWidth',1.8, 'DisplayName', 'Filtered pressure');
hold off
xlabel('Crank angle [deg]');
ylabel('Pressure [bar]');
title('Smoothed p–CA Diagram');
legend('Location', 'best');
grid on;

% %% =============================================================
% % ======= PERFORMANCE KPIs (IMEP CALCULATED FROM DATA) =======
% 
% 
% cycles_per_sec = n_test / 120;  % 4-stroke: n/120
% 
% % Indicated power from THIS measured IMEP
% P_i = W_ind * cycles_per_sec;   % [W]
% 
% % Fuel energy input 
% E_fuel = mfuel_mean * LHV;      % [W]
% 
% % Efficiencies & BSFC at this operating point
% ITE  = P_i / E_fuel;                 % indicated thermal efficiency
% BSFC = (mfuel_mean / P_i) * 3.6e9;   % [g/kWh]
% 
% fprintf('\n========== PERFORMANCE KPIs (from calculated IMEP) ==========\n');
% fprintf('IMEP (calculated)   = %6.2f bar\n', IMEP/bara);
% fprintf('Indicated power P_i = %6.2f kW\n', P_i/1e3);
% fprintf('Fuel power E_fuel   = %6.2f kW\n', E_fuel/1e3);
% fprintf('ITE                 = %6.3f [-]\n', ITE);
% fprintf('BSFC                = %6.1f g/kWh\n', BSFC);



%% =============================================================
%  SLOW SENSOR PLOTS 

figure;
plot(mfuel,'LineWidth',1.5);
xlabel('Index');
ylabel('Fuel flow [kg/s]');
title('Fuel Mass Flow (slow sensor)');
grid on;

figure;
plot(Tint,'LineWidth',1.5); hold on;
plot(Texh,'LineWidth',1.5);
legend('Intake','Exhaust');
xlabel('Index');
ylabel('Temperature [°C]');
title('Intake & Exhaust Temperatures (slow sensors)');
grid on;

figure;
plot(Pint/bara,'LineWidth',1.5);
xlabel('Index');
ylabel('Pressure [bar]');
title('Intake Pressure (slow sensor)');
grid on;

% %% =============================================================
% %  EMISSIONS KPIs  (single operating point / single load)
% %% =============================================================
% 
% % ====== 2. CONSTANTS & BASIC ASSUMPTIONS ======
% M_exhaust    = 29e-3;      % [kg/mol] approx. exhaust molar mass
% M_CO2        = 44.01e-3;   % [kg/mol]
% M_CO         = 28.01e-3;   % [kg/mol]
% M_HC         = 44.1e-3;    % [kg/mol] Propane equivalent
% M_NOx        = 46.0e-3;    % [kg/mol] NO2-equivalent
% V_molar_STP  = 22.414e-3;  % [m^3/mol] molar volume at STP
% Texh_mean_K  = Texh_mean_C + 273.15;  % [K]
% GWP20_CH4    = 79;         % GWP values
% GWP100_CH4   = 27.2;       % GWP values
% 
% % Temperature correction factor (STP -> exhaust conditions)
% T_correction = 273.15 / Texh_mean_K;
% 
% % ====== 3. AIR & EXHAUST MASS FLOWS ======
% m_air     = mfuel_mean * lambda * AFR_stoich;   % [kg/s]
% m_exhaust = m_air + mfuel_mean;                 % [kg/s]
% 
% % Exhaust density and volumetric flow
% R_univ      = 8.314;                                     % [J/mol/K]
% rho_exhaust = p_ref * M_exhaust / (R_univ * Texh_mean_K);% [kg/m^3]
% V_exhaust   = m_exhaust / rho_exhaust;                   % [m^3/s]
% C_soot_mg_m3 = (4.95/0.405) * FSN * exp(0.38 * FSN);
% 
% 
% % ====== 4. EMISSION MASS FLOWS (from analyser readings) ======
% mass_CO  = (CO_vol  / 100) * V_exhaust * (M_CO  / V_molar_STP) * T_correction;   % [kg/s]
% mass_CO2 = (CO2_vol / 100) * V_exhaust * (M_CO2 / V_molar_STP) * T_correction;   % [kg/s]
% mass_HC  = (HC_vol / 100) * V_exhaust * (M_HC / V_molar_STP) * T_correction;     % [kg/s]
% mass_NOx = (NOx_ppm / 1e6) * V_exhaust * (M_NOx / V_molar_STP) * T_correction;   % [kg/s]
% mass_soot = C_soot_mg_m3 * V_exhaust;
% norm_nox = mass_NOx * ((21-15)/21 - O2_vol);
% 
% 
% % ====== 5. GHG-20 & DHD 100 calculations ======
% CO2eq_from_CO = mass_CO * (M_CO2 / M_CO);
% GHG20 = mass_CO2 + CO2eq_from_CO + mass_HC * GWP20_CH4;
% GHG100 = mass_CO2 + CO2eq_from_CO + mass_HC * GWP100_CH4;
% 
% % ====== 5. BRAKE-SPECIFIC EMISSIONS (per kWh) ======
% BSCO2 = (mass_CO2 / P_i) * 3.6e9;   % [g/kWh]
% BSCO2nonren = (mass_CO2*nonrenfactor / P_i) * 3.6e9;   % [g/kWh]
% BSCO  = (mass_CO  / P_i) * 3.6e9;   % [g/kWh]
% BSNOx = (mass_NOx / P_i) * 3.6e9;   % [g/kWh]
% BSGHG20 = (GHG20 / P_i) * 3.6e9;    % [g/kWh]
% BSGHG100 = (GHG100 / P_i) * 3.6e9;  % [g/kWh]
% BSSOOT = (mass_soot / P_i) * 3.6e6; % [mg/kWh]
% 
% 
% fprintf('\n========== EMISSIONS KPIs ==========\n');
% fprintf('Air mass flow       = %8.4f kg/s\n', m_air);
% fprintf('Exhaust mass flow   = %8.4f kg/s\n', m_exhaust);
% fprintf('Lambda              = %8.3f [-]\n', lambda);
% fprintf('N NOx               = %8.3f kg/s\n', norm_nox);
% fprintf('\nEmission mass flows:\n');
% fprintf('  CO2               = %8.5f kg/s\n', mass_CO2);
% fprintf('  CO                = %8.5f kg/s\n', mass_CO);
% fprintf('  HC                = %8.5f kg/s\n', mass_HC);
% fprintf('  NOx               = %8.5f kg/s\n', mass_NOx);
% fprintf('  GHG20             = %8.5f kg/s\n', GHG20);
% fprintf('  GHG100            = %8.5f kg/s\n', GHG100);
% fprintf('\nBrake-specific emissions:\n');
% fprintf('  BSCO2             = %8.1f g/kWh\n', BSCO2);
% fprintf('  BSCO2nonren       = %8.1f g/kWh\n', BSCO2);
% fprintf('  BSCO              = %8.2f g/kWh\n', BSCO);
% fprintf('  BSNOx             = %8.2f g/kWh\n', BSNOx);
% fprintf('  BSGHG20           = %8.2f g/kWh\n', BSGHG20);
% fprintf('  BSGHG100          = %8.2f g/kWh\n', BSGHG100);
% fprintf('  BSSOOT            = %8.2f mg/kWh\n', BSSOOT);

%% =============================================================
% CYLINDER VOLUME FUNCTION
function V = CylinderVolume(Ca, Cyl)
    B  = Cyl.Bore;
    S  = Cyl.Stroke;
    CR = Cyl.CompressionRatio;
    L  = Cyl.ConRod;

    Ap = pi*(B^2)/4;
    Vd = Ap*S;
    Vc = Vd/(CR-1);

    theta = deg2rad(Ca - Cyl.TDCangle);
    r = S/2;

    under = L^2 - (r*sin(theta)).^2;
    under(under<0) = 0;

    x = r*(1 - cos(theta)) + (L - sqrt(under));
    V = Vc + Ap*x;
end
