%% ===============================================================
%  ENGINE COMBUSTION + DATA ANALYSIS (FINAL VERSION + KPIs)
%% ===============================================================

clear; clc; close all;
addpath("Functions","Nasa");

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

%% ======================
% LOAD FAST DATA (fdaq)
fdaq_file = fullfile('Data','20251120_0000001_example_fdaq.txt');
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
sdaq_file = fullfile('Data','20251120_0000001_example_sdaq.txt');
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
p_corr  = p_corr_all(:, iselect);      % corrected (single cycle) [not directly used]

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

ROHR = (gamma/(gamma-1))*p_smooth.*dV_dtheta + ...
       (1/(gamma-1))*V_sel.*dp_dtheta;       % [J/rad]

% Convert to J/deg:
ROHR_deg = ROHR * (pi/180);                  % [J/deg]

% Restrict to combustion window and smooth again:
comb_mask   = (Ca_sel >= -40) & (Ca_sel <= 120);
Ca_c        = Ca_sel(comb_mask);
ROHR_c      = ROHR_deg(comb_mask);
ROHR_smooth = sgolayfilt(ROHR_c, 3, 31);

figure;
plot(Ca_c, ROHR_smooth,'LineWidth',1.7);
xlabel('Crank angle [deg]');
ylabel('ROHR [J/deg]');
title('Corrected Smoothed Rate of Heat Release');
grid on;

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
plot(Ca_sel, p_smooth/bara,'LineWidth',1.8);
xlabel('Crank angle [deg]');
ylabel('Pressure [bar]');
title('Smoothed p–CA Diagram');
grid on;

%% =============================================================
% ======= PERFORMANCE KPIs (IMEP CALCULATED FROM DATA) =======

n_test = 1500;                  % [rpm] measurement speed
cycles_per_sec = n_test / 120;  % 4-stroke: n/120

% Indicated power from THIS measured IMEP
P_i = W_ind * cycles_per_sec;   % [W]

% Fuel energy input (using diesel LHV)
LHV    = 42.7e6;                % [J/kg]
E_fuel = mfuel_mean * LHV;      % [W]

% Friction MEP (assumed) and mechanical efficiency from IMEP
pf_bar = 0.5;                   % [bar] typical friction mean effective pressure
pf     = pf_bar * bara;         % [Pa]

eta_mech = 1 - pf/IMEP;         % derived from IMEP and FMEP

% Brake power from indicated power and eta_mech
P_b    = eta_mech * P_i;        % [W]
P_b_kW = P_b / 1e3;             % [kW]

% Efficiencies & BSFC at this operating point
ITE  = P_i / E_fuel;                 % indicated thermal efficiency
BTE  = P_b / E_fuel;                 % brake thermal efficiency
BSFC = (mfuel_mean / P_b) * 3.6e9;   % [g/kWh]

fprintf('\n========== PERFORMANCE KPIs (from calculated IMEP) ==========\n');
fprintf('IMEP (calculated)   = %6.2f bar\n', IMEP/bara);
fprintf('Indicated power P_i = %6.2f kW\n', P_i/1e3);
fprintf('Brake power P_b     = %6.2f kW\n', P_b/1e3);
fprintf('Fuel power E_fuel   = %6.2f kW\n', E_fuel/1e3);
fprintf('ITE                 = %6.3f [-]\n', ITE);
fprintf('BTE                 = %6.3f [-]\n', BTE);
fprintf('Mechanical eff.     = %6.3f [-]\n', eta_mech);
fprintf('BSFC                = %6.1f g/kWh\n', BSFC);



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

%% =============================================================
%  EMISSIONS KPIs  (single operating point / single load)
%  -> Set the analyser readings below for each run (30%, 50%, 70% load)
%% =============================================================

% ====== 1. INPUTS FROM GAS ANALYSER (CHANGE PER LOAD) ======
lambda   = 1.6;     % [-] excess-air ratio at this load (EXAMPLE)
CO2_vol  = 9.0;     % [vol-%] CO2 dry gas
CO_vol   = 0.05;    % [vol-%] CO  dry gas
NOx_ppm  = 400;     % [ppm] NOx

% ====== 2. CONSTANTS & BASIC ASSUMPTIONS ======
AFR_stoich   = 14.5;       % [-] stoich AFR for diesel
M_exhaust    = 0.029;      % [kg/mol] approx. exhaust molar mass
M_CO2        = 0.044;      % [kg/mol]
M_CO         = 0.028;      % [kg/mol]
M_NOx        = 0.046;      % [kg/mol] NO2-equivalent
V_molar_STP  = 22.414e-3;  % [m^3/mol] molar volume at STP
Texh_mean_K  = Texh_mean_C + 273.15;  % [K]


% Temperature correction factor (STP -> exhaust conditions)
T_correction = 273.15 / Texh_mean_K;

% ====== 3. AIR & EXHAUST MASS FLOWS ======
m_air     = mfuel_mean * lambda * AFR_stoich;   % [kg/s]
m_exhaust = m_air + mfuel_mean;                 % [kg/s]

% Exhaust density and volumetric flow
R_univ      = 8.314;                                     % [J/mol/K]
rho_exhaust = p_ref * M_exhaust / (R_univ * Texh_mean_K);% [kg/m^3]
V_exhaust   = m_exhaust / rho_exhaust;                   % [m^3/s]

% ====== 4. EMISSION MASS FLOWS (from analyser readings) ======
mass_CO  = (CO_vol  / 100) * V_exhaust * (M_CO  / V_molar_STP) * T_correction;   % [kg/s]
mass_CO2 = (CO2_vol / 100) * V_exhaust * (M_CO2 / V_molar_STP) * T_correction;   % [kg/s]
mass_NOx = (NOx_ppm / 1e6) * V_exhaust * (M_NOx / V_molar_STP) * T_correction;   % [kg/s]

% ====== 5. BRAKE-SPECIFIC EMISSIONS (per kWh) ======
BSCO2 = (mass_CO2 / P_b) * 3.6e9;   % [g/kWh]
BSCO  = (mass_CO  / P_b) * 3.6e9;   % [g/kWh]
BSNOx = (mass_NOx / P_b) * 3.6e9;   % [g/kWh]

fprintf('\n========== EMISSIONS KPIs ==========\n');
fprintf('Air mass flow       = %8.4f kg/s\n', m_air);
fprintf('Exhaust mass flow   = %8.4f kg/s\n', m_exhaust);
fprintf('Lambda              = %8.3f [-]\n', lambda);
fprintf('\nEmission mass flows:\n');
fprintf('  CO2               = %8.5f kg/s\n', mass_CO2);
fprintf('  CO                = %8.5f kg/s\n', mass_CO);
fprintf('  NOx               = %8.5f kg/s\n', mass_NOx);
fprintf('\nBrake-specific emissions:\n');
fprintf('  BSCO2             = %8.1f g/kWh\n', BSCO2);
fprintf('  BSCO              = %8.2f g/kWh\n', BSCO);
fprintf('  BSNOx             = %8.2f g/kWh\n', BSNOx);

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
