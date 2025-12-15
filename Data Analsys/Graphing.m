%% ===============================================================
%  ENGINE COMBUSTION + DATA ANALYSIS (FINAL VERSION + KPIs)
%  FULL SCRIPT (PASTE INTO MATLAB)
%  Update included: ROBUST CA10/50/90 (fixes CA90 instability)
%   - Wider window (-60 to 180 deg)
%   - Threshold-based burn start/end (ROHR > 2% of peak)
%   - Enforced monotonic cumHR using cummax()
%% ===============================================================

clear; clc; close all;
addpath("Functions","Nasa");

%% ========================================================================
%  NASA THERMODYNAMIC DATABASE SETUP
%% ========================================================================

global Runiv
Runiv = 8.314;

TdataBase = fullfile('Nasa', 'NasaThermalDatabase');
load(TdataBase);

iSp = myfind({Sp.Name}, {'Diesel', 'O2', 'N2', 'CO2', 'H2O'});
SpS = Sp(iSp);
NSp = length(SpS);
Mi = [SpS.Mass];

%% ===============================================================
%  Starting Parameters
%% ===============================================================

% Import emission data
fdaq_file = fullfile('Data\Diesel','20251125_0000001_CA12_IMEP1.5_fdaq.txt');
sdaq_file = fullfile('Data\Diesel','20251125_0000001_CA12_IMEP1.5_sdaq.txt');
csv_file  = fullfile('Data', 'Emissions_csv\Diesel_emissions.csv');
emission_data = readtable(csv_file);

% User inputs
CA_value = 12;
IMEP_value = 1.5;
use_variable_gamma = false;

% Static parameters
AFR_stoich   = 14.5;
LHV          = 42.7e6;
nonrenfactor = 1;
n_test       = 1500;

% Construct name to search in CSV
imep_str = strrep(num2str(IMEP_value), '.', '_');
search_name = sprintf('CA%d_IMEP%s', CA_value, imep_str);

row_idx = find(strcmp(emission_data.Name, search_name));
if isempty(row_idx)
    error('No matching data found for CA=%d and IMEP=%g', CA_value, IMEP_value);
end

% Extract emissions
lambda   = emission_data.lambda(row_idx);
CO_vol   = emission_data.CO_volpct(row_idx);
CO2_vol  = emission_data.CO2_volpct(row_idx);
HC_ppm   = emission_data.HC_ppm(row_idx);
NOx_ppm  = emission_data.NOx_ppm(row_idx);
O2_vol   = emission_data.O2_volpct(row_idx);
FSN      = emission_data.FSN(row_idx);

HC_vol = HC_ppm / 10000;

fprintf('Loaded data for %s:\n', search_name);
fprintf('  Lambda:  %.2f\n', lambda);
fprintf('  CO:      %.2f vol-%%\n', CO_vol);
fprintf('  CO2:     %.2f vol-%%\n', CO2_vol);
fprintf('  HC:      %.2f vol-%% (%d ppm)\n', HC_vol, HC_ppm);
fprintf('  NOx:     %d ppm\n', NOx_ppm);
fprintf('  O2:      %.2f vol-%%\n', O2_vol);
fprintf('  FSN:     %.2f\n', FSN);

%% ========================================================================
%  BASE OPERATING CONDITIONS
%% ========================================================================

Tamb = 21 + 273.15;

FuelWeightPercent = struct();
FuelWeightPercent.C = 0.86;
FuelWeightPercent.H = 0.14;
FuelWeightPercent.O = 0.00;

rho_air = 1.200012;
eta_v = 0.95;

M_atom = struct('C', 12.01e-3, 'H', 1.008e-3, 'O', 16.00e-3);

Xair = [0 0.21 0.79 0 0];
MAir = Xair*Mi';
Yair = Xair.*Mi/MAir;

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
Cyl.TDCangle         = 0;
n_engine             = 1500;

%% ======================
% LOAD FAST DATA (fdaq)
fdaq = table2array(readtable(fdaq_file));

Ca_raw   = fdaq(:,1);
p_raw    = fdaq(:,2) * bara;
Iinj_raw = fdaq(:,3);

dCa = 0.2;
pts_per_cycle = 720/dCa;
Ncycles = floor(length(Ca_raw)/pts_per_cycle);

Ca   = reshape(Ca_raw ,[],Ncycles);
p    = reshape(p_raw  ,[],Ncycles);
Iinj = reshape(Iinj_raw,[],Ncycles);

iselect = 10;

%% ======================
% LOAD SLOW DATA (sdaq)
sdaq = table2array(readtable(sdaq_file));

mfuel = sdaq(:,1) * 1e-3;
Tint  = sdaq(:,2);
Texh  = sdaq(:,3);
Pint  = sdaq(:,4) * bara;

mfuel_mean  = mean(mfuel);
Texh_mean_C = mean(Texh);

%% =============================================================
% DRIFT CORRECTION (PEGGING at BDC) - all cycles
[~, idxBDC] = min(abs(Ca(:,1) + 180));
p_ref = mean(Pint);

p_corr_all = zeros(size(p));
for i = 1:Ncycles
    p_meas = p(idxBDC, i);
    p_corr_all(:, i) = p(:, i) + (p_ref - p_meas);
end

%% =============================================================
% PRESSURE SMOOTHING
p_avg = mean(p_corr_all, 2);
p_smooth = sgolayfilt(p_avg, 3, 81);

Ca_sel  = Ca(:, iselect);
p_sel   = p(:, iselect);
p_corr  = p_corr_all(:, iselect);

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
% CYLINDER VOLUME
V_sel = CylinderVolume(Ca_sel, Cyl);

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
% IMEP (smoothed avg pressure)
W_ind = trapz(V_sel, p_smooth);
IMEP  = W_ind / (max(V_sel)-min(V_sel));

fprintf('Indicated work: %.2f J/cycle\n', W_ind);
fprintf('IMEP (from p-V): %.2f bar\n', IMEP/bara);

%% =============================================================
% ROHR (using smoothed pressure)
gamma_const  = 1.35;
dtheta = deg2rad(mean(diff(Ca_sel)));

dp_dtheta = gradient(p_smooth, dtheta);
dV_dtheta = gradient(V_sel    , dtheta);

ROHRnonvar = (gamma_const/(gamma_const-1))*p_smooth.*dV_dtheta + ...
           (1/(gamma_const-1))*V_sel.*dp_dtheta;
ROHR_degnonvar = ROHRnonvar * (pi/180);

%% =============================================================
% VARIABLE GAMMA ROHR CALCULATION (unchanged from your version)
% NOTE: This section is left intact, so your variable-gamma option still works
%% =============================================================

T_intake = mean(Tint) + 273.15;

cycles_per_sec = n_test / 120;
mfuelinj = mfuel_mean / cycles_per_sec;

R_mixpre = (Runiv / MAir);

P_ivc = p_ref;
T_ivc = T_intake;
V_ivc = CylinderVolume(-180, Cyl);
m_trapped = (P_ivc * V_ivc) / (R_mixpre * T_ivc);

m_C_fuel = mfuelinj * FuelWeightPercent.C;
m_H_fuel = mfuelinj * FuelWeightPercent.H;
m_O_fuel = mfuelinj * FuelWeightPercent.O;
n_C_fuel = m_C_fuel / M_atom.C;
n_H_fuel = m_H_fuel / M_atom.H;
n_O_fuel = m_O_fuel / M_atom.O;

O2_for_C = n_C_fuel;
O2_for_H = n_H_fuel / 4;
O2_from_fuel = n_O_fuel / 2;
O2_consumed = O2_for_C + O2_for_H - O2_from_fuel;
CO2_produced = n_C_fuel;
H2O_produced = n_H_fuel / 2;

m_species = Yair * m_trapped;
n_species = m_species ./ Mi;
n_post = n_species;
n_post(1) = 0;
n_post(2) = n_post(2) - O2_consumed;
n_post(4) = n_post(4) + CO2_produced;
n_post(5) = n_post(5) + H2O_produced;

m_post = n_post .* Mi;
Ypost = m_post / sum(m_post);
M_mean_post = sum(n_post.*Mi) / sum(n_post);
R_mixpost = Runiv / M_mean_post;

m_total_pre  = m_trapped;
m_total_post = m_trapped + mfuelinj;

gamma_var = zeros(size(Ca_sel));
T = zeros(size(Ca_sel));

for idx = 1:length(Ca_sel)
    if Ca_sel(idx) < 0
        Y_current = Yair;
        R_current = R_mixpre;
        m_current = m_total_pre;
    else
        Y_current = Ypost;
        R_current = R_mixpost;
        m_current = m_total_post;
    end

    if (Ca_sel(idx)) > -180 && (Ca_sel(idx)) < 180
        T(idx) = (p_smooth(idx) * V_sel(idx)) / (m_current * R_current);
    else
        T(idx) = T_intake;
    end

    for i = 1:NSp
        Cp_mass(i) = CpNasa(T(idx), SpS(i)) / Mi(i);
        Cv_mass(i) = CvNasa(T(idx), SpS(i)) / Mi(i);
    end

    cp = Y_current * Cp_mass';
    cv = Y_current * Cv_mass';
    gamma_var(idx) = cp / cv;
end

ROHRvar = gamma_var./(gamma_var-1) .* p_smooth .* dV_dtheta + ...
          1./(gamma_var-1) .* V_sel .* dp_dtheta;
ROHR_degvar = ROHRvar * (pi/180);

% Diagnostics (kept)
figure;
subplot(2,1,1)
plot(Ca_sel, gamma_var, 'LineWidth', 1.5);
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

%% =============================================================
% SELECT ROHR VERSION (variable or constant)
%% =============================================================
if use_variable_gamma == true
    ROHR_for_CA = ROHR_degvar;       % full-length ROHR [J/deg]
else
    ROHR_for_CA = ROHR_degnonvar;    % full-length ROHR [J/deg]
end

%% =============================================================
% ROBUST CA10/50/90 FIX (REPLACES YOUR ORIGINAL BLOCK)
%% =============================================================

% Wider window to capture high-load burn tails
comb_mask = (Ca_sel >= -60) & (Ca_sel <= 180);
Ca_c = Ca_sel(comb_mask);
ROHR_c = ROHR_for_CA(comb_mask);

% Smooth ROHR in the wider window
ROHR_c_smooth = sgolayfilt(ROHR_c, 3, 31);

% Threshold-based burn start/end (2% peak)
thr = 0.02 * max(ROHR_c_smooth);
idx_start = find(ROHR_c_smooth > thr, 1, 'first');
idx_end   = find(ROHR_c_smooth > thr, 1, 'last');

if isempty(idx_start) || isempty(idx_end) || idx_end <= idx_start
    warning('CA10/50/90 thresholding failed. Falling back to full window.');
    idx_start = 1;
    idx_end = numel(Ca_c);
end

Ca_b   = Ca_c(idx_start:idx_end);
ROHR_b = ROHR_c_smooth(idx_start:idx_end);

% Cumulative heat release (monotone enforced)
cumHR = cumtrapz(Ca_b, ROHR_b);
cumHR = cumHR - min(cumHR);
cumHR = cummax(cumHR);

Q_total = max(cumHR);

CA10 = interp1(cumHR, Ca_b, 0.10 * Q_total, 'linear', 'extrap');
CA50 = interp1(cumHR, Ca_b, 0.50 * Q_total, 'linear', 'extrap');
CA90 = interp1(cumHR, Ca_b, 0.90 * Q_total, 'linear', 'extrap');

% Plot cumulative HR with markers
figure;
plot(Ca_b, cumHR,'LineWidth',1.8); hold on;
yl = ylim;
y_label_pos = yl(1) + 0.05*(yl(2)-yl(1));

xline(CA10,'--r'); xline(CA50,'--g'); xline(CA90,'--b');
text(CA10, y_label_pos, 'CA10', 'Color','r','FontWeight','bold','HorizontalAlignment','center');
text(CA50, y_label_pos, 'CA50', 'Color','g','FontWeight','bold','HorizontalAlignment','center');
text(CA90, y_label_pos, 'CA90', 'Color','b','FontWeight','bold','HorizontalAlignment','center');

xlabel('Crank angle [deg]');
ylabel('Cumulative Heat Release [J]');
title('Robust Cumulative Heat Release (CA10/50/90)');
grid on;

% Plot ROHR with CA markers (diagnostic)
figure;
plot(Ca_c, ROHR_c_smooth, 'LineWidth', 1.7); hold on;
xline(CA10,'--r'); xline(CA50,'--g'); xline(CA90,'--b');
xlabel('Crank angle [deg]');
ylabel('ROHR [J/deg]');
title('ROHR (smoothed) with CA10/50/90 markers');
grid on;
legend('ROHR','CA10','CA50','CA90','Location','best');

%% =============================================================
% --------- COMBUSTION KPIs ----------
[p_max, idx_pmax] = max(p_smooth);
CA_pmax = Ca_sel(idx_pmax);

[ROHR_max, idx_rmax] = max(ROHR_b);
CA_ROHR_max = Ca_b(idx_rmax);

comb_duration = CA90 - CA10;

fprintf('\n========== COMBUSTION KPIs (ROBUST) ==========\n');
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
