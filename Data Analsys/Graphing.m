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

% For this example file we assume a single steady point:
mfuel_mean = mean(mfuel);      % average fuel mass flow [kg/s]

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
% IMEP (using smoothed avg pressure)
W_ind = trapz(V_sel, p_smooth);              % [J/cycle]
IMEP  = W_ind / (max(V_sel)-min(V_sel));     % [Pa] = [J/m^3]

fprintf('Indicated work: %.2f J/cycle\n', W_ind);
fprintf('IMEP: %.2f bar\n', IMEP/bara);

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
% ======= PERFORMANCE KPIs =======

n_test = 1500;        % [rpm] measurement speed

% Indicated power: work per cycle * cycles per second
cycles_per_sec = n_test / 120;        % 4-stroke: n/120
P_i = W_ind * cycles_per_sec;         % [W]

% Fuel energy input (using diesel LHV)
LHV    = 42.7e6;                      % [J/kg]
E_fuel = mfuel_mean * LHV;            % [W]

% ---- Estimate brake power from assumed mechanical efficiency ----
eta_mech_assumed = 0.85;              % typical value for small diesel
P_b   = eta_mech_assumed * P_i;       % [W] estimated brake power
P_b_kW = P_b / 1e3;                   % [kW]

% Efficiencies & BSFC
ITE      = P_i / E_fuel;              % indicated thermal efficiency
BTE      = P_b / E_fuel;              % brake thermal efficiency
eta_mech = P_b / P_i;                 % (= eta_mech_assumed, by construction)
BSFC     = (mfuel_mean / P_b) * 3.6e9; % [g/kWh]

fprintf('\n========== PERFORMANCE KPIs ==========\n');
fprintf('Indicated work  W_ind = %6.2f J/cycle\n', W_ind);
fprintf('IMEP                = %6.2f bar\n', IMEP/bara);
fprintf('Indicated power P_i = %6.2f kW\n', P_i/1e3);
fprintf('Brake power P_b     = %6.2f kW (estimated, η_mech = %.2f)\n', P_b/1e3, eta_mech_assumed);
fprintf('Fuel power E_fuel   = %6.2f kW\n', E_fuel/1e3);
fprintf('ITE                 = %6.3f [-]\n', ITE);
fprintf('BTE                 = %6.3f [-]\n', BTE);
fprintf('Mechanical eff.     = %6.3f [-]\n', eta_mech);
fprintf('BSFC                = %6.1f g/kWh\n', BSFC);

%% =============================================================
% KPI PLOTS (single operating point) (most of these plots are redundant)

% 1) Efficiency bar plot
figure;
eff_names = categorical({'ITE','BTE','\eta_{mech}'});
eff_names = reordercats(eff_names, {'ITE','BTE','\eta_{mech}'});
eff_vals  = [ITE, BTE, eta_mech];

bar(eff_names, eff_vals);
ylim([0 1]);
ylabel('Efficiency [-]');
title('Engine Efficiencies at Operating Point');
grid on;

% 2) BSFC plot (single operating point)
figure;
bar(1, BSFC);
set(gca,'XTick',1,'XTickLabel',{sprintf('P_b = %.1f kW', P_b/1000)});
ylabel('BSFC [g/kWh]');
title('Brake Specific Fuel Consumption');
grid on;

% 3) Combustion phasing bar plot
figure;
phase_names = categorical({'CA10','CA50','CA90'});
phase_names = reordercats(phase_names, {'CA10','CA50','CA90'});
phase_vals  = [CA10, CA50, CA90];

bar(phase_names, phase_vals);
ylabel('Crank angle [deg]');
title('Combustion Phasing Angles');
grid on;

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
% EMISSIONS KPIs 
% Fill these from your lab emissions table (make sure lengths match load_perc)

load_perc = [30, 50, 70];      % [%] engine/genset load

NOx_perc = [];   % [%] NOx  (fill with 3 values)
CO_perc  = [];   % [%] CO   (fill with 3 values)
CO2_perc = [];   % [vol.%] CO2 (3 values)
THC_perc = [];   % [%] THC  (3 values)

% -------- NOx vs load --------
if numel(NOx_perc) == numel(load_perc)
    figure;
    plot(load_perc, NOx_perc, '-o', 'LineWidth', 1.8, 'MarkerSize', 8);
    xlabel('Engine load [%]');
    ylabel('NO_x [%]');
    title('NO_x emissions vs load');
    grid on;
else
    warning('Fill NOx_perc with data (same length as load_perc) to enable NOx plot.');
end

% -------- CO vs load --------
if numel(CO_perc) == numel(load_perc)
    figure;
    plot(load_perc, CO_perc, '-s', 'LineWidth', 1.8, 'MarkerSize', 8);
    xlabel('Engine load [%]');
    ylabel('CO [%]');
    title('CO emissions vs load');
    grid on;
else
    warning('Fill CO_perc with data (same length as load_perc) to enable CO plot.');
end

% -------- CO2 vs load --------
if numel(CO2_perc) == numel(load_perc)
    figure;
    plot(load_perc, CO2_perc, '-^', 'LineWidth', 1.8, 'MarkerSize', 8);
    xlabel('Engine load [%]');
    ylabel('CO_2 [vol.%]');
    title('CO_2 emissions vs load');
    grid on;
else
    warning('Fill CO2_perc with data (same length as load_perc) to enable CO2 plot.');
end

% -------- THC vs load (optional) --------
if numel(THC_perc) == numel(load_perc)
    figure;
    plot(load_perc, THC_perc, '-d', 'LineWidth', 1.8, 'MarkerSize', 8);
    xlabel('Engine load [%]');
    ylabel('THC [%]');
    title('Unburned HC vs load');
    grid on;
else
    warning('Fill THC_perc with data (same length as load_perc) to enable THC plot.');
end

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
