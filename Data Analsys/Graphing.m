%% ===============================================================
%  ENGINE COMBUSTION + DATA ANALYSIS (FINAL – FIXED VERSION)
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
Cyl.Bore   = 104*mm;
Cyl.Stroke = 85*mm;
Cyl.CompressionRatio = 21.5;
Cyl.ConRod = 136.5*mm;
Cyl.TDCangle = 0;

%% ======================
% LOAD FAST DATA (fdaq)
fdaq_file = fullfile('Data','20251120_0000001_example_fdaq.txt');
fdaq = table2array(readtable(fdaq_file));

Ca_raw   = fdaq(:,1);
p_raw    = fdaq(:,2) * bara;
Iinj_raw = fdaq(:,3);

%% Reshape into cycles
dCa = 0.2;
pts_per_cycle = 720/dCa;
Ncycles = floor(length(Ca_raw)/pts_per_cycle);

Ca   = reshape(Ca_raw ,[],Ncycles);
p    = reshape(p_raw  ,[],Ncycles);
Iinj = reshape(Iinj_raw,[],Ncycles);

iselect = 10;

%% ======================
% LOAD SLOW DATA (sdaq)
sdaq_file = fullfile('Data','20251120_0000001_example_sdaq.txt');
sdaq = table2array(readtable(sdaq_file));

mfuel = sdaq(:,1) * 1e-3;    % g/s → kg/s
Tint  = sdaq(:,2);
Texh  = sdaq(:,3);
Pint  = sdaq(:,4) * bara;    % bar → Pa

%% =============================================================
% DRIFT CORRECTION (PEGGING at BDC) - Applied to all cycles
[~, idxBDC] = min(abs(Ca(:,1) + 180));

p_ref = mean(Pint);           % intake pressure reference

% Apply drift correction to all cycles
p_corr_all = zeros(size(p));
for i = 1:Ncycles
    p_meas = p(idxBDC, i);       % measured value at BDC for each cycle
    p_corr_all(:, i) = p(:, i) + (p_ref - p_meas);
end

%% =============================================================
% PRESSURE SMOOTHING (noise reduction)
% First take average of all cycles to reduce noise
p_avg = mean(p_corr_all, 2);  % Average across all cycles

% Then apply Savitzky-Golay filter to the averaged signal
p_smooth = sgolayfilt(p_avg, 3, 81);

% Extract the selected cycle for consistency with rest of code
Ca_sel = Ca(:, iselect);
p_sel = p(:, iselect);
p_corr = p_corr_all(:, iselect);  % Corrected pressure for selected cycle


%% ======================
% PLOT: p–CA
figure;
plot(Ca, p/bara,'LineWidth',1);
hold on;
plot(Ca_sel, p_sel/bara,'r','LineWidth',2);
xlabel('Crank angle [deg]');
ylabel('Pressure [bar]');
title('p–CA Diagram');
grid on;

%% ======================
% CYLINDER VOLUME
V_sel = CylinderVolume(Ca_sel, Cyl);

%% ======================
% PLOT: p–V
figure;

subplot(2,1,1)
plot(V_sel/dm^3, p_sel/bara,'LineWidth',1.5);
xlabel('V [dm^3]');
ylabel('p [bar]');
title('p–V (Linear)');
grid on;


subplot(2,1,2)
mask = (V_sel>0) & (p_sel>0);
loglog(V_sel(mask)/dm^3, p_sel(mask)/bara,'LineWidth',1.5);
xlabel('V [dm^3]');
ylabel('p [bar]');
title('p–V (Log-Log)');
grid on;


%% =============================================================
% IMEP
W_ind = trapz(V_sel, p_smooth);
IMEP  = W_ind / (max(V_sel)-min(V_sel));

fprintf('Indicated work: %.2f J\n', W_ind);
fprintf('IMEP: %.2f bar\n', IMEP/bara);

%% =============================================================
% === FIXED ROHR CALCULATION ===
gamma  = 1.35;
dtheta = deg2rad(mean(diff(Ca_sel)));

dp_dtheta = gradient(p_smooth, dtheta);
dV_dtheta = gradient(V_sel    , dtheta);

ROHR = (gamma/(gamma-1))*p_smooth.*dV_dtheta + ...
       (1/(gamma-1))*V_sel.*dp_dtheta;

ROHR_deg = ROHR * (pi/180);

% === COMBUSTION WINDOW FIX ===
comb_mask = (Ca_sel >= -40) & (Ca_sel <= 120);

Ca_c = Ca_sel(comb_mask);
ROHR_c = ROHR_deg(comb_mask);

% Final smoothing (smaller window!)
ROHR_smooth = sgolayfilt(ROHR_c, 3, 31);

figure;
plot(Ca_c, ROHR_smooth,'LineWidth',1.7);
xlabel('Crank angle [deg]');
ylabel('ROHR [J/deg]');
title('Corrected Smoothed Rate of Heat Release');
grid on;
legend 

%% =============================================================
% === FIXED CUMULATIVE HEAT RELEASE ===
cumHR = cumtrapz(Ca_c, ROHR_smooth);

% Find the minimum point (lowest value) in cumHR
[min_val, min_idx] = min(cumHR);

% Shift cumHR so the minimum point is at zero
cumHR = cumHR - min_val;

% Use only the data from the minimum point onwards for CA10/50/90
cumHR_from_min = cumHR(min_idx:end);
Ca_c_from_min = Ca_c(min_idx:end);

Q_total = max(cumHR_from_min);

CA10 = interp1(cumHR_from_min, Ca_c_from_min, 0.10 * Q_total);
CA50 = interp1(cumHR_from_min, Ca_c_from_min, 0.50 * Q_total);
CA90 = interp1(cumHR_from_min, Ca_c_from_min, 0.90 * Q_total);

figure;
plot(Ca_c, cumHR,'LineWidth',1.8); hold on;

% Get y-axis limits to position labels at bottom
yl = ylim;
y_label_pos = yl(1) + 0.05*(yl(2)-yl(1));  % 5% from bottom

xline(CA10,'--r');
xline(CA50,'--g');
xline(CA90,'--b');
text(CA10, y_label_pos, 'CA10', 'Color', 'r', 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(CA50, y_label_pos, 'CA50', 'Color', 'g', 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(CA90, y_label_pos, 'CA90', 'Color', 'b', 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

xlabel('Crank angle [deg]');
ylabel('Cumulative Heat Release [J]');
title('Corrected Cumulative Heat Release');
grid on;

%% =============================================================
% INJECTOR CURRENT
figure;
plot(Ca_sel, Iinj(:,iselect),'LineWidth',1.5);
xlabel('Crank angle [deg]');
ylabel('Injector Current [A]');
title('Injector Current Profile');
grid on;

%% ======================
% PLOT: p–V (Smoothed)
figure;

subplot(2,1,1)
plot(V_sel/dm^3, p_smooth/bara,'LineWidth',1.8);
xlabel('V [dm^3]');
ylabel('p [bar]');
title('Smoothed p–V (Linear)');
grid on;

subplot(2,1,2)
mask = (V_sel>0) & (p_smooth>0);
loglog(V_sel(mask)/dm^3, p_smooth(mask)/bara,'LineWidth',1.8);
xlabel('V [dm^3]');
ylabel('p [bar]');
title('Smoothed p–V (Log-Log)');
grid on;


%% ======================
% PLOT: p–CA (Smoothed)
figure;
plot(Ca_sel, p_smooth/bara,'LineWidth',1.8);
xlabel('Crank angle [deg]');
ylabel('Pressure [bar]');
title('Smoothed p–CA Diagram');
grid on;

%% =============================================================
% SLOW SENSOR PLOTS

figure;
plot(mfuel,'LineWidth',1.5);
xlabel('Index');
ylabel('Fuel flow [kg/s]');
title('Fuel Mass Flow');
grid on;


figure;
plot(Tint,'LineWidth',1.5); hold on;
plot(Texh,'LineWidth',1.5);
legend('Intake','Exhaust');
xlabel('Index');
ylabel('Temperature [°C]');
title('Temperatures');
grid on;

figure;
plot(Pint/bara,'LineWidth',1.5);
xlabel('Index');
ylabel('Pressure [bar]');
title('Intake Pressure');
grid on;
 

%% =============================================================
% CYLINDER VOLUME FUNCTION
function V = CylinderVolume(Ca, Cyl)
    B = Cyl.Bore; 
    S = Cyl.Stroke;
    CR = Cyl.CompressionRatio;
    L = Cyl.ConRod;

    Ap = pi*(B^2)/4;
    Vd = Ap*S;
    Vc = Vd/(CR-1);

    theta = deg2rad(Ca - Cyl.TDCangle);
    r = S/2;

    under = L^2 - (r*sin(theta)).^2;
    under(under<0)=0;

    x = r*(1-cos(theta)) + (L - sqrt(under));
    V = Vc + Ap*x;
end