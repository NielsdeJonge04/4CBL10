clear; clc;

% ============================================================
% Directories
% ============================================================
dieselDir     = 'Data/Diesel';
hvoDir        = 'Data/HVO';
gtlDir        = 'Data/GTL';
hvo_dieselDir = 'Data/HVO_Diesel';
gtl_dieselDir = 'Data/GTL_Diesel';   

dieselEmisFile     = 'Data/Emissions_csv/Diesel_emissions.csv';
hvoEmisFile        = 'Data/Emissions_csv/HVO_emissions.csv';
GTLEmisFile        = 'Data/Emissions_csv/GTL_emissions.csv';
hvo_dieselEmisFile = 'Data/Emissions_csv/HVO&diesel_emissions.csv';
GTL_dieselEmisFile = 'Data/Emissions_csv/GTL&diesel_emissions.csv';

% ============================================================
% Fuel properties  (LHV in MJ/kg, AFRst in kg air / kg fuel)
% ============================================================
LHV_diesel   = 43.00;   % MJ/kg
AFRst_diesel = 14.50;

LHV_hvo      = 44.10;   % MJ/kg
AFRst_hvo    = 14.97;

LHV_GTL      = 44.00;   % MJ/kg
AFRst_GTL    = 14.98;

LHV_HVO_diesel   = 43.55;   % MJ/kg
AFRst_HVO_diesel = 14.74;

LHV_GTL_diesel   = 43.50;   % MJ/kg
AFRst_GTL_diesel = 14.74;   

% Non-renewable carbon factors (0 = fully renewable, 1 = fossil)
nonren_diesel      = 1.0;
nonren_hvo         = 0.0;   % pure renewable (adjust if needed)
nonren_GTL         = 1.0;
nonren_HVO_diesel  = 0.5;   % 50/50 blend example (adjust if needed)
nonren_GTL_diesel = 1.0;

% ============================================================
% Engine & operating conditions for brake-specific metrics
% ============================================================
mm   = 1e-3;
Cyl.Bore             = 104*mm;
Cyl.Stroke           = 85*mm;
Cyl.CompressionRatio = 21.5;
Cyl.ConRod           = 136.5*mm;
Cyl.TDCangle         = 0;          % TDC reference
v_min = CylinderVolume(0, Cyl);
v_max = CylinderVolume(180, Cyl);
Vd = v_max - v_min;
n_engine = 1500;   % [rpm] engine speed (constant in tests)

% ============================================================
% Gas properties & constants
% ============================================================
M_mix = 29e-3;  % kg/mol (average exhaust gas)
M_exhaust = 29e-3;  % kg/mol (average exhaust gas)
M_CO  = 28.01e-3;   % kg/mol
M_CO2 = 44.01e-3;   % kg/mol
M_HC  = 44.1e-3;    % kg/mol (propane-equivalent)
M_NOx = 46.0e-3;    % kg/mol

R_u   = 8.314;      % J/(mol K)

% ============================================================
% Compute results (g/MJ + brake-specific emissions)
% ============================================================
resultsDiesel = compute_g_per_MJ_with_soot_and_BS( ...
    dieselDir, dieselEmisFile, LHV_diesel, AFRst_diesel, ...
    M_mix, M_CO, M_CO2, M_NOx, M_HC, M_exhaust, R_u, Cyl, ...
    Vd, n_engine, nonren_diesel);

resultsHVO = compute_g_per_MJ_with_soot_and_BS( ...
    hvoDir, hvoEmisFile, LHV_hvo, AFRst_hvo, ...
    M_mix, M_CO, M_CO2, M_NOx, M_HC, M_exhaust, R_u, Cyl, ...
    Vd, n_engine, nonren_hvo);

resultsGTL = compute_g_per_MJ_with_soot_and_BS( ...
    gtlDir, GTLEmisFile, LHV_GTL, AFRst_GTL, ...
    M_mix, M_CO, M_CO2, M_NOx, M_HC, M_exhaust, R_u, Cyl, ...
    Vd, n_engine, nonren_GTL);

resultsHVODiesel = compute_g_per_MJ_with_soot_and_BS( ...
    hvo_dieselDir, hvo_dieselEmisFile, LHV_HVO_diesel, AFRst_HVO_diesel, ...
    M_mix, M_CO, M_CO2, M_NOx, M_HC, M_exhaust, R_u, Cyl, ...
    Vd, n_engine, nonren_HVO_diesel);

resultsGTLDiesel = compute_g_per_MJ_with_soot_and_BS( ...
    gtl_dieselDir, GTL_dieselEmisFile, LHV_GTL_diesel, AFRst_GTL_diesel, ...
    M_mix, M_CO, M_CO2, M_NOx, M_HC, M_exhaust, R_u, Cyl, ...
    Vd, n_engine, nonren_GTL_diesel);

% ============================================================
% Display in command window
% ============================================================
disp('=== Diesel emission factors [g/MJ fuel] + BS emissions ===');
disp(resultsDiesel);

disp('=== HVO emission factors [g/MJ fuel] + BS emissions ===');
disp(resultsHVO);

disp('=== GTL emission factors [g/MJ fuel] + BS emissions ===');
disp(resultsGTL);

disp('=== HVO/Diesel emission factors [g/MJ fuel] + BS emissions ===');
disp(resultsHVODiesel);

% disp('=== GTL/Diesel emission factors [g/MJ fuel] + BS emissions ===');
% disp(resultsGTLDiesel);

% ============================================================
% Write to CSV (for LaTeX tables etc.)
% ============================================================
writetable(resultsDiesel,    'Diesel_g_per_MJ.csv');
writetable(resultsHVO,       'HVO_g_per_MJ.csv');
writetable(resultsGTL,       'GTL_g_per_MJ.csv');
writetable(resultsHVODiesel, 'HVO_Diesel_g_per_MJ.csv');
writetable(resultsGTLDiesel, 'GTL_Diesel_g_per_MJ.csv');











% ============================================================
% FUNCTION: compute g/MJ + soot + brake-specific emissions
% ============================================================

function results = compute_g_per_MJ_with_soot_and_BS( ...
    sdaqDir, emisFile, LHV_MJkg, AFRst, ...
    M_mix, M_CO, M_CO2, M_NOx, M_HC, M_exhaust, R_u, Cyl, ...
    Vd, n_engine, nonren_factor)

    % Read emissions table (one row per operating point)
    emis = readtable(emisFile);

    % Slow data files (for mfuel, Tint, Texh, Pint)
    sdaqFiles = dir(fullfile(sdaqDir, '*sdaq*.txt'));
    fdaqFiles = dir(fullfile(sdaqDir,'*fdaq*.txt')); % fast data file
    nOP = height(emis);

    % Preallocate
    opName   = strings(nOP,1);
    BMEP_BA = zeros(nOP,1);
    mfuel_gs = zeros(nOP,1);
    lambda   = zeros(nOP,1);

    EI_CO   = zeros(nOP,1);
    EI_CO2  = zeros(nOP,1);
    EI_HC   = zeros(nOP,1);
    EI_NOx  = zeros(nOP,1);
    EI_soot = zeros(nOP,1);  

    GHG20  = zeros(nOP,1);
    GHG100 = zeros(nOP,1);

    % Brake-specific emissions
    P_b_kW                 = zeros(nOP,1);
    BS_CO2_g_per_kWh       = zeros(nOP,1);
    BS_CO2nonren_g_per_kWh = zeros(nOP,1);
    BS_CO_g_per_kWh        = zeros(nOP,1);
    BS_NOx_g_per_kWh       = zeros(nOP,1);
    BS_GHG20_g_per_kWh     = zeros(nOP,1);
    BS_GHG100_g_per_kWh    = zeros(nOP,1);
    BS_soot_mg_per_kWh     = zeros(nOP,1);
    BSFC_g_per_kWh         = zeros(nOP,1);
    BS_NNOx_g_per_kWh      = zeros(nOP,1);

    % GWP values (same as in single-point KPI script)
    GWP20_CH4  = 79;
    GWP100_CH4 = 27.2;

    bara = 1e5; % [Pa/bar]

    V_molar_STP  = 22.414e-3;  % [m^3/mol] molar volume at STP


    for k = 1:nOP

        % --------------------------------------------------------
        % Load slow data for this operating point
        % --------------------------------------------------------
        sdata = readmatrix(fullfile(sdaqDir, sdaqFiles(k).name));
        
        mfuel_gs(k) = mean(sdata(:,1));          % g/s
        Texh_mean_K      = mean(sdata(:,3)) + 273.15; % exhaust T [K]
        Pint  = sdata(:,4) * bara;      % [Pa]     (bar(a) -> Pa)
        % --------------------------------------------------------
        % Load fast data for this operating point
        % --------------------------------------------------------

        Fdata = readmatrix(fullfile(sdaqDir, fdaqFiles(k).name));

        Ca_raw   = Fdata(:,1);          % crank angle [deg]
        p_raw    = Fdata(:,2) * bara;   % in-cylinder pressure [Pa] (NOT pegged)
        
        % Reshape into cycles
        dCa = 0.2;                     % deg/step
        pts_per_cycle = 720/dCa;       % 4-stroke -> 720 deg/cycle
        Ncycles = floor(length(Ca_raw)/pts_per_cycle);
        
        Ca   = reshape(Ca_raw ,[],Ncycles);
        p    = reshape(p_raw  ,[],Ncycles);

        % DRIFT CORRECTION (PEGGING at BDC) - Applied to all cycles
        [~, idxBDC] = min(abs(Ca(:,1) + 180));  % BDC ~ -180 deg CA
        
        p_ref = mean(Pint);                     % intake pressure reference [Pa]
        
        p_corr_all = zeros(size(p));
        for i = 1:Ncycles
            p_meas = p(idxBDC, i);             % measured p at BDC (cycle i)
            p_corr_all(:, i) = p(:, i) + (p_ref - p_meas);
        end
        
        % PRESSURE SMOOTHING (noise reduction)
        % First average all pegged cycles:
        p_avg = mean(p_corr_all, 2);           % [Pa]
        
        % Then apply Savitzky-Golay filter:
        p_smooth = sgolayfilt(p_avg, 3, 81);   % (order 3, window 81 samples)

        % --------------------------------------------------------
        % Emission inputs from CSV
        % --------------------------------------------------------
        opName(k) = string(emis.Name{k});
        lambda(k) = emis.lambda(k);

        CO_volpct  = emis.CO_volpct(k);
        CO2_volpct = emis.CO2_volpct(k);
        HC_ppm     = emis.HC_ppm(k);
        NOx_ppm    = emis.NOx_ppm(k);
        FSN        = emis.FSN(k);
        O2_volpct  = emis.O2_volpct(k);

        % --------------------------------------------------------
        % Mass & energy flows
        % --------------------------------------------------------
        mfuel_kg_s = mfuel_gs(k) / 1000;     % [kg/s]
        AFR_actual = lambda(k) * AFRst;      % [-]

        mair_kg_s = AFR_actual * mfuel_kg_s; % [kg/s]
        mexh_kg_s = mair_kg_s + mfuel_kg_s;  % [kg/s]

        % Fuel energy rate [MJ/s]
        Edot_MJ_s = mfuel_kg_s * LHV_MJkg;   % [MJ/s]

        % --------------------------------------------------------
        % Dry-gas concentrations
        % --------------------------------------------------------
        x_CO  = CO_volpct  / 100;   % vol fraction
        x_CO2 = CO2_volpct / 100;
        x_HC  = HC_ppm  / 1e5; % ppm*10
        x_NOx = NOx_ppm / 1e6;
        
        % --------------------------------------------------------
        % temperature corrected exhaust density and volume
        % --------------------------------------------------------
        T_correction = 273.15 / Texh_mean_K;

        % Exhaust density and volumetric flow
        rho_exhaust = p_ref * M_exhaust / (R_u * Texh_mean_K);% [kg/m^3]
        V_exhaust   = mexh_kg_s / rho_exhaust;                   % [m^3/s]

        % --------------------------------------------------------
        % Gas species mass flows [kg/s]
        % --------------------------------------------------------
        mdot_CO   = x_CO  * V_exhaust * (M_CO / V_molar_STP) * T_correction;     % [kg/s]
        mdot_CO2  = x_CO2 * V_exhaust * (M_CO2 / V_molar_STP) * T_correction;    % [kg/s]
        mdot_HC   = x_HC  * V_exhaust * (M_HC / V_molar_STP) * T_correction;     % [kg/s]
        mdot_NOx  = x_NOx * V_exhaust * (M_NOx / V_molar_STP) * T_correction;    % [kg/s]
        mdot_NNox = mdot_NOx * ((21-15)/21 - O2_volpct);

        % Emission indices [g/MJ fuel]
        EI_CO(k)  = (mdot_CO  * 1000) / Edot_MJ_s;
        EI_CO2(k) = (mdot_CO2 * 1000) / Edot_MJ_s;
        EI_HC(k)  = (mdot_HC  * 1000) / Edot_MJ_s;
        EI_NOx(k) = (mdot_NOx * 1000) / Edot_MJ_s;

        % --------------------------------------------------------
        % Soot from FSN  (AVL/ISO-like correlation)
        % --------------------------------------------------------
        C_soot_mg_m3 = (4.95/0.405) * FSN * exp(0.38 * FSN);  % [mg/m^3]

        % Exhaust density & volume flow
        mdot_soot_mg_s = C_soot_mg_m3 * V_exhaust;
        EI_soot(k)     = (mdot_soot_mg_s / 1000) / Edot_MJ_s;  % [g/MJ]

        % --------------------------------------------------------
        % GHG (gas-phase only) in g/MJ (same as before)
        % --------------------------------------------------------
        CO2eq_from_CO_g_per_MJ = EI_CO(k) * (M_CO2 / M_CO);
        GHG20(k)  = EI_CO2(k) + CO2eq_from_CO_g_per_MJ + EI_HC(k) * GWP20_CH4;
        GHG100(k) = EI_CO2(k) + CO2eq_from_CO_g_per_MJ + EI_HC(k) * GWP100_CH4;

        % --------------------------------------------------------
        % Brake power from nominal IMEP in operating point name
        % --------------------------------------------------------
        % IMEP (using smoothed average pressure)
        Ca_rad = deg2rad(Ca(:,1));  % Convert crank angle to radians
        V = CylinderVolume(Ca(:,1), Cyl);  % Calculate volume at each crank angle
        W_ind = trapz(V, p_smooth);  % [J/cycle] - integrate p dV
        
        BMEP_Pa  = W_ind / Vd;     % [Pa]
        BMEP_BA(k) = BMEP_Pa / bara;% [Bar]

        % 4-stroke single-cylinder: P_b = BMEP * Vd * n/120
        P_b_W = BMEP_Pa * Vd * (n_engine / 120); % [W]
        P_b_kW(k) = P_b_W / 1e3;

        % --------------------------------------------------------
        % Brake-specific emissions (per kWh)
        % --------------------------------------------------------
        % Species mass flows [kg/s]
        mass_CO   = mdot_CO;
        mass_CO2  = mdot_CO2;
        mass_HC   = mdot_HC;
        mass_NOx  = mdot_NOx;
        mass_soot = mdot_soot_mg_s;  % [mg/s]

        % GHG mass flows [kg/s]
        GHG20_mass  = mass_CO2 + mass_CO*(M_CO2/M_CO) + mass_HC*GWP20_CH4;
        GHG100_mass = mass_CO2 + mass_CO*(M_CO2/M_CO) + mass_HC*GWP100_CH4;

        % Conversion factors:
        %   - species kg/s → g/kWh: * (3.6e9 / P_b_W)
        %   - soot mg/s   → mg/kWh: * (3.6e6 / P_b_W)
        fact_g  = 3.6e9 / P_b_W;
        fact_mg = 3.6e6 / P_b_W;

        BS_CO2_g_per_kWh(k)       = mass_CO2       * fact_g;
        BS_CO2nonren_g_per_kWh(k) = mass_CO2*nonren_factor * fact_g;
        BS_CO_g_per_kWh(k)        = mass_CO        * fact_g;
        BS_NOx_g_per_kWh(k)       = mass_NOx       * fact_g;
        BS_GHG20_g_per_kWh(k)     = GHG20_mass     * fact_g;
        BS_GHG100_g_per_kWh(k)    = GHG100_mass    * fact_g;
        BS_soot_mg_per_kWh(k)     = mass_soot      * fact_mg;
        BS_NNOx_g_per_kWh(k)      = mdot_NNox      * fact_g;

        % Brake-specific fuel consumption [g/kWh]
        BSFC_g_per_kWh(k) = (mfuel_kg_s * fact_g);  % mfuel_kg_s * (3.6e9/P_b_W)

    end

    % ------------------------------------------------------------
    % Build output table
    % ------------------------------------------------------------
    results = table(...
            opName, BMEP_BA ,mfuel_gs, lambda, ...
            EI_CO, EI_CO2, EI_HC, EI_NOx, EI_soot, ...
            GHG20, GHG100, ...
            P_b_kW, ...
            BSFC_g_per_kWh, ...
            BS_CO2_g_per_kWh, BS_CO2nonren_g_per_kWh, ...
            BS_CO_g_per_kWh, BS_NOx_g_per_kWh, ...
            BS_GHG20_g_per_kWh, BS_GHG100_g_per_kWh, ...
            BS_soot_mg_per_kWh, BS_NNOx_g_per_kWh, ...
    'VariableNames', { ...
    'OperatingPoint','BMEP_BA' ,'mFuel_g_per_s','lambda', ...
    'EI_CO_g_per_MJ','EI_CO2_g_per_MJ', ...
    'EI_HC_g_per_MJ','EI_NOx_g_per_MJ', ...
    'EI_soot_g_per_MJ', ...
    'GHG20_gCO2eq_per_MJ','GHG100_gCO2eq_per_MJ', ...
    'P_b_kW', ...
    'BSFC_g_per_kWh', ...
    'BS_CO2_g_per_kWh','BS_CO2nonren_g_per_kWh', ...
    'BS_CO_g_per_kWh','BS_NOx_g_per_kWh', ...
    'BS_GHG20_gCO2eq_per_kWh','BS_GHG100_gCO2eq_per_kWh', ...
    'BS_soot_mg_per_kWh', 'BS_NNOx_g_per_kWh' ...
            });

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
