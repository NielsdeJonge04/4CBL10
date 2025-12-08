clear; clc;

% ============================================================
% Directories
% ============================================================
dieselDir     = 'Diesel';
hvoDir        = 'HVO';
gtlDir        = 'GTL';
hvo_dieselDir = 'HVO_Diesel';
gtl_dieselDir = 'GTL_Diesel';   %#ok<NASGU>  % not used yet

dieselEmisFile     = 'Diesel_emissions.csv';
hvoEmisFile        = 'HVO_emissions.csv';
GTLEmisFile        = 'GTL_emissions.csv';
hvo_dieselEmisFile = 'HVO&diesel_emissions.csv';
%GTL_dieselEmisFile = 'GTL&diesel_emissions.csv';

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
AFRst_GTL_diesel = 14.74;   %#ok<NASGU>

% Non-renewable carbon factors (0 = fully renewable, 1 = fossil)
nonren_diesel      = 1.0;
nonren_hvo         = 0.0;   % pure renewable (adjust if needed)
nonren_GTL         = 1.0;
nonren_HVO_diesel  = 0.5;   % 50/50 blend example (adjust if needed)
%nonren_GTL_diesel = 1.0;

% ============================================================
% Engine & operating conditions for brake-specific metrics
% ============================================================
mm   = 1e-3;
B    = 104*mm;     % bore [m]
S    = 85*mm;      % stroke [m]
Vd   = pi*(B^2)/4 * S;   % displaced volume per cycle [m^3] (single cyl)
n_engine = 1500;   % [rpm] engine speed (constant in tests)

% ============================================================
% Gas properties & constants
% ============================================================
M_mix = 29e-3;      % kg/mol (average exhaust gas)
M_CO  = 28.01e-3;   % kg/mol
M_CO2 = 44.01e-3;   % kg/mol
M_HC  = 44.1e-3;    % kg/mol (propane-equivalent)
M_NOx = 46.0e-3;    % kg/mol

R_u   = 8.314;      % J/(mol K)
p_exh = 1e5;        % Pa (assumed atmospheric)

% ============================================================
% Compute results (g/MJ + brake-specific emissions)
% ============================================================
resultsDiesel = compute_g_per_MJ_with_soot_and_BS( ...
    dieselDir, dieselEmisFile, LHV_diesel, AFRst_diesel, ...
    M_mix, M_CO, M_CO2, M_NOx, M_HC, R_u, p_exh, ...
    Vd, n_engine, nonren_diesel);

resultsHVO = compute_g_per_MJ_with_soot_and_BS( ...
    hvoDir, hvoEmisFile, LHV_hvo, AFRst_hvo, ...
    M_mix, M_CO, M_CO2, M_NOx, M_HC, R_u, p_exh, ...
    Vd, n_engine, nonren_hvo);

resultsGTL = compute_g_per_MJ_with_soot_and_BS( ...
    gtlDir, GTLEmisFile, LHV_GTL, AFRst_GTL, ...
    M_mix, M_CO, M_CO2, M_NOx, M_HC, R_u, p_exh, ...
    Vd, n_engine, nonren_GTL);

resultsHVODiesel = compute_g_per_MJ_with_soot_and_BS( ...
    hvo_dieselDir, hvo_dieselEmisFile, LHV_HVO_diesel, AFRst_HVO_diesel, ...
    M_mix, M_CO, M_CO2, M_NOx, M_HC, R_u, p_exh, ...
    Vd, n_engine, nonren_HVO_diesel);

% resultsGTLDiesel = compute_g_per_MJ_with_soot_and_BS( ...
%     gtl_dieselDir, gtl_dieselEmisFile, LHV_GTL_diesel, AFRst_GTL_diesel, ...
%     M_mix, M_CO, M_CO2, M_NOx, M_HC, R_u, p_exh, ...
%     Vd, n_engine, nonren_GTL_diesel);

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
%writetable(resultsGTLDiesel, 'GTL_Diesel_g_per_MJ.csv');











% ============================================================
% FUNCTION: compute g/MJ + soot + brake-specific emissions
% ============================================================

function results = compute_g_per_MJ_with_soot_and_BS( ...
    sdaqDir, emisFile, LHV_MJkg, AFRst, ...
    M_mix, M_CO, M_CO2, M_NOx, M_HC, R_u, p_exh, ...
    Vd, n_engine, nonren_factor)

    % Read emissions table (one row per operating point)
    emis = readtable(emisFile);

    % Slow data files (for mfuel, Tint, Texh, Pint)
    sdaqFiles = dir(fullfile(sdaqDir, '*sdaq*.txt'));
    nOP = height(emis);

    % Preallocate
    opName   = strings(nOP,1);
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

    % GWP values (same as in single-point KPI script)
    GWP20_CH4  = 79;
    GWP100_CH4 = 27.2;

    bara = 1e5; % [Pa/bar]

    for k = 1:nOP

        % --------------------------------------------------------
        % Load slow data for this operating point
        % --------------------------------------------------------
        data = readmatrix(fullfile(sdaqDir, sdaqFiles(k).name));

        mfuel_gs(k) = mean(data(:,1));          % g/s
        Texh_K      = mean(data(:,3)) + 273.15; % exhaust T [K]

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
        x_HC  = HC_ppm  / 1e6;
        x_NOx = NOx_ppm / 1e6;

        % --------------------------------------------------------
        % Gas species mass flows [kg/s]
        % --------------------------------------------------------
        mdot_CO   = mexh_kg_s * x_CO  * (M_CO  / M_mix);
        mdot_CO2  = mexh_kg_s * x_CO2 * (M_CO2 / M_mix);
        mdot_HC   = mexh_kg_s * x_HC  * (M_HC  / M_mix);
        mdot_NOx  = mexh_kg_s * x_NOx * (M_NOx / M_mix);
        mdot_NNox = mdot_NOx * ((21-15)/21 - O2_vol);

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
        rho_exh   = p_exh * M_mix / (R_u * Texh_K);   % [kg/m^3]
        Vdot_exh  = mexh_kg_s / rho_exh;              % [m^3/s]

        mdot_soot_mg_s = C_soot_mg_m3 * Vdot_exh;     % [mg/s]
        EI_soot(k)     = (mdot_soot_mg_s / 1000) / Edot_MJ_s;  % [g/MJ]

        % --------------------------------------------------------
        % GHG (gas-phase only) in g/MJ (same as before)
        % --------------------------------------------------------
        CO2eq_from_CO_g_per_MJ = EI_CO(k) * (M_CO2 / M_CO);
        GHG20(k)  = EI_CO2(k) + CO2eq_from_CO_g_per_MJ + EI_HC(k) * GWP20_CH4;
        GHG100(k) = EI_CO2(k) + CO2eq_from_CO_g_per_MJ + EI_HC(k) * GWP100_CH4;

        % --------------------------------------------------------
        % Brake power from nominal IMEP in operating point name
        % (assume "CA12_IMEP1_5" etc → BMEP = 1.5 bar)
        % --------------------------------------------------------
        nameStr = emis.Name{k};  % e.g. 'CA12_IMEP1_5'
        % Extract substring after 'IMEP'
        idxIMEP = strfind(nameStr, 'IMEP');
        if isempty(idxIMEP)
            error('Operating point name "%s" does not contain "IMEP".', nameStr);
        end
        imep_sub = nameStr(idxIMEP+4:end);        % e.g. '1_5'
        imep_bar = str2double(strrep(imep_sub, '_', '.'));  % [bar]
        BMEP_Pa  = imep_bar * bara;              % [Pa]

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
        fact_g = 3.6e9 / P_b_W;
        fact_mg = 3.6e6 / P_b_W;

        BS_CO2_g_per_kWh(k)       = mass_CO2       * fact_g;
        BS_CO2nonren_g_per_kWh(k) = mass_CO2*nonren_factor * fact_g;
        BS_CO_g_per_kWh(k)        = mass_CO        * fact_g;
        BS_NOx_g_per_kWh(k)       = mass_NOx       * fact_g;
        BS_GHG20_g_per_kWh(k)     = GHG20_mass     * fact_g;
        BS_GHG100_g_per_kWh(k)    = GHG100_mass    * fact_g;
        BS_soot_mg_per_kWh(k)     = mass_soot      * fact_mg;

        % Brake-specific fuel consumption [g/kWh]
        BSFC_g_per_kWh(k) = (mfuel_kg_s * fact_g);  % mfuel_kg_s * (3.6e9/P_b_W)

    end

    % ------------------------------------------------------------
    % Build output table
    % ------------------------------------------------------------
    results = table(...
        opName, mfuel_gs, lambda, ...
        EI_CO, EI_CO2, EI_HC, EI_NOx, EI_soot, ...
        GHG20, GHG100, ...
        P_b_kW, ...
        BSFC_g_per_kWh, ...
        BS_CO2_g_per_kWh, BS_CO2nonren_g_per_kWh, ...
        BS_CO_g_per_kWh, BS_NOx_g_per_kWh, ...
        BS_GHG20_g_per_kWh, BS_GHG100_g_per_kWh, ...
        BS_soot_mg_per_kWh, mdot_NNox, ...
        'VariableNames', { ...
            'OperatingPoint','mFuel_g_per_s','lambda', ...
            'EI_CO_g_per_MJ','EI_CO2_g_per_MJ', ...
            'EI_HC_g_per_MJ','EI_NOx_g_per_MJ', ...
            'EI_soot_g_per_MJ', ...
            'GHG20_gCO2eq_per_MJ','GHG100_gCO2eq_per_MJ', ...
            'P_b_kW', ...
            'BSFC_g_per_kWh', ...
            'BS_CO2_g_per_kWh','BS_CO2nonren_g_per_kWh', ...
            'BS_CO_g_per_kWh','BS_NOx_g_per_kWh', ...
            'BS_GHG20_gCO2eq_per_kWh','BS_GHG100_gCO2eq_per_kWh', ...
            'BS_soot_mg_per_kWh', 'mdot_NNox' ...
        } ...
    );

end

