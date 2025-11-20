%% Info
% Toerental: 1500 RPM
% SOA van 4.2º voor TDC
% Resolutie van 0.2º CA
% Data voor 69 cycles (maximale van de Smetec, de OGO gensets kunnen in principe “onbeperkt” aan)
% 
%% init
clear; clc;close;
addpath( "Functions","Nasa");
%% Units
mm      = 1e-3;dm=0.1;
bara    = 1e5;   %bar to pascal
MJ      = 1e6;
kWhr    = 1000*3600;
volperc = 0.01; % Emissions are in volume percentages
ppm     = 1e-6; % Some are in ppm (also a volume- not a mass-fraction)
g       = 1e-3;
s       = 1;
kelvin  = 273.15;
%% Load NASA maybe you need it at some point?
% Global (for the Nasa database in case you wish to use it).
global Runiv 
Runiv = 8.314;
TdataBase=fullfile('Nasa','NasaThermalDatabase');
load(TdataBase);

% Select species for the case at hand
iSp = myfind({Sp.Name},{'Diesel','O2','N2','CO2','H2O'});                      % Find indexes of these species
SpS=Sp(iSp);                                                                % Subselection of the database in the order
NSp = length(SpS);
Mi = [SpS.Mass];

%% Engine geom data (check if these are correct) (get this data from the rest)
Cyl.Bore                = 104*mm;
Cyl.Stroke              = 85*mm;
Cyl.CompressionRatio    = 21.5;
Cyl.ConRod              = 136.5*mm;
Cyl.TDCangle            = 0;
% -- Valve closing events can sometimes be seen in fast oscillations in the pressure signal (due
% to the impact when the Valve hits its seat).  (obsolute I think)
CaIVO = -355; % Intake valve open angle
CaIVC = -135; % Intake valve closing angle
CaEVO = 149; % Exhaust valve opening angle
CaEVC = -344; % Exhaust valve closing angle
CaSOI = -3.2; % Start of Injection Angle
% Maximal and minimal values of the Cylinder sizing
Vmin = CylinderVolume(0,Cyl);
Vmax = CylinderVolume(180,Cyl);

%% Base Variables
% Base conditions (please update these with correct values later on)
pamb = 1*bara; %ambient pressure of 1 bar 
Tamb = 21+kelvin; %ambient temperature converted to Kelvin
Fuelcomposition = 'CH4';
Q_LHV = 1; %NEEDS A UPDATE
mfuelinj = 0.001;                                                               %fuel injected per cycle in kg (NEEDS TO BE FOUND STILL)
p_exhaust = pamb+1*bara; % pressure of exhaust in bar

% air composition
Xair = [0 0.21 0.79 0 0];                                                   % Molar fraction air composition
MAir = Xair*Mi';                                                            
Yair = Xair.*Mi/MAir;                                                       % mass fraction air

% Initialise the volume vector
% volume is taken as a function of crank angle and then used throughout the
% rest of the file to compute the pressures
a_start = -360; %starting angle
a_fin = 360; %finishing angle
a_res = 0.2; %º CA
a_steps = (abs(a_start)+abs(a_fin))/a_res; %amount of steps of Ca
Ca = linspace(a_start,a_fin,a_steps); % vector with all crank angles that are sweeped through.

% volume vector
V_theory = CylinderVolume(Ca , Cyl);

% Create pressure vector
p_theory = zeros(size(Ca)); % Initialize pressure vector




%% Thermodynamical Analysis of the engine Cycle
% Reserved for comments

%% Intake cycle
v_6 = Vmax;
p_6 = pamb;

% Intake stroke: -360° to -180°
intake_indices = (Ca > -360) & (Ca <= -180);
p_theory(intake_indices) = pamb;

%% Compression Cycle
p1 = pamb;
V1 = Vmax;
V2 = Vmin;
T1 = Tamb;
R1 = Runiv/MAir;

% mass of air in cylinder
m_air = (p1*V1)/(R1*T1);

% Use NASA tables for cp and cv
for i = 1:NSp
    Cpi1(i) = CpNasa(T1, SpS(i));
    Cvi1(i) = CvNasa(T1, SpS(i));
end
cp1 = Yair * Cpi1';
cv1 = Yair * Cvi1';
gamma1 = cp1/cv1;

% Compression stroke: -180° to 0° (or -180° to TDC)
compression_indices = (Ca > -180) & (Ca <= 0);
V_compression = V_theory(compression_indices); % Get volumes during compression

% Calculate pressure for these specific volumes
p_compression = @(Vcom) p1 * (V1 ./ Vcom).^gamma1;
p_theory(compression_indices) = p_compression(V_compression);

% Calculate end pressure
p2 = p_compression(V2);

%end temperature
T2 = (p2*V2)/(m_air*R1);

%% Combustion
% AF ratio for fuel injected for now can be changed later to just plain
% fuel mass
mabsair = Yair .* m_air;
mtot = m_air + mfuelinj;
mabsair(1) = mfuelinj;

%calculate new mass and molar fractions
Yprecomb = mabsair ./ mtot;                                                 % Mass fraction of pre combustion mix
Mmixpre = mabsair ./ Mi;                                                    % Moles pre burn
Xprecomb = Mmixpre ./ sum(Mmixpre);                                         % molar fractions of pre combustion mix

% Post combustion molar and mass fractions
Cr = combustion_coefficients(Fuelcomposition);                              % combustion ratios 
Mmixpost = Mmixpre+(Mmixpre(1)*Cr);                                         % Molar Mix post burn
Xpostcomb = Mmixpost ./ sum(Mmixpost);
Mpostcomb = Xpostcomb*Mi'; 
Ypostcomb = Xpostcomb.*Mi/Mpostcomb;

% Temperature increase due to fuel burn
Q_in = mfuelinj*Q_LHV;
T3 = T2 + (Q_in/(cp1*mtot));

% Use NASA tables for new cp and cv
for i = 1:NSp
    Cpi2(i) = CpNasa(T3, SpS(i));
    Cvi2(i) = CvNasa(T3, SpS(i));
end
cp2 = Yair * Cpi2';
cv2 = Yair * Cvi2';
gamma2 = cp2/cv2;
R3 = Runiv/Mpostcomb;

% post combustion pressure
p3 = p2;

% post combustion Volume
V3 = (mtot*R3*T3)/p3;

% Add debugging
disp(['V3 = ', num2str(V3)])
disp(['Vmin = ', num2str(Vmin)])
disp(['Vmax = ', num2str(Vmax)])

% Find the angle where V reaches V3
% Only search in the expansion region (Ca > 0 and Ca < 180)
expansion_region = (Ca > 0) & (Ca <= 180);
combustion_end_index = find((V_theory >= V3) & expansion_region, 1, 'first');

% Check if combustion end point was found
if isempty(combustion_end_index)
    warning('V3 exceeds maximum volume in expansion stroke. Setting combustion end to 180°');
    combustion_end_angle = 180; % Use end of expansion stroke
else
    combustion_end_angle = Ca(combustion_end_index);
end

disp(['Combustion ends at: ', num2str(combustion_end_angle), '° CA'])

% Set pressure for combustion range
combustion_indices = (Ca > 0) & (Ca <= combustion_end_angle);
p_theory(combustion_indices) = p3; % Constant pressure

%% Expansion Cycle
V4=Vmax;

% Expansion stroke: endofcombustion to 180°
expansion_indices = (Ca > combustion_end_angle) & (Ca <= 180);
V_expansion = V_theory(expansion_indices); % Get volumes during compression

% Calculate pressure for these specific volumes
p_expansion = @(Vexp) p3 * (V3 ./ Vexp).^gamma2;
p_theory(expansion_indices) = p_expansion(V_expansion);

% Calculate end pressure
p4 = p_expansion(V4);

%end temperature
T4 = (p4*V4)/(mtot*R3);

% heat lost
Q_out = cv2*mtot*(T4-T1);

%% Exhaust cycle
v_5 = Vmin;
p_5 = p_exhaust;

% exhaust stroke: 180° to 360°
intake_indices = (Ca > 180) & (Ca <= 360);
p_theory(intake_indices) = p_exhaust;

% set point 1 as exhaust pressure to draw nice line
p_theory(1) = p_exhaust;

%% Data Analysis
% Below is given the code for the data analysis and viewing (NOT YET
% IMPLEMENTED)

% work done by system
W_theoretical = Q_in - Q_out - (p_exhaust-pamb)*(Vmax-Vmin)
% Calculation of work goes here

% reserved for bugfixing can be removed later


%% Load data (if txt file)
FullName        = fullfile('Data','ExampleDataSet.txt');
dataIn          = table2array(readtable(FullName));
[Nrows,Ncols]   = size(dataIn);                    % Determine size of array
NdatapointsperCycle = 720/0.2;                     % Nrows is a multitude of NdatapointsperCycle
Ncycles         = Nrows/NdatapointsperCycle;       % This must be an integer. If not checkwhat is going on
Ca              = reshape(dataIn(:,1),[],Ncycles); % Both p and Ca are now matrices of size (NCa,Ncycles)
p_measured               = reshape(dataIn(:,2),[],Ncycles)*bara; % type 'help reshape' in the command window if you want to know what it does (reshape is a Matlab buit-in command
%% Plotting 
f1=figure(1);
set(f1,'Position',[ 200 800 1200 400]);             % Just a size I like. Your choice
pp = plot(Ca,p_measured/bara,'LineWidth',1);                 % Plots the whole matrix
xlabel('Ca');ylabel('p [bar]');                     % Always add axis labels
xlim([-360 360]);ylim([0 50]);                      % Matter of taste
iselect = 10;                                    % Plot cycle 10 again in the same plot to emphasize it. Just to show how to access individual cycles.
line(Ca(:,iselect),p_measured(:,iselect)/bara,'LineWidth',2,'Color','r');
YLIM = ylim;
% Add some extras to the plot
line([CaIVC CaIVC],YLIM,'LineWidth',1,'Color','b'); % Plot a vertical line at IVC. Just for reference not a particular reason.
line([CaEVO CaEVO],YLIM,'LineWidth',1,'Color','r'); % Plot a vertical line at EVO. Just for reference not a particular reason.
set(gca,'XTick',-360:60:360,'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title('All cycles in one plot.')


%% pV-diagram
V_measured = CylinderVolume(Ca(:,iselect),Cyl);
f2 = figure(2);
set(f2,'Position',[ 200 400 600 800]);              % Just a size I like. Your choice
subplot(2,1,1)
plot(V_measured/dm^3,p_measured(:,iselect)/bara);
xlabel('V [dm^3]');ylabel('p [bar]');               % Always add axis labels
xlim([0 0.8]);ylim([0.5 50]);                      % Matter of taste
set(gca,'XTick',0:0.1:0.8,'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title({'pV-diagram','(with wrong Volume function btw)'})
subplot(2,1,2)
loglog(V_measured/dm^3,p_measured(:,iselect)/bara);
xlabel('V [dm^3]');ylabel('p [bar]');               % Always add axis labels
xlim([0.02 0.8]);ylim([0 50]);                      % Matter of taste
set(gca,'XTick',[0.02 0.05 0.1 0.2 0.5 0.8],...
    'YTick',[0.5 1 2 5 10 20 50],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title({'pV-diagram','(with wrong Volume function btw)'})


%% Plot p-V diagram with different colors for each stroke

figure;
hold on;

% Intake stroke: -360° to -180° (e.g., blue)
intake_indices = (Ca > -360) & (Ca <= -180);
intake_indices = intake_indices(1:length(V_theory));
plot(V_theory(intake_indices), p_theory(intake_indices), 'b-', 'LineWidth', 2, 'DisplayName', 'Intake');

% Compression stroke: -180° to 0° (e.g., red)
compression_indices = (Ca > -180) & (Ca <= 0);
compression_indices = compression_indices(1:length(V_theory));
plot(V_theory(compression_indices), p_theory(compression_indices), 'r-', 'LineWidth', 2, 'DisplayName', 'Compression');

% Expansion/Power stroke: 0° to 180° (e.g., green)
expansion_indices = (Ca > 0) & (Ca <= 180);
expansion_indices = expansion_indices(1:length(V_theory));
plot(V_theory(expansion_indices), p_theory(expansion_indices), 'g-', 'LineWidth', 2, 'DisplayName', 'Expansion');

% Exhaust stroke: 180° to 360° (e.g., magenta)
exhaust_indices = (Ca > 180) & (Ca <= 360) & 1;
exhaust_indices = exhaust_indices(1:length(V_theory));
plot(V_theory(exhaust_indices), p_theory(exhaust_indices), 'm-', 'LineWidth', 2, 'DisplayName', 'Exhaust');

xlabel('Volume (m³)');
ylabel('Pressure (Pa)');
title('p-V Diagram - Engine Cycle');
legend('Location', 'best');
grid on;
hold off;