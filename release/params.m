% These parameters are used across models. 

function a = params()
% Gas Tanks
a.volume = 20/1000; % m^3

% Turbine (Pelton Wheel) Parameters
a.mass_turb = 10; % kg
a.radius_turb = 0.21; % m
a.J_turb = 0.5 * a.mass_turb * a.radius_turb^2; % approximate as a disc
a.b_turb = 0.2; % coefficient of friction

% Generator
a.J_gen = a.J_turb/2; % less than the turbine
a.k_emf = 1.4; % 
a.L_gen = 0.315; % H
a.R_gen = 0.27;  % ?
a.b_gen = 0.35; % coefficient of friction
a.R_load = 40.7; % ? - from electrolyzer model
a.C_gen = 0.005; % F

% Shaft parameters
a.G = 80E9; % stainless steel, Pa
a.d = 1/100; % diameter, metres
a.l = 0.10; % length, metres
calculate_k = @(G, d, l) -(((pi * d^4)/32)*G)/l;
a.K = calculate_k(a.G,a.d,a.l); 

% Penstock Pipe 
a.head = 50; % metres
a.pipe_diam = 0.0736; % metres
a.pipe_length = a.head/sin(45); % metres
a.pipe_area = (a.pipe_diam/2)^2 * pi; % metres
a.Rf = (128 * 8.9E-4 * a.pipe_length) / (pi * a.pipe_diam^4); % fluid resistance
a.If = 2 * 1000 * a.pipe_length / a.pipe_area; % fluid inertence
a.Q_static = 5.4/1000; % Static flow rate used as input to linear model

% Water jet
a.jet_diam = 0.015; % metres
a.jet_area = pi*(a.jet_diam/2)^2;
a.beta = a.jet_diam/a.pipe_diam; 
a.jet_coefficient = 0.97; 

%% Consumption
a.T_A = 20;               % ambient temperature [C]

% Can
a.canD = 0.0682625;               % [m]
a.canH = 0.123825;                % [m]
a.canT = 0.00025;                 % [m]
a.canArea = pi * a.canD * a.canH;
a.canHeaterAreaRatio = 0.5;
a.canHeaterArea = a.canArea * a.canHeaterAreaRatio;

% Resistive heater
a.R = 1;                 % ?
a.resThickness = 0.01;            % [m]
a.resDensity = 2400;              % [kg/m3]
a.resHeatCapacity = 1085;         % [kJ/kgC]
a.resBaseArea = ((a.canD/2 + a.resThickness)^2 - (a.canD/2)^2)*pi;
a.resVolume = a.resBaseArea * a.canHeaterAreaRatio * a.canH;
a.resSurfaceArea = 2 * a.resBaseArea + a.canHeaterAreaRatio * a.canH;

a.C_TR = a.resVolume * a.resDensity * a.resHeatCapacity;        % [J/C]

a.k_RB = 50;                  % [W/mK]

% Beans
a.beanVolume = 400/1000000;       % [m^3] - 400mL
a.beanDensity = 1082.05;             % [kg/m^3] - Water
a.beanSpecificHeat = 4184;        % [J/(kg*K)] - Water
a.beanSurfaceArea = pi*(a.canD/2)^2;
a.C_TB = a.beanVolume * a.beanDensity * a.beanSpecificHeat;
a.hc_A = 10.45;               % [Estimate]





end