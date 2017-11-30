function params = params()
% Tanks
params.volume = 20/1000;

% Turbine (Pelton Wheel) Parameters
params.mass_turb = 10;
params.radius_turb = 0.21;
params.J_turb = 0.5 * params.mass_turb * params.radius_turb^2; % approximate as a disc
params.b_turb = 0.2;

% Generator
params.J_gen = params.J_turb/2; % less than the turbine
params.k_emf = 1.4;
params.L_gen = 0.315;
params.R_gen = 0.27; 
params.b_gen = 0.35;
params.R_load = 40.7;
params.C_gen = 0.005;

% Shaft parameters
params.G = 80E9; % stainless steel
params.d = 1/100; % 1cm
params.l = 0.1; % 10cm
calculate_k = @(G, d, l) -(((pi * d^4)/32)*G)/l;
params.K = calculate_k(params.G,params.d,params.l);

% Penstock Pipe 
params.head = 50;
params.pipe_diam = 0.0736;
params.pipe_length = 50/sin(45);
params.pipe_area = (params.pipe_diam/2)^2 * pi;
params.Rf = (128 * 8.9E-4 * params.pipe_length) / (pi * params.pipe_diam^4); 
params.If = 2 * 1000 * params.pipe_length / params.pipe_area;
params.Q_static = 5.4/1000;

% Water jet
params.jet_diam = 0.015;
params.jet_area = pi*(params.jet_diam/2)^2;
params.beta = params.jet_diam/params.pipe_diam;
params.jet_coefficient = 0.97;
end