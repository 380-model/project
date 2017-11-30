function a = params()
% Tanks
a.volume = 20/1000;

% Turbine (Pelton Wheel) Parameters
a.mass_turb = 10;
a.radius_turb = 0.21;
a.J_turb = 0.5 * a.mass_turb * a.radius_turb^2; % approximate as a disc
a.b_turb = 0.2;

% Generator
a.J_gen = a.J_turb/2; % less than the turbine
a.k_emf = 1.4;
a.L_gen = 0.315;
a.R_gen = 0.27; 
a.b_gen = 0.35;
a.R_load = 40.7;
a.C_gen = 0.005;

% Shaft parameters
a.G = 80E9; % stainless steel
a.d = 1/100; % 1cm
a.l = 0.1; % 10cm
calculate_k = @(G, d, l) -(((pi * d^4)/32)*G)/l;
a.K = calculate_k(a.G,a.d,a.l);

% Penstock Pipe 
a.head = 50;
a.pipe_diam = 0.0736;
a.pipe_length = 50/sin(45);
a.pipe_area = (a.pipe_diam/2)^2 * pi;
a.Rf = (128 * 8.9E-4 * a.pipe_length) / (pi * a.pipe_diam^4); 
a.If = 2 * 1000 * a.pipe_length / a.pipe_area;
a.Q_static = 5.4/1000;

% Water jet
a.jet_diam = 0.015;
a.jet_area = pi*(a.jet_diam/2)^2;
a.beta = a.jet_diam/a.pipe_diam;
a.jet_coefficient = 0.97;
end