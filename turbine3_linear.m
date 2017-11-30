clear all
%% Setup Variables
% Tanks
volume = 20/1000; % 20 L tanks (each gas)

% Turbine (Pelton Wheel) Parameters
mass_turb = 10;
radius_turb = 0.21;
J_turb = 0.5 * mass_turb * radius_turb^2; % approximate as a disc
b_turb = 0.2;

% Generator
J_gen = 0.005; % less than the turbine
k_emf = 1.4; % emf constant from DC-540 generator
L_gen = 0.315;
R_gen = 0.27; %http://forums.pelicanparts.com/porsche-911-technical-forum/199928-alternator-stator-coil-resistance.html
b_gen = 0.35;
R_load = 40.7;
C_gen = 0.005;

% Shaft parameters
G = 80E9; % stainless steel
d = 1/100; % 1cm
l = 0.1; % 10cm
% calculate_k = @(G, d, l) -(((pi * d^4)/32)*G)/l;
K = -785.3982;%calculate_k(G,d,l);

% Penstock Pipe 
head = 50;
pipe_diam = 0.0736;
pipe_length = 50/sin(45);
pipe_area = (pipe_diam/2)^2 * pi;
Rf = (128 * 8.9E-4 * pipe_length) / (pi * pipe_diam^4); 
If = 2 * 1000 * pipe_length / pipe_area;

% Water jet
jet_diam = 0.015;
jet_area = pi*(jet_diam/2)^2;
beta = jet_diam/pipe_diam;
jet_coefficient = 0.97;

mass_turb = 10;
radius_turb = 0.15;
J_turb = 0.5 * mass_turb * radius_turb^2; % approximate as a disc
b_turb = 0.2;

%% Setup Simulation
Q = 5.4/1000; % flow rate. From steady state flow rate determined using nonlinear model. 

a = [((-2*radius_turb^2*Q*1000)/J_turb - b_turb/J_turb) 0 -1/J_turb 0 0;
    0 -b_gen/J_gen 1/J_gen 0 0;
    -K K 0 0 0;
    0 k_emf/sqrt(2)/L_gen 0 -1/(L_gen) 1/(L_gen * R_load)
    0 0 0 1/C_gen -1/(R_gen * C_gen)];

b = [2000*radius_turb*Q*jet_coefficient*sqrt(2*9.8*head)/J_turb; 
    0 
    0 
    -1.4/L_gen/R_load
    0];


c = [1 0 0 0 0
     0 1 0 0 0
     0 0 1 0 0
     0 0 0 1 0
     0 0 0 0 1];
 
d = 0;

%% Simulate the System
sys = ss(a, b, c, d);
dt = 0.001;
t = 0:dt:15;
u = ones(length(t), 1);

simresult = lsim(sys, u, t);
plot(t,simresult)
legend('Turbine Shaft Rad/s', 'Generator Shaft Rad/s', 'Shaft Torque (N*m)','Output Current', 'Output Voltage');
xlabel('Time (s)')

%% Calculate volume of H2 produced at each time step
moles_H2 = simresult(:,5) .* simresult(:,4) * (dt /285000 * 8.314 * 298.15 / 101300 * 0.0899 / 0.00201588); % moles of H2
moles_O2 = simresult(:,5) .* simresult(:,4) * (dt /285000 * 8.314 * 298.15 /2 /101300 * 1.429 / 0.03199); % moles of H2

moles_H2 = cumsum(moles_H2); % integrate
moles_O2 = cumsum(moles_O2);

p_H2 = moles_H2(end) * 8.314 * 298.15 / volume; % pascals
p_O2 = moles_O2(end) * 8.314 * 298.15 / volume; % pascals

fprintf('Total moles of H2 produced in %d seconds: %d \r\n',t(end),moles_H2(end))
fprintf('Total moles of O2 produced in %d seconds: %d \r\n',t(end),moles_O2(end)/2)