x0 = [0]
fun = @nonlin
x = fsolve(fun,x0)
function xprime = nonlin(p)
% Turbine (Pelton Wheel) Parameters
mass_turb = 10;
radius_turb = 0.21;
J_turb = 0.5 * mass_turb * radius_turb^2; % approximate as a disc
b_turb = 0.2;

% Generator
J_gen = J_turb/2; % less than the turbine
k_emf = 1.4;
L_gen = 0.315;
R_gen = 0.27; 
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
y = 5.4/1000
x = 72.22
z = 57.1429

jet_velocity = sqrt(2*9.8*(head - (Rf * y / 1000 / 9.8)));

% xprime = (1000 * 9.8 * head - Rf*y - 0.5*1000*(1-beta^4)*(y/(jet_coefficient*jet_area))^2)/If; % penstock flow rate
% xprime = (2000 * y * radius_turb / J_turb)*(jet_coefficient * jet_velocity - x*radius_turb) - x*b_turb/J_turb - 20/J_turb; % turbine speed
% xprime = (20 - z*b_gen)/J_gen; % Generator speed

end