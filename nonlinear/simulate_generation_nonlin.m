function [ xprime ] = simulate_generation_nonlin(t,x,a)
%SIMULATE_GENERATION_NONLIN Summary of this function goes here
%   Detailed explanation goes here
radius_turb = a.radius_turb;
J_turb = a.J_turb; % approximate as a disc
b_turb = a.b_turb;
J_gen = a.J_gen; % less than the turbine
k_emf = a.k_emf; % emf constant from DC-540 generator
L_gen = a.L_gen;
R_gen = a.R_gen; %http://forums.pelicanparts.com/porsche-911-technical-forum/199928-alternator-stator-coil-resistance.html
b_gen = a.b_gen;
R_load = a.R_load;
C_gen = a.C_gen;
K = a.K;
head = a.head;
Rf = a.Rf; 
If = a.If;
jet_area = a.jet_area;
beta = a.beta;
jet_coefficient = a.jet_coefficient;

% simulation
jet_velocity = sqrt(2*9.8*(head - (Rf * x(1) / 1000 / 9.8)));

xprime(1,1) = (1000 * 9.8 * head - Rf*x(1) - 0.5*1000*(1-beta^4)*(x(1)/(jet_coefficient*jet_area))^2)/If; % penstock flow rate
xprime(2,1) = (2000 * x(1) * radius_turb / J_turb)*(jet_coefficient * jet_velocity - x(2)*radius_turb) - x(2)*b_turb/J_turb - x(4)/J_turb; % turbine speed
xprime(3,1) = (x(4) - x(3)*b_gen)/J_gen; % Generator speed
xprime(4,1) = K * (-x(2) + x(3));  % shaft torque
xprime(5,1) = (k_emf*x(3)/sqrt(2) - (x(5) - (x(6)-1.4)/R_load)) / L_gen; % inductor current
xprime(6,1) = (x(5) - x(6)/R_gen)/C_gen; % output voltage
xprime(7,1) = x(5) * x(6) * 8.314 * 298.15 /285000 /101300 * 0.0899 /0.00201588; %vol * density / molar mass =  moles of H2
xprime(8,1) = x(5) * x(6) * 8.314 * 298.15 /285000 /101300 /2 * 1.429 / .03199; %vol * density / molar mass =  moles of O2
end




