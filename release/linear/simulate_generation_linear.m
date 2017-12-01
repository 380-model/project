function [ t, y, moles_H2, moles_O2, p_H2, p_O2 ] = simulate_generation_linear(simulation_time, a )
%SIMULATE_GENERATION_LINEAR Summary of this function goes here
%   Detailed explanation goes here

%% Setup Variables
volume = a.volume; % 20 L tanks (each gas)
radius_turb = a.radius_turb;
J_turb = a.J_turb; 
b_turb = a.b_turb;
J_gen = a.J_gen;
k_emf = a.k_emf; 
L_gen = a.L_gen;
R_gen = a.R_gen; 
b_gen = a.b_gen;
R_load = a.R_load;
C_gen = a.C_gen;
K = a.K;
head = a.head;
jet_coefficient = a.jet_coefficient;
Q = a.Q_static; % flow rate. From steady state flow rate determined using nonlinear model. 


%% Setup Simulation
A = [((-2*radius_turb^2*Q*1000)/J_turb - b_turb/J_turb) 0 -1/J_turb 0 0;
    0 -b_gen/J_gen 1/J_gen 0 0;
    -K K 0 0 0;
    0 k_emf/sqrt(2)/L_gen 0 -1/(L_gen) 1/(L_gen * R_load)
    0 0 0 1/C_gen -1/(R_gen * C_gen)];

B = [2000*radius_turb*Q*jet_coefficient*sqrt(2*9.8*head)/J_turb; 
    0 
    0 
    -1.4/L_gen/R_load
    0];


C = [1 0 0 0 0
     0 1 0 0 0
     0 0 1 0 0
     0 0 0 1 0
     0 0 0 0 1];
 
D = 0;

%% Simulate the System
sys = ss(A, B, C, D);
dt = 0.001;
t = 0:dt:simulation_time;
u = ones(length(t), 1);

y = lsim(sys, u, t);

figure('NumberTitle', 'off', 'Name', 'Linear Generation Results')
plot(t,y)
title('Data - Linear Generation Model')
legend('Turbine Shaft Rad/s', 'Generator Shaft Rad/s', 'Shaft Torque (N*m)','Output Current', 'Output Voltage');
xlabel('Time (s)')

%% Calculate volume of H2 produced at each time step
moles_H2 = y(:,5) .* y(:,4) * (dt /285000 * 8.314 * 298.15 / 101300 * 0.0899 / 0.00201588); % moles of H2
moles_O2 = y(:,5) .* y(:,4) * (dt /285000 * 8.314 * 298.15 /2 /101300 * 1.429 / 0.03199); % moles of H2

moles_H2 = cumsum(moles_H2); % integrate
moles_O2 = cumsum(moles_O2);

p_H2 = moles_H2(end) * 8.314 * 298.15 / volume /101300; % bar
p_O2 = moles_O2(end) * 8.314 * 298.15 / volume /101300; % bar


%% Plot the result
% Plot of gas production
subplot(2,2,1)
plot(t,moles_H2,t,moles_O2)
title('Gas Production Over Time (moles)')
xlabel('Simulation time (s)')
ylabel('Moles of Gas')
legend('Quantity of H2 (moles)','Quantity of O2 (moles)')

% Plot of generator variables
subplot(2,2,3)
plot(t,y(:,4),t,y(:,5))
title('Generator')
xlabel('Simulation time (s)')
ylabel('Value')
legend('Generator Current (A)','Generator Voltage (V)')

% Plot of mechanical things
subplot(2,2,4)
plot(t,y(:,1),t,y(:,2),t,y(:,3))
title('Mechanical Outputs')
xlabel('Simulation time (s)')
ylabel('Value')
legend('Turbine Speed (RPM)','Generator Speed (RPM)','Shaft Torque (N*m)')

%% Prints
shaft_power = y(:,2) .* y(:,3);
electrical_power = y(:,4) .* y(:,5);
disp('------------ LINEAR MODEL -----------------')
fprintf('Average shaft power: %f W\r\n',mean(shaft_power))
fprintf('Average electrical power: %f W\r\n',mean(electrical_power))
fprintf('Generator efficiency: %f percent\r\n',mean(electrical_power)/mean(shaft_power)*100)
fprintf('Total moles of H2 produced in %d seconds: %d \r\n',t(end),moles_H2(end))
fprintf('Total moles of O2 produced in %d seconds: %d \r\n',t(end),moles_O2(end))


end

