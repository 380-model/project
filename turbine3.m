%% ---------- Added Rectifier Circuit --------------
close all
% penstock flow rate, turbine speed, generator speed, shaft torque,
% generator current, generator voltage, moles of H2, moles of O2

y0 = [0,0,0,0,0,0,0,0]; % initial angular velocity = 0
[t,y] = ode45(@turbine4, [0 150], y0);

%% Calculate Pressure of H2 and O2 for fuel cell
volume = 20/1000; % 20 L tanks (each gas)
p_H2 = y(:,7);
p_O2 = y(:,8);
p_H2 = p_H2(end) * 8.314 * 298.15 / volume; % pascals
p_O2 = p_O2(end) * 8.314 * 298.15 / volume; % pascals

%% Do various conversions
shaft_power = y(:,3) .* y(:,4);
electrical_power = y(:,5) .* y(:,6);

y(:,2) = y(:,2) / (2*pi) * 60;  % convert to RPM
y(:,3) = y(:,3) / (2*pi) * 60;  % convert to RPM

%% Plotting 
figure
plot(t,y(:,2),t,y(:,3),t,y(:,4),t,y(:,5),t,y(:,6))
hold on
plot(t, shaft_power)
plot(t, electrical_power)
title('Pelton Wheel Turbine')
legend('Turbine RPM','Generator RPM','Shaft Torque (N*m)', 'Inductor current (A)', 'Rectifier Voltage (V)','Shaft Power (W)','Electrical Power (W)')
xlabel('Simulation time (s)')
ylabel('Value')

% Plot of Flow rate
figure
subplot(2,2,1)
plot(t,y(:,1))
title('Penstock Flow Rate')
xlabel('Simulation time (s)')
ylabel('Penstock Flow Rate (m^3/s)')

% Plot of gas production
subplot(2,2,2)
plot(t,y(:,7),t,y(:,8))
title('Gas Production Over Time (moles)')
xlabel('Simulation time (s)')
ylabel('Moles of Gas')
legend('Quantity of H2 (moles)','Quantity of O2 (moles)')

% Plot of generator variables
subplot(2,2,3)
plot(t,y(:,5),t,y(:,6))
title('Generator')
xlabel('Simulation time (s)')
ylabel('Value')
legend('Generator Current (A)','Generator Voltage (V)')

% Prints
fprintf('Average shaft power: %f W\r\n',mean(shaft_power))
fprintf('Average electrical power: %f W\r\n',mean(electrical_power))
fprintf('Generator efficiency: %f percent\r\n',mean(electrical_power)/mean(shaft_power)*100)

%% Nonlinear Input Model
function xprime = turbine4(t,y)
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

jet_velocity = sqrt(2*9.8*(head - (Rf * y(1) / 1000 / 9.8)));

xprime(1,1) = (1000 * 9.8 * head - Rf*y(1) - 0.5*1000*(1-beta^4)*(y(1)/(jet_coefficient*jet_area))^2)/If; % penstock flow rate
xprime(2,1) = (2000 * y(1) * radius_turb / J_turb)*(jet_coefficient * jet_velocity - y(2)*radius_turb) - y(2)*b_turb/J_turb - y(4)/J_turb; % turbine speed
xprime(3,1) = (y(4) - y(3)*b_gen)/J_gen; % Generator speed
xprime(4,1) = K * (-y(2) + y(3));  % shaft torque
xprime(5,1) = (k_emf*y(3)/sqrt(2) - (y(5) - (y(6)-1.4)/R_load)) / L_gen; % inductor current
xprime(6,1) = (y(5) - y(6)/R_gen)/C_gen; % output voltage
xprime(7,1) = y(5) * y(6) * 8.314 * 298.15 /285000 /101300 * 0.0899 /0.00201588; %vol * density / molar mass =  moles of H2
xprime(8,1) = y(5) * y(6) * 8.314 * 298.15 /285000 /101300 /2 * 1.429 / .03199; %vol * density / molar mass =  moles of O2
end
