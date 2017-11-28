y0 = [0,0,0]; % initial angular velocity = 0
[t,y] = ode45(@turbine3, [0 10], y0);
power = y(:,2) .* y(:,3);
y = y / (2*pi) * 60; % convert to RPM
figure
plot(t,y)
hold on
plot(t, power)
title('Pelton Wheel Turbine')
legend('Turbine RPM','Generator RPM','Shaft Torque (N*m)', 'Shaft Power (W)')
xlabel('Simulation time (s)')
ylabel('Value')

function xprime = turbine3(t,y)
Q = 5/1000; % (5 lps) this is the input, assume constant flow for now
head = 50;
mass_turb = 10;
radius_turb = 0.15;
a_jet = 0.002^2 * pi;
J_turb = 0.5 * mass_turb * radius_turb^2; % approximate as a disc
b_turb = 0.2;
J_gen = 0.005; % less than the turbine
jet_factor = 0.97; % to account for losses and such

% Shaft parameters
G = 80E9; % stainless steel
d = 1/100; % 1cm
l = 0.1; % 10cm
% calculate_k = @(G, d, l) -(((pi * d^4)/32)*G)/l;
K = -785.3982;%calculate_k(G,d,l);
b_gen = 0.35;

xprime(1,1) = (2000 * Q * radius_turb / J_turb)*(jet_factor * sqrt(2*9.8*head) - y(1)*radius_turb) - y(1)*b_turb/J_turb - y(3)/J_turb;
xprime(2,1) = (y(3) - y(2)*b_gen)/J_gen;
xprime(3,1) = K * (-y(1) + y(2)); % spring
end