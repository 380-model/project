%% ---------------- Initial model, just inertia of turbine
% y0 = [0];
% [t,y]  =ode45(@turbine, [0 10], y0);
%
% function xprime = turbine(t,y)
%     Q = 2; % this is the input
%     mass_turb = 20;
%     radius_turb = 15;
%     a_jet = 0.01^2 * pi;
%     J = 0.5 * mass_turb * radius_turb^2; % approximate as a disc
%
%
%     xprime = (1/J) * (-radius_turb * 2 * 1000 * Q * y(1) + ((2 * 1000 * Q^2)/a_jet));
% end

%% ---------- Added viscous damping --------------
%
% y0 = [0]; % initial angular velocity = 0
% [t,y] = ode45(@turbine, [0:0.01:1], y0);
% y = y / (2*pi) * 60; % convert to RPM
% figure
% plot(t,y)
% title('Initial plot of Pelton Wheel Turbine')
% legend('Angular Velocity (RPM)')
% xlabel('Simulation time (s)')
% ylabel('Angular Velocity (RPM)')
%
% function xprime = turbine(t,y)
%     Q = 1; % this is the input, assume constant flow for now
%     mass_turb = 5;
%     radius_turb = 0.1;
%     a_jet = 0.005^2 * pi;
%     J = 0.5 * mass_turb * radius_turb^2; % approximate as a disc
%     b_turb = 5000;
%
%     xprime = (radius_turb/J) * (-2 * 1000 * Q * y(1) + (( 2 * 1000 * Q^2)/a_jet)) - y(1)*b_turb/J;
% end

%% ---------- Added spring and other damper --------------

% y0 = [0,0,0]; % initial angular velocity = 0
% [t,y] = ode45(@turbine3, [0 10], y0);
% power = y(:,2) .* y(:,3);
% y = y / (2*pi) * 60; % convert to RPM
% figure
% plot(t,y)
% hold on
% plot(t, power)
% title('Initial plot of Pelton Wheel Turbine')
% legend('Turbine RPM','Generator RPM','Shaft Torque (N*m)', 'Shaft Power (W)')
% xlabel('Simulation time (s)')
% ylabel('Angular Velocity (RPM)')
% 
% function xprime = turbine3(t,y)
% Q = 5/1000; % (5 lps) this is the input, assume constant flow for now
% head = 50;
% mass_turb = 10;
% radius_turb = 0.15;
% a_jet = 0.002^2 * pi;
% J_turb = 0.5 * mass_turb * radius_turb^2; % approximate as a disc
% b_turb = 0.2;
% J_gen = 0.005; % less than the turbine
% jet_factor = 0.97; % to account for losses and such
% 
% % Shaft parameters
% G = 80E9; % stainless steel
% d = 1/100; % 1cm
% l = 0.1; % 10cm
% % calculate_k = @(G, d, l) -(((pi * d^4)/32)*G)/l;
% K = -785.3982;%calculate_k(G,d,l);
% b_gen = 0.35;
% 
% xprime(1,1) = (2000 * Q * radius_turb / J_turb)*(jet_factor * sqrt(2*9.8*head) - y(1)*radius_turb) - y(1)*b_turb/J_turb - y(3)/J_turb;
% xprime(2,1) = (y(3) - y(2)*b_gen)/J_gen;
% xprime(3,1) = K * (-y(1) + y(2)); % spring
% end

%% ---------- Added generator --------------
% 
% y0 = [0,0,0,0]; % initial angular velocity = 0
% [t,y] = ode45(@turbine3, [0 10], y0);
% power = y(:,2) .* y(:,3);
% y = y / (2*pi) * 60; % convert to RPM
% figure
% plot(t,y)
% hold on
% plot(t, power)
% title('Initial plot of Pelton Wheel Turbine')
% legend('Turbine RPM','Generator RPM','Shaft Torque (N*m)', 'Output current (A)','Shaft Power (W)')
% xlabel('Simulation time (s)')
% ylabel('Value')
% 
% function xprime = turbine3(t,y)
% Q = 5/1000; % (5 lps) this is the input, assume constant flow for now
% head = 50;
% mass_turb = 10;
% radius_turb = 0.15;
% a_jet = 0.002^2 * pi;
% J_turb = 0.5 * mass_turb * radius_turb^2; % approximate as a disc
% b_turb = 0.2;
% jet_factor = 0.97; % to account for losses and such
% 
% % Generator
% J_gen = 0.005; % less than the turbine
% k_emf = 0.699; % emf constant from DC-540 generator
% k_gen = 0.999;% motor (generator constant for DC-540 generator)
% L_gen = 0.315;
% R_gen = 0.27; %http://forums.pelicanparts.com/porsche-911-technical-forum/199928-alternator-stator-coil-resistance.html
% b_gen = 0.35;
% R_load = 30;
% 
% % Shaft parameters
% G = 80E9; % stainless steel
% d = 1/100; % 1cm
% l = 0.1; % 10cm
% % calculate_k = @(G, d, l) -(((pi * d^4)/32)*G)/l;
% K = -785.3982;%calculate_k(G,d,l);
% 
% % Pipe 
% pipe_diam = 0.15;
% pipe_length = 10;
% pipe_area = (pipe_diam/2)^2 * pi;
% Rf = 128 * 8.9E-4 * pipe_length / (pi * pipe_diam^4);
% If = 2 * 1000 * pipe_length / pipe_area;
% 
% 
% xprime(1,1) = (2000 * Q * radius_turb / J_turb)*(jet_factor * sqrt(2*9.8*head) - y(1)*radius_turb) - y(1)*b_turb/J_turb - y(3)/J_turb;
% xprime(2,1) = (y(3) - y(2)*b_gen)/J_gen;
% xprime(3,1) = K * (-y(1) + y(2)); 
% xprime(4,1) = (1/L_gen) * (k_emf * y(2) - R_gen*y(4) - R_load*y(4));
% end

%% ---------- Added Penstock Pipe --------------

y0 = [0,0,0,0,0]; % initial angular velocity = 0
[t,y] = ode45(@turbine3, [0 100], y0);
power = y(:,2) .* y(:,3);
y = y / (2*pi) * 60; % convert to RPM
figure
plot(t,y)
hold on
plot(t, power)
title('Initial plot of Pelton Wheel Turbine')
legend('Turbine RPM','Generator RPM','Shaft Torque (N*m)', 'Output current (A)','Shaft Power (W)','Penstock Flow Rate')
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
jet_factor = 0.97; % to account for losses and such

% Generator
J_gen = 0.005; % less than the turbine
k_emf = 0.699; % emf constant from DC-540 generator
k_gen = 0.999;% motor (generator constant for DC-540 generator)
L_gen = 0.315;
R_gen = 0.27; %http://forums.pelicanparts.com/porsche-911-technical-forum/199928-alternator-stator-coil-resistance.html
b_gen = 0.35;
R_load = 30;

% Shaft parameters
G = 80E9; % stainless steel
d = 1/100; % 1cm
l = 0.1; % 10cm
% calculate_k = @(G, d, l) -(((pi * d^4)/32)*G)/l;
K = -785.3982;%calculate_k(G,d,l);

% Penstock Pipe 
pipe_diam = 0.0736;
pipe_length = head*sin(45);
pipe_area = (pipe_diam/2)^2 * pi;
Rf = (128 * 8.9E-4 * pipe_length) / (pi * pipe_diam^4);
If = 2 * 1000 * pipe_length / pipe_area;


xprime(1,1) = (2000 * Q * radius_turb / J_turb)*(jet_factor * sqrt(2*9.8*(head - (Rf * Q)/(1000 *9.8))) - y(1)*radius_turb) - y(1)*b_turb/J_turb - y(3)/J_turb; % turbine speed
xprime(2,1) = (y(3) - y(2)*b_gen)/J_gen; % Generator speed
xprime(3,1) = K * (-y(1) + y(2));  % shaft torque
xprime(4,1) = (1/L_gen) * (k_emf * y(2) - R_gen*y(4) - R_load*y(4)); % Output current
xprime(5,1) = (1000 * 9.8 * head - Rf*y(5)) / If; % flow rate through penstock ... values seem weird but the power output is right ??
end
