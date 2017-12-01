function [  ] = plot_generation_nonlin( t,y )
%PLOT_GENERATION_NONLIN Summary of this function goes here
%   Detailed explanation goes here
%% Unit conversion
shaft_power = y(:,3) .* y(:,4);
electrical_power = y(:,5) .* y(:,6);

y(:,2) = y(:,2) / (2*pi) * 60;  % convert to RPM
y(:,3) = y(:,3) / (2*pi) * 60;  % convert to RPM

%% Plotting 
% figure('NumberTitle', 'off', 'Name', 'Nonlinear Generation Results')
% plot(t,y(:,2),t,y(:,3),t,y(:,4),t,y(:,5),t,y(:,6))
% hold on
% plot(t, shaft_power)
% plot(t, electrical_power)
% title('Pelton Wheel Turbine')
% legend('Turbine RPM','Generator RPM','Shaft Torque (N*m)', 'Inductor current (A)', 'Rectifier Voltage (V)','Shaft Power (W)','Electrical Power (W)')
% xlabel('Simulation time (s)')
% ylabel('Value')

% Plot of Flow rate
figure('NumberTitle', 'off', 'Name', 'Nonlinear Generation Results')

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

% Plot of mechanical things
subplot(2,2,4)
plot(t,y(:,2),t,y(:,3),t,y(:,4))
title('Mechanical Outputs')
xlabel('Simulation time (s)')
ylabel('Value')
legend('Turbine Speed (RPM)','Generator Speed (RPM)','Shaft Torque (N*m)')

%% Prints
disp('------------ NONLINEAR MODEL -----------------')
fprintf('Average shaft power: %f W\r\n',mean(shaft_power))
fprintf('Average electrical power: %f W\r\n',mean(electrical_power))
fprintf('Generator efficiency: %f percent\r\n',mean(electrical_power)/mean(shaft_power)*100)
moles_H2 = y(:,7);
moles_O2 = y(:,8);
fprintf('Total moles of H2 produced in %d seconds: %d \r\n',t(end),moles_H2(end))
fprintf('Total moles of O2 produced in %d seconds: %d \r\n',t(end),moles_O2(end))

end

