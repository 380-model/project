clc;
clear all;
close all;

% === Electrical === Order: 0
% Voltage Source from FUELCELL
% Resistor dissipates conducts into CAN

% === Thermal === Order: 2
% Capacitor: HEATER -> Q = C_TH * dT/dt
% Convection from HEATER to AMBIENT - > deltaT = R_HA * Q
% Conduction from HEATER to BEANS -> deltaT = Rk * Q
% Capacitor: BEANS -> Q = CT * dT/dt
% Convection from BEANS to AMBIENT -> deltaT = Rc * Q

simulationTime = 24*3600; %s
deltaTime = 0.01;
iterations = simulationTime/deltaTime;

beanVolume = 400/1000;      % [m^3] - 400mL
beanDensity = 1000;         % [kg/m^3] - Water
beanSpecificHeat = 4184;    % [J/(kg*K)] - Water

canDiameter = 0.0682625;
canHeight = 0.123825;
canThickness = 0.00025;

canArea = pi * canDiameter * canHeight;

canHeaterAreaRatio = 0.5; % Amount of the can covered by the heater
canHeaterArea = canArea * canHeaterAreaRatio;

canAirArea = (pi * canDiameter^2)/4; % + canArea * (1 - resistorCanAreaRatio) This is insulated so no effect

% Ambient Temperature
T_A = 20;

%Q_Dissipate = 1;
V2 = 10000;
R_H = 1;

% === Specific Capacities === 
% volume [m^3] * density [kg/m^3] * specific heat [J/(kg*K)]

% Heater: Guessing?
C_TH = 100

% Beans
C_TB = beanVolume * beanDensity * beanSpecificHeat


% === Conductive Resistance === 
% length (m) / (conductivity * area)

% Heater-Beans: 0.00025m / (50W/mK * PI * 0.0682625m * 0.123825/2m)
R_HB = canThickness / (50 * canHeaterArea)

% === Convective Resistance ===
% 1 / (convection coefficient * area)


R_HA = 0.01;

% h_c estimated to be 10
R_BA = 1 / (10 * canAirArea);

A = [[-(1/R_HB + 1/R_HA)/C_TH,   1/(R_HB * C_TH)      ]
     [  1/(R_HB * C_TB),       -(1/R_HB + 1/R_BA)/C_TB]]

B = [[1/(R_HA * C_TH), 1/(R_H * C_TH)]
     [1/(R_BA * C_TB), 0             ]] 

T_H = T_A;
T_B = T_A;

T = zeros(2,iterations);
t = 0:deltaTime:iterations*deltaTime;
T(:,1) = [T_H, T_B]';
input = [T_A, V2]';

for i=1:iterations
    derivative = A * T(:,i) + B * input;
    
    T(:,i+1) = T(:,i) + derivative * deltaTime;
end

figure;
plot(t,T(1,:));
title('Heater Temperature')
figure;
plot(t,T(2,:));
title('Bean Temperature')