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

%%  Define Simulation Parameters
simulationTime = 24*3600; %s
deltaTime = 0.01;
iterations = simulationTime/deltaTime;

%%  Define material and geometric properties 
%   of the system for the Thermal Capacitors

%   Beans
    beanVolume = 400/1000;          % [m^3] - 400mL
    beanDensity = 1000;             % [kg/m^3] - Water
    beanSpecificHeat = 4184;        % [J/(kg*K)] - Water
    C_TB = beanVolume * beanDensity * beanSpecificHeat;

%   Can
    canD = 0.0682625;               % [m]
    canH = 0.123825;                % [m]
    canT = 0.00025;                 % [m]
    %canVolume = canD * canH * canT; % [m^3]

%   Resistor
    C_TR = 100; %IDK
    
%%  Thermal Resistors

%   Resistor to Beans [Conduction]
    canArea = pi * canD * canH;
    canHeaterAreaRatio = 0.5; % Amount of the can covered by the heater
    canHeaterArea = canArea * canHeaterAreaRatio;

    k_RB = 50;                  % [W/mK]
    A_RB = canHeaterArea;       % [m^2]
    dx_RB = canT;               % [m]
    Rk_RB = dx_RB/(k_RB * A_RB);
    
%   Resistor to Amb [Convection]
    canAirArea = (pi * canD^2)/4; % + canArea * (1 - resistorCanAreaRatio) This is insulated so no effect

%     hc_RA = 1;                         % [Estimate]
%     A_RA = canAirArea;
%     Rc_RA = 1/(hc_RA * A_RA);
    Rc_RA = 0.01;
    
%   Beans to Amb    [Convection]
    hc_BA = 10;                         % [Estimate]
    A_BA = pi * ((canD/2) - canT)^2;
    Rc_BA = 1/(hc_BA * A_BA);
    
%%  Define System Inputs

% Ambient Temperature
T_A = 20;

%Q_Dissipate = 1;
V = 10;
R = 1;
Qin = (V^2)/R;



A = [-(1/Rk_RB + 1/Rc_RA)/C_TR   1/(Rk_RB * C_TR)      
     1/(Rk_RB * C_TB)           -(1/Rk_RB + 1/Rc_BA)/C_TB];

B = [1/(Rc_RA * C_TR)    1/(C_TR)
     1/(Rc_BA * C_TB)    0       ];

C = [1  0
     0  1];
 
D = [0  0
     0  0];
 
sys = ss(A, B, C, D);

T_R = T_A;
T_B = T_A;

T = zeros(2,iterations);
t = 0:deltaTime:iterations*deltaTime;
T(:,1) = [T_R, T_B]';
input = [T_A, Qin]';

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