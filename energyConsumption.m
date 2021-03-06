function [yH, yO, yB] = energyConsumptionLinear(P_O0, P_H0, dt, maxtime)
% P_O0 = 2.8;       % [atm]
% P_H0 = 3;       % [atm]

%%  Define Run Time
time = 0:dt:maxtime;
uH = zeros(length(time), 1);
uO = zeros(length(time), 1);

%%  Model Parameter Constants
%%  Tank Gas Flow Model
%   Hydrogen
RfH = 5002.51;   % [atm s/m^3]
IfH = 107856;    % [kg/m^4]
VH = 0.04;       % [m^3]
kH = 1.405;      % unitless
PH = P_H0/4;     % [atm]
CfH = VH/(PH*kH);% [m^3/atm]
QH = 0;          % [J]

ICH = [QH
       P_H0];  
  
%   Oxygen
RfO = 11596.4;   % [atm s/m^3]
IfO = 1714800;   % [kg/m^4]
VO = 0.04;       % [m^3]
kO = 1.395;      % unitless
PO = P_O0/4;     % [atm]
CfO = VO/(PO*kO);% [m^3/atm]
QO = 0;          % [J]

ICO = [QO
       P_O0]; 
  
%%  Fuel Cell Constants
T = 343;                % 70C
Rc = 0.0003;            % Ohms
A = 430;                % [cm^2]
t = 0.000178;           % [m]
n = 100;                % number of cells
psi = 14;               % saturated
i_max = 2;              % max current density [A/cm2]
i = 1.5;                % Current Pull [A/cm2]
%   Can
canD = 0.0682625;               % [m]
canH = 0.123825;                % [m]
canT = 0.00025;                 % [m]
%%  Heating Beans State Model Constants

%   Heater Resistance
R = 1;                  % [IDK]

%   Beans
beanVolume = 400/1000000;       % [m^3] - 400mL
beanDensity = 1082.05;             % [kg/m^3] - Water
beanSpecificHeat = 4184;        % [J/(kg*K)] - Water
beanSurfaceArea = pi*(canD/2)^2;
C_TB = beanVolume * beanDensity * beanSpecificHeat;



%   Resistor to Beans [Conduction]
canArea = pi * canD * canH;
canHeaterAreaRatio = 0.5; % Amount of the can covered by the heater
canHeaterArea = canArea * canHeaterAreaRatio;

%   Resistor
resThickness = 0.01;            % [m]
resDensity = 2400;              % [kg/m3]
resHeatCapacity = 1085;          % [kJ/kgC]
resBaseArea = ((canD/2 + resThickness)^2 - (canD/2)^2)*pi;
resVolume = resBaseArea * canHeaterAreaRatio * canH;
resSurfaceArea = 2 * resBaseArea + canHeaterAreaRatio * canH;

C_TR = resVolume * resDensity * resHeatCapacity;        % [J/C]


k_RB = 50;                  % [W/mK]
A_RB = canHeaterArea;       % [m^2]
dx_RB = canT;               % [m]
Rk_RB = dx_RB/(k_RB * A_RB);

hc_A = 10.45;               % [Estimate]

%   Resistor to Amb [Convection]
Rc_RA = 1/(hc_A * resSurfaceArea);

%   Beans to Amb    [Convection]
Rc_BA = 1/(hc_A * beanSurfaceArea);

% Ambient Temperature
T_A = 20;               % [C]

ICB = [T_A
       T_A];

%%  Define Hydrogen State Space Matrix
AH = [-RfH/IfH     1/IfH
     -1/CfH      0   ];
 
BH = [0
      0];
  
CH = [1 0
      0 1];
 
DH = [0
      0];
  
sysH = ss(AH, BH, CH, DH);
 
%%  Define Hydrogen State Space Matrix
AO = [-RfO/IfO     1/IfO
      -1/CfO       0    ];
 
BO = [0
      0];
  
CO = [1 0
      0 1];
 
DO = [0
      0];
  
 sysO = ss(AO, BO, CO, DO);
 
%%  Define Bean Heater State Space Matrix
AB = [-(1/Rk_RB + 1/Rc_RA)/C_TR   1/(Rk_RB * C_TR)      
      1/(Rk_RB * C_TB)           -(1/Rk_RB + 1/Rc_BA)/C_TB];

BB = [1/(Rc_RA * C_TR)    1/(C_TR)
      1/(Rc_BA * C_TB)    0       ];
  
CB = [1 0
      0 1];
  
DB = [0 0
      0 0];

sysB = ss(AB, BB, CB, DB);

%%  Run Simulation

yH = lsim(sysH, uH, time, ICH);
yO = lsim(sysO, uO, time, ICO);

for w = 1:length(time)
    PO2 = yO(w,2);
    PH2 = yH(w,2);
    
    
    CH2 = P_H0/((5.08*10^6)*exp(-498/T));
    zeta1 = -0.948;
    zeta2 = (0.00286+0.0002*log(A)+(4.38*10^-5)*log(CH2));
    zeta3 = 7.6*(10^-5);
    zeta4 = -1.93*(10^-4);
    %%  Kinetic Voltage
    Co2 = PO2/((5.08*10^6)*exp(-498/T));
    Vk = -(zeta1 + zeta2*T + zeta3*T*log(Co2)+ zeta4*T*log(i*A));
    
    %%  Ohmic Voltage
    Rm = (t/A)*(181.6*(1+0.03*i+0.062*((T/303)^2)*(i^2.5)))/...
         (psi-0.634-3*i*exp(4.18*((T-303)/T))); 

    Vo = i*A*(Rm+Rc);

    %%  Transportation Voltage
    Vt = -0.016*log(1-(i/i_max));

    %%  Nernst Voltage
    Vn = 1.229-(0.85*10^-3)*(T-298.15)+(4.3085*10^-5)*T*...
         (log(PH2)+0.5*log(PO2));

    %%  Stack Voltage
    z(w,1) = Vn;
    z(w,2) = Vk;
    z(w,3) = Vo;
    z(w,4) = Vt;
    if(n*(Vn - Vk - Vo - Vt) > 0)
        uB(w,2) = dt*((n*(Vn - Vk - Vo - Vt))^2)/R;
    else
        uB(w,2) = 0;
    end
    uB(w,1) = T_A;

end

%%  Bean Simulator
yB = lsim(sysB, uB, time, ICB);

figure('NumberTitle', 'off', 'Name', 'Linear Consumption Results')

subplot(2,2,1)
plot(time(1,:),yH(:,2))
title('Hydrogen Tank Pressure')
xlabel('Time [s]');
ylabel('Pressure [atm]');
grid on

subplot(2,2,2)
plot(time(1,:),yO(:,2))
title('Oxygen Tank Pressure')
xlabel('Time [s]');
ylabel('Pressure [atm]');
grid on

subplot(2,2,3)
plot(time(1,:),yB(:,2))
title('Bean Temperature')
xlabel('Time [s]');
ylabel('Bean Temperature [C]');
grid on

subplot(2,2,4)
plot(time(1,:),yB(:,1))
title('Resistor Temperature')
xlabel('Time [s]');
ylabel('Temperature [C]');
grid on
end
