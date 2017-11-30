function [yH, yO, yB] = energyConsumptionODE(P_O0, P_H0)
% P_O0 = 2.8;       % [atm]
% P_H0 = 3;       % [atm]

%%  Model Parameter Constants
%%  Tank Gas Flow Model
%   Hydrogen
QH = 0;          % [J]

ICH = [QH
       P_H0];  
  
%   Oxygen

QO = 0;          % [J]

ICO = [QO
       P_O0]; 
  
%%  Fuel Cell Constants
T = 343;                % 70C
Rc = 0.0003;            % Ohms
A = 430;                % [cm^2]
t = 0.000178;           % [m]
n = 800;                % number of cells
psi = 13;               % saturated
i_max = 2;              % max current density [A/cm2]
i = 1.5;                % Current Pull [A/cm2]

%%  Heating Beans State Model Constants

%   Heater Resistance
R = 1;                  % [IDK]

%   Beans
beanVolume = 400/1000000;       % [m^3] - 400mL
beanDensity = 1000;             % [kg/m^3] - Water
beanSpecificHeat = 4184;        % [J/(kg*K)] - Water
C_TB = beanVolume * beanDensity * beanSpecificHeat;

%   Can
canD = 0.0682625;               % [m]
canH = 0.123825;                % [m]
canT = 0.00025;                 % [m]

%   Resistor
C_TR = 100;                     % [IDK]

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
Rc_RA = 0.01;

%   Beans to Amb    [Convection]
hc_BA = 10;                         % [Estimate]
A_BA = pi * ((canD/2) - canT)^2;
Rc_BA = 1/(hc_BA * A_BA);

% Ambient Temperature
T_A = 20;               % [C]

ICB = [T_A
       T_A];


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

%%  Define Run Time
dt = 0.01;
time = 0:dt:360;
uH = zeros(length(time), 1);
uO = zeros(length(time), 1);

%%  Define Hydrogen State Space Matrix

 
%%  Define Hydrogen State Space Matrix

[tH, yH] = ode45(@ODEfunctionH, time, ICH); 
[tO, yO] = ode45(@ODEfunctionO, time, ICO);
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

figure(1)

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
function dvariables = ODEfunctionH(t, y)
   RfH = 7000;      % [Ns/m^5]
    IfH = 50;        % [kg/m^4]
    VH = 0.04;       % [m^3]
    kH = 1.405;      % unitless

    dvariables = zeros(2,1);
    dvariables(1) = -(RfH*y(1)/IfH) + y(2)/IfH;
    dvariables(2) = -y(1)*kH*y(2)/VH;
end 

function dvariables = ODEfunctionO(t, y)
    RfO = 7000;      % [Ns/m^5]
    IfO = 50;        % [kg/m^4]
    VO = 0.04;       % [m^3]
    kO = 1.405;      % unitless

    dvariables = zeros(2,1);
    dvariables(1) = -(RfO*y(1)/IfO) + y(2)/IfO;
    dvariables(2) = -y(1)*kO*y(2)/VO;
end 
