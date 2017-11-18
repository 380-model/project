% === Electrical === Order: 0
% Voltage Source from FUELCELL
% Resistor dissipates conducts into CAN

% === Thermal === Order: 3
% Capacitor: RESISTOR -> Q = CT * dT/dt
% Convection from RESISTOR to AMBIENT - > deltaT = Rc * Q
% Conduction from RESISTOR to CAN -> deltaT = Rk * Q
% Capacitor: CAN -> Q = CT * dT/dt
% Convection from CAN to AMBIENT -> deltaT = Rc * Q
% Conduction from CAN to BEANS -> deltaT = Rk * Q
% Capacitor: BEANS -> Q = CT * dT/dt
% Convection from BEANS to AMBIENT -> deltaT = Rc * Q

T_A = 20;

%Q_Dissipate = 1;
R_H = 1;

C_TH = 1;
C_TC = 1;
C_TB = 1;


R_HC = 1;
R_HA = 1;
R_CB = 1;
R_CA = 1;
R_BA = 1;

A = [[-(1/R_HC + 1/R_HA)/C_TH,  1/(R_HC * C_TH),                0                      ]
     [ 1/(C_TC * R_HC),        (1/R_HC - 1/R_CB - 1/R_CA)/C_TC, 1/(C_TC * R_CB)        ]
     [ 0,                      1/(C_TB * R_CB),                 -(1/R_BA + 1/R_CB)/C_TB]]

B = [[1/(C_TH * R_HA), 1/(C_TH * R_H)]
     [1/(C_TC * R_CA), 0             ]
     [1/(C_TB * R_BA), 0             ]] 

T_H = T_A;
T_C = T_A;
T_B = T_A;

T = [T_H, T_C, T_B]';