close all

%%  Tank Flow Out
%%  Initial Conditions
P_1a = 3;   % [atm]
Q = 0;

Rf = 7;
If = 5;
V = 0.5;
k = 1.405;
P = P_1a/2;
Cf = V/(P*k)

IC = [Q
      P_1a];  
  
dt = 0.001;

%%  Define State Space Matrix
A = [-Rf/If     1/If
     -1/Cf      0   ];
 
 B = [0
      0];
  
 C = [1 0
      0 1];
 
 D = [0
      0];
  
 sys = ss(A, B, C, D);
 
 t = 0:dt:15;
 u = zeros(length(t), 1);
 
 %% Simulation
 y = lsim(sys, u, t, IC);
 
 %% Plot Results

figure(1)
subplot(2,1,1)
plot (t,y(:,1),'k','LineWidth',2)
hold on
grid on
myxlabel=xlabel ('time [s]')
myylabel=ylabel ('Flow Rate')
set (myxlabel,'FontSize',12);
set (myylabel,'FontSize',12);
title('MSE 380 Assignment 1 - Q1 - y1')

subplot(2,1,2)
plot (t, y(:,2),'k', 'LineWidth', 2)
hold on
grid on
myxlabel=xlabel ('time [s]');
myylabel=ylabel ('Pressure');
set (myxlabel,'FontSize',12)
set (myylabel,'FontSize',12)
title('MSE 380 Assignment 1 - Q1 - y2')

