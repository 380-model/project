% A turbine simulation requires the input pressure and flow rate

% River Parameters http://www.homepower.ca/data_tables.htm
head = 50; % metres, to calculate pressure
P_pa = head * 1000 * 9.8;  % should we consider dynamic pressure (1/2 * rho * v^2)?
Q_lps = 5; %litres per second

% model a known micro hydro thing to get the pump parameters
% in the main model, use those parameters and also use 
