%%  Define Model Inputs
PH2 = 0.2;              %[atm]
PO2 = 0.18;             %[atm]

%%  Define Constants
T = 343;                %70C
Rc = 0.0003;            %Ohms
A = 430;                %[cm^2]
t = 0.000178;           %[m]
n = 1;                  %number of cells
psi = 13;               %saturated
i_max = 2;              %max current density [A/cm2]

CH2 = PH2/((5.08*10^6)*exp(-498/T));
zeta1 = -0.948;
zeta2 = (0.00286+0.0002*log(A)+(4.38*10^-5)*log(CH2));
zeta3 = 7.6*(10^-5);
zeta4 = -1.93*(10^-4);


%   Preallocate Array
y = zeros(2001,6);

for i = 0:0.001:2
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
    Vstack = n*(Vn - Vk - Vo - Vt);
    
    %   Acts weird
    y(round(i*1000+1), 1) = i;
    y(round(i*1000+1), 2) = Vstack;
    y(round(i*1000+1), 3) = Vn;
    y(round(i*1000+1), 4) = Vk;
    y(round(i*1000+1), 5) = Vo;
    y(round(i*1000+1), 6) = Vt;
end

figure(1)

subplot(2,3,1)
plot(y(:,1),y(:,2))
title('Stack Voltage Pol Curve at 70{\circ}C')
xlabel('Current Density [A/cm2]');
ylabel('Cell Voltage [V]');
grid on

subplot(2,3,2)
plot(y(:,1),y(:,3))
title('Nernest Voltage Pol Curve at 70{\circ}C')
xlabel('Current Density [A/cm2]');
ylabel('Cell Voltage [V]');
grid on

subplot(2,3,3)
plot(y(:,1),y(:,4))
title('Kinetic Voltage Loss Pol Curve at 70{\circ}C')
xlabel('Current Density [A/cm2]');
ylabel('Cell Voltage [V]');
grid on

subplot(2,3,4)
plot(y(:,1),y(:,5))
title('Ohmic Voltage Loss Pol Curve at 70{\circ}C')
xlabel('Current Density [A/cm2]');
ylabel('Cell Voltage [V]');
grid on

subplot(2,3,5)
plot(y(:,1),y(:,6))
title('Mass Transport Voltage Loss Pol Curve at 70{\circ}C')
xlabel('Current Density [A/cm2]');
ylabel('Cell Voltage [V]');
grid on

subplot(2,3,6)
plot(y(:,1),y(:,2).*y(:,1))
title('Power Pol Curve at 70{\circ}C')
xlabel('Current Density [A/cm2]');
ylabel('Power [Watts]');
grid on