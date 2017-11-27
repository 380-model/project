close all
y0 = [0]
% options = odeset('RelTol',1e-2,'AbsTol',1e-4);
[t,y] = ode45(@turbine, [0 1000], y0)

velocity = sqrt(2*9.8*(50 - (5.2576E3 * y/(1000*9.8))))

plot(t,y)



function xprime = turbine(t,y)
pipe_diam = 0.0736;
pipe_length = 50*sin(45);
pipe_area = (pipe_diam/2)^2 * pi;
Rf = (128 * 8.9E-4 * pipe_length) / (pi * pipe_diam^4);
If = 2 * 1000 * pipe_length / pipe_area;
head = 50;

xprime(1,1) =  (1000 * 9.8 * head - Rf*y(1)) / If;

end