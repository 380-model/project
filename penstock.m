close all
y0 = [5/1000];
[t,y] = ode45(@turbine, [0 5], y0);

pipe_diam = 0.0736;
jet_diam = 0.015;
jet_area = pi*(jet_diam/2)^2;
beta = jet_diam/pipe_diam;
head = 50;
jet_coefficient = 0.97;
pipe_length = 50*sin(45);
pipe_area = (pipe_diam/2)^2 * pi;
Rf = (128 * 8.9E-4 * pipe_length) / (pi * pipe_diam^4); 
If = 2 * 1000 * pipe_length / pipe_area;

plot(t,y)

function xprime = turbine(t,y)
pipe_diam = 0.0736;
jet_diam = 0.015;
jet_area = pi*(jet_diam/2)^2;
beta = jet_diam/pipe_diam;
jet_coefficient = 0.97;

head = 50;
pipe_length = 50*sin(45);
pipe_area = (pipe_diam/2)^2 * pi;
Rf = (128 * 8.9E-4 * pipe_length) / (pi * pipe_diam^4); 
If = 2 * 1000 * pipe_length / pipe_area;

xprime(1,1) = (1000 * 9.8 * head - Rf*y(1) - 0.5*1000*(1-beta^4)*(y(1)/(jet_coefficient*jet_area))^2)/If;
jet_velocity = sqrt(2*9.8*(head - (Rf * y(1) / 1000 / 9.8)));

end