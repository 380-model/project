function [ p_H2, p_O2 ] = gaspressure_nonlin( y,a )
%GASPRESSURE_NONLIN Summary of this function goes here
%   Detailed explanation goes here
volume = a.volume; % 20 L tanks (each gas)
p_H2 = y(:,7);
p_O2 = y(:,8);
p_H2 = p_H2(end) * 8.314 * 298.15 / volume; % pascals
p_O2 = p_O2(end) * 8.314 * 298.15 / volume; % pascals

end

