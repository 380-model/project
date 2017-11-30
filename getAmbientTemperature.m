function [ ambientTemperature ] = getAmbientTemperature( month )
%getAmbientTemperature returns the ambient temperature for the month
%   The month is 0 to 11 with decimals inbetween (ex 3.5 is 15th of March)
    
    % Location: Prince George
    temperatureAverage = 3.68;
    temperatureRange = 25.2;
    
    ambientTemperature = ((sin((month/12) * 2 * pi - pi/2))/2)*temperatureRange + temperatureAverage;
end

