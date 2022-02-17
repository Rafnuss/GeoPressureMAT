function P = altitude2pressure(dh,varargin)
% PRESSUREALT  Approximate altitude in troposphere for a given pressure.
% 
%     P in hectoPascal

% standard temperature lapse rate [K/m] = -0.0065 [K/m]
Lb = -0.0065;
% universal gas constant = 8.31432 [N * m / mol /K]
R = 8.31432;
% gravitational acceleration constant = 9.80665 [m/s^2]
g0 = 9.80665;
% molar mass of Earth’s air = 0.0289644 [kg/mol]
M = 0.0289644;
%  standard temperature (temperature at sea level) [K]
T0 = 273.15+15;

%  static pressure (pressure at sea level) [Pa]
if nargin > 1
    P0 = varargin{1};
else
    P0 = 1013.25;
end

P = P0 .* (dh.*Lb/T0 +1).^(-g0*M/Lb/R);

end