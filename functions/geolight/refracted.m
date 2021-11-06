function out = refracted(zenith)
%REFRACTED Adjust the solar zenith angle for atmospheric refraction.
%   Given a vector of solar zeniths computed by zenith, refracted
%   calculates the solar zeniths adjusted for the effect of atmospheric
%   refraction.

elev = 90 - zenith;
te = tand(elev);
r = zeros(size(elev));
id = elev<=85&elev>5;
r(id) = 58.1/te(id) - 0.07/te(id).^3 + 8.6e-05/te(id).^5;
id = elev<=5& elev > -0.575;
r(id) = 58.1/te(id) - 0.07/te(id).^3 + 8.6e-05/te(id).^5;
id = elev<= -0.575;
r(id) = -20.772/te(id);
out = zenith - r/3600;
end