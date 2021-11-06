function out = twilightSolartime(solar, lon, lat, rise, varagin) 
%TWILIGHTSOLARTIME estimate time of sunrise or sunset for a given location
%given the approximate solar time of twilight
%   Solar declination and equation of time vary slowly over the day, and so
%   the values of the Solar declination and equation of time at
%   sunrise/sunset can be caclulated approximately if an approximate time
%   of sunrise/sunset is known. The sun's hour angle and hence
%   sunrise/sunset for the required zenith can then be calculated from
%   these approximations. Note this function returns the time of twilight
%   in solar time.

if nargin>4
    zenith=varagin(1);
else
    zenith=96;
end
rad = pi/180;
cosz = cosd(zenith);
cosHA = (cosz - sind(lat) .* solar.sinSolarDec)./(cosd(lat) .* solar.cosSolarDec);
hourAngle = rise*360 + (1-2*rise) .* acos(cosHA)/rad;
solarTime = mod((hourAngle + 180 - lon),360);
out = mod(solarTime - solar.solarTime + 180,360) - 180 + solar.solarTime;
end