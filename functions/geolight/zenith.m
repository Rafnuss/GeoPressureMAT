function out = zenith(Twilight, known_coord)
%ZENITH calculate the solar zenith angle for given times and locations
%   zenith uses the solar time and declination calculated by solar to
%   compute the solar zenith angle for given times and locations, using the
%   same methods as www.esrl.noaa.gov/gmd/grad/solcalc/. This function does
%   not adjust for atmospheric refraction see refracted.

s = solar(Twilight);
hourAngle = s.solarTime + known_coord(:,1)' - 180;
cosZenith = sind(known_coord(:,2)') .* s.sinSolarDec + cosd(known_coord(:,2)') .* s.cosSolarDec .* cosd(hourAngle);
cosZenith(cosZenith > 1) = 1;
cosZenith(cosZenith < -1) = -1;
out = acosd(cosZenith);
end
