function out = solar(tm)
%SOLAR Calculate solar time, the equation of time and solar declination
%   The solar time, the equation of time and the sine and cosine of the
%   solar  declination are calculted for the times specified by tm using
%   the same methods as www.esrl.noaa.gov/gmd/grad/solcalc/.

rad = pi/180;
if all(isdatetime(tm))
    Jd = posixtime(tm)/86400 + 2440587.5;
else
    tm2=(tm-719529)*24*60*60; % datenum('1-Jan-1970')
    Jd = tm2/86400 + 2440587.5;
end
Jc = (Jd - 2451545)/36525;
L0 = mod(280.46646 + Jc .* (36000.76983 + 0.0003032 * Jc),360);
M = 357.52911 + Jc .* (35999.05029 - 0.0001537 * Jc);
e = 0.016708634 - Jc .* (4.2037e-05 + 1.267e-07 * Jc);
eqctr = sind(M) .* (1.914602 - Jc .* (0.004817 + 1.4e-05 .* Jc)) + sind(2 .* M) .* (0.019993 - 0.000101 .* Jc) +  sind(3 .* M) .* 0.000289;
lambda0 = L0 + eqctr;
omega = 125.04 - 1934.136 * Jc;
lambda = lambda0 - 0.00569 - 0.00478 * sind(omega);
seconds = 21.448 - Jc .* (46.815 + Jc .* (0.00059 - Jc *(0.001813)));
obliq0 = 23 + (26 + (seconds/60))/60;
omega = 125.04 - 1934.136 * Jc;
obliq = obliq0 + 0.00256 * cosd(omega);
y = tand(obliq/2).^2;
eqnTime = 4/rad .* (y .* sind(2 .* L0) - 2 .* e .* sind(M) + 4 .* e .* y .* sind(M) .* cosd(2 * L0) - 0.5 .* y.^2 .* sind(4 .* L0) - 1.25 .* e.^2 .* sind(2 .* M));
solarDec = asin(sind(obliq) .* sind(lambda));
sinSolarDec = sin(solarDec);
cosSolarDec = cos(solarDec);
solarTime = (mod((Jd - 0.5),1) * 1440 + eqnTime)/4;
out = table(solarTime, eqnTime, sinSolarDec, cosSolarDec,'VariableNames',{'solarTime', 'eqnTime', 'sinSolarDec', 'cosSolarDec'});
end
