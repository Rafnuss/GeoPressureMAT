function twl = twilight(tm, lon, lat, rise, varagin) 
%TWILIGHT estimate time of sunrise or sunset for a given day and location
%   Twilight uses an iterative algorithm to estimate times of sunrise and
%   sunset. Note that these functions return the twilight that occurs on
%   the same date GMT as tm, and so sunset may occur before sunrise,
%   depending upon latitude. Solar declination and equation of time vary
%   slowly over the day, and so the values of the Solar declination and
%   equation of time at sunrise/sunset are well approximated by their
%   values at 6AM/6PM local time. The sun's hour angle and hence
%   sunrise/sunset for the required zenith can then be caclulates from
%   these approximations. The calculation is then repeated using the
%   approximate sunrise/sunset times to derive more accurate values of the
%   Solar declination and equation of time and hence better approximations
%   of sunrise/sunset. The process is repreated and is accurate to less
%   than 2 seconds within 2 or 3 iterations.

if nargin>4
    zenith = varagin(1);
else
    zenith = 96;
end

if nargin>5
    iters = varagin(2);
else
    iters = 10;
end

if nargin>6
    closest = varagin(3);
else
    closest = false;
end

lat=lat(:)';
lon=lon(:)';


date = dateshift(tm,'start','day');
lon = mod(lon + 180,360) - 180;
twl = date + 240 .* (90+180*~rise - lon)/60/60/24;

for k =1:iters
    s = solar(twl);
    s.solarTime = mod(s.solarTime,360);
    solarTime = 4 * twilightSolartime(s, lon, lat, rise, zenith) - s.eqnTime;
    twl = date + 60 * solarTime / 60/60/24;
end

if (closest)
    error('nor corrected')
    delta = (as.numeric(tm) - as.numeric(twl))/3600;
    off = double(length(delta));
    off(delta > 12) = 86400;
    off(delta < -12) = -86400;
    twl = twilight(tm + off, lon, lat, rise, zenith, iters, false);
end

end