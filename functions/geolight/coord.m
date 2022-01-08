function out = coord(twl, varargin)
%COORD estimate location from consecutive twilights
%   This function estimates the location given the times at which the
%   observer sees two successive twilights. Longitude is estimated by
%   computing apparent time of local noon from sunrise and sunset, and
%   determining the longitude for which this is noon. Latitude is estimated
%   from the required zenith and the sun's hour angle for both sunrise and
%   sunset, and averaged. When the solar declination is near zero (at the
%   equinoxes) latitude estimates are extremely sensitive to errors. Where
%   the sine of the solar declination is less than tol, the latitude
%   estimates are returned as NA. The format (date and time) of tFirst and
%   tSecond has to be "yyyy-mm-dd hh:mm" corresponding to Universal Time
%   Zone UTC (see: as.POSIXct, time zones)
%   
%   Modified with ThresholdEstimate for consistency

% if nargin>1
%     tEMed = varargin{1};
% else
%     tEMed = 0;
% end
if nargin>1
    z = varargin{1}(:);
    if numel(z)==2
        zRise = z(1);
        zSet = z(2);
    elseif (numel(z)==height(twl))
        zRise = z(twl.Rise);
        zSet = z(~twl.Rise);
    elseif numel(z)>2
        zRise = z(:)';
        zSet = z(:)';
    else
        zRise = z;
        zSet = z;
    end
else
    zRise = 96;
    zSet =  96;
end
if nargin>2
    tol = varargin{2};
else
    tol = 0.05;
end

% sr = solar(twl.Twilight(twl.Rise)-tEMed/24/60);
% ss = solar(twl.Twilight(~twl.Rise)+tEMed/24/60);

% Take only pair number of twilight
if height(twl.Twilight)>1 && mod(height(twl.Twilight),2)~=0
    twl(end,:)=[];
end

sr = solar(twl.Twilight(twl.Rise));
ss = solar(twl.Twilight(~twl.Rise));

lon = -(sr.solarTime + ss.solarTime + 360*(sr.solarTime < ss.solarTime))/2;
lon = mod(lon + 180, 360) - 180;
hourAngle = sr.solarTime + lon - 180;

a = sr.sinSolarDec;
b = sr.cosSolarDec .* cosd(hourAngle);
x = (a .* cosd(zRise) - sign(a) .* b .* realsqrtnan(a.^2 + b.^2 - cosd(zRise).^2, tol))./(a.^2 + b.^2);
latRise = asind(x);
% latRise(abs(a) <= tol)=nan;

hourAngle = ss.solarTime + lon - 180;
a = ss.sinSolarDec;
b = ss.cosSolarDec .* cosd(hourAngle);
x = (a .* cosd(zSet) - sign(a) .* b .* realsqrtnan(a.^2 + b.^2 - cosd(zSet).^2, tol))./(a.^2 + b.^2);
latSet = asind(x);
% latSet(abs(a) <= tol)=nan;

dim = numel(size(latRise))+1;
lat = mean(cat(dim,latRise,latSet),dim);
% lon=repmat(lon,1,numel(zernith));

out = table(lon, lat, latRise, latSet,'VariableNames',{'lon','lat','latRise','latSet'});

if any(strcmp('isOutliar',twl.Properties.VariableNames))
    out.isOutliar = twl.isOutliar(twl.Rise) | twl.isOutliar(~twl.Rise);
end

end

function x = realsqrtnan(X,tol)
x=real(sqrt(X)); 
x(x<tol)=nan;
end
