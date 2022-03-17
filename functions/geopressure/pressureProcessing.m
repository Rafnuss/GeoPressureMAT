function data = pressureProcessing(pressure,sta)

% Convert to table
data = struct2table(pressure);

% create label
Ml = (sta.start<pressure.date' & pressure.date' < sta.end) .* sta.staID;
% assert(sum(sum(Ml>0)>1)==0,'pressure date belong to two stationary period')
data.label = sum(Ml,1)';

% remove flight period
data.obs(data.label==0)=nan;

% smooth, and convert to pascal, round to unit
data.obs = movmean(data.obs,3);
data.obs = round(data.obs*100);

% remove outliar
data = data(~isnan(data.obs),:);

% Downscale to 1hr
t_1 = dateshift(pressure.date(1),'start','hour'):1/24:dateshift(pressure.date(end),'end','hour');
data =  data(ismember(data.date,t_1),:);

% convert to UNIX time in int
data.date = uint32(posixtime(data.date));

%
% data = data(data.label==1,:);

% Assert time
assert(min(data.date)>=posixtime(datetime('1-jan-1981')) & max(data.date)<=posixtime(datetime()-3*31),'Time should be between 1-janv-1981 and three month from now')
assert(all(diff(data.date)/60/60>=1),'time resoution of pressure should be greater than 1hour')


end