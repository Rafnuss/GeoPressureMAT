function out = getPressueTimeseries(sta,lon,lat, varargin)
%GETPRESSUEMAP Summary of this function goes here
%   Detailed explanation goes here

base_url = "http://glp.mgravey.com:24853/GeoPressure/v1/timeseries";

clear f;

if nargin>3
   data = pressureProcessing(varargin{1},sta);
   query_staID=unique(data.label);
else
    query_staID=1:height(sta);
end

for i_qs=1:numel(query_staID)
    i_s=query_staID(i_qs);

    if nargin>3
        f(i_qs) = parfeval(gcp, @(lon, lat, date, obs) webwrite(base_url,...
            "lon",lon,"lat",lat,...
            "time", jsonencode(date),"pressure",jsonencode(obs),... 
            weboptions('Timeout',60*5))...
            ,1,lon(i_s),lat(i_s),data.date(data.label==i_s),data.obs(data.label==i_s));
    else
        f(i_qs) = parfeval(gcp, @(lon,lat,startTime, endTime) webwrite(base_url,...
            "lon",lon,"lat",lat,...
            "startTime",startTime,"endTime",endTime,... 
            weboptions('Timeout',60*5))...
            ,1,lon(i_s),lat(i_s),uint32(posixtime(sta.start(i_s))),uint32(posixtime(sta.end(i_s))));
    end
end

wait(f)

clear f2;
for i_qs=1:numel(query_staID)
    if (f(i_qs).OutputArguments{1}.status~="success")
        error('Error in the request')
        return
    end
    f2(i_qs) = parfeval(gcp, @(x) webread(x,weboptions('Timeout',60*5)),1,f(i_qs).OutputArguments{1}.data.url);
end

wait(f2)

c=cellfun(@(x) x{1},{f2.OutputArguments},'UniformOutput',false);
for i_qs=1:numel(query_staID)
    if ~isempty(c{i_qs})
        c{i_qs}.label(:) = query_staID(i_qs);
    end
end
out = vertcat(c{:});

out.time=datetime(out.time,'ConvertFrom','posixtime');
out.pressure=out.pressure/100; % pa to hpa;

end