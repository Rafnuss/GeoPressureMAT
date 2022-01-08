function out = getPressueTimeseries(map,sta,lon,lat)
%GETPRESSUEMAP Summary of this function goes here
%   Detailed explanation goes here

[gLON,gLAT] = meshgrid(lon,lat);

base_url = "http://glp.mgravey.com/GeoPressure/v1/timeseries";

clear f;
for i_s=1:height(sta)
    
    [~,id]=max(reshape(map(:,:,i_s),1,[]));
    % plot(gLON(id), gLAT(id),'or','linewidth',2,'MarkerSize',12)

    f(i_s) = parfeval(gcp, @(lon,lat,startTime, endTime) webwrite(base_url,...
        "lon",lon,"lat",lat,...
        "startTime",startTime,"endTime",endTime,... 
        weboptions('Timeout',60*5))...
        ,1,gLON(id),gLAT(id),uint32(posixtime(sta.start(i_s))),uint32(posixtime(sta.end(i_s))));
end

wait(f)

clear f2;
for i_s=1:height(sta)
    if (f(i_s).OutputArguments{1}.status~="success")
        error('Error in the request')
        return
    end
    f2(i_s) = parfeval(gcp, @(x) webread(x,weboptions('Timeout',60*5)),1,f(i_s).OutputArguments{1}.data.url);
end

wait(f2)

c=cellfun(@(x) x{1},{f2.OutputArguments},'UniformOutput',false);
out = vertcat(c{:});

out.time=datetime(out.time,'ConvertFrom','posixtime');
out.pressure=out.pressure/100; % pa to hpa;

end