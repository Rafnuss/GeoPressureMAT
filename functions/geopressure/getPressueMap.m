function [pres_prob, pres_thr, pres_n] = getPressueMap(pressure,sta,lon,lat)
%GETPRESSUEMAP Summary of this function goes here
%   Detailed explanation goes here

%% Prepare data 

% Convert to table
data = struct2table(pressure);

% create label
Ml = (sta.start<pressure.date' & pressure.date' < sta.end) .* sta.staID;
assert(sum(sum(Ml>0)>1)==0,'pressure date belong to two stationary period')
data.label = sum(Ml)';

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

%%
% Assert grid and define scale
assert(std([diff(lon(:)); diff(lat(:))])<1e-5,'lat lon need to have the a unique and constant resolution')
res = diff(lon(1:2));
assert(res>=0.1,'A resolution smaller than 10 will resulting in interpolating ERA5 data.')
scale = 1/res; % 10 -> 0.1° | 4 -> 0.25°

% Assert time
assert(min(data.date)>=posixtime(datetime('1-jan-1981')) & max(data.date)<=posixtime(datetime()-3*31),'Time should be between 1-janv-1981 and three month from now')
assert(all(diff(data.date)/60/60>=1),'time resoution of pressure should be greater than 1hour')


base_url = "http://glp.mgravey.com/GeoPressure/v1/map";
%base_url = "http://glp.mgravey.com/test.py?";

queries = webwrite(base_url,...
    "W",min(lon)-res/2,"S",min(lat)-res/2,"E",max(lon)+res/2,"N",max(lat)+res/2,...
    "scale",scale,... 
    "time", jsonencode(data.date),"pressure",jsonencode(data.obs),"label",jsonencode(data.label),...
    "maxSample",250,"margin",30,...
    weboptions('Timeout',60*5));

if (queries.status~="success")
    error('Error in the request')
    return
else
    disp('Queries received successfully')
end
%display(urls)
%labelName = fieldnames(queries);

pres_n = nan(1,height(sta));

pres_prob = nan(numel(lat),numel(lon),height(sta));
pres_thr = false(numel(lat),numel(lon),height(sta));

% sort by longer first
sta.dur=sta.end-sta.start;
sta = sortrows(sta,'dur','descend');

clear f;
for i_l=1:numel(queries.data.urls)
    f(i_l) = parfeval(gcp, @(x) webread(x,weboptions('Timeout',60*5)),1,queries.data.urls{i_l});
end

wait(f);

for i_s=1:height(sta)
    
    i_l = find(sta.staID(i_s) == queries.data.labels);
    
    if (~isempty(i_l))
        
        A = f(i_l).OutputArguments{1};

        % Get pressure threashold
        pres_thr(:,:,sta.staID(i_s)) = flipud(A(:,:,2)>.1);

        % number of element
        pres_n(sta.staID(i_s)) = sum(data.label==queries.data.labels(i_l));

        % Compute the weight
        w = log(pres_n(sta.staID(i_s)))-1;

        % replace water by nan
        A(A==0)=nan;
        %pres_rmse = w .* mean((et/s).^2,3,'omitnan');
        %pres_prob(:,:,i_s) = exp(-pres_rmse);

        pres_prob(:,:,sta.staID(i_s)) = flipud(exp(-w.*A(:,:,1)/1e4./(sta.s(i_s).^2)));
        
%         figure; hold on;
%         mask = 0.3+0.7*double(pres_thr(:,:,sta.staID(i_s)));
%         img_tmp = real2rgb(pres_prob(:,:,sta.staID(i_s)),colormap);
%         img_tmp = img_tmp.*mask;
%         imagesc(lon,lat,img_tmp)
%         borders('countries','w')
%         axis equal; axis([min(lon) max(lon) min(lat) max(lat) ]);
    end
end

end