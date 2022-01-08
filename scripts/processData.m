
% Setting-up
tic

project="StudyPressure";

addpath(genpath('../functions'))
load coastlines

% Default figure
set(0,'DefaultAxesBox','on')
% set(0,'defaultLineLineWidth',2)
% set(0,'DefaultFigureColormap',crameri('batlow'));
% colorder = [17, 138, 178;82, 43, 41; 204, 164, 59; 67, 100, 54; 228, 87, 46]/255;
% set(0,'defaultAxesColorOrder',colorder)
colorder = get(0,'defaultAxesColorOrder');
colorder = repmat(colorder,10,1);

% Read metadata
tblLog = readtable('../data/tracks/rawData.xlsx');
tblLog = tblLog(tblLog.Usable,:);
tblLog = tblLog(tblLog.(project),:);
tblLog.Color(cellfun(@isempty,tblLog.Color))={'000000'};
tblLog.Color = hex2rgb(tblLog.Color);

%tblLog=tblLog(1,:);

%% Load DEM and ERA5 data
load('../data/DEM/MatlabData.mat')

% Create grid
lat=cell(height(tblLog),1);
lon=cell(height(tblLog),1);
for lt=1:height(tblLog)
    lon{lt} = glon(glon>=tblLog.bndy_W(lt) & glon<=tblLog.bndy_E(lt));
    lat{lt} = glat(glat>=tblLog.bndy_S(lt) & glat<=tblLog.bndy_N(lt));
end

%% Load and process geolocator data
raw=cell(height(tblLog),1);
twl=cell(height(tblLog),1);
sta=cell(height(tblLog),1);
activityMig = cell(height(tblLog),1);
gE=cell(height(tblLog),1);
for lt=1:height(tblLog)
    raw{lt} = importPAM(['../data/tracks/' tblLog.DataFileName{lt}], tblLog.SubsetStart(lt), tblLog.SubsetEnd(lt));
    if any(raw{lt}.light.obs<0)
        warning('Removing negative light value')
        raw{lt}.light.obs(raw{lt}.light.obs<0) = nan;
    end
    raw{lt}.light.obstrans = log(double(raw{lt}.light.obs)+0.0001) + abs(min(log(double(raw{lt}.light.obs)+0.0001)));
    
    % Add calibration period
    raw{lt}.calib.first_period = [tblLog.CalibFirstStart(lt) tblLog.CalibFirstEnd(lt)];
    raw{lt}.calib.second_period = [tblLog.CalibSecondStart(lt) tblLog.CalibSecondEnd(lt)];
    raw{lt}.calib.lon = tblLog.LongitudeAttached(lt); % Inverted in the spreadsheet
    raw{lt}.calib.lat = tblLog.LatitudeAttached(lt);
    
    twl{lt} = twilightEdit(findTwilightsRaf(raw{lt}.light),false);
    twl{lt} = twilightEditTrainset(raw{lt},twl{lt},false);
    
    id_calib = raw{lt}.calib.first_period(1)<twl{lt}.Twilight&twl{lt}.Twilight<raw{lt}.calib.first_period(2) | raw{lt}.calib.second_period(1)<twl{lt}.Twilight&twl{lt}.Twilight<raw{lt}.calib.second_period(2);
    gE{lt} = getElevation(twl{lt}(id_calib,:),[raw{lt}.calib.lon raw{lt}.calib.lat],false,80:100);
    gE{lt}.id_calib = id_calib';
    bd=1.2;
    xi = 60:.01:120;%floor(min(gE{lt}.z(~twl{lt}.isOutliar(id_calib)))):0.01:ceil(max(gE{lt}.z(~twl{lt}.isOutliar(id_calib))));
    [gE_f,gE_xi]=ksdensity(gE{lt}.z(~twl{lt}.isOutliar(gE{lt}.id_calib)),xi,'Bandwidth',bd);
    gE{lt}.pdf = @(x) interp1(gE_xi,gE_f,x,'nearest','extrap');
    
    % Compute Activity Data
    activity = classifyActivityTrainsetRaf(raw{lt},[]);
    activityMig{lt} = activity(activity.km==3,:);

    % Compute stationary period
    sta{lt} = table();
    sta{lt}.start = [raw{lt}.pressure.date(1);activityMig{lt}.date_max];
    sta{lt}.end = [activityMig{lt}.date_min;raw{lt}.pressure.date(end)];
    %delete stationary period below 7 hours
    % sta{lt}(hours(sta{lt}.end-sta{lt}.start)<7,:)=[];
    for i_s=1:height(sta{lt})
        if i_s~=height(sta{lt})
            idAct = sta{lt}.end(i_s) <= activityMig{lt}.date_min &  activityMig{lt}.date_max <= sta{lt}.start(i_s+1);
            sta{lt}.actNb(i_s) = sum(idAct);
            sta{lt}.actSum(i_s) = sum(activityMig{lt}.act_sum(idAct));
            sta{lt}.actDuration(i_s) = sum(activityMig{lt}.duration(idAct));
            activityMig{lt}.staID(idAct)=i_s;
        end
        idtwl = sta{lt}.start(i_s) <= twl{lt}.Twilight & twl{lt}.Twilight <= sta{lt}.end(i_s);
        sta{lt}.twlNb(i_s) = sum(idtwl & ~twl{lt}.isOutliar);
        twl{lt}.staID(idtwl) = i_s;
    end
    sta{lt}.staID = (1:height(sta{lt}))';
    if any(isnan(sta{lt}.actSum)) || all(sta{lt}.actSum==0)
        sta{lt}.actEffort =  sta{lt}.actDuration;
        warning(['Cannot compute effort of flight for' raw{lt}.GDL_ID])
    else
        sta{lt}.actEffort =  sta{lt}.actSum./sum(sta{lt}.actSum)*sum(sta{lt}.actDuration);
    end


    % label stationary period & Migration
    sta{lt}.status=repmat("",height(sta{lt}),1);
    [~,id]=maxk(sta{lt}.end-sta{lt}.start,3);
    if any(tblLog.GDL_ID{lt}==["20OA"])
        iws=id(2);
    else
        iws = id(1);
    end
    
    if ~any(strcmp(tblLog.GDL_ID{lt},["16IT","16JB","TZ","DK","24EA","24EP","24IS"]))
        sta{lt}.status(iws)='wintering';
    end
    sta{lt}.status(1:height(sta{lt})<iws)='migration_post_equipement';
    sta{lt}.status(1:height(sta{lt})>iws)='migration_pre_retrival';
    sta{lt}.status(1)='equipment';
    if ~isnat(tblLog.CalibSecondStart(lt))
        sta{lt}.status(end)='retrieval';
    end

    
    % Add pressure outliar
    filename = ['../data/labels/activity_label/' raw{lt}.GDL_ID '_act_pres-labeled.csv'];
    if exist(filename)>0
        T2 = readtable(filename);
        % tmp = T2.timestamp(strcmp(T2.series,'pres'));
        % tmp = datetime(cellfun(@(x) [x(1:10) ' ' x(12:16)],tmp,'UniformOutput',false));
        raw{lt}.pressure.isOutliar = ~isnan(T2.label(strcmp(T2.series,'pres')));
        raw{lt}.pressure.obsWithOutliars=raw{lt}.pressure.obs;
        raw{lt}.pressure.obs(raw{lt}.pressure.isOutliar)=nan;
    end
    % Add day/night information
    % raw{lt}.pressure.isDay = isDay(raw{lt}.pressure.date,twl{lt});
end
clear activity idtwl idAct bd i_s km xi bd id_calib

%% Assert stopover duration
% sta{lt}.start(hours(sta{lt}.end-sta{lt}.start)<5)

%% Light Map

light_prob = cell(height(tblLog),1);
for lt=1:height(tblLog)
    % Compute probability map
    tmp = coordMapProb(twl{lt}, gE{lt}, lon{lt}, lat{lt});
    w=0.1;
    light_prob{lt} = nan(numel(lat{lt}),numel(lon{lt}),height(sta{lt}));
    for i_s = 1:height(sta{lt})
        light_prob{lt}(:,:,i_s) = exp(w*sum(log(tmp(:,:,twl{lt}.staID==i_s&~twl{lt}.isOutliar)),3));
    end
    % This normalization is still not clear to me! We are in one way ignoring the relative number of twl per stationary.
    % light_prob = light_prob ./sum(sum(light_prob,1),2);
end

%% Test on prob aggr
% lt=7
% A=tmp(:,:,twl{lt}.staID==i_s&~twl{lt}.isOutliar);
% 
% figure; hold on;
% % imagesc(lon{lt},lat{lt},A(:,:,17));
% % imagesc(lon{lt},lat{lt},A(:,:,17)>=max(max(A(:,:,17))-.0005));
% n=2;
% alpha=0.1;
% w=alpha+(1-alpha)*1/n;
% imagesc(lon{lt},lat{lt},prod(A(:,:,1:n).^w,3));
% borders('countries','w');
% plot(raw{lt}.calib.lon,raw{lt}.calib.lat,'xr','linewidth',2,'MarkerSize',12)
% axis equal; axis([13 40 -33 -10 ]);
% xticks('');yticks('')



%% Mask water
tmp = find(isnan(coastlat));
mask_water_g = true(numel(glat),numel(glon));
[gLON,gLAT] = meshgrid(glon,glat);
for id=1:numel(tmp)-1
    mask_water_g(inpolygon(gLAT,gLON,coastlat(tmp(id):tmp(id+1)),coastlon(tmp(id):tmp(id+1))))=false;
end
mask_water = cell(height(tblLog),1);
for lt=1:height(tblLog)
    mask_water{lt}=mask_water_g(glat>=tblLog.bndy_S(lt) & glat<=tblLog.bndy_N(lt), glon>=tblLog.bndy_W(lt) & glon<=tblLog.bndy_E(lt));
end

% light_prob{lt}(repmat(mask_water,1,1,height(sta{lt})))=0;


%% Pressure Map

% Load Pressure data
file='../data/ECMWF/surface_pressure.nc'; %ncdisp(file);
% glon=double(ncread(file,'longitude'));
% glat=flip(double(ncread(file,'latitude')));
spttime = datetime(double(ncread(file,'time'))/24 + datenum('1900-01-01 00:00:00'),'convertFrom','datenum');

% Allow for some error for the threahodl of pressure
dp_margin = [3 -3];

pres_rmse = cell(height(tblLog),1);
pres_thr = cell(height(tblLog),1);
pres_prob = cell(height(tblLog),1);
pres_n = cell(height(tblLog),1);

for lt=1:height(tblLog)
    
    id_lon = find(glon>=tblLog.bndy_W(lt) & glon<=tblLog.bndy_E(lt));
    id_lat = find(glat>=tblLog.bndy_S(lt) & glat<=tblLog.bndy_N(lt));
    id_t = find(raw{lt}.pressure.date(1)<=spttime & spttime <= raw{lt}.pressure.date(end));
    splt = permute(flip(ncread(file,'sp',[id_lon(1) numel(glat)-id_lat(end)+1 1 id_t(1)],[numel(id_lon) id_lat(end)-id_lat(1)+1 1 id_t(end)-id_t(1)+1]),2)/100,[2 1 4 3]);
    if raw{lt}.pressure.date(end)>datetime('1-jul-2021')
        splt_0 = permute(flip(ncread(file,'sp',[id_lon(1) numel(glat)-id_lat(end)+1 2 id_t(1)],[numel(id_lon) id_lat(end)-id_lat(1)+1 1 id_t(end)-id_t(1)+1]),2)/100,[2 1 4 3]);
        splt(:,:,spttime(id_t)>=datetime('1-jul-2021'))=splt_0(:,:,spttime(id_t)>=datetime('1-jul-2021'));
    end
    
    % Calibration location
    [~,id_lat_tmp]=min(abs(raw{lt}.calib.lat-glat(id_lat)));
    [~,id_lon_tmp]=min(abs(raw{lt}.calib.lon-glon(id_lon)));
    
    % Compute Pressue error map for each stationary period
    pres_rmse{lt} = nan(numel(lat{lt}),numel(lon{lt}),height(sta{lt}));
    pres_prob{lt} = nan(numel(lat{lt}),numel(lon{lt}),height(sta{lt}));
    pres_thr{lt} = false(numel(lat{lt}),numel(lon{lt}),height(sta{lt}));
    pres_n{lt} = nan(height(sta{lt}),1);

    for i_s = 1:height(sta{lt})
        id_tgr = find(sta{lt}.start(i_s)<spttime(id_t) & spttime(id_t) < sta{lt}.end(i_s));
        id_tge = sta{lt}.start(i_s)<raw{lt}.pressure.date & raw{lt}.pressure.date < sta{lt}.end(i_s);
        % pres_rmse(:,:,i_s) = sqrt(mean(((  pres(:,:,id_tgr)-reshape(press_movmean(id_tgr),1,1,[])  ).^2),3));
        
        % Apply a movinn median to remove outliars and smooth with 1hr
        dt = 1./hours(diff(raw{lt}.pressure.date(1:2)));
        % pres_ge = movmean(movmedian(raw{lt}.pressure.obs(id_tge),3),dt);
        pres_ge = movmean(raw{lt}.pressure.obs(id_tge),3);

        % Downscale to same time scale
        id_t_gre = ismember(raw{lt}.pressure.date(id_tge),spttime(id_t(id_tgr)));
        pres_gr = pres_ge(id_t_gre);

        pres_n{lt}(i_s)=sum(~isnan(pres_gr));
        
        % id_tgr=id_tgr(~isnan(pres_gr));
        % pres_gr=pres_gr(~isnan(pres_gr));
        
        if any(~isnan(pres_gr))
            
            % PART I: Offset / Absolute pressure error
            
            % Perform a second smoothing to remove daily fluctuation
            % (mainly for pressure tide which can trigger the threadhold
            % too easily)
            pres_gr_s = movmean(pres_gr,24,'omitnan');
            pres_gr_s(isnan(pres_gr))=nan;
            splt_s = movmean(splt(:,:,id_tgr),24,3);
            
            dh = pressure2altitude(reshape(pres_gr_s,1,1,[]),splt_s);
            % use mink and maxk to not be sensitive to outliar
%             n_trim = ceil(size(dh,3).*.01)+1;
%             dh_min = mink(dh,n_trim,3); dh_min=dh_min(:,:,end);
%             dh_max = maxk(dh,n_trim,3);dh_max=dh_max(:,:,end);
            dh_min = min(dh,[],3,'omitnan');
            dh_max = max(dh,[],3,'omitnan');

            % convert the pressure margin to altitude margin
            dh_margin = pressure2altitude(mean(splt_s,3,'omitnan')+reshape(dp_margin,1,1,[]),mean(splt_s,3,'omitnan'));

            pres_thr{lt}(:,:,i_s) = (dh_min-dh_margin(:,:,1)+geoDEM(id_lat,id_lon))>=DEM_min(id_lat,id_lon) & ...
                (dh_max-dh_margin(:,:,2)+geoDEM(id_lat,id_lon))<=DEM_max(id_lat,id_lon);

            if all(~pres_thr{lt}(:,:,i_s),'all')
               warning(['no pressure threashold match for this stationary period: ' num2str(i_s) '-' tblLog.GDL_ID{lt}])
            end
            
            % PART II: temporal change
            x = splt(:,:,id_tgr)-mean(splt(:,:,id_tgr),3,'omitnan');
            y = reshape(pres_gr-mean(pres_gr,'omitnan'),1,1,[]);
            e = x-y;
 
            if tblLog.Truncated(lt) && i_s~=0 && i_s~=height(height(sta{lt})) && min(min(std(e,[],3)))>1
                e(:)=0;
            end
            
            % autocovariance
            % C = @(h) .9*exp(-h/4) + .1*exp(-h/24/5).*cos(h/24*2*pi);
            % figure; plot(0:.1:(10*24),C(0:.1:(10*24)))
            % Cr = squareform(C(hours(median(diff(spttime))).*pdist(id_tgr)));
            % Cr(eye(size(Cr,1))==1)=1;
            % imagesc(Cr)

            % Assess the match
            et=e(:,:,:);
            if sta{lt}.status(i_s)=="equipment" || sta{lt}.status(i_s)=="retrieval"
                s=tblLog.std_pres_calib(lt);
            else
                s=tblLog.std_pres(lt);
            end
            if isnan(s)
                warning('s cannot be zero. Use default 1')
                s=1;
            end

            w = log(pres_n{lt}(i_s))-1;
            % pres_rmse{lt}(:,:,i_s) = w .* mean((et/s).^2,3,'omitnan');
            pres_rmse{lt}(:,:,i_s) = w .* mean((et/s).^2,3,'omitnan');
            
            % RMSE to prob
            pres_prob{lt}(:,:,i_s) = exp(-pres_rmse{lt}(:,:,i_s));
            
            % figure; imagesc(pres_prob{lt}(:,:,i_s))
       
        else
            warning(['no pressure available for this stationary period: ' num2str(i_s) '-' tblLog.GDL_ID{lt}])
            pres_prob{lt}(:,:,i_s)=1;
            pres_thr{lt}(:,:,i_s)=true;
        end

    end
end


% RMSE to prob
% s=1;
% 
% for lt=1:height(tblLog)
%     pres_prob{lt} = exp(-pres_rmse{lt}/s);
% end

toc

%% 
clear splt T2 err prese_rmse dh e x x_r sp* pres_ge id* pres_rmse pres_gr* y tmp file* i* s w dt et ans
% clear gE % need to remove it to save
save('../data/processedData'+ project +'.mat')