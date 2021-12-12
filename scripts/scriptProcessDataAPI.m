% Setting-up
tic
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
tblLog = tblLog(tblLog.StudyPressure,:);
%tblLog = tblLog(tblLog.StudyTraquet,:);
%tblLog = tblLog(tblLog.Kenya,:);
tblLog.Color(cellfun(@isempty,tblLog.Color))={'000000'};
tblLog.Color = hex2rgb(tblLog.Color);

tblLog=tblLog(1,:);

%% Define grid
% Use area-cell (). 

res = 0.1; % 0.1° = 10km
glon=(-19+res/2):res:(43-res/2);
glat=(-33+res/2):res:(53-res/2);

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
    xi = 0:.01:360;%floor(min(gE{lt}.z(~twl{lt}.isOutliar(id_calib)))):0.01:ceil(max(gE{lt}.z(~twl{lt}.isOutliar(id_calib))));
    [gE{lt}.f,gE{lt}.xi]=ksdensity(gE{lt}.z(~twl{lt}.isOutliar(gE{lt}.id_calib)),xi,'Bandwidth',bd);
    gE{lt}.pdf = @(x) interp1(gE{lt}.xi,gE{lt}.f,x,'nearest','extrap');
    
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
    [~,id]=maxk(sta{lt}.twlNb,3);
    if tblLog.GDL_ID{lt}=="24FF"
        iws=id(2);
    elseif tblLog.GDL_ID{lt}=="20OA"
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


%% Pressure Map New

pres_thr = cell(height(tblLog),1);
pres_prob = cell(height(tblLog),1);
pres_n = cell(height(tblLog),1);

for lt=1:height(tblLog)
    sta{lt}.s = repmat(tblLog.std_pres(lt),height(sta{lt}),1);
    sta{lt}.s(sta{lt}.status=="equipment" | sta{lt}.status=="retrieval") = tblLog.std_pres_calib(lt);
    tic
    [pres_prob{lt}, pres_thr{lt}, pres_n{lt}] = getPressueMap(raw{lt}.pressure,sta{lt},lon{lt},lat{lt});
    toc
end

%% Check
for lt=1:height(tblLog)
    map = pres_thr{lt} .* pres_prob{lt};
    out = getPressueTimeseries(map,sta{lt},lon{lt},lat{lt});
    
    figure; hold on;
    plot(raw{lt}.pressure.date,raw{lt}.pressure.obs,'-k')
    plot(raw{lt}.pressure.date,raw{lt}.pressure.obsWithOutliars,'-r','linewidth',2)
    plot(out.time,out.pressure)
end

%% 
clear splt T2 err prese_rmse dh e x x_r sp* pres_ge id* pres_rmse pres_gr* y tmp file* i* s w dt et ans
clear gE % need to remove it to save
save('../data/processedDataTraquet.mat')