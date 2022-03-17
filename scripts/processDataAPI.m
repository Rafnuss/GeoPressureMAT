
% Setting-up
tic

% project="StudyKenya";
% project="StudyPressure";
% project="StudyTraquet";
project="StudyStarling";


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
tblLog = tblLog(tblLog.multisensor,:);
tblLog = tblLog(tblLog.(project),:);
tblLog.Color(cellfun(@isempty,tblLog.Color))={'000000'};
tblLog.Color = hex2rgb(tblLog.Color);

%tblLog=tblLog(1,:);

%% Load DEM and ERA5 data
dll = 0.25;

% Create grid
lat=cell(height(tblLog),1);
lon=cell(height(tblLog),1);
for lt=1:height(tblLog)
    lon{lt} = floor(tblLog.bndy_W(lt)/dll)*dll:dll:ceil(tblLog.bndy_E(lt)/dll)*dll;
    lat{lt} = floor(tblLog.bndy_S(lt)/dll)*dll:dll:ceil(tblLog.bndy_N(lt)/dll)*dll;
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
    % delete stationary period below 7 hours
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
    % assert(any(hours(sta{lt}.actEffort(1:end-1))>.5))


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

mask_water = cell(height(tblLog),1);
for lt=1:height(tblLog)
    tmp = find(isnan(coastlat));
    mask_water_g = true(numel(lat{lt}),numel(lon{lt}));
    [gLON,gLAT] = meshgrid(lon{lt},lat{lt});
    for id=1:numel(tmp)-1
        mask_water_g(inpolygon(gLAT,gLON,coastlat(tmp(id):tmp(id+1)),coastlon(tmp(id):tmp(id+1))))=false;
    end
    mask_water{lt}=mask_water_g;
end

% light_prob{lt}(repmat(mask_water,1,1,height(sta{lt})))=0;


%% Pressure Map

pres_thr = cell(height(tblLog),1);
pres_prob = cell(height(tblLog),1);
pres_n = cell(height(tblLog),1);

for lt=1:height(tblLog)
    sta{lt}.s(:) = tblLog.std_pres(lt);
    sta{lt}.s(sta{lt}.status=="equipment" | sta{lt}.status=="retrieval") = tblLog.std_pres_calib(lt);
    tic
    [pres_prob{lt}, pres_thr{lt}, pres_n{lt}] = getPressueMap(raw{lt}.pressure,sta{lt},lon{lt},lat{lt});
    toc
end


%%
ts = cell(height(tblLog),1);
ts_calib = cell(height(tblLog),1);

% Add best location
for lt=1:height(tblLog)
    map =pres_prob{lt}.*pres_thr{lt};
    [~,id]=max(map,[],[1 2],'linear');
    [id_lat,id_lon,~]=ind2sub(size(map),squeeze(id));

    ts{lt} = getPressueTimeseries(sta{lt}, lon{lt}(id_lon), lat{lt}(id_lat), raw{lt}.pressure);

    sta_calib = sta{lt}(sta{lt}.status=="equipment" | sta{lt}.status=="retrieval",:);
    ts_calib{lt} = getPressueTimeseries(sta_calib, repmat(raw{lt}.calib.lon,1,height(sta_calib)), repmat(raw{lt}.calib.lat,1,height(sta_calib)), raw{lt}.pressure);
end



%% save
clear splt T2 err prese_rmse dh e x x_r sp* pres_ge id* pres_rmse pres_gr* y tmp file* i* s w dt et ans map tmp*
save('../data/processedData'+ project +'.mat')
















%% figure timeseries


lt=1;

col = [42,71,94;126 71 149;38 38 38;5 102 47;108 49 14]/255;

figure; hold on;
plot(raw{lt}.pressure.date,raw{lt}.pressure.obsWithOutliars,'color',[.4 .4 .4]);  hold on
data = pressureProcessing(raw{lt}.pressure,sta{lt});
for i_s = 1:height(sta{lt})
    plot(datetime(data.date(data.label==i_s),'ConvertFrom','posixtime'),  data.obs(data.label==i_s)/100,'LineWidth',2,'color',col(mod(i_s+1,2)+1,:))

    ids = ts{lt}.label==i_s;
    plot(ts{lt}.time(ids), ts{lt}.pressure(ids)-mean(ts{lt}.pressure(ids))+mean(data.obs(data.label==i_s)/100),'r')
    
    ids = ts_calib{lt}.label==i_s;
    plot(ts_calib{lt}.time(ids), ts_calib{lt}.pressure(ids)-mean(ts_calib{lt}.pressure(ids))+mean(data.obs(data.label==i_s)/100),'--r')

end
ylabel({raw{lt}.GDL_ID ,'Pressure(hPa)'})
grid on; box on; axis tight;

%%
sta{lt}(minutes(sta{lt}.actDuration)<30,:)

%% Map

sta_sm=cell(1,height(tblLog));
for lt=1:height(tblLog)
    grp_id = hours(sta{lt}.end-sta{lt}.start)>0;%sta{lt}.twlNb>=4;
    grp_id(1) = true;
    if ~isnat(tblLog.CalibSecondStart(lt))
        grp_id(end) = true;
    end
    sta_sm{lt} = sta{lt}(grp_id,:);
    sta_sm{lt}.actNb =  splitapply(@sum, sta{lt}.actNb,cumsum(grp_id));
    sta_sm{lt}.actEffort =  splitapply(@sum, sta{lt}.actEffort,cumsum(grp_id));
    sta_sm{lt}.actDuration =  splitapply(@sum, sta{lt}.actDuration,cumsum(grp_id));
    sta_sm{lt}.twlNbStopover =  splitapply(@sum, sta{lt}.twlNb,cumsum(grp_id))-sta_sm{lt}.twlNb;
    sta_sm{lt}.staID = find(grp_id);
end


% subp_row=2*ones(1,height(tblLog));
% subp_row([11 12])=3;
% subp_row([15 16])=1;
% for lt=1:height(tblLog)
 
    figure('position',[0 0 1200 750], 'Name', [raw{lt}.GDL_ID ' | ' tblLog.CommonName{lt}] );
    tiledlayout('flow','TileSpacing','tight','Padding','tight')
    
    % mvt_pdf = movementModel('energy',tblLog.CommonName{1});
    mvt_pdf = movementModel('gam');
    
    [gLON,gLAT] = meshgrid(lon{lt},lat{lt});
    for i_s = 1:height(sta_sm{lt})
        nexttile; hold on;
        
        figure('position',[0 0 1600 900]); hold on ; xticks([]); yticks([])
        
        %set(gca,'Color','k')
        tmp_p = pres_prob{lt}(:,:,sta_sm{lt}.staID(i_s));
        tmp_pt = pres_thr{lt}(:,:,sta_sm{lt}.staID(i_s));
        tmp_l = light_prob{lt}(:,:,sta_sm{lt}.staID(i_s));
        
        
        % Option 1
        mask = 0.3+0.7*double(tmp_pt);
        mask(mask_water{lt}) = 0;
        img_tmp = real2rgb(tmp_p,colormap);
        img_tmp = img_tmp.*mask;
        
        imagesc(lon{lt},lat{lt},img_tmp)
        borders('countries','w')
        % imagesc(lon{lt},lat{lt},tmp_p,'AlphaData',mask); 
        c_axis=caxis();
        contour(lon{lt},lat{lt},tmp_l,3,'color',[255, 252, 49]/255,'linewidth',2);
        

        
        tt = datestr(sta_sm{lt}.start(i_s),'dd-mmm HH:MM') + " to " + datestr(sta_sm{lt}.end(i_s),'dd-mmm HH:MM');
        tt = tt+ " (pres.: "+ num2str(pres_n{lt}(i_s)) + " hrs, light: "+ num2str(sta_sm{lt}.twlNb(i_s) + " twls)");
        
        if i_s>1
            tt = tt+ " | flight: "+num2str(round(hours(sta_sm{lt}.actEffort(max(1,i_s-1))))) + " hrs"  ;
            % tt = {tt{:} num2str(sta_sm{lt}.actNb(max(1,i_s-1))) + "act";
            if exist('gr','var')
                tmpd = reshape(lldistkm([gr{lt}.sp.lat(i_s-1) gr{lt}.sp.lon(i_s-1)],[gLAT(:) gLON(:)]),size(gLAT))./hours(sta_sm{lt}.actEffort(i_s-1));
                plot(gr{lt}.sp.lon(i_s-1), gr{lt}.sp.lat(i_s-1),'.w','linewidth',2,'MarkerSize',40)
                plot(gr{lt}.sp.lon(i_s-1), gr{lt}.sp.lat(i_s-1),'.g','linewidth',2,'MarkerSize',30)
                plot(gr{lt}.sp.lon(i_s), gr{lt}.sp.lat(i_s),'or','linewidth',2,'MarkerSize',12)
            else
                tmpd = reshape(lldistkm([gLAT(id) gLON(id)],[gLAT(:) gLON(:)]),size(gLAT))./hours(sta_sm{lt}.actEffort(i_s-1));
                plot(gLON(id), gLAT(id),'.w','linewidth',2,'MarkerSize',40)
                plot(gLON(id), gLAT(id),'.g','linewidth',2,'MarkerSize',30)
                [~,id]=max(tmp_p(:).*tmp_pt(:).*tmp_l(:).*mvt_pdf(tmpd(:)));
                plot(gLON(id), gLAT(id),'or','linewidth',2,'MarkerSize',12)
            end
            contour(lon{lt},lat{lt},tmpd,[7 17]/1000*60*60,'-g','linewidth',2); % ,'ShowText','on'
            
        else
            [~,id]=max(tmp_p(:).*tmp_pt(:).*tmp_l(:));
            plot(gLON(id), gLAT(id),'.w','MarkerSize',30)
        end

        % t=text(.01,.8,tt,'Units','normalized','FontSize',40,'FontWeight','bold','Color','w','VerticalAlignment','top');

        title(tt)

        caxis(c_axis); xticklabels('');yticklabels('');
        
        plot(raw{lt}.calib.lon,raw{lt}.calib.lat,'.w','MarkerSize',40)
        plot(raw{lt}.calib.lon,raw{lt}.calib.lat,'.r','MarkerSize',30)
        
        axis equal; axis([min(lon{lt}) max(lon{lt}) min(lat{lt}) max(lat{lt}) ]);
        
%         x0=[raw{lt}.calib.lon-4 raw{lt}.calib.lat-4];
%         dx = 1./[lldistkm(x0,x0+[1 0]) lldistkm(x0,x0+[0 1])]*100;
%         h=plot([x0(1) x0(1)], [x0(2) x0(2)+dx(2)], '-r', 'LineWidth', 2);%label(h,'100km','location','middle','slope')
%         h=plot([x0(1) x0(1)+dx(1)], [x0(2) x0(2)],'-r', 'LineWidth', 2); label(h,'100km')
%         
        exportgraphics(gcf,['combined_map_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '_' num2str(i_s) '.png'],'Resolution',300)
        %keyboard
        close all
        
    end
     %keyboard
    % exportgraphics(gcf,['combined_map_48h_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '.png'],'Resolution',300)
% end


%% Altitude
for lt=1:height(tblLog)
    dt = 1./hours(diff(raw{lt}.pressure.date(1:2)));

    figure('position',[0 0 1200 675], 'Name', [tblLog.CommonName{lt} ' (' raw{lt}.GDL_ID ')' ]);
    hold on;


    for i_s = 1:height(sta{lt})
        % plot(datetime(data.date(data.label==i_s),'ConvertFrom','posixtime'),  data.obs(data.label==i_s)/100,'LineWidth',2,'color',col(mod(i_s+1,2)+1,:))
    
        ids = ts{lt}.label==i_s;
        plot(ts{lt}.time(ids), ts{lt}.altitude(ids),'color',col(mod(i_s+1,2)+1,:))
        
        ids = ts_calib{lt}.label==i_s;
        plot(ts_calib{lt}.time(ids), ts_calib{lt}.altitude(ids),'--r')
    
    end

    % p3=yline(tblLog.DEMAttached(lt),'--r'); hold on;
    ylabel({raw{lt}.GDL_ID ,'Altitude (m asl)'})
    grid on; box on; axis tight;
    %legend([p1,p2,p3],'Alitude of equipement', 'Altitude according to the best matched location','Altitude according to the equipement location','Position', [.51 .97 0 0],'orientation','horizontal')

    % exportgraphics(gcf,['altitude_ts_' tblLog.CommonName{lt} '_' raw{lt-1}.GDL_ID '_' raw{lt}.GDL_ID '.png'],'Resolution',300)

end