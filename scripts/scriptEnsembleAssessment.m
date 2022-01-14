%% Assess ensemble spread for two species

addpath(genpath('../functions'))
load("../data/processedDataStudyPressure.mat")
scriptAltPres()

%% Choose track
ltt=[1 7];%1 7

%% Write and make requested

tblLog.CalibFirstStart(lt)
tblLog.CalibFirstEnd(lt)
tblLog.LatitudeAttached(lt)+[-10 10]
tblLog.LongitudeAttached(lt)+[-10 10]

%% 
for i_lt=1:numel(ltt)
    lt=ltt(i_lt);
    file="../data/ECMWF/surface_pressure_"+raw{lt}.GDL_ID+".nc"; 
    spttime_ens{i_lt} = datetime(double(ncread(file,'time'))/24 + datenum('1900-01-01 00:00:00'),'convertFrom','datenum');
    splt_ens = permute(flip(ncread(file,'sp'),2)/100,[2 1 4 3]);
    lat_ens = flip(ncread(file,'latitude'),2);
    lon_ens = ncread(file,'longitude');

    % imagesc(lon_ens,lat_ens,splt_ens(:,:,1,1))

    % Calibration site DEM
    [~,id_lat_calib_ens]=min(abs(raw{lt}.calib.lat-lat_ens));
    [~,id_lon_calib_ens]=min(abs(raw{lt}.calib.lon-lon_ens));

    sp_ens_lt_presCalib{i_lt} = reshape(ncread(file,'sp',[id_lon_calib_ens numel(lat_ens)-id_lat_calib_ens+1 1 1],[1 1 10 numel(spttime_ens{i_lt})])/100,10,[]);
end

%%
figure; tiledlayout('flow')
for i_lt=1:numel(ltt)
    lt=ltt(i_lt); nexttile; hold on;
    plot(spttime_ens{i_lt}, sp_ens_lt_presCalib{i_lt}')
    plot(sp{lt}{1}.time,sp{lt}{1}.presCalib,'k','LineWidth',2)
    plot(raw{lt}.pressure.date,raw{lt}.pressure.obsWithOutliars,'r','LineWidth',2);
    xlim([sp{lt}{1}.time(1) sp{lt}{1}.time(end)])
end

figure; tiledlayout('flow')
for i_lt=1:numel(ltt)
    lt=ltt(i_lt); nexttile; hold on;
    n = sp{lt}{1}.presCalib-mean(sp{lt}{1}.presCalib);
    t=sp{lt}{1}.time;
    tmp=interp1(raw{lt}.pressure.date,raw{lt}.pressure.obsWithOutliars,t);
    plot(t, tmp'-n-mean(tmp),'r','LineWidth',2);
    tmp=interp1(spttime_ens{i_lt}, sp_ens_lt_presCalib{i_lt}',t);
    plot(t, tmp'-n-mean(tmp',2))
    xlim([sp{lt}{1}.time(1) sp{lt}{1}.time(end)])
end

%%

splt = permute(flip(ncread(file,'sp',[id_lon(1) numel(glat)-id_lat(end)+1 1 id_t(1)],[numel(id_lon) id_lat(end)-id_lat(1)+1 1 id_t(end)-id_t(1)+1]),2)/100,[2 1 4 3]);
    