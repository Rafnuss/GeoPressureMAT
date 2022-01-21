function gr = windSpeedGraph(gr,raw,sta,sta_sm,activityMig)


[Slat,Slon,St]=ind2sub(gr.snds,gr.s);
[Tlat,Tlon,~]=ind2sub(gr.snds,gr.t);


% Initialized zero windspeed for all edges.
WS = nan(height(gr.s),2);


for i_sm = 1:height(sta_sm)-1

    % Edge i_s
    st_id = find(St==i_sm);

    i_sm_staID = sta_sm.staID(i_sm):(sta_sm.staID(i_sm+1)-1);

    assert(sum(sta.actDuration(i_sm_staID))==sta_sm.actDuration(i_sm),'issue with duration')
    ratio_sta = cumsum(sta.actDuration(i_sm_staID)')./sta_sm.actDuration(i_sm);
    ratio_sta = [0 ratio_sta];
    
    lat_sta = gr.lat(Slat(st_id)) + ratio_sta .* (gr.lat(Tlat(st_id))-gr.lat(Slat(st_id)));
    lon_sta = gr.lon(Slon(st_id)) + ratio_sta .* (gr.lon(Tlon(st_id))-gr.lon(Slon(st_id)));

    Usta=nan(numel(st_id),numel(i_sm_staID)); Vsta=Usta;

    for i_s = 1:numel(i_sm_staID)

        i_s_staID = i_sm_staID(i_s);

        filename = "../data/ECMWF/" + raw.GDL_ID + "/" + "sta_" + num2str(i_s_staID) + ".nc";

        wpres=double(ncread(filename,'level'));
        wlat=double(ncread(filename,'latitude'));
        wlon=double(ncread(filename,'longitude'));
        wtime = datetime(double(ncread(filename,'time'))/24 + datenum('1900-01-01 00:00:00'),'convertFrom','datenum');

        wtime_regular_id = wtime>=(dateshift(sta.end(i_s_staID),'start','hour')-1/24) & wtime <= (dateshift(sta.start(i_s_staID+1),'end','hour')+1/24);

        u = flip(ncread(filename,'u'),2)/1000*60*60;
        Fu = griddedInterpolant({wlon flipud(wlat) wpres datenum(wtime(wtime_regular_id))},u(:,:,:,wtime_regular_id),'linear','nearest');
        v = flip(ncread(filename,'v'),2)/1000*60*60;
        Fv = griddedInterpolant({wlon flipud(wlat) wpres datenum(wtime(wtime_regular_id))},v(:,:,:,wtime_regular_id),'linear','nearest');

        % mean(sqrt(u.^2+v.^2),'all')

        % Find activity
        id_activityMig = find(activityMig.staID==i_s_staID);
        assert(sum(activityMig.duration(id_activityMig))==sta.actDuration(i_s_staID),'issue with duration')
        ratio_act = cumsum(activityMig.duration(id_activityMig)')./sta.actDuration(i_s_staID);
        ratio_act = [0 ratio_act];

        lat_act = lat_sta(:,i_s) + ratio_act .* (lat_sta(:,i_s+1)-lat_sta(:,i_s));
        lon_act = lon_sta(:,i_s) + ratio_act .* (lon_sta(:,i_s+1)-lon_sta(:,i_s));

        Uact=nan(numel(st_id),numel(id_activityMig)); Vact=Uact;
        for i_m = 1:numel(id_activityMig)

            q_time = linspace(datenum(activityMig.date_min(id_activityMig(i_m))), datenum(activityMig.date_max(id_activityMig(i_m))),2+hours(activityMig.duration(id_activityMig(i_m))));
            q_pres = interp1(datenum(raw.pressure.date),raw.pressure.obsWithOutliars,q_time);
            
            % First approx of linear spatial spacing (constant groundspeed)
            ratio_hr=linspace(0,1,numel(q_time));
            lat_diff = (lat_act(:,i_m+1)-lat_act(:,i_m));
            lon_diff = (lon_act(:,i_m+1)-lon_act(:,i_m));
            % orientation = atan(lat_diff./lon_diff);
            q_lat = lat_act(:,i_m) + ratio_hr .* lat_diff;
            q_lon = lon_act(:,i_m) + ratio_hr .* lon_diff;
            Uact_appox = Fu(q_lon, q_lat, repmat(q_pres,numel(st_id),1), repmat(q_time,numel(st_id),1));
            Vact_appox = Fv(q_lon, q_lat, repmat(q_pres,numel(st_id),1), repmat(q_time,numel(st_id),1));
            
            % Correct the spacing for constant airspeed 
            % first, compute the average airspeed for each eadge
            aslon = real(gr.gs(st_id)) - mean(Uact_appox,2);
            % now using this average (and assumed constant airspeed, we can
            % compute the ratio of the position of the bird along its 
            % trajectory assuming this airspeed and the winspeed discretize 
            w = movmean(Uact_appox,[0 1],2,'endpoints','discard')+aslon; 
            % we can then estimate the ratio of distance covered by
            % normalizing w to 1 and computing the cumulative. 
            ratio_hr_corr_lon = [zeros(height(w),1) cumsum(w./sum(w,2),2)];
            % one issue appears when lon_diff is zero. the ratio fof
            % distance is getting to infinitity. In this case the ratio
            % don't really matter, it will give the same lat lon
            ratio_hr_corr_lon(lon_diff==0,:)=1;
            % for diff_lon close to zero, we get the same effect, we can
            % simply block them to 0 - 1.
            ratio_hr_corr_lon(ratio_hr_corr_lon<0)=0;
            ratio_hr_corr_lon(ratio_hr_corr_lon>1)=1;
            % Another issue appear if gs=1. Then sum(w) is zero. We fix it
            % by forcing all ratio to 1. 
            ratio_hr_corr_lon(isnan(ratio_hr_corr_lon))=1;
            % do the same for lat
            aslat = imag(gr.gs(st_id)) - mean(Vact_appox,2);
            w = movmean(Vact_appox,[0 1],2,'endpoints','discard')+aslat; 
            ratio_hr_corr_lat = [zeros(height(w),1) cumsum(w./sum(w,2),2)];
            ratio_hr_corr_lat(ratio_hr_corr_lat<0)=0;
            ratio_hr_corr_lat(ratio_hr_corr_lat>1)=1;
            ratio_hr_corr_lat(isnan(ratio_hr_corr_lat))=1;

            q_lat = lat_act(:,i_m) + ratio_hr_corr_lat .* (lat_act(:,i_m+1)-lat_act(:,i_m));
            q_lon = lon_act(:,i_m) + ratio_hr_corr_lon .* (lon_act(:,i_m+1)-lon_act(:,i_m));
            Uact(:,i_m) = mean(Fu(q_lon, q_lat, repmat(q_pres,numel(st_id),1), repmat(q_time,numel(st_id),1)),2);
            Vact(:,i_m) = mean(Fv(q_lon, q_lat, repmat(q_pres,numel(st_id),1), repmat(q_time,numel(st_id),1)),2);

            assert(~any(isnan(Uact(:,i_m))))
            assert(~any(isnan(Vact(:,i_m))))
        end
        Usta(:,i_s) = sum(activityMig.duration(id_activityMig)'./sta.actEffort(i_s_staID) .* Uact,2);
        Vsta(:,i_s) = sum(activityMig.duration(id_activityMig)'./sta.actEffort(i_s_staID) .* Vact,2);
    end

    WS(st_id,1) = sum(sta.actEffort(i_sm_staID)'./sta_sm.actEffort(i_sm) .* Usta,2);
    WS(st_id,2) = sum(sta.actEffort(i_sm_staID)'./sta_sm.actEffort(i_sm) .* Vsta,2);
    % mean(sqrt(sum(WS(st_id,:).^2,2)))
end

gr.ws = sum(WS.*[1 1i],2);
gr.as=gr.gs-gr.ws;


end