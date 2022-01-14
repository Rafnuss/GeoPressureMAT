% Script compute pressure at most likely position and altitude

sp = cell(height(tblLog),1);
alt = cell(height(tblLog),1);
file='../data/ECMWF/surface_pressure.nc';
spttime = datetime(double(ncread(file,'time'))/24 + datenum('1900-01-01 00:00:00'),'convertFrom','datenum');
%glon=double(ncread(file,'longitude'));
%glat=flip(double(ncread(file,'latitude')));
calibGeoDEM=nan(height(tblLog),3);
for lt=1:height(tblLog)
    sp{lt}=cell(height(sta{lt}),1);
    alt{lt}=cell(height(sta{lt}),1);
    % Calibration site DEM
    [~,id_lat_calib]=min(abs(raw{lt}.calib.lat-glat));
    [~,id_lon_calib]=min(abs(raw{lt}.calib.lon-glon));
    calibGeoDEM(lt,1) = geoDEM(id_lat_calib,id_lon_calib);
    calibGeoDEM(lt,2) = DEM_min(id_lat_calib,id_lon_calib);
    calibGeoDEM(lt,3) = DEM_max(id_lat_calib,id_lon_calib);

    for i_s = 1:height(sta{lt})
        id_t = find(sta{lt}.start(i_s)<=spttime & spttime <= sta{lt}.end(i_s));
        % id_tge = sta{lt}.start(i_s)<raw{lt}.pressure.date & raw{lt}.pressure.date < sta{lt}.end(i_s);

        % Get pressure at best match
        map_prob = pres_thr{lt}(:,:,i_s).*pres_prob{lt}(:,:,i_s);
        [~,id_max] = max(map_prob(:)); [i_y,i_x]=ind2sub(size(map_prob),id_max);
        [~,id_lat_tmp]=min(abs(lat{lt}(i_y)-glat));
        [~,id_lon_tmp]=min(abs(lon{lt}(i_x)-glon));
        sp{lt}{i_s}.time = spttime(id_t);
        sp{lt}{i_s}.pres = reshape(ncread(file,'sp',[id_lon_tmp numel(glat)-id_lat_tmp+1 1 id_t(1)],[1 1 1 id_t(end)-id_t(1)+1])/100,1,[]);
        if raw{lt}.pressure.date(end)>datetime('1-jul-2021')
            splt_0 = reshape(ncread(file,'sp',[id_lon_tmp numel(glat)-id_lat_tmp+1 2 id_t(1)],[1 1 1 id_t(end)-id_t(1)+1])/100,1,[]);
            sp{lt}{i_s}.pres(spttime(id_t)>=datetime('1-jul-2021'))=splt_0(spttime(id_t)>=datetime('1-jul-2021'));
        end
        % Known location
        if sta{lt}.status(i_s)=="equipment" || sta{lt}.status(i_s)=="retrieval"
            sp{lt}{i_s}.presCalib = reshape(ncread(file,'sp',[id_lon_calib numel(glat)-id_lat_calib+1 1 id_t(1)],[1 1 1 id_t(end)-id_t(1)+1])/100,1,[]);
            if raw{lt}.pressure.date(end)>datetime('1-jul-2021')
                splt_0 = reshape(ncread(file,'sp',[id_lon_calib numel(glat)-id_lat_calib+1 2 id_t(1)],[1 1 1 id_t(end)-id_t(1)+1])/100,1,[]);
                sp{lt}{i_s}.presCalib(spttime(id_t)>=datetime('1-jul-2021'))=splt_0(spttime(id_t)>=datetime('1-jul-2021'));
            end
        end


        % Get altitudinal profile 
        alt{lt}{i_s}.dem = geoDEM(id_lat_tmp,id_lon_tmp);
        % get pressure for flight before and after also
        if i_s==1
            id_tm = find(sta{lt}.start(i_s)<=spttime & spttime <= sta{lt}.start(i_s+1));
            id_tgem = sta{lt}.start(i_s)<raw{lt}.pressure.date & raw{lt}.pressure.date < sta{lt}.start(i_s+1);
        elseif i_s==height(sta{lt})
            id_tm = find(sta{lt}.end(i_s-1)<=spttime & spttime <= sta{lt}.end(i_s));
            id_tgem = sta{lt}.end(i_s-1)<raw{lt}.pressure.date & raw{lt}.pressure.date < sta{lt}.end(i_s);
        else
            id_tm = find(sta{lt}.end(i_s-1)<=spttime & spttime <= sta{lt}.start(i_s+1));
            id_tgem = sta{lt}.end(i_s-1)<raw{lt}.pressure.date & raw{lt}.pressure.date < sta{lt}.start(i_s+1);
        end
        tmp_pres = reshape(ncread(file,'sp',[id_lon_tmp numel(glat)-id_lat_tmp+1 1 id_tm(1)],[1 1 1 id_tm(end)-id_tm(1)+1])/100,1,[]);
        if raw{lt}.pressure.date(end)>datetime('1-jul-2021')
            splt_0 = reshape(ncread(file,'sp',[id_lon_tmp numel(glat)-id_lat_tmp+1 2 id_tm(1)],[1 1 1 id_tm(end)-id_tm(1)+1])/100,1,[]);
            tmp_pres(spttime(id_tm)>=datetime('1-jul-2021'))=splt_0(spttime(id_tm)>=datetime('1-jul-2021'));
        end
        alt{lt}{i_s}.time = raw{lt}.pressure.date(id_tgem);
        
        % downscale to geolocator resolution
        if numel(tmp_pres)>1
            tmp_pres_interp=interp1(spttime(id_tm),tmp_pres',alt{lt}{i_s}.time);
        else
            tmp_pres_interp = tmp_pres;
        end

        % convert pressure to altitude based on the raw data
        alt{lt}{i_s}.alt = pressure2altitude(raw{lt}.pressure.obsWithOutliars(id_tgem), tmp_pres_interp)+alt{lt}{i_s}.dem;
        alt{lt}{i_s}.alt_basic = pressure2altitude(raw{lt}.pressure.obsWithOutliars(id_tgem));
        
        % Known location
        if sta{lt}.status(i_s)=="equipment" || sta{lt}.status(i_s)=="retrieval"
            tmp_pres = reshape(ncread(file,'sp',[id_lon_calib numel(glat)-id_lat_calib+1 1 id_tm(1)],[1 1 1 id_tm(end)-id_tm(1)+1])/100,1,[]);
            if raw{lt}.pressure.date(end)>datetime('1-jul-2021')
                splt_0 = reshape(ncread(file,'sp',[id_lon_tmp numel(glat)-id_lat_tmp+1 2 id_tm(1)],[1 1 1 id_tm(end)-id_tm(1)+1])/100,1,[]);
                tmp_pres(spttime(id_tm)>=datetime('1-jul-2021'))=splt_0(spttime(id_tm)>=datetime('1-jul-2021'));
            end
            tmp_pres_interp=interp1(spttime(id_tm),tmp_pres',alt{lt}{i_s}.time);
            alt{lt}{i_s}.altCalib = pressure2altitude(raw{lt}.pressure.obsWithOutliars(id_tgem), tmp_pres_interp)+calibGeoDEM(lt,1);
        end
    end
    %lt
end
