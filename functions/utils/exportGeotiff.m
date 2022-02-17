function exportGeotiff(lat,lon,M,raw,tblLog,sta,alt,sp)

assert(size(M,3)==height(sta))

% define the R
R = georefpostings([min(lat) max(lat)],[min(lon) max(lon)],size(M,[1 2]));

% define the manifest
manifest0={};
%manifest0.tilesets(1).id="id";
%manifest0.tilesets(1).crs="crs";
manifest0.bands(1).id="pressure_threshold";
manifest0.bands(1).tileset_band_index=0;
manifest0.bands(2).id="pressure_mismatch";
manifest0.bands(2).tileset_band_index=1;
manifest0.bands(3).id="light";
manifest0.bands(3).tileset_band_index=2;
manifest0.bands(4).id="graph";
manifest0.bands(4).tileset_band_index=3;
manifest0.pyramiding_policy="MEAN";

for i_s=1:size(M,3)
    
    folder="../data/export/geotiff/"+raw.GDL_ID+"/"+i_s;
    mkdir(folder)
    
    Ms = squeeze(M(:,:,i_s,:));
    geotiffwrite(folder+"/"+raw.GDL_ID+"_"+i_s,Ms,R)
    
    manifest = manifest0;
    manifest.tilesets={struct("sources",{{struct("uris",{{raw.GDL_ID+"_"+i_s+".tif"}})}})};
    manifest.name= "projects/earthengine-legacy/assets/projects/earthimages4unil/PostDocProjects/rafnuss/Geolocator/"+raw.GDL_ID+"_"+i_s;
    manifest.start_time = datestr(sta.start(i_s),'yyyy-mm-ddTHH:MM:SSZ');
    manifest.end_time = datestr(sta.end(i_s),'yyyy-mm-ddTHH:MM:SSZ');

    manifest.properties.actDuration = minutes(sta.actDuration(i_s));
    manifest.properties.actEffort = minutes(sta.actEffort(i_s));
    manifest.properties.twlNb = sta.twlNb(i_s);
    manifest.properties.staID = sta.staID(i_s);
    manifest.properties.CommonName = tblLog.CommonName{1};
    manifest.properties.status = sta.status(i_s);
    manifest.properties.GDL_ID = raw.GDL_ID;
    manifest.properties.calib_lon = raw.calib.lon;
    manifest.properties.calib_lat = raw.calib.lat;
    manifest.properties.alt_mean = nanmean(alt{i_s}.alt);
    manifest.properties.alt_min = nanmin(alt{i_s}.alt);    
    manifest.properties.alt_max = nanmax(alt{i_s}.alt);
    manifest.properties.alt_q05 = quantile(alt{i_s}.alt,.05);  
    manifest.properties.alt_q95 = quantile(alt{i_s}.alt,.95);  
    manifest.properties.alt_max = nanmax(alt{i_s}.alt);
    manifest.properties.sp_lon = sp(i_s,1);
    manifest.properties.sp_lat = sp(i_s,2);
    
    fileID = fopen(folder+"/manifest.json","w");
    fprintf(fileID,jsonencode(manifest));
    fclose(fileID);
end

end