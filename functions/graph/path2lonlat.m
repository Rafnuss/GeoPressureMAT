function [lon, lat, i_s] = path2lonlat(path,gr)

[lat_id,lon_id,i_s] = ind2sub(gr.snds,path);
lon = gr.lon(lon_id);
lat = gr.lat(lat_id);
end