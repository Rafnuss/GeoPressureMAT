file="../data/DEM/geopotential_land.nc";

glon=double(ncread(file,'longitude'));
glat=flip(double(ncread(file,'latitude')));

geopotential = permute(flip(ncread(file,'z'),2),[2 1]); %ncdisp(file);

% The geopotential can be converted to elevation with: 
% http://stcorp.github.io/harp/doc/html/algorithms/derivations/geopotential_height.html
% http://wind101.net/geopotential-altitude/
%https://www.springer.com/gp/book/9783211335444
% WMO-defined gravity constant
% g0 = 9.80665;
% Surface geopotential to surface geopotential height
% geopotentialheight = geopotential/g0;
% The gravity constant and earth radius are adjusted based on latittude: nominal gravity at sea level [m/s]
g = 9.7803253359 * (1 + 0.00193185265241*sind(glat).^2) ./ sqrt(1-0.00669437999013*sind(glat).^2);
% local earth curvature radius [m]
R = sqrt( ( (6378137^2 * cosd(glat)).^2 + (6356752^2 * sind(glat)).^2 ) ./ ( (6378137 * cosd(glat)).^2 + (6356752 * sind(glat)).^2 ));
% surface altitude (z)
% geopotentialheight = g/g0 * geoDEM.*R ./ (R+geoDEM);
% geopotential = g * geoDEM.*R ./ (R+geoDEM);
% geoDEM = geopotential.*R ./ (g.*R-geopotential);
% Compute the other way
% geopotential = g .* R.*z/(R+z);

%% Illustration
% figure('position',[0 0 1000 500]);
% subplot(1,3,1); hold on; title('ERA5 geopotential')
% imagesc(geolon,geolat,geoDEM)
% scatter(raw.calib.lon, raw.calib.lat,'or','filled')
% demcmap([-1 3000]); axis tight equal; colorbar;
% borders('countries','k')
% axis equal; axis([min(dem_30_lon) max(dem_30_lon) min(dem_30_lat) max(dem_30_lat) ]);
% subplot(1,3,2); hold on; title('DEM updascaled')
% imagesc(dem_30_lon,dem_30_lat,DEM_30)
% scatter(raw.calib.lon, raw.calib.lat,'or','filled')
% demcmap([-1 3000]); axis tight equal; colorbar;
% borders('countries','k')
% axis equal; axis([min(dem_30_lon) max(dem_30_lon) min(dem_30_lat) max(dem_30_lat) ]);
% subplot(1,3,3); hold on; title('difference')
% imagesc(dem_30_lon,dem_30_lat,geoDEM-DEM_30)
% axis tight equal; colorbar;
% borders('countries','k')
% axis equal; axis([min(dem_30_lon) max(dem_30_lon) min(dem_30_lat) max(dem_30_lat)]);
Export
% a=[geoDEM(:,end/2:end) geoDEM(:,1:end/2)];
% b=[glon(end/2:end)-360; glon(1:end/2)];
% R2 = georefpostings([min(glat) max(glat)],[min(b) max(b)],size(a));
% geotiffwrite("../data/DEM/geoDEM",a,R2)

High resoluation DEM STRM90
files = {'export_DEM-0000000000-0000000000','export_DEM-0000000000-0000056320';'export_DEM-0000076800-0000000000','export_DEM-0000076800-0000056320'};
DEM=zeros(106373,76687,'int16'); R=cell(size(files));
[DEM(1:76800,1:56320),R{1,1}] = readgeoraster(['../data/DEM/' files{1,1} '.tif']); 
[DEM(1:76800,56321:end),R{1,2}] = readgeoraster(['../data/DEM/' files{1,2} '.tif']); 
[DEM(76801:end,1:56320),R{2,1}] = readgeoraster(['../data/DEM/' files{2,1} '.tif']); 
[DEM(76801:end,56321:end),R{2,2}] = readgeoraster(['../data/DEM/' files{2,2} '.tif']); 

dem_lat = (R{2,2}.LatitudeLimits(1)+R{1,1}.CellExtentInLatitude/2):R{1,1}.CellExtentInLatitude:(R{1,1}.LatitudeLimits(2)-R{1,1}.CellExtentInLatitude/2);
dem_lon = (R{1,1}.LongitudeLimits(1)+R{1,1}.CellExtentInLongitude/2):R{1,1}.CellExtentInLongitude:(R{2,2}.LongitudeLimits(2)-R{1,1}.CellExtentInLongitude/2);
Compute max and min elevation
Method 1: Blockproc. Working but not accurate as the grid is not perfectly regular.
%     dres = round(0.25/R{1,1}.CellExtentInLongitude);
%     npadlat = (dres*numel(glat)-numel(dem_lat));
%     npadlon = (dres*numel(glon)-numel(dem_lon))/2;
%     DEM_max =  blockproc(padarray(DEM,[npadlat npadlon],nan),[dres dres],@(x) nanmax(x.data(:)));
%     DEM_min =  blockproc(padarray(DEM,[npadlat npadlon],nan),[dres dres],@(x) nanmin(x.data(:)));
%     DEM_max = double(flipud(DEM_max));
%     DEM_min = double(flipud(DEM_min));
Method 2: looping. More time possibly but more accurate. <5 min.
DEM_max = nan(numel(glat),numel(glon)); DEM_min = DEM_max;
for i_lat=1:numel(glat)
    id_lat = (glat(i_lat)-.25/2) <= dem_lat & dem_lat <= (glat(i_lat)+.25/2);
    for i_lon=1:numel(glon)
        id_lon = (glon(i_lon)-.25/2) <= dem_lon & dem_lon <= (glon(i_lon)+.25/2);
        DEM_max(i_lat,i_lon) =  max(DEM(id_lat,id_lon),[],'all');
        DEM_min(i_lat,i_lon) =  min(DEM(id_lat,id_lon),[],'all'); 
        if DEM_max>0 & DEM_max==DEM_min
            keyboard
        end
    end
    i_lat
end
DEM_max = double(flipud(DEM_max));
DEM_min = double(flipud(DEM_min));
Compute the DEM at the calibration site
for lt=1:height(tblLog)
    [~,id_lat] = min(abs(tblLog.LatitudeAttached(lt)-fliplr(dem_lat)));
    [~,id_lon] = min(abs(tblLog.LongitudeAttached(lt)-dem_lon));
    tblLog.DEMAttached(lt) = DEM(id_lat,id_lon);
end
Save and Load variable
save('../data/DEM/MatlabData.mat','DEM_max','DEM_min','glon','glat','geoDEM');