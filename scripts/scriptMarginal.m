
% Script margiinal compute the marginal prob dist of position at each
% stationary period + sthortest path from the prob matrix only.

% It is working but too slow compare to the graph. It does not implement
% windspeed.

%%
addpath(genpath('../functions'))
load("../data/processedDataStudyPressure.mat")

%%
lt=1;

prob_map = pres_prob{lt} .* pres_thr{lt} .* light_prob{lt} .* ~mask_water{lt};
mvt_pdf = movementModel('energy',tblLog.mass(lt),tblLog.wingSpan(lt));
thr_speed=150;
flight_duration = hours(sta{lt}.actEffort);
shortestPath = true;


tic
[map,path] = probMapMarginal(prob_map, mvt_pdf, lon{lt}, lat{lt}, thr_speed, flight_duration, shortestPath);
toc % 120sec without shortes path 148 with


%% Figure
figure; hold on

imagesc(lon{lt}, lat{lt}, map(:,:,6)); axis tight equal; a_xis= axis();
borders('countries','w')
[LON,LAT] = meshgrid(lon, lat);
plot(LON(path),LAT(path),'.-r');
axis(a_xis);