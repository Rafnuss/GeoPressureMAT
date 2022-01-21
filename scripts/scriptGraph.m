
%startup
addpath(genpath('../functions'))
addpath('~/Documents/GitHub/Flight-Matlab/functions/')
addpath('~/Documents/GitHub/Flight-Matlab/data/')
load("../data/processedDataStudyPressure.mat")

% filter s
% {'18IC','18LX','22BK','22BN','22KT','24FF','24TA','24UL','16LP','20IK','22QL','22QO','20OA','20OE','16AQ','16DM'}
skip_gdl = {'22KT','24FF', '20IK'};
% raw
% tblLog
% sta
% lat
% light_prob
% lon


%% Build the graph

gr=cell(height(tblLog),1);
% lt=find(tblLog.GDL_ID == "22QO"); % 22QL

thr_prob_percentile = .99;
thr_gs = 150;
thr_as_percentile = .95;

for lt=1:height(tblLog)
    if any(strcmp(skip_gdl, tblLog.GDL_ID{lt})), continue; end
    tic
    disp(['Sarting: ' tblLog.GDL_ID{lt} ' ' num2str(lt)])
    prob_map = pres_prob{lt}(:,:,sta{lt}.staID) .* pres_thr{lt}(:,:,sta{lt}.staID) .* light_prob{lt}(:,:,sta{lt}.staID) .* ~mask_water{lt};
    % prob_map = pres_prob{lt}(:,:,sta{lt}.staID) .* pres_thr{lt}(:,:,sta{lt}.staID) .* ~mask_water{lt};
    [grt,nds] = createGraph(prob_map,lat{lt},lon{lt},raw{lt}.calib,sta{lt}.actEffort,thr_prob_percentile,thr_gs);
    t=toc; disp(['Creating graph in ' num2str(t,3) ' sec'])
    grt = filterGraph(grt,'gs',thr_gs);
    t=toc; disp(['Filter groundspeed at ' num2str(thr_gs) 'km/h in ' num2str(t,3) ' sec'])
    grt = windSpeedGraph(grt,raw{lt},sta{lt},sta{lt},activityMig{lt});
    t=toc; disp(['Adding windspeed ' num2str(t,3) ' sec'])
    % grt.mvt_pdf = movementModel('step');
    grt.mvt_pdf = movementModel('energy',tblLog.CommonName{lt});
    % grt.thr_as = find(cumsum(grt.mvt_pdf(1:1000))./sum(grt.mvt_pdf(1:1000))>thr_as_percentile,1);
    grt.thr_as = 100;
    grt = filterGraph(grt,'as',grt.thr_as);
    t=toc; disp(['Filter airspeed at ' num2str(grt.thr_as) 'km/h in ' num2str(t,3) ' sec'])

    grt.p = grt.ps .* grt.mvt_pdf(abs(grt.as));
    grt.M = probMapGraph(grt);
    grt.sp = shortestPathGraph(grt);
    gr{lt} = grt;
    t=toc; disp(['Finished in ' num2str(t,3) ' sec'])
    disp('----')
end


% Simulation path with gibbs
nj=1000;

for lt=1:height(tblLog)
    if any(strcmp(skip_gdl, tblLog.GDL_ID{lt})), continue; end
    disp(raw{lt}.GDL_ID)
    
    % Define intial path as the shortest path
    path0=gr{lt}.sp.path;
    % Define the fixed node/sta for first and possibly last
    fixPath = false(size(path0));
    fixPath(1)=true;
    if ~isnat(tblLog.CalibSecondStart(lt))
        fixPath(end)=true;
    end
    tic
    gr{lt}.psim = gibbsGraph(nj,path0,fixPath,gr{lt});
    toc
end



%% Save
save('../data/graph'+project+'.mat','gr','-v7.3')

%% Export to geotiff
scriptAltPres()

for lt=1:height(tblLog)
    if any(strcmp(skip_gdl, tblLog.GDL_ID{lt})), continue; end
    M = cat(4, pres_thr{lt}(:,:,sta{lt}.staID),pres_prob{lt}(:,:,sta{lt}.staID), light_prob{lt}(:,:,sta{lt}.staID), gr{lt}.M);
    exportGeotiff(lat{lt},lon{lt},M,raw{lt},tblLog(lt,:),sta{lt},alt{lt},[gr{lt}.sp.lon gr{lt}.sp.lat])
end












