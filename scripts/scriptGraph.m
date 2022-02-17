
%startup
addpath(genpath('../functions'))
addpath('~/Documents/GitHub/Flight-Matlab/functions/')
addpath('~/Documents/GitHub/Flight-Matlab/data/')
load("../data/processedDataStudyPressure.mat")

% filter s
% {'18IC','18LX','22BK','22BN','22KT','24FF','24TA','24UL','16LP','20IK','22QL','22QO','20OA','20OE','16AQ','16DM'}
skip_gdl = {'24FF'};
% raw
% tblLog
% sta
% lat
% light_prob
% lon

sta_sm=cell(1,height(tblLog));
for lt=1:height(tblLog)
    if any(strcmp(skip_gdl, tblLog.GDL_ID{lt}))
        grp_id = hours(sta{lt}.end-sta{lt}.start)>12;%sta{lt}.twlNb>=4;
    else
        grp_id =true(height(sta{lt}),1);
    end
    
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


%% Build the graph

gr=cell(height(tblLog),1);
% lt=find(tblLog.GDL_ID == "22QO"); % 22QL
trun = cell(height(tblLog),1);

thr_prob_percentile = .99;
thr_gs = 150;
thr_as_percentile = .95;

for lt=7:9%:height(tblLog)
    % if any(strcmp(skip_gdl, tblLog.GDL_ID{lt})), continue; end
    tic
    disp(['Sarting: ' tblLog.GDL_ID{lt} ' ' num2str(lt)])
    prob_map = pres_prob{lt}(:,:,sta_sm{lt}.staID) .* pres_thr{lt}(:,:,sta_sm{lt}.staID) .* light_prob{lt}(:,:,sta_sm{lt}.staID) .* ~mask_water{lt};
    % prob_map = pres_prob{lt}(:,:,sta_sm{lt}.staID) .* pres_thr{lt}(:,:,sta_sm{lt}.staID) .* ~mask_water{lt};
    [grt,nds] = createGraph(prob_map,lat{lt},lon{lt},raw{lt}.calib,sta_sm{lt}.actEffort,thr_prob_percentile,thr_gs);
    t1=toc; disp(['Creating graph in ' num2str(t1,3) ' sec'])
    grt = filterGraph(grt,'gs',thr_gs);
    t2=toc; disp(['Filter groundspeed at ' num2str(thr_gs) 'km/h in ' num2str(t2,3) ' sec'])
    grt = windSpeedGraph(grt,raw{lt},sta_sm{lt},sta_sm{lt},activityMig{lt});
    t3=toc; disp(['Adding windspeed ' num2str(t2,3) ' sec'])
    % grt.mvt_pdf = movementModel('step');
    grt.mvt_pdf = movementModel('energy',tblLog.CommonName{lt});
    % grt.thr_as = find(cumsum(grt.mvt_pdf(1:1000))./sum(grt.mvt_pdf(1:1000))>thr_as_percentile,1);
    grt.thr_as = 100;
    grt = filterGraph(grt,'as',grt.thr_as);
    t4=toc; disp(['Filter airspeed at ' num2str(grt.thr_as) 'km/h in ' num2str(t4,3) ' sec'])
    trun{lt}.create_graph=[t1 t2 t3 t4]; 

    grt.p = grt.ps .* grt.mvt_pdf(abs(grt.as));
    tic
    grt.M = probMapGraph(grt);
    trun{lt}.prob_map=toc;
    tic
    grt.sp = shortestPathGraph(grt);
    trun{lt}.shortestpath=toc;
    gr{lt} = grt;
    t=toc; disp(['Finished in ' num2str(t,3) ' sec'])
    trun{lt}.prob_map=toc; 
    disp('----')
end


% Simulation path with gibbs
nj=1000;
for lt=7:height(tblLog)

    % if any(strcmp(skip_gdl, tblLog.GDL_ID{lt})), continue; end
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
    trun{lt}.run_gibbs = toc;
end


%% Save
save('../data/graph'+project+'.mat','gr','trun','sta_sm','-v7.3')

%% Export to geotiff
scriptAltPres()

for lt=1:height(tblLog)
    if any(strcmp(skip_gdl, tblLog.GDL_ID{lt})), continue; end
    M = cat(4, pres_thr{lt}(:,:,sta_sm{lt}.staID),pres_prob{lt}(:,:,sta_sm{lt}.staID), light_prob{lt}(:,:,sta_sm{lt}.staID), gr{lt}.M);
    exportGeotiff(lat{lt},lon{lt},M,raw{lt},tblLog(lt,:),sta_sm{lt},alt{lt},[gr{lt}.sp.lon gr{lt}.sp.lat])
end












