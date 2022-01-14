
%startup
addpath(genpath('../functions'))
load("../data/processedDataStudyPressure.mat")

%% Define the stationary period for the graphs

sta_sm=cell(1,height(tblLog));
for lt=1:height(tblLog)
    if strcmp(tblLog.CommonName{lt},'Eurasian Nightjar')
        grp_id = hours(sta{lt}.end-sta{lt}.start)>48;%sta{lt}.twlNb>=4;
    else
        grp_id = hours(sta{lt}.end-sta{lt}.start)>0;%sta{lt}.twlNb>=4;
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

thr_prob_percentile = .99;
thr_gs = 100;

for lt=1:height(tblLog)
    try
        tic
        disp(['Sarting: ' tblLog.GDL_ID{lt} ' ' num2str(lt)])
        prob_map = pres_prob{lt}(:,:,sta_sm{lt}.staID) .* pres_thr{lt}(:,:,sta_sm{lt}.staID) .* light_prob{lt}(:,:,sta_sm{lt}.staID) .* ~mask_water{lt};
        % prob_map = pres_prob{lt}(:,:,sta_sm{lt}.staID) .* pres_thr{lt}(:,:,sta_sm{lt}.staID) .* ~mask_water{lt};
        [grt,nds] = createGraph(prob_map,lat{lt},lon{lt},raw{lt}.calib,sta_sm{lt}.actEffort,thr_prob_percentile,thr_gs);
        t=toc; disp(['Creating graph in ' num2str(t,3) ' sec'])
        grt = filterGraph(grt,'gs',thr_gs);
        t=toc; disp(['Filter groundspeed ' num2str(t,3) ' sec'])
        grt = windSpeedGraph(grt,raw{lt},sta{lt},sta{lt},activityMig{lt});
        t=toc; disp(['Adding windspeed ' num2str(t,3) ' sec'])
        grt = filterGraph(grt,'as',100);
        t=toc; disp(['Filter airspeed ' num2str(t,3) ' sec'])
        mvt_pdf = movementModel('energy',tblLog.mass(lt),tblLog.wingSpan(lt));
        grt.p = grt.ps .* mvt_pdf(abs(grt.as));
        grt.M = probMapGraph(grt);
        grt.sp = shortestPathGraph(grt);
        gr{lt} = grt;
        t=toc; disp(['Finished in ' num2str(t,3) ' sec'])
        disp('----')
    catch ME
        disp(['Erro with ' tblLog.GDL_ID{lt} ' ' ME.message])
    end
end

%% Simulation path with gibbs
nj=1000;

for lt=6:height(tblLog)
    disp(raw{lt}.GDL_ID)
    
    % Define initial and fixed path
    % intial path as the shortest path
    path0=sub2ind(gr{lt}.snds,gr{lt}.sp.lat,gr{lt}.sp.lon,1:gr{lt}.snds(3));
    % Fixed path for first and possibly last
    fixPath = false(size(path0));
    if isnat(tblLog.CalibSecondStart(lt))
        fixPath(1)=true;
    else
        fixPath([1 end])=true;
    end
    tic
    gr{lt}.psim = gibbsGraph(nj,path0,fixPath,gr{lt});
toc
end

%% Save
save('../data/graph'+project+'.mat','gr','-v7.3')

%% Export to geotiff
for lt=1:height(tblLog)
    M = cat(4, pres_thr{lt}(:,:,sta_sm{lt}.staID),pres_prob{lt}(:,:,sta_sm{lt}.staID), light_prob{lt}(:,:,sta_sm{lt}.staID),gr{lt}.M);
    exportGeotiff(lat{lt},lon{lt},M,raw{lt},tblLog(lt,:),sta_sm{lt},alt{lt},[gr{lt}.lon(gr{lt}.sp.lon), gr{lt}.lat(gr{lt}.sp.lat)])
end












