
project="StudyPressure";
% project="StudyKenya";

addpath(genpath('../functions'))
addpath('~/Documents/GitHub/Flight-Matlab/functions/')
addpath('~/Documents/GitHub/Flight-Matlab/data/')
load("../data/processedData"+project+".mat")


sta_sm=cell(1,height(tblLog));
for lt=1:height(tblLog)
    if any(strcmp({'24FF','24SZ'}, tblLog.GDL_ID{lt}))
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

for lt=1:height(tblLog)
   
    tic; disp(['Sarting: ' tblLog.GDL_ID{lt} ' ' num2str(lt)])
    prob_map = pres_prob{lt}(:,:,sta_sm{lt}.staID) .* pres_thr{lt}(:,:,sta_sm{lt}.staID) .* light_prob{lt}(:,:,sta_sm{lt}.staID) .* ~mask_water{lt};
    % prob_map = pres_prob{lt}(:,:,sta_sm{lt}.staID) .* pres_thr{lt}(:,:,sta_sm{lt}.staID) .* ~mask_water{lt};
    [grt,~] = createGraph(prob_map,lat{lt},lon{lt},raw{lt}.calib,sta_sm{lt}.actEffort,thr_prob_percentile,thr_gs);
    t1=toc; disp(['Creating graph in ' num2str(t1,3) ' sec'])
    
    grt = filterGraph(grt,'gs',thr_gs);
    t2=toc; disp(['Filter groundspeed at ' num2str(thr_gs) 'km/h in ' num2str(t2-t1,3) ' sec'])
    
    grt = windSpeedGraph(grt,raw{lt},sta{lt},sta_sm{lt},activityMig{lt});
    t3=toc; disp(['Adding windspeed ' num2str(t3-t2,3) ' sec'])
    
    % grt.mvt_pdf = movementModel('step');
    grt.mvt_pdf = movementModel('energy',tblLog.CommonName{lt});
    % grt.thr_as = find(cumsum(grt.mvt_pdf(1:1000))./sum(grt.mvt_pdf(1:1000))>thr_as_percentile,1);
    grt.thr_as = 100;
    grt = filterGraph(grt,'as',grt.thr_as);
    t4=toc; disp(['Filter airspeed at ' num2str(grt.thr_as) 'km/h in ' num2str(t4-t3,3) ' sec'])
    grt.pt = grt.mvt_pdf(abs(grt.as));
    trunt.create_graph=[t1 t2 t3 t4]; 

    
    tic; [grt.M, grt.E ] = probMapGraph(grt);
    trunt.prob_map=toc;
    disp(['Prob map in ' num2str(trunt.prob_map,3) ' sec'])

    tic; grt.sp = shortestPathGraph(grt);toc
    trunt.shortestpath=toc;
    disp(['Shortest path in ' num2str(trunt.shortestpath,3) ' sec'])

    % Simulation path
    tic
    nj=1000;
    grt.psim = sequantialSimulationGraph(nj,grt);
    trunt.sim=toc;
    disp(['Simulation of ' num2str(nj) ' paths in ' num2str(trunt.sim,3) ' sec']) 

    t=toc; disp(['Finished in ' num2str(t4+trunt.prob_map+trunt.shortestpath+trunt.sim,3) ' sec'])
    disp('----')

    sta_smt=sta_sm{lt};
    % save("../data/graph/"+raw{lt}.GDL_ID,"grt","trunt","sta_smt","thr_prob_percentile","thr_gs")
    gr{lt} = grt;
    trun{lt} = trunt;
end


% Simulation path with gibbs
% nj=1000;
% for lt=7:height(tblLog)
% 
%     % if any(strcmp(skip_gdl, tblLog.GDL_ID{lt})), continue; end
%     disp(raw{lt}.GDL_ID)
%     
%     % Define intial path as the shortest path
%     path0=gr{lt}.sp.path;
%     % Define the fixed node/sta for first and possibly last
%     fixPath = false(size(path0));
%     fixPath(1)=true;
%     if ~isnat(tblLog.CalibSecondStart(lt))
%         fixPath(end)=true;
%     end
%     tic
%     gr{lt}.psim = gibbsGraph(nj,path0,fixPath,gr{lt});
%     trunt.run_gibbs = toc;
% end


%% Save
save('../data/graph'+project+'.mat','gr','trun','sta_sm','-v7.3')
load('../data/graph'+project+'.mat')
%% Export to geotiff
scriptAltPres()

for lt=1:height(tblLog)
    if any(strcmp(skip_gdl, tblLog.GDL_ID{lt})), continue; end
    M = cat(4, pres_thr{lt}(:,:,sta_sm{lt}.staID),pres_prob{lt}(:,:,sta_sm{lt}.staID), light_prob{lt}(:,:,sta_sm{lt}.staID), gr{lt}.M);
    exportGeotiff(lat{lt},lon{lt},M,raw{lt},tblLog(lt,:),sta_sm{lt},alt{lt},[gr{lt}.sp.lon gr{lt}.sp.lat])
end












