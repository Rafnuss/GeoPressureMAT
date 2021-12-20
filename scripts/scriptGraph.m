
%startup
addpath(genpath('../functions'))
load("../data/processedDataStudyPressure.mat")

%% 
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

%% 

gr=cell(height(tblLog),1);
% lt=find(tblLog.GDL_ID == "22QO"); % 22QL

thr_prob_percentile = .99;


for lt=1:height(tblLog)
    prob_map = pres_prob{lt}(:,:,sta_sm{lt}.staID) .* pres_thr{lt}(:,:,sta_sm{lt}.staID) .* light_prob{lt}(:,:,sta_sm{lt}.staID) .* ~mask_water{lt};
    tic
    [grt,nds] = createGraph(prob_map,lat{lt},lon{lt},raw{lt}.calib,sta_sm{lt}.actEffort,thr_prob_percentile);
    grt = filterGraph(grt,'gs',100);
    grt = windSpeedGraph(grt,raw{lt},sta{lt},sta{lt},activityMig{lt});
    grt = filterGraph(grt,'as',70);
    mvt_pdf = movementModel('energy',tblLog.mass(lt),tblLog.wingSpan(lt));
    grt.p = grt.ps .* mvt_pdf(abs(grt.as));
    grt.M = probMapGraph(grt);
    grt.sp = shortestPathGraph(grt);
    gr{lt} = grt;
    t=toc;
    disp([tblLog.GDL_ID{lt} ' ' num2str(t,3) ' sec'])
end

save('../data/graph'+project+'.mat','gr','-v7.3')




%% Export to geotiff
for lt=1:height(tblLog)
    M = cat(4, pres_thr{lt}(:,:,sta_sm{lt}.staID),pres_prob{lt}(:,:,sta_sm{lt}.staID), light_prob{lt}(:,:,sta_sm{lt}.staID),gr{lt}.M);
    exportGeotiff(lat{lt},lon{lt},M,raw{lt},tblLog(lt,:),sta_sm{lt},alt{lt},[gr{lt}.lon(gr{lt}.sp.lon), gr{lt}.lat(gr{lt}.sp.lat)])
end












