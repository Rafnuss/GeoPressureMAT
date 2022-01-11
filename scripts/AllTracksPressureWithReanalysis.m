load('../data/processedDataStudyPressure.mat')
addpath(genpath('../functions'))
scriptAltPres()


tlt4=[];
for lt=1:height(tblLog)
    
    tlt1 = struct2table(raw{lt}.pressure);
    tlt1.obs=tlt1.obsWithOutliars;
    tlt1 = removevars(tlt1,'obsWithOutliars');
    tlt1.series(:) = "presGL";
    tlt1.sta_id(:) = "";
    
    tmp1=[];  tmp2=[]; tmp3=[];
    for i_s = 1:height(sta{lt})
        id_tgr = find(sta{lt}.start(i_s)<spttime & spttime < sta{lt}.end(i_s));
        id_tge = sta{lt}.start(i_s)<raw{lt}.pressure.date & raw{lt}.pressure.date < sta{lt}.end(i_s);
        pres_ge = movmean(raw{lt}.pressure.obs(id_tge),3);
        id_t_gre = ismember(raw{lt}.pressure.date(id_tge),spttime(id_tgr));
        pres_gr = pres_ge(id_t_gre);
        
        tmp1 = [tmp1 ; sp{lt}{i_s}.time ];
        tmp2 = [tmp2; sp{lt}{i_s}.pres'-nanmean(sp{lt}{i_s}.pres)+nanmean(pres_gr)];
        tmp3 = [tmp3; i_s*ones(size(sp{lt}{i_s}.time))];
    end
    
    tlt2 = table(tmp1,tmp2,tmp3,'VariableNames',{'date','obs','sta_id'});
    tlt2.isOutliar(:) = false;
    tlt2.series(:) = "presERA5";
    
    tlt3= [tlt1;tlt2];
    % tlt3.GDL(:) = string(raw{lt}.GDL_ID);
    tlt3.date.Format="yyyy-MM-dd HH:mm:SS";
    
% % %     figure; hold on
% % %     plot(tlt3.date(tlt3.series=="presERA5"),tlt3.obs(tlt3.series=="presERA5"))
% % %     plot(tlt3.date(tlt3.series~="presERA5"),tlt3.obs(tlt3.series~="presERA5"))


    writetable(tlt3,"../data/export/label/AllTracksPressureWithReanalysis_"+string(raw{lt}.GDL_ID)+".csv")
end

str = "";
for lt=1:height(tblLog)
    str= str+"<option value='"+string(raw{lt}.GDL_ID)+"'>"+string(tblLog.CommonName{lt})+" ("+string(raw{lt}.GDL_ID)+")</option>";
end