load('processedDataStudyPressure.mat')
addpath(genpath('../functions'))
scriptAltPres()


tlt4=[];
for lt=1:height(tblLog)
    
    tlt1 = struct2table(raw{lt}.pressure);
    tlt1.obs=tlt1.obsWithOutliars;
    tlt1 = removevars(tlt1,'obsWithOutliars');
    tlt1.series(:) = "presGL";
    
    tmp1=[];  tmp2=[];
    for i_s = 1:height(sta{lt})
        id_tgr = find(sta{lt}.start(i_s)<spttime & spttime < sta{lt}.end(i_s));
        id_tge = sta{lt}.start(i_s)<raw{lt}.pressure.date & raw{lt}.pressure.date < sta{lt}.end(i_s);
        pres_ge = movmean(raw{lt}.pressure.obs(id_tge),3);
        id_t_gre = ismember(raw{lt}.pressure.date(id_tge),spttime(id_tgr));
        pres_gr = pres_ge(id_t_gre);
        
        tmp1 = [tmp1 ; sp{lt}{i_s}.time ];
        tmp2 = [tmp2; sp{lt}{i_s}.pres'-nanmean(sp{lt}{i_s}.pres)+nanmean(pres_gr)];
    end
    
    tlt2 = table(tmp1,tmp2,'VariableNames',{'date','obs'});
    tlt2.isOutliar(:) = false;
    tlt2.series(:) = "presERA5";
    
    tlt3= [tlt;tlt2];
    tlt3.GDL(:) = string(raw{lt}.GDL_ID);
    
    tlt4 = [tlt4;tlt3];
end

writetable(tlt4,'../data/labels/AllTracksPressureWithReanalysis.csv')
