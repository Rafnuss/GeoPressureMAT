
% Pressure Map New

pres_thr = cell(height(tblLog),1);
pres_prob = cell(height(tblLog),1);
pres_n = cell(height(tblLog),1);

for lt=1:height(tblLog)
    sta{lt}.s = repmat(tblLog.std_pres(lt),height(sta{lt}),1);
    sta{lt}.s(sta{lt}.status=="equipment" | sta{lt}.status=="retrieval") = tblLog.std_pres_calib(lt);
    tic
    [pres_prob{lt}, pres_thr{lt}, pres_n{lt}] = getPressueMap(raw{lt}.pressure,sta{lt},lon{lt},lat{lt});
    toc
end

% Check
for lt=1:height(tblLog)
    map = pres_thr{lt} .* pres_prob{lt};
    out = getPressueTimeseries(map,sta{lt},lon{lt},lat{lt});
    
    figure; hold on;
    plot(raw{lt}.pressure.date,raw{lt}.pressure.obs,'-k')
    plot(raw{lt}.pressure.date,raw{lt}.pressure.obsWithOutliars,'-r','linewidth',2)
    plot(out.time,out.pressure)
end