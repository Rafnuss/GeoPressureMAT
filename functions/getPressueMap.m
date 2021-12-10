function [pres_prob, pres_thr, pres_n] = getPressueMap(pressure,sta,tblLog)
%GETPRESSUEMAP Summary of this function goes here
%   Detailed explanation goes here

%base_url = "https://glp.mgravey.com/test.py?";
base_url = "https://glp.mgravey.com/GeoLocGEEserver.py";

data = struct2table(pressure);
data.date = posixtime(pressure.date);
Ml = (sta.start<pressure.date' & pressure.date' < sta.end) .* sta.staID;
assert(sum(sum(Ml>0)>1)==0,'pressure date belong to two stationary period')
data.label = sum(Ml)';

% remove flight period
data.obs(data.label==0)=nan;

% smoothing
data.obs = movmean(data.obs,3);

% remove outliar
data = data(~isnan(data.obs),:);

% shorten the timeserie for a test
% data = data(1:100,:);
% data = data(data.label==3,:);

urls=webwrite(base_url,...
    "W",tblLog.bndy_W,"S",tblLog.bndy_S,"E",tblLog.bndy_E,"N",tblLog.bndy_N,...
    "time", jsonencode(uint32(data.date)),"pressure",jsonencode(data.obs),"label",jsonencode(data.label),...
    weboptions('Timeout',20));



for i_url=1:numel(urls)
    %A = imread(urls{i_url});
    %[A,R] = readgeoraster(urls{i_url});

    A = webread(urls{i_url},weboptions('Timeout',60*5));

    A(:,:,2)

    % Get std
    s=repmat(tblLog.std_pres(lt),1,height(sta));
    s(sta.status=="equipment" || sta.status=="retrieval") = tblLog.std_pres_calib(lt);

    if isnan(s)
        warning('s cannot be zero. Use default 1')
        s=ones(1,height(sta));
    end

    w = log(pres_n)-1;

    %pres_rmse = w .* mean((et/s).^2,3,'omitnan');
    %pres_prob(:,:,i_s) = exp(-pres_rmse);

    pres_prob(:,:,i_s) = exp(-w.*A./(s.^2));

end

pres_n = splitapply(@numel,data.obs,findgroups(data.label));


end

