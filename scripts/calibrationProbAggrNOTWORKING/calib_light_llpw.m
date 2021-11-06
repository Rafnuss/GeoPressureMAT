
% weighting scheme
w = @(n,alpha) alpha + (1-alpha)/n;
w = @(n,alpha) alpha;

% alpha=[0 .1 0.5 .9 1];
alpha=0:.1:1;

nn=[1 4 10 20 50];
Ntest = 30;
q=nan(Ntest,numel(nn),numel(alpha),2,height(tblLog));

for lt=1:2:height(tblLog)
    lon = round(raw{lt}.calib.lon)+(-20:.25:20);
    lat = round(raw{lt}.calib.lat)+(-20:.25:20);

    % Compute probability map
    tmp = coordMapProb(twl{lt}, gE{lt}, lon, lat);

    % calib location
    [gLON,gLAT] = meshgrid(lon, lat);
    [~,id_calib] = min((gLAT(:) - raw{lt}.calib.lat).^2 + (gLON(:) - raw{lt}.calib.lon).^2);

%     w=0.1;
%     light_prob{lt} = nan(numel(lat{lt}),numel(lon{lt}),height(sta{lt}));
%     for i_s = 1:height(sta{lt})
%         light_prob{lt}(:,:,i_s) = exp(w*sum(log(tmp(:,:,twl{lt}.staID==i_s&~twl{lt}.isOutliar)),3));
%     end
    
    for i_sta=1:2

        if i_sta==1
            i_s = find(sta{lt}.status=="equipment");
        else
            i_s = find(sta{lt}.status=="retrieval");
            if isempty(i_s)
                break
            end
        end
        id_possible = find(twl{lt}.staID==i_s&~twl{lt}.isOutliar);
    
        for i_n=1:numel(nn)
            for ii=1:Ntest
                ii_t = round(rand*(numel(id_possible)-nn(i_n)))+(1:nn(i_n));
                if numel(id_possible)>max(ii_t)
                    for i_alpha=1:numel(alpha)
                        prob = exp(w(nn(i_n),alpha(i_alpha)).*sum(log(tmp(:,:,id_possible(ii_t))),3));
                        q(ii,i_n,i_alpha,i_sta,lt) = sum(prob(prob<prob(id_calib)))./sum(prob(:));

%                         figure; hold on;
%                         imagesc(lon, lat,prob)
%                         plot(gLON(id_calib),gLAT((id_calib)),'rx')


                    end
                end
            end
        end
    end
    lt
end


%%

figure; 
for lt=1:height(tblLog)
ha=tight_subplot(numel(nn),numel(alpha)); u=1;
for i_n=1:numel(nn)
    for i_alpha=1:numel(alpha)
        axes(ha(u)); u=u+1; hold on;
        
        histogram(q(:,i_n,i_alpha,:,lt)); 
        xlim([0 1])
        if (i_alpha==1); ylabel(nn(i_n));end
        if (i_n==numel(nn)); xlabel(alpha(i_alpha));end
    end
end
end

figure;
ha=tight_subplot(numel(nn),numel(alpha)); u=1;
for i_n=1:numel(nn)
    for i_alpha=1:numel(alpha)
        axes(ha(u)); u=u+1; hold on;
        histogram(nanmedian(q(:,i_n,i_alpha,:,:),1),0:.1:1); 
      xlim([0 1])
        if (i_alpha==1); ylabel(nn(i_n));end
        if (i_n==numel(nn)); xlabel(alpha(i_alpha));end
    end
end
