% Figure Map
err=nan(height(tblLog),3,2);
sz=nan(height(tblLog),2,2);
q=nan(height(tblLog),2,2);

block_time = [.5 1 2 5 10 50];
Ntest = 100;
q=nan(height(tblLog),2,numel(block_time),Ntest);

alpha=0.1;

for lt=1:height(tblLog)

    id_lon = find(glon>=tblLog.bndy_W(lt) & glon<=tblLog.bndy_E(lt));
    id_lat = find(glat>=tblLog.bndy_S(lt) & glat<=tblLog.bndy_N(lt));
    id_t = find(raw{lt}.pressure.date(1)<=spttime & spttime <= raw{lt}.pressure.date(end));
    splt = permute(flip(ncread(file,'sp',[id_lon(1) numel(glat)-id_lat(end)+1 1 id_t(1)],[numel(id_lon) id_lat(end)-id_lat(1)+1 1 id_t(end)-id_t(1)+1]),2)/100,[2 1 4 3]);
    if lt==16
        splt_0 = permute(flip(ncread(file,'sp',[id_lon(1) numel(glat)-id_lat(end)+1 2 id_t(1)],[numel(id_lon) id_lat(end)-id_lat(1)+1 1 id_t(end)-id_t(1)+1]),2)/100,[2 1 4 3]);
        splt(:,:,spttime(id_t)>=datetime('1-jul-2021'))=splt_0(:,:,spttime(id_t)>=datetime('1-jul-2021'));
    end
    
    [gLON,gLAT] = meshgrid(lon{lt},lat{lt});
    [~,id_calib] = min((gLAT(:) - raw{lt}.calib.lat).^2 + (gLON(:) - raw{lt}.calib.lon).^2);


    for i=1:(1+~isnat(tblLog.CalibSecondStart(lt)))

        i_s=max(1,(i-1)*(height(sta{lt})));
        
        id_tgr = find(sta{lt}.start(i_s)<spttime(id_t) & spttime(id_t) < sta{lt}.end(i_s));
        id_tge = sta{lt}.start(i_s)<raw{lt}.pressure.date & raw{lt}.pressure.date < sta{lt}.end(i_s);
        % pres_rmse(:,:,i_s) = sqrt(mean(((  pres(:,:,id_tgr)-reshape(press_movmean(id_tgr),1,1,[])  ).^2),3));
        
        % Apply a movinn median to remove outliars and smooth with 1hr
        dt = 1./hours(diff(raw{lt}.pressure.date(1:2)));
        pres_ge = movmean(movmedian(raw{lt}.pressure.obs(id_tge),3),dt);

        % Downscale to same time scale
        id_t_gre = ismember(raw{lt}.pressure.date(id_tge),spttime(id_t(id_tgr)));
        pres_gr = pres_ge(id_t_gre);
        
        id_tgr=id_tgr(~isnan(pres_gr));
        pres_gr=pres_gr(~isnan(pres_gr));
       
        x = splt(:,:,id_tgr);
        y = reshape(pres_gr,1,1,[]);
        N=numel(pres_gr);

        for i_b=3:numel(block_time)
            fun = @(wi) (0.5-nanmean(calibWi(x,y,N,Ntest,block_time(i_b),id_calib,wi))).^2;
            wi(i_b)=fminsearch(fun,1);
        end
    end
end

figure; ha=tight_subplot(numel(block_time),1);
for i_b=1:numel(block_time)
    axes(ha(i_b))
    figure; histogram(q(:,:,i_b,:)); xlim([0 1])
    ylabel(num2str(block_time(i_b)));
    
end