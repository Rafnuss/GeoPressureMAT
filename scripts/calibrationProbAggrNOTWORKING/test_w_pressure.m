addpath(genpath('../functions'))
load("../data/processedDataStudyPressure.mat")
scriptAltPres()
%%
tic

% Load Pressure data
file='../data/ECMWF/surface_pressure.nc'; 
spttime = datetime(double(ncread(file,'time'))/24 + datenum('1900-01-01 00:00:00'),'convertFrom','datenum');

pres_mse = cell(height(tblLog),1);
pres_n = cell(height(tblLog),1);
pres_probC{1} = cell(height(tblLog),1);
pres_probC{2} = cell(height(tblLog),1);

for lt=1:height(tblLog)
    
    id_lon = find(glon>=tblLog.bndy_W(lt) & glon<=tblLog.bndy_E(lt));
    id_lat = find(glat>=tblLog.bndy_S(lt) & glat<=tblLog.bndy_N(lt));
    id_t = find(raw{lt}.pressure.date(1)<=spttime & spttime <= raw{lt}.pressure.date(end));
    splt = permute(flip(ncread(file,'sp',[id_lon(1) numel(glat)-id_lat(end)+1 1 id_t(1)],[numel(id_lon) id_lat(end)-id_lat(1)+1 1 id_t(end)-id_t(1)+1]),2)/100,[2 1 4 3]);
    if raw{lt}.pressure.date(end)>datetime('1-jul-2021')
        splt_0 = permute(flip(ncread(file,'sp',[id_lon(1) numel(glat)-id_lat(end)+1 2 id_t(1)],[numel(id_lon) id_lat(end)-id_lat(1)+1 1 id_t(end)-id_t(1)+1]),2)/100,[2 1 4 3]);
        splt(:,:,spttime(id_t)>=datetime('1-jul-2021'))=splt_0(:,:,spttime(id_t)>=datetime('1-jul-2021'));
    end
    
    % Calibration location
    [~,id_lat_tmp]=min(abs(raw{lt}.calib.lat-glat(id_lat)));
    [~,id_lon_tmp]=min(abs(raw{lt}.calib.lon-glon(id_lon)));
    
    % Compute Pressue error map for each stationary period
    pres_mse{lt} = nan(numel(lat{lt}),numel(lon{lt}),height(sta{lt}));
    pres_probC{1}{lt} = pres_mse{lt};
    pres_probC{2}{lt} = pres_mse{lt};
    pres_n{lt} = nan(height(sta{lt}),1);
    
    ss=find(sta{lt}.status=="equipment" | sta{lt}.status=="retrieval");
    for i=1:numel(ss)
        i_s=ss(i);
        
        id_tgr = find(sta{lt}.start(i_s)<spttime(id_t) & spttime(id_t) < sta{lt}.end(i_s));
        id_tge = sta{lt}.start(i_s)<raw{lt}.pressure.date & raw{lt}.pressure.date < sta{lt}.end(i_s);

        pres_ge = movmean(raw{lt}.pressure.obs(id_tge),3);
        
        % Downscale to same time scale
        id_t_gre = ismember(raw{lt}.pressure.date(id_tge),spttime(id_t(id_tgr)));
        pres_gr = pres_ge(id_t_gre);
        
        pres_n{lt}(i_s)=sum(~isnan(pres_gr));
        
        if any(~isnan(pres_gr))
            
            x = splt(:,:,id_tgr)-mean(splt(:,:,id_tgr),3,'omitnan');
            y = reshape(pres_gr-mean(pres_gr,'omitnan'),1,1,[]);
            e = x-y;
            
            if tblLog.Truncated(lt) && i_s~=0 && i_s~=height(height(sta{lt})) && min(min(std(e,[],3)))>1
                e(:)=0;
            end

            % Covariance
            et=pres_gr(:)-sp{lt}{i_s}.presCalib(1:numel(pres_gr))';
            et = et(:)-mean(et,'omitnan');
            xi=0:24*2;
            covxi = autocorr(et(:), min(sum(~isnan(et))-1,24*5)).*var(et,'omitnan');
            Cline = zeros(size(et));Cline(1:numel(covxi))=covxi;
            C = toeplitz(Cline);
            
            A = squeeze(pagemtimes(permute(e,[4 3 1 2]),pagemtimes(inv(C),permute(e,[3 4 1 2]))));
            pres_probC{1}{lt}(:,:,i_s) = A;
            % imagesc(exp(-1/2*A))
            % 1/((2*pi)^2*sqrt(det(C))) * exp(-1/2*A)
 
             ff = fit((1:numel(covxi))',covxi,@(c,x) var(et,'omitnan').*exp(-x/c));
             cs(lt,1,(i_s==1)+1) = ff.c;
             cs(lt,2,(i_s==1)+1) = var(et,'omitnan');
%             % plot(ff); hold on; plot(covxi)
             Ch = feval(ff,1:numel(et));
             Ch(1)=var(et,'omitnan');
             C = toeplitz(Ch);
%             (2*pi)^numel(et).*det(C)
            A = squeeze(pagemtimes(permute(e,[4 3 1 2]),pagemtimes(inv(C),permute(e,[3 4 1 2]))));
            pres_probC{2}{lt}(:,:,i_s) = A;
            
            % Assess the match
            pres_mse{lt}(:,:,i_s) = mean(e(:,:,:).^2,3,'omitnan');
        end
        
    end
end
toc % 150s

%%
nn= 12:3500;

scheme = {@(n) ones(size(n)), @(n) 1./n, @(n) log(n)./n, @(n) 3*log(n)./n, @(n) 5*log(n)./n};

figure; hold on
for i_sc = 1:numel(scheme)
    plot(nn,scheme{i_sc}(nn),'DisplayName',func2str(scheme{i_sc}));
end
set(gca,'xscale','log'); legend

%%
pres_prob = cell(numel(scheme)+1,1);

for i_sc = 1:numel(scheme)
    pres_prob{i_sc} = cell(height(tblLog),1);
    for lt=1:height(tblLog)
        ss=find(sta{lt}.status=="equipment" | sta{lt}.status=="retrieval");
        for i=1:numel(ss)
            i_s=ss(i);

            if sta{lt}.status(i_s)=="equipment" || sta{lt}.status(i_s)=="retrieval"
                s=tblLog.std_pres_calib(lt);
            else
                s=tblLog.std_pres(lt);
            end
            if isnan(s)
                warning('s cannot be zero. Use default 1')
                s=1;
            end

            n = pres_n{lt}(i_s);
            w = scheme{i_sc}(n);

            f_prob  = @(x) (1/(2*pi*s^2))^(n*w/2)*exp(-w*n/2/(s^2)*x);

            pres_prob{i_sc}{lt}(:,:,i_s) = f_prob(pres_mse{lt}(:,:,i_s));
        end
    end
end

pres_prob{i_sc+1}=pres_probC{1};
pres_prob{i_sc+2}=pres_probC{2};
%%
% figure('position',[0 0 675 675]); ha = tight_subplot(16,2,0,[0.03 0],[0.03 0]);
sz=nan(height(tblLog), 2, 2);
err = sz;
for lt=1:height(tblLog)
    [gLON,gLAT] = meshgrid(lon{lt},lat{lt});

    ss=find(sta{lt}.status=="equipment" | sta{lt}.status=="retrieval");
    for i=1:numel(ss)
        i_s=ss(i);
        
        for i_sc = 1:numel(pres_prob)
            tmp_p = pres_prob{i_sc}{lt}(:,:,i_s);
            [~,id_ppt]=max(tmp_p(:));

            err(lt,2,(i_s==1)+1)=lldistkm([gLAT(id_ppt) gLON(id_ppt)], [raw{lt}.calib.lat raw{lt}.calib.lon]);
            sz(lt,2,(i_s==1)+1) = days(sum(~isnan(raw{lt}.pressure.obs(raw{lt}.pressure.date > sta{lt}.start(i_s) & raw{lt}.pressure.date < sta{lt}.end(i_s))))*diff(raw{lt}.pressure.date(1:2)));

            [~,id_min] = min((gLAT(:) - raw{lt}.calib.lat).^2 + (gLON(:) - raw{lt}.calib.lon).^2);

            q(lt,i_sc,(i_s==1)+1) = sum(tmp_p(tmp_p<tmp_p(id_min)))./sum(tmp_p(:));
        end
    end
end

figure('position',[0 0 500 500]); hold on;
for i_sc = 1:numel(scheme)
    [f,x] = ecdf(reshape(q(:,i_sc,:),1,[]));
    stairs(x,f,'Linewidth',2,'DisplayName',func2str(scheme{i_sc}))
end
[f,x] = ecdf(reshape(q(:,i_sc+1,:),1,[]));
stairs(x,f,'Linewidth',2,'DisplayName','Covariance')
[f,x] = ecdf(reshape(q(:,i_sc+2,:),1,[]));
stairs(x,f,'Linewidth',2,'DisplayName','Covariance fit')
plot([0 1],[0 1],'--k')
xlabel('Quantile of the equipement site in the distribution of prob. map')
ylabel('CDF(x)'); legend('Location','northoutside','Orientation','horizontal')
box on

% figure;
% scatter(24*reshape(sz(:,2,:),1,[]),reshape(q(:,2,:),1,[]),'o','filled')
% xlabel('Number of datapoint');ylabel('quantile')

%% Assess
err=nan(height(tblLog),3,2);
sz=nan(height(tblLog),2,2);
q=nan(height(tblLog),2,2);

for lt=1:height(tblLog)

    [gLON,gLAT] = meshgrid(lon{lt},lat{lt});
        
    if mod(lt,2)==1
        figure('position',[0 0 675 675], 'Name', [tblLog.CommonName{lt} ' (' raw{lt}.GDL_ID ',' raw{lt+1}.GDL_ID ')' ]);
        ha = tight_subplot(2,2,[.03 .01],[.02 .03],[.03 .01]);
        u=1;
    end
    ss=find(sta{lt}.status=="equipment" | sta{lt}.status=="retrieval");
    for i=1:numel(ss)
        i_s=ss(i);
        axes(ha(u)); hold on;set(gca,'Color','k')
        
        tmp_p = pres_prob{lt}(:,:,i_s);
        tmp_pt = ones(size(pres_thr{lt}(:,:,i_s)));
        tmp_l = light_prob{lt}(:,:,i_s);
        tmp_w = ~mask_water{lt};
        
        [~,id_l]=max(tmp_l(:) .* tmp_w(:));
        [~,id_ppt]=max(tmp_p(:).*tmp_pt(:).* tmp_w(:));
        
        mask = 0.3+0.7*double(tmp_pt);
        mask(mask_water{lt}) = 0;
        img_tmp = real2rgb(tmp_p,colormap);
        img_tmp = img_tmp.*mask;
        
        h0=imagesc(lon{lt},lat{lt},img_tmp);
        c_axis=caxis();
        [~,h1]=contour(lon{lt},lat{lt},tmp_l,3,'color',[255, 252, 49]/255,'linewidth',2);
        caxis(c_axis);
        borders('countries','w')
        
        tt = datestr(sta{lt}.start(i_s),'dd-mmm');
        tt = tt + "-" + datestr(sta{lt}.end(i_s),'dd-mmm')+ " | ";
        tt = tt + num2str(sta{lt}.twlNb(i_s))+" twls"+ " | ";
        tt = tt + num2str(pres_n{lt}(i_s) ) + " pres. pts";
        title(tt)
        
        h2=scatter(gLON(id_l), gLAT(id_l),200,'MarkerFaceColor','y','MarkerEdgecolor','k');
        h3=scatter(gLON(id_ppt), gLAT(id_ppt),200,'MarkerFaceColor','b','MarkerEdgecolor','k');
        h4=plot(raw{lt}.calib.lon,raw{lt}.calib.lat,'xr','markersize',15,'linewidth',3);
        
        axis equal; axis([raw{lt}.calib.lon-5 raw{lt}.calib.lon+5 raw{lt}.calib.lat-5 raw{lt}.calib.lat+5]);
        
        x0=[raw{lt}.calib.lon-4 raw{lt}.calib.lat-4];
        dx = 1./[lldistkm(x0,x0+[1 0]) lldistkm(x0,x0+[0 1])]*100;
        h=plot([x0(1) x0(1)], [x0(2) x0(2)+dx(2)], '-r', 'LineWidth', 2);%label(h,'100km','location','middle','slope')
        h=plot([x0(1) x0(1)+dx(1)], [x0(2) x0(2)],'-r', 'LineWidth', 2); label(h,'100km');
        
        if i_s==1
            ylabel(raw{lt}.GDL_ID)
        else
            ax_tmp = [min([get(ha(u-1),'ylim') get(ha(u),'ylim')]) max([get(ha(u-1),'ylim') get(ha(u),'ylim')])];
            set(ha(u),'ylim',ax_tmp)
            set(ha(u-1),'ylim',ax_tmp)
            set(ha(u),'yticklabels','');
        end
        
        if mod(lt,2)==1 && i_s==1
            legend([h2 h3 h4],'Best light','Best pressure','Equipment site','color','w')
        end
        if mod(lt,2)==0
            if i_s==1
                xlabel('First Calibration period (before leaving)')
            else
                xlabel('Second Calibration period (after returning)')
            end
        end
        u=u+1;
    end
    u=3;
    if mod(lt,2)==0
        % exportgraphics(gcf,['pressure_light_calib_map_' tblLog.CommonName{lt} '_' raw{lt-1}.GDL_ID '_' raw{lt}.GDL_ID '.png'],'Resolution',300)
    end
    % keyboard
end
