
addpath(genpath('../functions'))
project= "StudyPressure";
load("../data/processedData"+project+".mat")
load('../data/graph'+project+'.mat')


%% Figure illustration graph

lt=1;
% Build the graph and save the graph status at all the steps
thr_prob_percentile = .99;
prob_map = pres_prob{lt}(:,:,sta{lt}.staID) .* pres_thr{lt}(:,:,sta{lt}.staID) .* light_prob{lt}(:,:,sta{lt}.staID) .* ~mask_water{lt};
[gr1,nds] = createGraph(prob_map,lat{lt},lon{lt},raw{lt}.calib,sta{lt}.actEffort,thr_prob_percentile);
id = find(abs(gr1.gs)>100);
gr2=gr1; gr2.s(id)=[]; gr2.t(id)=[]; gr2.gs(id)=[];
gr3 = filterGraph(gr2,'gs',100);
gr4 = windSpeedGraph(gr3,raw{lt},sta{lt},sta{lt},activityMig{lt});
gr4 = filterGraph(gr4,'as',70);

% Illustrate with a figure
mat = nan(gr1.snds);
mat(unique(gr1.s))=1;
mat(unique(gr2.s))=2;
mat(unique(gr3.s))=3;
mat(unique(gr4.s))=4;

figure('position',[0 0 1400 900]); 
tiledlayout('flow','TileSpacing','none','Padding','none')
col = brewermap(4,'YlOrRd'); colormap(col)
for i_s = 1:find(sta{lt}.status=="wintering") % height(sta{lt})
    nexttile(i_s);hold on;
    set(gca,'Color','k')   
    axis equal; axis([3 max(lon{lt}) min(lat{lt}) max(lat{lt}) ]);
    xticklabels(''); yticklabels('')
    
    img_bkg = real2rgb(~mask_water{lt},[0 0 0;19 12 51]/255);
    img_tmp = real2rgb(mat(:,:,i_s),col);
    img_bkg(repmat(~isnan(mat(:,:,i_s)),1,1,3)) = img_tmp(repmat(~isnan(mat(:,:,i_s)),1,1,3));

    imagesc(lon{lt},lat{lt},img_bkg)
    
    %imagesc(gr1.lon,gr1.lat,prob_map(:,:,i_s))
    %imagesc(gr1.lon,gr1.lat,mat(:,:,i_s),'alphadata',~isnan(mat(:,:,i_s)))
    %surf(gr1.lon,gr1.lat,mat(:,:,i_s),'facealpha',0,'edgecolor', 'interp','alphadata',~isnan(mat(:,:,i_s)))
    % s= pcolor(gr1.lon,gr1.lat,mat(:,:,i_s))
    borders('countries','color',[.6 .6 .6]);
end
% exportgraphics(gcf,'test.eps','BackgroundColor','none','ContentType','vector')



%% Map of transition
% grf = gr{lt};
% [Slat,Slon,St]=ind2sub(grf.snds,grf.s);
% [Tlat,Tlon,Tt]=ind2sub(grf.snds,grf.t);
% nrow=4;
% figure('position',[0 0 1400 900]); colormap(col)
% tiledlayout('flow','TileSpacing','none','Padding','none')
% for i_s = 1:find(sta{lt}.status=="wintering")
%     nexttile; hold on;
%     st_id = find(St==i_s);
%     borders('countries','color',[.6 .6 .6]);
%     X = [grf.lon(Slon(st_id)) grf.lat(Slat(st_id))]; X=unique(X ,"rows");
%     plot(X(:,1),X(:,2),'.k')
%     X = [grf.lon(Tlon(st_id)) grf.lat(Tlat(st_id))]; X=unique(X ,"rows");
%     plot(X(:,1),X(:,2),'.r')
%     axis equal; axis([3 max(grf.lon) min(grf.lat) max(grf.lat) ]);
%     xticklabels(''); yticklabels('')
% end



%% 
for lt=1:height(tblLog)
    if isempty(gr{lt}), continue; end
    figure('position',[0 0 900 1200], 'Name', [raw{lt}.GDL_ID ' | ' tblLog.CommonName{lt}] ); hold on
    
    set(gca,'Color',[.5 .5 .5])
    path_mean=nan(height(sta{lt}),2);
    imagesc(lon{lt},lat{lt},ones(numel(lat{lt}),numel(lon{lt}),3).*reshape([0 0 0],1,1,3),'AlphaData',~mask_water{lt});
    
    borders('countries','w')
    plot(raw{lt}.calib.lon,raw{lt}.calib.lat,'xr','linewidth',2)
    axis equal; 
    % axis([min(lon{lt}) max(lon{lt}) min(lat{lt}) max(lat{lt}) ]);
    [lonext, latext]=path2lonlat(gr{lt}.s,gr{lt});
    axis([min(lonext) max(lonext) min(latext) max(latext) ]);
    set(gca,'ydir','normal');xticks([]);yticks([])
    set(gca,'LooseInset',get(gca,'TightInset'));
    
    colordergrad=crameri('batlow',height(sta{lt}));
    t = sta{lt}.start + (sta{lt}.end-sta{lt}.start)/2;
    
    colorinterp = interp1(datenum(t),colordergrad,datenum(sta{lt}.start(1):sta{lt}.end(end)));
    colormap(colorinterp)
    c=colorbar('south'); c.Color='w'; c.FontSize=12;
    c.Ticks=datenum(unique(dateshift(t(1):t(end),'start','month'))-t(1))/datenum(t(end)-t(1));
    c.TickLabels=datestr(unique(dateshift(t(1):t(end),'start','month')),'mmm');
       
    
    gif(['graph_probMap_sp_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '.gif'],'overwrite',true,'DelayTime',(5+height(sta{lt})/10)/height(sta{lt}),'frame',gca)
    
    
    % plot(lon{lt}(gr{lt}.psim.lon(:,1:100:end)), lat{lt}(gr{lt}.psim.lat(:,1:100:end)),'color',[.7 .7 .7]);
    % idrd = randsample( size(gr{lt}.psim.lon,2) , 30 );
    % plot(gr{lt}.psim.lon(:,idrd), gr{lt}.psim.lat(:,idrd),'color',[.7 .7 .7]);
    
    % plot(mean(lon{lt}(gr{lt}.psim.lon),2), mean(lat{lt}(gr{lt}.psim.lat),2),'linewidth',10,'color',[.7 .7 .7]);

    for i_s = 1:height(sta{lt})
        f=gr{lt}.M(:,:,i_s);
        imagesc(lon{lt},lat{lt},ones(size(f,1),size(f,2),3).*reshape(colordergrad(i_s,:),1,1,3),'AlphaData',f./max(f(:)));
        
       
        p=plot(gr{lt}.sp.lon(1:i_s),gr{lt}.sp.lat(1:i_s),'w','linewidth',3);
        
        for i_ss = 1:i_s
            sz = 20+(hours(sta{lt}.end(i_ss)-sta{lt}.start(i_ss)))/100;
            p2(i_ss,1) = plot(gr{lt}.sp.lon(i_ss),gr{lt}.sp.lat(i_ss),'.w', 'MarkerSize',sz+15);
            p2(i_ss,2) = plot(gr{lt}.sp.lon(i_ss),gr{lt}.sp.lat(i_ss),'.', 'MarkerSize',sz,'color',colordergrad(i_ss,:));
        end
        % keyboard
        
        gif
        % tiadd = datenum([sta{lt}.start(i_s) sta{lt}.end(i_s)]-t(1))/datenum(t(end)-t(1));
        %pause(1)
        if i_s<height(sta{lt})
            delete(p);   delete(p2)
        end
    end
    
 
    % exportgraphics(gca,['graph_probMap_sp_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '.png'],'Resolution',600)
    % keyboard
    close all
end























%% Wind

for lt=1:height(tblLog)
    if isempty(gr{lt}), continue; end
    ss = find(hours(sta{lt}.actEffort)>3);

    figure('position',[0 0 900 1200], 'Name', [raw{lt}.GDL_ID ' | ' tblLog.CommonName{lt}] ); hold on
    set(gca,'LooseInset',get(gca,'TightInset'));
    borders('countries','facecolor',[.7 .7 .7],'edgecolor','w')
    plot(raw{lt}.calib.lon,raw{lt}.calib.lat,'xr','linewidth',2)

    plot(gr{lt}.sp.lon,gr{lt}.sp.lat,'color',[.3 .3 .3],'Linewidth',1)
    vwinds = round(linspace(0,ceil(max(abs(gr{lt}.psim.ws(:)))/10)*10,8));
    for i_ss= 1:numel(ss)
        i_s=ss(i_ss);
        plot(gr{lt}.sp.lon(i_s:i_s+1),gr{lt}.sp.lat(i_s:i_s+1),'color','k','Linewidth',2)
        WindRose(rad2deg(angle(gr{lt}.psim.ws(i_s,:))), abs(gr{lt}.psim.ws(i_s,:)),'X',mean(gr{lt}.sp.lon(i_s:i_s+1)),...
            'Y',mean(gr{lt}.sp.lat(i_s:i_s+1)), 'LegendType', 0, 'labels','','TitleString','',...
        'scalefactor',2,'nfreq',0,'vWinds', vwinds(1:end-1),'nDirections',24);
        text(mean(gr{lt}.sp.lon(i_s:i_s+1)),1+mean(gr{lt}.sp.lat(i_s:i_s+1)),string(i_s),'HorizontalAlignment','center','Color','white','FontWeight','bold')
        text(mean(gr{lt}.sp.lon(i_s:i_s+1)),1+mean(gr{lt}.sp.lat(i_s:i_s+1)),string(i_s),'HorizontalAlignment','center')
     % axis equal; axis([min(lon{lt}) max(lon{lt}) min(lat{lt}) max(lat{lt}) ]);
        % keyboard
    end
    
    colormap(crameri('batlow'));
    scatter(gr{lt}.sp.lon,gr{lt}.sp.lat,100,1:height(sta{lt}),'filled','MarkerEdgeColor','k')
    
    axis equal; % axis([min(lon{lt}) max(lon{lt}) min(lat{lt}) max(lat{lt}) ]);
    dll=2; xticks([]); yticks([]); box on;
    axis([min(gr{lt}.sp.lon)-dll max(gr{lt}.sp.lon)+dll min(gr{lt}.sp.lat)-2*dll max(gr{lt}.sp.lat)+dll ]);
    
    ax2=axes;ax2.Visible='off';
    colormap(ax2,parula(numel(vwinds)-1)); c=colorbar('Location','south'); c.Label.String='Windspeed [km/h]'; caxis([0 vwinds(end)]);
    axis equal; set(gca,'LooseInset',get(gca,'TightInset')); title(' '); c.Ticks=vwinds;
    axis([min(gr{lt}.sp.lon)-dll max(gr{lt}.sp.lon)+dll min(gr{lt}.sp.lat)-3*dll max(gr{lt}.sp.lat)+dll ]);
    
    exportgraphics(gcf,['graph_windrose_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '.png'],'resolution',600)
    close all
end


%% Histogram of speed
for lt=1:height(tblLog)
    if isempty(gr{lt}), continue; end
ss = find(hours(sta{lt}.actEffort)>0);

x_axis = (0:.1:100)';
gs = nan(numel(x_axis), numel(ss));
as = gs; ws=gs;
for i_ss= 1:numel(ss)
    i_s=ss(i_ss);
    gs(:,i_s) = ksdensity(abs(gr{lt}.psim.gs(i_s,:)),x_axis);
    ws(:,i_s) = ksdensity(abs(gr{lt}.psim.ws(i_s,:)),x_axis);
    as(:,i_s) = ksdensity(abs(gr{lt}.psim.as(i_s,:)),x_axis);
end

GS = [gs; -flipud(gs)];
AS = [as; -flipud(as)];
WS = [ws; -flipud(ws)];
X_axis = [x_axis;flipud(x_axis)];
scale_width=15;

colordergrad=crameri('batlow',height(sta{lt}));

figure('position',[0 0 1600 900], 'Name', [raw{lt}.GDL_ID ' | ' tblLog.CommonName{lt}]);
tiledlayout(4,1,'TileSpacing','tight','Padding','tight')
i_sl = (1:numel(ss))';
nexttile; hold on
for i_ss= 1:numel(ss)
    i_s=ss(i_ss);
    fill(i_ss+GS(:,i_s)*15,X_axis,colordergrad(i_s,:),'FaceAlpha',0.5)
end
plot([i_sl i_sl]', [quantile(abs(gr{lt}.psim.gs(ss,:)),.25,2) quantile(abs(gr{lt}.psim.gs(ss,:)),.75,2)]', '-w','linewidth',2)
plot(i_sl, mean(abs(gr{lt}.psim.gs(ss,:)),2), '.w', 'MarkerSize',30)
plot(i_sl, abs(gr{lt}.sp.gs(ss,:)), '.r', 'MarkerSize',20)
ylabel('Groundspeed'); xlim([0 numel(ss)+1]); xticks([])
nexttile; hold on
for i_ss= 1:numel(ss)
    i_s=ss(i_ss);
    fill(i_ss+WS(:,i_s)*3,X_axis,colordergrad(i_s,:),'FaceAlpha',0.5)
end
plot([i_sl i_sl]', [quantile(abs(gr{lt}.psim.ws(ss,:)),.25,2) quantile(abs(gr{lt}.psim.ws(ss,:)),.75,2)]', '-w','linewidth',2)
plot(i_sl, mean(abs(gr{lt}.psim.ws(ss,:)),2), '.w', 'MarkerSize',30)
plot(i_sl, abs(gr{lt}.sp.ws(ss,:)), '.r', 'MarkerSize',20)
ylabel('Windspeed'); xlim([0 numel(ss)+1]); xticks([])
nexttile; hold on
for i_ss= 1:numel(ss)
    i_s=ss(i_ss);
    fill(i_ss+AS(:,i_s)*scale_width,X_axis,colordergrad(i_s,:),'FaceAlpha',0.5)
end
plot([i_sl i_sl]', [quantile(abs(gr{lt}.psim.as(ss,:)),.25,2) quantile(abs(gr{lt}.psim.as(ss,:)),.75,2)]', '-w','linewidth',2)
plot(i_sl, mean(abs(gr{lt}.psim.as(ss,:)),2), '.w', 'MarkerSize',30)
plot(i_sl, abs(gr{lt}.sp.as(ss,:)), '.r', 'MarkerSize',20)
ylabel('Airspeed'); xlim([0 numel(ss)+1]); xticks([])
nexttile; hold on
for i_ss= 1:numel(ss)
    i_s=ss(i_ss);
    bar(i_ss,hours(sta{lt}.actEffort(i_s)),'facecolor',colordergrad(i_s,:),'FaceAlpha',0.5)
    plot(i_ss,hours(sta{lt}.actDuration(i_s)),'.k')
end
ylabel('Flight duration [hr]'); xlim([0 numel(ss)+1])
xticks(1:numel(ss)); xlabel('Sationary period')
xticklabels(ss+" | "+datestr(sta{lt}.end(ss),'dd-mmm'));%+"->"+datestr(sta{lt}.start(ss+1),'dd-mmm'))
% exportgraphics(gcf,['speed_hist_all_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '.png'],'Resolution',600)
close all
end

%% Airspeed vs Groundspeed
figure('position',[0 0 1600 900]); hold on; box on; grid on; axis equal;
for lt=1:height(tblLog)
    if isempty(gr{lt}), continue; end

    % sum(abs( (gr{lt}.psim.gs)-abs(gr{lt}.psim.as)) .* gr{lt}.actEffort(1:end-1)) % km/h *h -> km
    gsd=sum(abs(gr{lt}.psim.gs) .* gr{lt}.actEffort(1:end-1));
    asd=sum(abs(gr{lt}.psim.as) .* gr{lt}.actEffort(1:end-1));
    wsd=sum(abs(gr{lt}.psim.ws) .* gr{lt}.actEffort(1:end-1));
    
    gam = acos((wsd.^2+gsd.^2-asd.^2)./(2.*wsd.*gsd));
    
    awdx = wsd.*cos(gam);
    awdy = wsd.*sin(gam);
    
%     plot([zeros(1000,1) awdx' ]',[zeros(1000,1) awdy']','--','LineWidth',2,'color',tblLog.Color(lt,:)+.5*(1-tblLog.Color(lt,:)))
%     plot([awdx' gsd']',[awdy' zeros(1000,1)]','--','LineWidth',2,'color',tblLog.Color(lt,:)+.5*(1-tblLog.Color(lt,:)))
    plot([0 mean(awdx) ]',[0 mean(awdy)]','--','LineWidth',2,'color',tblLog.Color(lt,:))
    plot([mean(awdx) mean(gsd) ]',[ mean(awdy) 0]','LineWidth',2,'color',tblLog.Color(lt,:))
    plot( mean(gsd), 0,'.','MarkerSize',30,'color',tblLog.Color(lt,:))
    
    text(mean(awdx), mean(awdy), raw{lt}.GDL_ID,'HorizontalAlignment','center','VerticalAlignment','bottom')
end
axis tight;
xlabel('Total Distance traveled [km]'); ylabel('Drift [km]'); 
legend('Windspeed','Airspeed','Groundspeed')

exportgraphics(gcf,'speed_distance_GSS.png','Resolution',600)


for lt=1:height(tblLog)
    if isempty(gr{lt}), continue; end
    
    gsd=mean(gr{lt}.psim.gs .* gr{lt}.actEffort(1:end-1),2);
    asd=mean(gr{lt}.psim.as .* gr{lt}.actEffort(1:end-1),2);
    wsd=mean(gr{lt}.psim.ws .* gr{lt}.actEffort(1:end-1),2);
    
%     gsd=gr{lt}.sp.gs .* gr{lt}.actEffort(1:end-1);
%     asd=gr{lt}.psim.as .* gr{lt}.actEffort(1:end-1);
%     wsd=gr{lt}.psim.ws .* gr{lt}.actEffort(1:end-1);
    
    colordergrad=crameri('batlow',height(sta{lt}));
    pax = polaraxes; hold on;
    for i_s=1:height(sta{lt})-1
        if gr{lt}.actEffort(i_s)>3
            % polarplot(gr{lt}.psim.gs(i_s,:) .* gr{lt}.actEffort(i_s),'.','color',colordergrad(i_s,:),'markersize',3)
            polarplot([0 wsd(i_s)],'--','LineWidth',2,'Color',colordergrad(i_s,:)); 
            polarplot([wsd(i_s) wsd(i_s)+asd(i_s)],'LineWidth',2,'color',colordergrad(i_s,:))
            polarplot(gsd(i_s),'.','color',colordergrad(i_s,:),'markersize',30)
           text(angle(gsd(i_s)),abs(gsd(i_s))+100,string(i_s),'HorizontalAlignment','center')
        end
    end
    axis tight;
    pax.ThetaTick = 0:45:360;
    pax.ThetaTickLabel = {'E','NE','N','NW','W','SW','S','SE'};
    exportgraphics(gcf,['speed_polar_sp_3h_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '.png'],'Resolution',600)
    close all;
end

% Energy with duration
clear res

for lt=1:height(tblLog)
    if isempty(gr{lt}), continue; end
    
    bird = Bird(tblLog.CommonName{lt});
    Pmech_ms = mechanicalPower(bird);
    Pmech_ms = matlabFunction(Pmech_ms);
    Pmech = @(x) Pmech_ms(max(x*1000/60/60,5));
    P = Pmech(abs(gr{lt}.psim.as)).*gr{lt}.actEffort(1:end-1)*60*60; %in Joules (W=J/s)

    %Pgs = Pmech(abs(gr{lt}.psim.gs)).*gr{lt}.actEffort(1:end-1)*60*60; %in Joules (W=J/s)
    r = abs(gr{lt}.psim.gs) ./ abs(gr{lt}.psim.as);
    Pgsc = r.*Pmech(abs(gr{lt}.psim.as)).*gr{lt}.actEffort(1:end-1)*60*60;

    res(lt)=mean(sum(P))/mean(sum(Pgsc));
end
figure('position',[0 0 1000 400]); hold on; box on; grid on;
bar(res)
xticks(1:height(tblLog));
xticklabels(tblLog.GDL_ID)
yline(1,'k','LineWidth',2);
ylabel('Ratio of Energy spend with/without wind')

exportgraphics(gcf,'energy_ratio_wind.png','Resolution',600)

















%% Movement model 
lt=1;
grf = gr{lt};

w = [4 1 1/4];
for i_w=1:numel(w)
    grf.p = grf.ps .* grf.mvt_pdf(abs(grf.as).^w(i_w));
    sp{i_w} = shortestPathGraph(grf);
end

figure; tiledlayout('flow','TileSpacing','none','Padding','none')
nexttile; hold on;
xi=0:70;
for i_w=1:numel(w)
    f = grf.mvt_pdf(xi).^w(i_w);
    plot(xi,f./sum(f),'linewidth',2)
end
xlabel('Airspeed [km/h]'); ylabel('Probability');

nexttile; hold on
for i_w=1:numel(w)
    plot(grf.lon(sp{i_w}.lon), grf.lat(sp{i_w}.lat),'-o','linewidth',2)
end
borders('countries','color',[.6 .6 .6]);
axis equal; axis([3 max(grf.lon) min(grf.lat) max(grf.lat) ]);
xticklabels(''); yticklabels('')










%% Aggregation of prob
grf =gr{lt};

% weight of pressure, light and movement
wa = [1 1 1; 2 1 1];
wa=wa./sum(wa,2);
clear sp
for i_w=1:height(wa)
    prob_map = pres_prob{lt}(:,:,sta{lt}.staID).^wa(i_w,1) .* pres_thr{lt}(:,:,sta{lt}.staID) .* light_prob{lt}(:,:,sta{lt}.staID).^wa(i_w,2) .* ~mask_water{lt};
    grf.p = prob_map(grf.t) .* mvt_pdf(abs(grf.as)).^wa(i_w,3);
    sp{i_w} = shortestPathGraph(grf);
end
 
figure; hold on
for i_w=1:height(wa)
    plot(grf.lon(sp{i_w}.lon), grf.lat(sp{i_w}.lat),'-o','linewidth',2)
end
borders('countries','color',[.6 .6 .6]);
axis equal; axis([3 max(grf.lon) min(grf.lat) max(grf.lat) ]);
xticklabels(''); yticklabels('')










%% Graph output without light
for lt=1:height(tblLog)
    if isempty(gr{lt}), continue; end
    tic
    % remove light prob
    grt = gr{lt};
    % prob_map = pres_prob{lt}(:,:,sta{lt}.staID) .* pres_thr{lt}(:,:,sta{lt}.staID) .* light_prob{lt}(:,:,sta{lt}.staID) .* ~mask_water{lt};
    prob_map = pres_prob{lt}(:,:,sta{lt}.staID) .* pres_thr{lt}(:,:,sta{lt}.staID) .* ~mask_water{lt};
    prob_map = prob_map ./ sum(prob_map,[1 2]);
    % Set first and last one to true only at the known location.
    [~, tmp1] = min(abs(grt.lat(:)-raw{lt}.calib.lat));
    [~, tmp2] = min(abs(grt.lon(:)-raw{lt}.calib.lon));
    prob_map(:,:,1)=0;
    prob_map(tmp1,tmp2,1)=1;
    if ~any(isnat(raw{lt}.calib.second_period))
        prob_map(:,:,end)=0;
        prob_map(tmp1,tmp2,end)=1;
    end
    grt.p = prob_map(grt.t) .* grt.mvt_pdf(abs(grt.as));
    
    gr{lt}.sp_withoutlight = shortestPathGraph(grt);
    gr{lt}.M_withoutlight = probMapGraph(grt);

    path0=gr{lt}.sp.path;
    fixPath = false(size(path0));
    fixPath(1)=true;
    if ~isnat(tblLog.CalibSecondStart(lt))
        fixPath(end)=true;
    end
    gr{lt}.psim_withoutlight = gibbsGraph(nj,path0,fixPath,grt);
end

% Figure
colordergrad=crameri('batlow',5);
figure;tiledlayout('flow','TileSpacing','none','Padding','none')
for lt=1:height(tblLog)
    if isempty(gr{lt}), continue; end
    % figure('position',[0 0 900 1200], 'Name', [raw{lt}.GDL_ID ' | ' tblLog.CommonName{lt}] ); hold on
    %nexttile;
    set(gca,'Color',[.5 .5 .5])
    path_mean=nan(height(sta{lt}),2);
    borders('countries','edgecolor','w','facecolor','k')
    plot(raw{lt}.calib.lon,raw{lt}.calib.lat,'xr','linewidth',2)
    axis equal; axis([min(lon{lt}) max(lon{lt}) min(lat{lt}) max(lat{lt}) ]);
    dll=2; 
    %axis([min(gr{lt}.sp.lon)-dll max(gr{lt}.sp.lon)+dll min(gr{lt}.sp.lat)-2*dll max(gr{lt}.sp.lat)+dll ]);
    xticks([]);yticks([]); set(gca,'LooseInset',get(gca,'TightInset')); box on;
  
    %plot(gr{lt}.psim.lon(:,1:20:end),gr{lt}.psim.lat(:,1:20:end),'color',[.7 .7 .7].*colordergrad(2,:));
    plot(gr{lt}.sp.lon, gr{lt}.sp.lat,'-o','linewidth',2,'color',colordergrad(2,:))
    
    %plot(gr{lt}.psim.lon(:,1:20:end), gr{lt}.psim.lat(:,1:20:end),'color',[.7 .7 .7].*colordergrad(3,:));
    plot(gr{lt}.sp_withoutlight.lon, gr{lt}.sp_withoutlight.lat,'-o','linewidth',2,'color',colordergrad(3,:))
end

% Measure of distance 
ddist=cell(height(tblLog),1);
for lt=1:height(tblLog)
    if isempty(gr{lt}), continue; end
    ddist{lt}=nan(height(sta{lt}),6);
    for i_s=1:height(sta{lt})
        A=pdist([gr{lt}.psim.lon(i_s,:)' gr{lt}.psim.lat(i_s,:)'],@lldistkm);
        B=pdist([gr{lt}.psim_withoutlight.lon(i_s,:)' gr{lt}.psim_withoutlight.lat(i_s,:)'],@lldistkm);
        AB=pdist2([gr{lt}.psim.lon(i_s,:)' gr{lt}.psim.lat(i_s,:)'],[gr{lt}.psim_withoutlight.lon(i_s,:)' gr{lt}.psim_withoutlight.lat(i_s,:)'],@lldistkm);
        ddist{lt}(i_s,1)=mean(A);
        ddist{lt}(i_s,2)=std(A);
        ddist{lt}(i_s,3)=mean(B) ;
        ddist{lt}(i_s,4)=std(B);
        ddist{lt}(i_s,5)=mean(AB(:));
        ddist{lt}(i_s,6)=std(AB(:));
    end
end

% DL divergense
kldiv=cell(height(tblLog),1);
for lt=1:height(tblLog)
    if isempty(gr{lt}), continue; end
    P = gr{lt}.M ./ sum(gr{lt}.M,[1 2]);
    Q = gr{lt}.M_withoutlight ./ sum(gr{lt}.M_withoutlight,[1 2]);
    kldiv{lt}=squeeze(sum(P .* log(P./Q),[1 2],'omitnan'));
end

