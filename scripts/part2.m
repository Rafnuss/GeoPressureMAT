
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

mat = nan(grt.snds);
mat(unique(gr1.s))=1;
mat(unique(gr2.s))=2;
mat(unique(gr3.s))=3;
mat(unique(gr4.s))=4;

figure('position',[0 0 1400 900]); colormap(col)
tiledlayout('flow','TileSpacing','none','Padding','none')
col = brewermap(4,'YlOrRd');
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
grf =gr{lt};
[Slat,Slon,St]=ind2sub(grf.snds,grf.s);
[Tlat,Tlon,Tt]=ind2sub(grf.snds,grf.t);
nrow=4;
figure('position',[0 0 1400 900]); colormap(col)
tiledlayout('flow','TileSpacing','none','Padding','none')
for i_s = 1:find(sta{lt}.status=="wintering")
    nexttile; hold on;
    st_id = find(St==i_s);
    borders('countries','color',[.6 .6 .6]);
    X = [grf.lon(Slon(st_id)) grf.lat(Slat(st_id))]; X=unique(X ,"rows");
    plot(X(:,1),X(:,2),'.k')
    X = [grf.lon(Tlon(st_id)) grf.lat(Tlat(st_id))]; X=unique(X ,"rows");
    plot(X(:,1),X(:,2),'.r')
    axis equal; axis([3 max(grf.lon) min(grf.lat) max(grf.lat) ]);
    xticklabels(''); yticklabels('')
end



%% Movement model 
grf =gr{lt};
mvt_pdf = movementModel('energy',tblLog.mass(lt),tblLog.wingSpan(lt));

w = [2 1 1/2 1/4];
for i_w=1:numel(w)
    grf.p = grf.ps .* mvt_pdf(abs(grf.as).^w(i_w));
    sp{i_w} = shortestPathGraph(grf);
end

figure; tiledlayout('flow','TileSpacing','none','Padding','none')
nexttile; hold on;
xi=0:70;
for i_w=1:numel(w)
    f = mvt_pdf(xi).^w(i_w);
    plot(xi,f./sum(f))
end
xlabel('Airspeed [km/h]'); ylabel('Probability');

nexttile; hold on
for i_w=1:numel(w)
    plot(grf.lon(sp{i_w}.lon), grf.lat(sp{i_w}.lat),'-o')
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
 
figure;nexttile; hold on
for i_w=1:numel(w)
    plot(grf.lon(sp{i_w}.lon), grf.lat(sp{i_w}.lat),'-o')
end
borders('countries','color',[.6 .6 .6]);
axis equal; axis([3 max(grf.lon) min(grf.lat) max(grf.lat) ]);
xticklabels(''); yticklabels('')



%% Histogram of speed
grf =gr{lt};
[Slat,Slon,St]=ind2sub(grf.snds,grf.s);
[Tlat,Tlon,Tt]=ind2sub(grf.snds,grf.t);

figure('position',[0 0 2400 1000]);
tiledlayout('flow','TileSpacing','none','Padding','none')
for i_s = 1:height(sta{lt})-1
    nexttile; hold on
    histogram(abs(grf.gs(St==i_s)),0:5:100);
    histogram(abs(grf.as(St==i_s)),0:5:50);
    histogram(abs(grf.ws(St==i_s)),0:5:50);
    xlim([0 100])
    legend('GS','AS','WS')
end






%% Map of posteriori pdf
grf = gr{lt};
figure;
tiledlayout('flow','TileSpacing','none','Padding','none')
for i_s = 1:grf.snds(3)
    nexttile(i_s);
    hold on; set(gca,'Color','k')
    borders('countries','w')
    imagesc(grf.lon,grf.lat,grf.M(:,:,i_s),'alphadata',~grf.M(:,:,i_s)==0)
    axis equal; axis([min(grf.lon) max(grf.lon) min(grf.lat) max(grf.lat) ]);
    scatter(grf.lon(grf.sp.lon(i_s)), grf.lat(grf.sp.lat(i_s)),'filled','r','MarkerEdgeColor','k')
end
colormap("pink") % copper, bone


%% 
for lt=1:height(tblLog)
    
    figure('position',[0 0 900 1200], 'Name', [raw{lt}.GDL_ID ' | ' tblLog.CommonName{lt}] ); hold on
    
   
    set(gca,'Color',[.5 .5 .5])
    [gLON,gLAT] = meshgrid(lon{lt},lat{lt}); 
    path_mean=nan(height(sta_sm{lt}),2);
    imagesc(lon{lt},lat{lt},ones(numel(lat{lt}),numel(lon{lt}),3).*reshape([0 0 0],1,1,3),'AlphaData',~mask_water{lt});
    
    borders('countries','w')
    plot(raw{lt}.calib.lon,raw{lt}.calib.lat,'xr','linewidth',2)
    axis equal; axis([min(lon{lt}) max(lon{lt}) min(lat{lt}) max(lat{lt}) ]);
    set(gca,'ydir','normal');xticks([]);yticks([])
    set(gca,'LooseInset',get(gca,'TightInset'));
    
    colordergrad=crameri('batlow',height(sta_sm{lt}));
    t = sta_sm{lt}.start + (sta_sm{lt}.end-sta_sm{lt}.start)/2;
    
    colorinterp = interp1(datenum(t),colordergrad,datenum(sta_sm{lt}.start(1):sta_sm{lt}.end(end)));
    colormap(colorinterp)
    c=colorbar('south'); c.Color='w'; c.FontSize=12;
    c.Ticks=datenum(unique(dateshift(t(1):t(end),'start','month'))-t(1))/datenum(t(end)-t(1));
    c.TickLabels=datestr(unique(dateshift(t(1):t(end),'start','month')),'mmm');
       
    
    % gif(['combined_map_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '.gif'],'overwrite',true,'DelayTime',(5+height(sta_sm{lt})/10)/height(sta_sm{lt}),'frame',gca)
    
    path_short = [lon{lt}(gr{lt}.sp.lon) lat{lt}(gr{lt}.sp.lat)];

    for i_s = 1:height(sta_sm{lt})
        f=gr{lt}.M(:,:,i_s);
        imagesc(lon{lt},lat{lt},ones(size(f,1),size(f,2),3).*reshape(colordergrad(i_s,:),1,1,3),'AlphaData',f./max(f(:)));
        
       
        p=plot(path_short(1:i_s,1),path_short(1:i_s,2),'w','linewidth',2);
        
        for i_ss = 1:i_s
            sz = 20+(hours(sta_sm{lt}.end(i_ss)-sta_sm{lt}.start(i_ss)))/100;
            p2(i_ss,1) = plot(path_short(i_ss,1),path_short(i_ss,2),'.w', 'MarkerSize',sz+15);
            p2(i_ss,2) = plot(path_short(i_ss,1),path_short(i_ss,2),'.', 'MarkerSize',sz,'color',colordergrad(i_ss,:));
        end
        % keyboard
        
        % gif
        % tiadd = datenum([sta_sm{lt}.start(i_s) sta_sm{lt}.end(i_s)]-t(1))/datenum(t(end)-t(1));
        %pause(1)
        if i_s<height(sta_sm{lt})
            delete(p);   delete(p2)
        end
    end
    
 
    exportgraphics(gca,['graph_probMap_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '.png'])
    % keyboard
    %close all
end




%% Wind
thr_dist=200; % km

for lt=1:height(tblLog)
    [Slat,Slon,St]=ind2sub(gr{lt}.snds,gr{lt}.s);
    [Tlat,Tlon,Tt]=ind2sub(gr{lt}.snds,gr{lt}.t);


    path_short = [lon{lt}(gr{lt}.sp.lon) lat{lt}(gr{lt}.sp.lat)];
    path_short_dist = lldistkm(path_short(1:end-1,:),path_short(2:end,:));

    figure('position',[0 0 900 1200], 'Name', [raw{lt}.GDL_ID ' | ' tblLog.CommonName{lt}] ); hold on
    borders('countries','k')
    plot(raw{lt}.calib.lon,raw{lt}.calib.lat,'xr','linewidth',2)
    
    plot(path_short(:,1),path_short(:,2),'color',[.3 .3 .3],'Linewidth',1)
   
    ss = find(path_short_dist>thr_dist);
    for i_ss= 1:numel(ss)
        i_s=ss(i_ss);
        plot(path_short(i_s:i_s+1,1),path_short(i_s:i_s+1,2),'color','k','Linewidth',2)
        st_id = find(St==i_s);
        WindRose(rad2deg(angle(gr{lt}.ws(st_id))), abs(gr{lt}.ws(st_id)),'X',mean(path_short(i_s:i_s+1,1)),'Y',mean(path_short(i_s:i_s+1,2)), 'LegendType', 0, 'labels','','scalefactor',2,'nfreq',1,'vWinds', [0 5 10 20 30 40]);
        % axis equal; axis([min(lon{lt}) max(lon{lt}) min(lat{lt}) max(lat{lt}) ]);
        % keyboard
    end
    
    axis equal; axis([min(lon{lt}) max(lon{lt}) min(lat{lt}) max(lat{lt}) ]);
    colormap(crameri('batlow'));
    scatter(path_short(:,1),path_short(:,2),100,1:height(sta_sm{lt}),'filled','MarkerEdgeColor','k')
    exportgraphics(gca,['graph_windspeed_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '.png'])
end
 










%% Gibbs Sampling

G = digraph(S2, T2, P);

n=10000;
% Initiate the simulated nodes (path
sim_node = nan(height(sta),n);
sim_node(1,:) = find(nds,1);
sim_node(end,:) = find(nds,1,'last');
% Set first path to the shortest
sim_node(:,1) = p;
% Initiaze the simulated edges
sim_edge = nan(height(sta)-1,n);
tic
for i_p=2:n
    for i_s=1:height(sta)-2
        % Edges and Nodes going OUT from the PREVIOUS nodes sampled (current path)
        [eidOut,nidOut] = outedges(G,sim_node(i_s,i_p));
        % Edge and Nodes going IN the NEXT nodes from the previous path
        [eidIn,nidIn] = inedges(G,sim_node(i_s+2,i_p-1));
        % Find nodes that are presents in both (going out from previous and coming in in the next)
        [iOut, iIn] = ismember(nidOut, nidIn);
        
        % Possible nodes
        nid = nidOut(iOut);
        eid = eidIn(iIn(iIn~=0));
        % Probability of the possible nodes (products of edge in and edge out
        prob = G.Edges.Weight(eidOut(iOut)) .* G.Edges.Weight(eid);
        
        % Radom sample the nodes
        ids = randsample(numel(nid),1,true,prob);
        sim_node(i_s+1,i_p) = nid(ids);
        sim_edge(i_s,i_p) = eid(ids);
        assert(~isnan(nid(ids)))
    end
    sim_edge(i_s+1,i_p) = eidOut(iOut(ids));
end
toc

    [Px,Py,Pt] = ind2sub(size(nds),sim_node);


    figure; hold on;
    borders('countries','k')
    axis equal; axis([[15 35],[-30 15]]);
    for i_p=1:numel(p)
        plot(glon(Py(:,i_p)), glat(Px(:,i_p)),'k')
        scatter(glon(Py(:,i_p)), glat(Px(:,i_p)),10*log(sta.twlNb)+10,Pt(:,i_p),'filled')
    end


figure; hold on;
borders('countries','k')
axis equal; axis([[15 35],[-30 15]]);
plot(quantile(glon(Py),.5,2), quantile(glat(Px),.5,2),'k')
for i_s=1:height(sta)
    plot([quantile(glon(Py(i_s,:)),.05) quantile(glon(Py(i_s,:)),.95)], [quantile(glat(Px(i_s,:)),.5) quantile(glat(Px(i_s,:)),.5)],'k','linewidth',2)
    plot([quantile(glon(Py(i_s,:)),.5) quantile(glon(Py(i_s,:)),.5)], [quantile(glat(Px(i_s,:)),.05) quantile(glat(Px(i_s,:)),.95)],'k','linewidth',2)
end
scatter(quantile(glon(Py),.5,2), quantile(glat(Px),.5,2),50*log(sta.twlNb)+10,Pt(:,1),'filled')


figure;  boxplot(AS2(sim_edge(:,2:end))');

%% 

figure('position',[0 0 2400 1000]);
for i_s = 1:height(sta)

    subplot(3,ceil(height(sta)/3),i_s); hold on;
    id_tmp = sub2ind(numel(glat),numel(glon),Px(i_s,:),Py(i_s,:));
    tmp = histcounts(id_tmp,.5:(numel(gLAT)+.5));
    tmp = reshape(tmp,numel(glat),numel(glon));

    
     imagesc(glon,glat,tmp)
    
    colorbar;
    
    borders('countries','w')
    tt = num2str(sta.staID(i_s));
    tt = tt + "|" + num2str(sta.twlNb(i_s))+"twls";
    tt = tt + "|" + datestr(sta.end(i_s),'dd-mmm-yyyy');
    tt = tt + "|" + num2str(round(hours(sta.actDuration(max(1,i_s))))) + "hr";
    title(tt)
     axis equal; axis([min(glon) max(glon) min(glat) max(glat) ]);
end


