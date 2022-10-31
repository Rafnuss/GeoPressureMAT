
addpath(genpath('../functions'))
project= "StudyPressure";
load("../data/processedData"+project+".mat")
load('../data/graph'+project+'.mat')
addpath(genpath("/Users/raphael/Documents/GitHub/Flight-Matlab/"))

for lt=1:height(tblLog)
    [gr{lt}.psim.lon, gr{lt}.psim.lat, ~] = path2lonlat(gr{lt}.psim.path,gr{lt});

%     G = digraph(gr{lt}.s,gr{lt}.t,gr{lt}.p);
%     path_edge = reshape(findedge(G,gr{lt}.psim.path(1:end-1,:),gr{lt}.psim.path(2:end,:)),gr{lt}.snds(3)-1,[]);
%     gr{lt}.psim.ws = gr{lt}.ws(path_edge);
%     gr{lt}.psim.gs = gr{lt}.gs(path_edge);
%     gr{lt}.psim.as = gr{lt}.as(path_edge);
end


%% Figure illustration graph

lt=1;
% Build the graph and save the graph status at all the steps
thr_prob_percentile = .99;
prob_map = pres_prob{lt}(:,:,sta{lt}.staID) .* pres_thr{lt}(:,:,sta{lt}.staID) .* light_prob{lt}(:,:,sta{lt}.staID) .* ~mask_water{lt};
[gr1,nds] = createGraph(prob_map,lat{lt},lon{lt},raw{lt}.calib,sta{lt}.actEffort,thr_prob_percentile,Inf);
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
col = brewermap(4,'YlOrRd'); colormap(col);
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




%% Main result figure

for lt=1:height(tblLog)
    if isempty(gr{lt}), continue; end
    figure('position',[0 0 900 1200], 'Name', [raw{lt}.GDL_ID ' | ' tblLog.CommonName{lt}] ); hold on
    
    set(gca,'Color',[.5 .5 .5])
    %imagesc(lon{lt},lat{lt},ones(numel(lat{lt}),numel(lon{lt}),3).*reshape([0 0 0],1,1,3),'AlphaData',~mask_water{lt});
    
    borders('countries','facecolor',"k",'edgecolor',[.5 .5 .5])
    plot(raw{lt}.calib.lon,raw{lt}.calib.lat,'xr','linewidth',2)
    axis equal; 
    % axis([min(lon{lt}) max(lon{lt}) min(lat{lt}) max(lat{lt}) ]);
    [lonext, latext]=path2lonlat(gr{lt}.s,gr{lt});
    axis([min(lonext) max(lonext) min(latext) max(latext) ]);
    set(gca,'ydir','normal');xticks([]);yticks([])
    set(gca,'LooseInset',get(gca,'TightInset'));
    
    colordergrad=crameri('batlow',height(sta_sm{lt}));
    t = sta_sm{lt}.start + (sta_sm{lt}.end-sta_sm{lt}.start)/2;
    
    colorinterp = interp1(datenum(t),colordergrad,datenum(sta_sm{lt}.start(1):sta_sm{lt}.end(end)));
    colormap(colorinterp)
    c=colorbar('south'); c.Color='w'; c.FontSize=12;
    c.Ticks=datenum(unique(dateshift(t(1):t(end),'start','month'))-t(1))/datenum(t(end)-t(1));
    c.TickLabels=datestr(unique(dateshift(t(1):t(end),'start','month')),'mmm');
       
    
    % gif(['graph_probMap_sp_sim_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '.gif'],'overwrite',true,'DelayTime',(5+height(sta{lt})/10)/height(sta{lt}),'frame',gca)
    
    idrd = randsample( size(gr{lt}.psim.lon,2) , 10 );
    % plot(gr{lt}.psim.lon(:,idrd), gr{lt}.psim.lat(:,idrd),'color',[.6 .6 .6]);
    
    % plot(mean(gr{lt}.psim.lon,2), mean(gr{lt}.psim.lat,2),'linewidth',10,'color',[.7 .7 .7]);

 
    for i_s = 1:height(sta_sm{lt})

         p3 = plot(gr{lt}.psim.lon(1:i_s,idrd), gr{lt}.psim.lat(1:i_s,idrd),'color',[.6 .6 .6]);

        f = gr{lt}.M(:,:,i_s);
        imagesc(lon{lt},lat{lt},ones(size(f,1),size(f,2),3).*reshape(colordergrad(i_s,:),1,1,3),'AlphaData',f./max(f(:)));%./max(f(:)));
        
       
        p=plot(gr{lt}.sp.lon(1:i_s),gr{lt}.sp.lat(1:i_s),'w','linewidth',3);
        
        for i_ss = 1:i_s
            sz = 20+(hours(sta_sm{lt}.end(i_ss)-sta_sm{lt}.start(i_ss)))/100;
            % sz = hours(sta{lt}.end(i_ss)-sta{lt}.start(i_ss))^(0.3) * 7;
            p2(i_ss,1) = plot(gr{lt}.sp.lon(i_ss),gr{lt}.sp.lat(i_ss),'.w', 'MarkerSize',sz+15);
            p2(i_ss,2) = plot(gr{lt}.sp.lon(i_ss),gr{lt}.sp.lat(i_ss),'.', 'MarkerSize',sz,'color',colordergrad(i_ss,:));
        end
        % keyboard
        
       
        % gif
        % tiadd = datenum([sta{lt}.start(i_s) sta{lt}.end(i_s)]-t(1))/datenum(t(end)-t(1));
        %pause(1)
        if i_s<height(sta_sm{lt})
            delete(p);   delete(p2);  delete(p3)
        end
    end

    dll=2; xticks([]); yticks([]); box on;
    axis([min(gr{lt}.sp.lon)-dll max(gr{lt}.sp.lon)+dll min(gr{lt}.sp.lat)-3*dll max(gr{lt}.sp.lat)+dll ]);


    exportgraphics(gca,['graph_probMap_sp_sim_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '.png'],'Resolution',300)
    % keyboard
    close all
end






















%% WindRose figure

vwinds = 0:10:100;

for lt=1:height(tblLog)
    if isempty(gr{lt}), continue; end
    ss = find(hours(sta_sm{lt}.actEffort)>3);

    figure('position',[0 0 900 1200], 'Name', [raw{lt}.GDL_ID ' | ' tblLog.CommonName{lt}] ); hold on
    set(gca,'LooseInset',get(gca,'TightInset'));
    borders('countries','facecolor',[.7 .7 .7],'edgecolor','w')
    plot(raw{lt}.calib.lon,raw{lt}.calib.lat,'xr','linewidth',2)

    plot(gr{lt}.sp.lon,gr{lt}.sp.lat,'color',[.3 .3 .3],'Linewidth',1)
    axis equal
    axis([min(gr{lt}.sp.lon)-dll max(gr{lt}.sp.lon)+dll min(gr{lt}.sp.lat)-2*dll max(gr{lt}.sp.lat)+dll ]);
   
    % vwinds = round(linspace(0,ceil(max(abs(gr{lt}.psim.ws(:)))/10)*10,9));
    for i_ss= 1:numel(ss)
        i_s=ss(i_ss);
        plot(gr{lt}.sp.lon(i_s:i_s+1),gr{lt}.sp.lat(i_s:i_s+1),'color','k','Linewidth',2)
        axis equal 
axis([min(gr{lt}.sp.lon)-dll max(gr{lt}.sp.lon)+dll min(gr{lt}.sp.lat)-2*dll max(gr{lt}.sp.lat)+dll ]);
        WindRose(rad2deg(angle(gr{lt}.psim.ws(i_s,:))), abs(gr{lt}.psim.ws(i_s,:)),'X',mean(gr{lt}.sp.lon(i_s:i_s+1)),...
            'Y',mean(gr{lt}.sp.lat(i_s:i_s+1)), 'LegendType', 0, 'labels','','TitleString','',...
        'scalefactor',2,'nfreq',0,'vWinds', vwinds(1:end-1),'nDirections',24, 'anglenorth',90,'angleeast',0);
        % text(mean(gr{lt}.sp.lon(i_s:i_s+1)),1+mean(gr{lt}.sp.lat(i_s:i_s+1)),string(i_s),'HorizontalAlignment','center','Color','white','FontWeight','bold')
        % text(mean(gr{lt}.sp.lon(i_s:i_s+1)),1+mean(gr{lt}.sp.lat(i_s:i_s+1)),string(i_s),'HorizontalAlignment','center')
     % axis equal; axis([min(lon{lt}) max(lon{lt}) min(lat{lt}) max(lat{lt}) ]);
        % keyboard
    end
    
    colormap(crameri('batlow'));
    scatter(gr{lt}.sp.lon,gr{lt}.sp.lat,100,1:height(sta_sm{lt}),'filled','MarkerEdgeColor','k')
   

    dll=2; xticks([]); yticks([]); box on;
     %  ax2=axes;ax2.Visible='off';
    %   colormap(ax2,parula(numel(vwinds)-2)); c=colorbar('Location','south'); c.Label.String='Windspeed [km/h]'; caxis([0 vwinds(end)]);
%      axis equal; set(gca,'LooseInset',get(gca,'TightInset')); title(' '); c.Ticks=vwinds;
%      axis([min(gr{lt}.sp.lon)-dll max(gr{lt}.sp.lon)+dll min(gr{lt}.sp.lat)-3*dll max(gr{lt}.sp.lat)+dll ]);
   axis([min(gr{lt}.sp.lon)-dll max(gr{lt}.sp.lon)+dll min(gr{lt}.sp.lat)-2*dll max(gr{lt}.sp.lat)+dll ]);
 
    % keyboard
    % exportgraphics(gcf,['graph_windrose_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '.eps'],'ContentType','vector')
    exportgraphics(gcf,['graph_windrose_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '.png'],'resolution',300)
    
    close all
end



%% Wind Triangle 
figure('position',[0 0 1600 900]); hold on; box on; grid on; axis equal;
for lt=1:height(tblLog)
    if isempty(gr{lt}), continue; end
         
        [~,~,St] = ind2sub(gr{lt}.snds,gr{lt}.s);

        Es = splitapply(@sum,gr{lt}.E,St);
        gsd = sum(splitapply(@sum,double(abs(gr{lt}.gs)) .* gr{lt}.E,St) ./ Es .* gr{lt}.actEffort(1:end-1));
        asd = sum(splitapply(@sum,double(abs(gr{lt}.as)) .* gr{lt}.E,St) ./ Es .* gr{lt}.actEffort(1:end-1));
        wsd = sum(splitapply(@sum,double(abs(gr{lt}.ws)) .* gr{lt}.E,St) ./ Es .* gr{lt}.actEffort(1:end-1));
        
        gam = acos((wsd.^2+gsd.^2-asd.^2)./(2.*wsd.*gsd));
        
        awdx = wsd.*cos(gam);
        awdy = (i*2-1)*wsd.*sin(gam);
        
    %     plot([zeros(1000,1) awdx' ]',[zeros(1000,1) awdy']','--','LineWidth',2,'color',tblLog.Color(lt,:)+.5*(1-tblLog.Color(lt,:)))
    %     plot([awdx' gsd']',[awdy' zeros(1000,1)]','--','LineWidth',2,'color',tblLog.Color(lt,:)+.5*(1-tblLog.Color(lt,:)))
        plot([0 mean(awdx) ]',[0 mean(awdy)]','--','LineWidth',2,'color',tblLog.Color(lt,:))
        plot([mean(awdx) mean(gsd) ]',[ mean(awdy) 0]','LineWidth',2,'color',tblLog.Color(lt,:))
        plot( mean(gsd), 0,'.','MarkerSize',30,'color',tblLog.Color(lt,:))
    text(mean(gsd), 0-150, raw{lt}.GDL_ID,'HorizontalAlignment','center','VerticalAlignment','bottom')
end
axis tight;
xlabel('Total Distance traveled [km]'); ylabel('Drift [km]'); 
legend('Windspeed','Airspeed','Groundspeed')
 
% exportgraphics(gcf,'speed_distance_GSS.eps','Resolution',600)



%% Appendix: Trajectory id
gr_sp = [];










%% Appendix: Histogram of speed
% Version simulation
for lt=1:height(tblLog)
    if isempty(gr{lt}), continue; end
    ss = find(hours(gr{lt}.actEffort)>0);
    
    x_axis = (0:.1:100)';
    gs = nan(numel(x_axis), numel(ss));
    as = gs; ws=gs;

    for i_ss= 1:numel(ss)
        i_s=ss(i_ss);
        gs(:,i_s) = ksdensity(abs(gr{lt}.psim.gs(i_s,:)),x_axis);
        ws(:,i_s) = ksdensity(abs(gr{lt}.psim.ws(i_s,:)),x_axis);
        as(:,i_s) = ksdensity(abs(gr{lt}.psim.as(i_s,:)),x_axis);
    end


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
    keyboard
    close all
end

%% Appendix: Histogram of speed
% Version marginal
x_axis = (0:.5:100)';
for lt=1:height(tblLog)

    ss = find(gr{lt}.actEffort>3);

    [~,~,St] = ind2sub(gr{lt}.snds,gr{lt}.s);

    x_axis = (0:.1:100)';
    gs = nan(numel(x_axis), numel(ss));
    as = gs; 
    ws=gs;

    for i_ss= 1:numel(ss)
        id = find(St==ss(i_ss));
        
        [tmp1, tmp2] = sort(abs(gr{lt}.gs(id)));

        stairs(tmp1,cumsum(gr{lt}.E(id(tmp2))))
    end
    
    speed = zeros(gr{lt}.snds(3),numel(x_axis)-1,3);

    

    Y=discretize(abs(gr{lt}.gs),x_axis);
    [G,Stg,Yg] = findgroups(St, Y);
    id = sub2ind(size(speed),Stg,Yg,ones(size(Yg)));
    speed(id)=splitapply(@sum,gr{lt}.E,G);

    Y=discretize(abs(gr{lt}.as),x_axis);
    [G,Stg,Yg] = findgroups(St, Y);
    id = sub2ind(size(speed),Stg,Yg,2*ones(size(Yg)));
    speed(id)=splitapply(@sum,gr{lt}.E,G);

    Y=discretize(abs(gr{lt}.ws),x_axis);
    [G,Stg,Yg] = findgroups(St, Y);
    id = sub2ind(size(speed),Stg,Yg,3*ones(size(Yg)));
    speed(id)=splitapply(@sum,gr{lt}.E,G);


    ss
   
end



%% Appendix: polarplot

for lt=1:height(tblLog)
    figure
    if isempty(gr{lt}), continue; end
    
    gsd=mean(gr{lt}.psim.gs .* gr{lt}.actEffort(1:end-1),2);
    asd=mean(gr{lt}.psim.as .* gr{lt}.actEffort(1:end-1),2);
    wsd=mean(gr{lt}.psim.ws .* gr{lt}.actEffort(1:end-1),2);
    
%     gsd=gr{lt}.sp.gs .* gr{lt}.actEffort(1:end-1);
%     asd=gr{lt}.psim.as .* gr{lt}.actEffort(1:end-1);
%     wsd=gr{lt}.psim.ws .* gr{lt}.actEffort(1:end-1);
    
    colordergrad=crameri('batlow',height(sta{lt}));
    pax = polaraxes; hold on;
    for i_s=1:height(sta_sm{lt})-1
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
   % exportgraphics(gcf,['speed_polar_sp_3h_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '.png'],'Resolution',600)
    close all;
end


%% Not working: Prefered direction
% for lt=1:height(tblLog)
%     if isempty(gr{lt}), continue; end
%     
%     im=find(sta_sm{lt}.status=="wintering");
%     for i=0:1
%         id = ((1:height(gr{lt}.psim.as))>=im) == i;
%         
%         if i==0
%             dest = mean(gr{lt}.psim.lon(im,:))-gr{lt}.psim.lon(id,:) + 1i*(mean(gr{lt}.psim.lat(im,:))-gr{lt}.psim.lat(id,:));
%         else
%             dest = mean(gr{lt}.psim.lon(end,:))-gr{lt}.psim.lon(id,:) + 1i*(mean(gr{lt}.psim.lat(end,:))-gr{lt}.psim.lat(id,:));
%         end
% 
%         S = gr{lt}.psim.ws(id,:);
%         a1=mean( abs(S).* (real(dest).*real(S)+imag(dest).*imag(S)) ./ abs(dest) ,2);
%         a2=mean( abs(S).* (real(dest).*imag(S)-imag(dest).*real(S)) ./ abs(dest) ,2);
%         
%         S = gr{lt}.psim.as(id,:);
%         b1=mean( abs(S).* (real(dest).*real(S)+imag(dest).*imag(S)) ./ abs(dest) ,2);
%         b2=mean( abs(S) .*(real(dest).*imag(S)-imag(dest).*real(S)) ./ abs(dest) ,2);
%         
%         figure;
%         hold on
%         % x= (1:sum(id))';
%         x=-mean(abs(dest),2)*100; 
%         plot((x+[zeros(sum(id),1) b1 ])',[zeros(sum(id),1) b2]','-k')
%         plot((x+[zeros(sum(id),1) a1 ])',[zeros(sum(id),1) a2]','-','color',[.7 .7 .7])
%         scatter(x,zeros(1,sum(id)),gr{lt}.actEffort(id)*10,'filled')
%         text(x,zeros(1,sum(id)), string(sta_sm{lt}.staID(id)))
%     end
% end



%% Not included: Energy with duration

figure('position',[0 0 1600 900]); 
tiledlayout('flow','TileSpacing','tight','Padding','tight')

for lt=1:height(tblLog)
    if isempty(gr{lt}), continue; end
    
    bird = Bird(tblLog.CommonName{lt});
    Pmech_ms = mechanicalPower(bird);
    Pmech_ms = matlabFunction(Pmech_ms);
    % Pmech = @(x) Pmech_ms(max(x*1000/60/60,5));
    Pmech = @(x) Pmech_ms(x*1000/60/60);

    nexttile; hold on; box on; grid on;
    as_x =0:100;
    [b,h,M,S] = histvdens(abs(gr{lt}.psim.as(:)), repmat(gr{lt}.actEffort(1:end-1),1000,1),as_x,true);
    bar(b,h)
    yyaxis right;
    plot(as_x,log(Pmech(as_x)),'linewidth',2)
end

E = nan(1000,height(tblLog),2);
fl_dur = nan(height(tblLog),2);
for lt=1:height(tblLog)
    if isempty(gr{lt}), continue; end
    
    bird = Bird(tblLog.CommonName{lt});
    Pmech_ms = mechanicalPower(bird);
    Pmech_ms = matlabFunction(Pmech_ms);
    Pmech = @(x) Pmech_ms(max(x*1000/60/60,5));
    % Pmech = @(x) Pmech_ms(x*1000/60/60);

    im=find(sta_sm{lt}.status=="wintering");
    for i=0:1
        id = ((1:height(gr{lt}.psim.gs))>=im) == i;
        E(:,lt,i+1) = sum(Pmech(abs(gr{lt}.psim.as(id,:))).*gr{lt}.actEffort(id)*60*60); %in Joules (W=J/s)
        fl_dur(lt,i+1) = sum(gr{lt}.actEffort(id));
    end
end

figure('position',[0 0 1000 400]); hold on; box on; grid on;
bar(squeeze(mean(E)))
xticks(1:height(tblLog));
xticklabels(tblLog.GDL_ID)
yline(1,'k','LineWidth',2);
ylabel('Ratio of Energy spend with/without wind')

%exportgraphics(gcf,'energy_ratio_wind.png','Resolution',600)

figure('position',[0 0 1000 400]); hold on; box on; grid on;
bar(fl_dur)
xticks(1:height(tblLog));
xticklabels(tblLog.GDL_ID)
yline(1,'k','LineWidth',2);
ylabel('Ratio of Energy spend with/without wind')









%% Appendix: Computational time table
trunt=cell2table(cellfun(@(x) struct2table(x),trun,'UniformOutput',false));
trunt = trunt.Var1;


trunT = table(tblLog.CommonName,'VariableNames',"CommonName");
trunT.GDL_ID = tblLog.GDL_ID;
trunT.n_grid = cellfun(@(x) numel(x.lat)*numel(x.lon), gr);
trunT.n_t = cellfun(@(x) x.snds(3), gr);
trunT.n_edge = cellfun(@(x) numel(x.s), gr);
trunT.step13=seconds(trunt.create_graph(:,2));
trunT.step13.Format="mm:ss";
trunT.step4=seconds(trunt.create_graph(:,3)-trunt.create_graph(:,2));
trunT.step4.Format="mm:ss";
trunT.create_graph=seconds(trunt.create_graph(:,4));
trunT.create_graph.Format="mm:ss";
trunT.shortestpath=seconds(trunt.shortestpath);
trunT.shortestpath.Format="mm:ss";
trunT.prob_map=seconds(trunt.prob_map);
trunT.prob_map.Format="mm:ss";
trunT.sim=seconds(trunt.sim);
trunT.sim.Format="mm:ss";


writetable(trunT,'runtime.csv')   



%% Appendix table of sta and simplified figure

for lt=1:height(tblLog)
    T = table(sta{lt}.staID, sta{lt}.start, sta{lt}.end, sta{lt}.actDuration, sta{lt}.twlNb,...
        'VariableNames',{'ID','start_time','end_time','next_flight_duration','twilight_count'});
    writetable(T,'sta_info.xls','Sheet',tblLog.GDL_ID{lt})
end

group{1} = 1:6;
group{2} = 7:10;
group{3} = 11:17;

for i =1:numel(group)

    figure('position',[0 0 900 1200]); hold on
        
    set(gca,'Color',[.5 .5 .5])
    %imagesc(lon{lt},lat{lt},ones(numel(lat{lt}),numel(lon{lt}),3).*reshape([0 0 0],1,1,3),'AlphaData',~mask_water{lt});
    
    borders('countries','facecolor',"k",'edgecolor',[.5 .5 .5])
    axis equal; 
    set(gca,'ydir','normal');xticks([]);yticks([])
    set(gca,'LooseInset',get(gca,'TightInset'));
    
    
    for lt=group{i}
        if isempty(gr{lt}), continue; end
    
        colordergrad=crameri('batlow', gr{lt}.snds(3));
        sz = 20+(hours(sta_sm{lt}.end-sta_sm{lt}.start))/100;
    
        % plot(gr{lt}.sp.lon,gr{lt}.sp.lat,'.w', 'MarkerSize',sz+15);
        plot(gr{lt}.sp.lon,gr{lt}.sp.lat,'color', tblLog.Color(lt,:))
        scatter(gr{lt}.sp.lon,gr{lt}.sp.lat,sz,colordergrad,'filled','MarkerEdgeColor','w');
            
        text(gr{lt}.sp.lon, gr{lt}.sp.lat,string(sta_sm{lt}.staID),'color','w',"HorizontalAlignment",'center','VerticalAlignment','top');
    end
end

























%% Sensitivity Analysis

%% Movement model 
lt=1;
grf = gr{lt};

w = [2 1 1/4];
for i_w=1:numel(w)
    grf.p = grf.ps .* grf.mvt_pdf(abs(grf.as).^w(i_w));
    sp{i_w} = shortestPathGraph(grf);
    M{i_w} = probMapGraph(grf);
end

figure;  hold on;
xi=0:70;
for i_w=1:numel(w)
    f = grf.mvt_pdf(xi).^w(i_w);
    plot(xi,f./sum(f),'linewidth',2)
end
xlabel('Airspeed [km/h]'); ylabel('Probability');

figure; hold on;
for i_w=1:numel(w)
    plot(sp{i_w}.lon, sp{i_w}.lat,'-o','linewidth',2)
end
borders('countries','color',[.6 .6 .6]);
axis equal; axis([3 max(grf.lon) min(grf.lat) max(grf.lat) ]);
xticklabels(''); yticklabels('')


figure('position',[0 0 900*3 1200], 'Name', [raw{lt}.GDL_ID ' | ' tblLog.CommonName{lt}] ); hold on
tiledlayout('flow','TileSpacing','tight','Padding','tight') 
for i_w=1:numel(w)
    nexttile;
    imagesc(lon{lt},lat{lt},ones(numel(lat{lt}),numel(lon{lt}),3).*reshape([0 0 0],1,1,3),'AlphaData',~mask_water{lt});
    set(gca,'Color',[.5 .5 .5])
    
    borders('countries','w')%'edgecolor','w','facecolor','k')
    plot(raw{lt}.calib.lon,raw{lt}.calib.lat,'xr','linewidth',2)
    axis equal; 
    % axis([min(lon{lt}) max(lon{lt}) min(lat{lt}) max(lat{lt}) ]);
    [lonext, latext]=path2lonlat(gr{lt}.s,gr{lt});
    axis([min(lonext) max(lonext) min(latext) max(latext) ]);
    set(gca,'ydir','normal');xticks([]);yticks([])
    % set(gca,'LooseInset',get(gca,'TightInset'));
    
    colordergrad=crameri('batlow',height(sta{lt}));
    t = sta{lt}.start + (sta{lt}.end-sta{lt}.start)/2;
    
    colorinterp = interp1(datenum(t),colordergrad,datenum(sta{lt}.start(1):sta{lt}.end(end)));
    colormap(colorinterp)
    c=colorbar('south'); c.Color='w'; c.FontSize=12;
    c.Ticks=datenum(unique(dateshift(t(1):t(end),'start','month'))-t(1))/datenum(t(end)-t(1));
    c.TickLabels=datestr(unique(dateshift(t(1):t(end),'start','month')),'mmm');
       
    
    
    %idrd = randsample( size(gr{lt}.psim.lon,2) , 30 );
    %plot(gr{lt}.psim.lon(:,idrd), gr{lt}.psim.lat(:,idrd),'color',[.7 .7 .7]);
    
    plot(sp{i_w}.lon,sp{i_w}.lat,'w','linewidth',3);

    for i_s = 1:height(sta{lt})
        f=M{i_w}(:,:,i_s);
        imagesc(lon{lt},lat{lt},ones(size(f,1),size(f,2),3).*reshape(colordergrad(i_s,:),1,1,3),'AlphaData',f./max(f(:)));%./max(f(:)));
    end
    for i_s = 1:height(sta{lt})
        sz = 20+(hours(sta{lt}.end(i_s)-sta{lt}.start(i_s)))/100;
        p2(i_s,1) = plot(sp{i_w}.lon(i_s), sp{i_w}.lat(i_s),'.w', 'MarkerSize',sz+15);
        p2(i_s,2) = plot(sp{i_w}.lon(i_s), sp{i_w}.lat(i_s),'.', 'MarkerSize',sz,'color',colordergrad(i_s,:));
    end
end









%% Light
grf = gr{lt};

% weight of pressure, light and movement
w = [0 1 2];
clear sp M
for i_w=1:numel(w)
    prob_map = pres_prob{lt}(:,:,sta_sm{lt}.staID) .* pres_thr{lt}(:,:,sta_sm{lt}.staID) .* light_prob{lt}(:,:,sta_sm{lt}.staID).^w(i_w) .* ~mask_water{lt};
    prob_map = prob_map ./ sum(prob_map,[1 2]);
    grf.p = prob_map(grf.t) .* grf.mvt_pdf(abs(grf.as));
    sp{i_w} = shortestPathGraph(grf);
    M{i_w} = probMapGraph(grf);
end
 
figure('position',[0 0 900*3 1200], 'Name', [raw{lt}.GDL_ID ' | ' tblLog.CommonName{lt}] ); hold on
tiledlayout('flow','TileSpacing','tight','Padding','tight') 
for i_w=1:numel(w)
    nexttile;
    imagesc(lon{lt},lat{lt},ones(numel(lat{lt}),numel(lon{lt}),3).*reshape([0 0 0],1,1,3),'AlphaData',~mask_water{lt});
    set(gca,'Color',[.5 .5 .5])
    
    borders('countries','w')%'edgecolor','w','facecolor','k')
    plot(raw{lt}.calib.lon,raw{lt}.calib.lat,'xr','linewidth',2)
    axis equal; 
    % axis([min(lon{lt}) max(lon{lt}) min(lat{lt}) max(lat{lt}) ]);
    [lonext, latext]=path2lonlat(gr{lt}.s,gr{lt});
    axis([min(lonext) max(lonext) min(latext) max(latext) ]);
    set(gca,'ydir','normal');xticks([]);yticks([])
    % set(gca,'LooseInset',get(gca,'TightInset'));
    
    colordergrad=crameri('batlow',height(sta{lt}));
    t = sta{lt}.start + (sta{lt}.end-sta{lt}.start)/2;
    
    colorinterp = interp1(datenum(t),colordergrad,datenum(sta{lt}.start(1):sta{lt}.end(end)));
    colormap(colorinterp)
    c=colorbar('south'); c.Color='w'; c.FontSize=12;
    c.Ticks=datenum(unique(dateshift(t(1):t(end),'start','month'))-t(1))/datenum(t(end)-t(1));
    c.TickLabels=datestr(unique(dateshift(t(1):t(end),'start','month')),'mmm');
       
    
    
    %idrd = randsample( size(gr{lt}.psim.lon,2) , 30 );
    %plot(gr{lt}.psim.lon(:,idrd), gr{lt}.psim.lat(:,idrd),'color',[.7 .7 .7]);
    
    plot(sp{i_w}.lon,sp{i_w}.lat,'w','linewidth',3);

    for i_s = 1:height(sta{lt})
        f=M{i_w}(:,:,i_s);
        imagesc(lon{lt},lat{lt},ones(size(f,1),size(f,2),3).*reshape(colordergrad(i_s,:),1,1,3),'AlphaData',f./max(f(:)));%./max(f(:)));
    end
    for i_s = 1:height(sta{lt})
        sz = 20+(hours(sta{lt}.end(i_s)-sta{lt}.start(i_s)))/100;
        p2(i_s,1) = plot(sp{i_w}.lon(i_s), sp{i_w}.lat(i_s),'.w', 'MarkerSize',sz+15);
        p2(i_s,2) = plot(sp{i_w}.lon(i_s), sp{i_w}.lat(i_s),'.', 'MarkerSize',sz,'color',colordergrad(i_s,:));
    end
end











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
