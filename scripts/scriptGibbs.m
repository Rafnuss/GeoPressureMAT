
addpath(genpath('../functions'))
addpath('~/Documents/GitHub/Flight-Matlab/functions/')
addpath('~/Documents/GitHub/Flight-Matlab/data/')
load("../data/processedDataStudyPressure.mat")

% skip_gdl = {'22KT','24FF', '20IK'};

sta_sm=cell(1,height(tblLog));
for lt=1:height(tblLog)
    grp_id = hours(sta{lt}.end-sta{lt}.start)>12;%sta{lt}.twlNb>=4;
    grp_id(1) = true;
    if ~isnat(tblLog.CalibSecondStart(lt))
        grp_id(end) = true;
    end
    sta_sm{lt} = sta{lt}(grp_id,:);
    sta_sm{lt}.actNb =  splitapply(@sum, sta{lt}.actNb,cumsum(grp_id));
    sta_sm{lt}.actEffort =  splitapply(@sum, sta{lt}.actEffort,cumsum(grp_id));
    sta_sm{lt}.actDuration =  splitapply(@sum, sta{lt}.actDuration,cumsum(grp_id));
    sta_sm{lt}.twlNbStopover =  splitapply(@sum, sta{lt}.twlNb,cumsum(grp_id))-sta_sm{lt}.twlNb;
    sta_sm{lt}.staID = find(grp_id);
end

%% Gibbs Sampling

nj = 5000; % number of iteration for each path
np = 1; %number of different starting path
path = cell(height(tblLog),1);
prob_path = cell(height(tblLog),1);
trun = cell(height(tblLog),1);

for lt=1:height(tblLog)
    % Display Name
    disp(raw{lt}.GDL_ID)
    tic;
    % Map of probability
    prob_map = pres_prob{lt}(:,:,sta_sm{lt}.staID) .* pres_thr{lt}(:,:,sta_sm{lt}.staID) .* light_prob{lt}(:,:,sta_sm{lt}.staID).^1/10;
    
    prob_map_r = reshape(prob_map,[], height(sta_sm{lt}));
    % set to 0 onver water
    prob_map_r(mask_water{lt}(:),:)=0;
    % Normalized as each 
    prob_map_r = prob_map_r ./ nansum(prob_map_r,1); 

    [gLON,gLAT] = meshgrid(lon{lt},lat{lt});
    % Ddist = pdist([gLAT(~mask_water{lt}) gLON(~mask_water{lt})],@lldistkm);
    mvt_pdf = movementModel('gam');
    prob_mvt = @(pt,i_s,id) mvt_pdf(lldistkm([gLAT(pt) gLON(pt)],[gLAT(id) gLON(id)]) ./ hours(sta_sm{lt}.actEffort(i_s)));    
    
    % Define initial and fixed path
    % intial path as the maximum of probability excluding water
    [~,path0] = max(prob_map_r',[],2);
    % Expect on first and last location 
    [~,ll_calib] = min((raw{lt}.calib.lon-gLON(:)).^2+(raw{lt}.calib.lat-gLAT(:)).^2);

    fixPath = false(size(path0));
    if isnat(tblLog.CalibSecondStart(lt))
        path0(1)=ll_calib;
        fixPath(1)=true;
    else
        path0([1 end])=ll_calib;
        fixPath([1 end])=true;
    end
    
    % interpolate position for short stationary period
    path0(hours(sta_sm{lt}.end-sta_sm{lt}.start)<=24)=nan;
    [path_lat,path_lon] = ind2sub(size(gLAT),path0);
    path_lat0 = max(1,min(numel(lat{lt}),round(fillmissing(path_lat,'linear'))));
    path_lon0 = max(1,min(numel(lon{lt}),round(fillmissing(path_lon,'linear'))));

    path{lt} = nan(height(sta_sm{lt}),nj*np);
    prob_path{lt} = nan(1,nj*np);
    trun{lt}.prepartion = toc;
    tic
    for i_p=1:np
        path_lat=path_lat0;
        path_lon=path_lon0;

        % Add some noise on the initial path (bounded by the grid size
        %path_lat(2:end-1) = min(max(path_lat0(2:end-1)+round(randn(numel(path_lat0)-2,1)*2),1),size(gLAT,1));
        %path_lon(2:end-1) = min(max(path_lon0(2:end-1)+round(randn(numel(path_lon0)-2,1)*2),1),size(gLAT,2));
    
        path0p = sub2ind(size(gLAT),path_lat,path_lon);

        % figure; hold on; borders('countries','k');  plot(lon{lt}(path_lon0), lat{lt}(path_lat0),'or'); plot(gLON(path0p), gLAT(path0p),'r')
        
        [path{lt}(:,(i_p-1)*nj + (1:nj)'), prob_path{lt}((i_p-1)*nj + (1:nj)')] = Gibbs(nj,path0p,fixPath,prob_map_r,prob_mvt);
    end
    trun{lt}.run = toc;
    trun{lt}.nj = nj;
end

%save('../data/gibbs'+project+'.mat','path','prob_path','trun','-v7.3')
load('../data/gibbs'+project+'.mat')

%% Figure
% 
% 
% figure('position',[0 0 2400 1000], 'Name', [raw{lt}.GDL_ID ' | ' tblLog.CommonName{lt}] );
% %[ha, pos] = tight_subplot(2,ceil(height(sta_sm{lt})/2),[.03 .01],[.01 .03],.01);
% [ha, pos] = tight_subplot(3,ceil(height(sta_sm{lt})/3),[.03 .01],[.01 .03],.01);
% %[ha, pos] = tight_subplot(4,ceil(height(sta_sm{lt})/4),[.03 .01],[.01 .03],.01);
% [gLON,gLAT] = meshgrid(lon{lt},lat{lt});
% for i_s = 1:height(sta_sm{lt})
%     axes(ha(i_s)); hold on; set(gca,'Color','k')
%     [B,BG] = groupcounts(path{lt}(i_s,:)');
%     f=zeros(size(gLON));
%     f(BG)=B;
%     imagesc(lon{lt},lat{lt},f,'AlphaData',~mask_water{lt}); 
%     [~,id_max]=max(B);
%     
%     borders('countries','w')
%     plot(raw{lt}.calib.lon,raw{lt}.calib.lat,'xr')
%     % plot(gLON(path0(i_s)),gLAT(path0(i_s)),'or')
%     plot(gLON(BG(id_max)), gLAT(BG(id_max)),'or')
% 
%     tt = num2str(sta_sm{lt}.twlNb(i_s))+"twls" + " | ";
%     tt = tt + datestr(sta_sm{lt}.end(i_s),'dd-mmm');
%     if i_s>1
%         tt = tt + " | "  + num2str(round(hours(sta_sm{lt}.actEffort(max(1,i_s-1))))) + "hr"+ " | " ;
%     end
%     title(tt)
%     axis equal; axis([min(lon{lt}) max(lon{lt}) min(lat{lt}) max(lat{lt}) ]);
%     colormap(pink)
% end

%% Figure illustraion

lt=2;

[gLON,gLAT] = meshgrid(lon{lt},lat{lt});
mvt_pdf = movementModel('gam');
    
figure('position',[0 0 1200 750] );
tiledlayout('flow','TileSpacing','tight','Padding','tight')

[~,id]=min(prob_path{lt});
tmp_path = path{lt}(:,id);

i_s = 6;

nexttile; hold on;

tmp_p = pres_prob{lt}(:,:,sta_sm{lt}.staID(i_s));
tmp_pt = pres_thr{lt}(:,:,sta_sm{lt}.staID(i_s));
tmp_l = light_prob{lt}(:,:,sta_sm{lt}.staID(i_s));

tmp_mvt_prev = reshape(mvt_pdf(lldistkm([gLAT(tmp_path(i_s-1)) gLON(tmp_path(i_s-1))],[gLAT(:) gLON(:)]) ./ hours(sta_sm{lt}.actEffort(i_s-1))),size(tmp_l));    
tmp_mvt_next = reshape(mvt_pdf(lldistkm([gLAT(tmp_path(i_s+1)) gLON(tmp_path(i_s+1))],[gLAT(:) gLON(:)]) ./ hours(sta_sm{lt}.actEffort(i_s))),size(tmp_l));    

tmp_c = tmp_p.*tmp_pt.*tmp_l.*tmp_mvt_prev.*tmp_mvt_next;
path_i_s_sim = randsample(numel(tmp_c),20,true,tmp_c(:));

mask = 0.3+0.7*double(tmp_pt);
mask(mask_water{lt}) = 0;
img_tmp = real2rgb(tmp_p,colormap);
img_tmp = img_tmp.*mask;

imagesc(lon{lt},lat{lt},img_tmp)
borders('countries','w')
c_axis=caxis();
contour(lon{lt},lat{lt},tmp_l,2,'color',[255, 252, 49]/255,'linewidth',1);
contour(lon{lt},lat{lt},tmp_mvt_prev,1,'color','g','linewidth',1);
contour(lon{lt},lat{lt},tmp_mvt_next,1,'color','g','linewidth',1);

for i=1:10
    plot([gLON(tmp_path(i_s-1)) gLON(path_i_s_sim(i)) gLON(tmp_path(i_s+1))],[gLAT(tmp_path(i_s-1)) gLAT(path_i_s_sim(i)) gLAT(tmp_path(i_s+1))],'-w','linewidth',.7)
    plot(gLON(path_i_s_sim(i)),gLAT(path_i_s_sim(i)),'ok','MarkerSize',5,'markerFaceColor','w')
end

x=gLON(tmp_path); y=gLAT(tmp_path); x(i_s)=nan; y(i_s)=nan;
plot(x,y,'-w','linewidth',2)
plot(x,y,'ok','MarkerSize',10,'markerFaceColor','w')

plot(gLON(tmp_path(i_s-1)),gLAT(tmp_path(i_s-1)),'ok','MarkerSize',10,'markerFaceColor','g')
plot(gLON(tmp_path(i_s+1)),gLAT(tmp_path(i_s+1)),'ok','MarkerSize',10,'markerFaceColor','g')

plot(raw{lt}.calib.lon,raw{lt}.calib.lat,'ow','MarkerSize',12,'markerFaceColor','r')


axis equal; axis([min(lon{lt}) max(lon{lt}) min(lat{lt}) max(lat{lt}) ]);
set(gca,'ydir','normal');
yticks([]);xticks([])



%% Convert Path to matrix
probMapGibbs=cell(height(tblLog),1);

for lt=1:height(tblLog)
    probMapGibbs{lt}=nan(numel(lon{lt}),numel(lat{lt}),height(sta_sm{lt}));
    for i_s = 1:height(sta_sm{lt})
        [B,BG] = groupcounts(path{lt}(i_s,:)');
        f=zeros(size(probMapGibbs{lt}(:,:,i_s)));
        f(BG)=B;
        probMapGibbs{lt}(:,:,i_s)=f;
    end
end


%% 

for lt=1:height(tblLog)
    
    figure('position',[0 0 900 1200], 'Name', [raw{lt}.GDL_ID ' | ' tblLog.CommonName{lt}] ); hold on
    
   
    set(gca,'Color',[.5 .5 .5])
    [gLON,gLAT] = meshgrid(lon{lt},lat{lt}); 
    path_mean=nan(height(sta_sm{lt}),2);
    [~,id]=maxk(prob_path{lt},3);
    path_min = [gLON(path{lt}(:,id(2))) gLAT(path{lt}(:,id(2)))];
    imagesc(lon{lt},lat{lt},ones(numel(lat{lt}),numel(lon{lt}),3).*reshape([0 0 0],1,1,3),'AlphaData',~mask_water{lt});
    
    borders('countries','w')
    plot(raw{lt}.calib.lon,raw{lt}.calib.lat,'xr','linewidth',2)
    axis equal; 
    % axis([min(lon{lt}) max(lon{lt}) min(lat{lt}) max(lat{lt}) ]);
    tmp = sum(probMapGibbs{lt},3)>0;
    lonext = sort(gLON(tmp));
    latext = sort(gLAT(tmp));
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
       
    idrd = randsample( size(path{lt},2) , 30 );
    plot(gLON(path{lt}(:,idrd)), gLAT(path{lt}(:,idrd)),'color',[.7 .7 .7]);
    
    
    %gif(['combined_map_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '.gif'],'overwrite',true,'DelayTime',(5+height(sta_sm{lt})/10)/height(sta_sm{lt}),'frame',gca)
    
    for i_s = 1:height(sta_sm{lt})
        [B,BG] = groupcounts(path{lt}(i_s,:)');
        %f(BG)=f(BG)+B;
        f=zeros(size(gLON));
        f(BG)=B;
        imagesc(lon{lt},lat{lt},ones(size(f,1),size(f,2),3).*reshape(colordergrad(i_s,:),1,1,3),'AlphaData',f./max(f(:)));
        [~,id_max]=max(B);
        path_mean(i_s,:) = mean([gLON(path{lt}(i_s,:)') gLAT(path{lt}(i_s,:)')]);

        p=plot(path_min(1:i_s,1),path_min(1:i_s,2),'w','linewidth',2);
        
        for i_ss = 1:i_s
            sz = 20+(hours(sta_sm{lt}.end(i_ss)-sta_sm{lt}.start(i_ss)))/100;
            p2(i_ss,1) = plot(path_min(i_ss,1),path_min(i_ss,2),'.w', 'MarkerSize',sz+15);
            p2(i_ss,2) = plot(path_min(i_ss,1),path_min(i_ss,2),'.', 'MarkerSize',sz,'color',colordergrad(i_ss,:));
        end
        % keyboard
        
        %gif
        % tiadd = datenum([sta_sm{lt}.start(i_s) sta_sm{lt}.end(i_s)]-t(1))/datenum(t(end)-t(1));
        %pause(1)
        if i_s<height(sta_sm{lt})
            delete(p);   delete(p2)
        end
    end
    
 
    exportgraphics(gca,['gibbs_5000_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '.png'])
    % keyboard
    close all
end
