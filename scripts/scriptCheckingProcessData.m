%%
addpath(genpath('../functions'))
load("../data/processedDataTraquet.mat")
scriptAltPres

%%
lt=10;

col = [42,71,94;126 71 149;38 38 38;5 102 47;108 49 14]/255;

figure; hold on;
plot(raw{lt}.pressure.date,raw{lt}.pressure.obsWithOutliars,'color',[.4 .4 .4]);  hold on
    
for i_s = 1:height(sta{lt})
    
    id_tgr = find(sta{lt}.start(i_s)<spttime & spttime < sta{lt}.end(i_s));
    id_tge = sta{lt}.start(i_s)<raw{lt}.pressure.date & raw{lt}.pressure.date < sta{lt}.end(i_s);
    % pres_ge = movmean(movmedian(raw{lt}.pressure.obs(id_tge),3),dt);
    pres_ge = movmean(raw{lt}.pressure.obs(id_tge),3);
    id_t_gre = ismember(raw{lt}.pressure.date(id_tge),spttime(id_tgr));
    pres_gr = pres_ge(id_t_gre);
    plot(spttime(id_tgr), pres_gr,'LineWidth',2,'color',col(mod(i_s+1,2)+1,:))

    if sta{lt}.status(i_s)=="equipment" || sta{lt}.status(i_s)=="retrieval"
        plot(sp{lt}{i_s}.time,sp{lt}{i_s}.presCalib-nanmean(sp{lt}{i_s}.presCalib)+nanmean(pres_gr),'--r')
    end

    plot(sp{lt}{i_s}.time,sp{lt}{i_s}.pres-nanmean(sp{lt}{i_s}.pres)+nanmean(pres_gr),'color','r')
end
ylabel({raw{lt}.GDL_ID ,'Pressure(hPa)'})
grid on; box on; axis tight;

%%
sta_sm=cell(1,height(tblLog));
for lt=1:height(tblLog)
    grp_id = hours(sta{lt}.end-sta{lt}.start)>0;%sta{lt}.twlNb>=4;
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


%%
% subp_row=2*ones(1,height(tblLog));
% subp_row([11 12])=3;
% subp_row([15 16])=1;
for lt=10:height(tblLog)
 
    %figure('position',[0 0 1200 750], 'Name', [raw{lt}.GDL_ID ' | ' tblLog.CommonName{lt}] );
    %tiledlayout('flow','TileSpacing','tight','Padding','tight')
    
    mvt_pdf = movementModel('energy',tblLog.mass(lt),tblLog.wingSpan(lt));

    
    [gLON,gLAT] = meshgrid(lon{lt},lat{lt});
    for i_s = 1:height(sta_sm{lt})
        %nexttile; hold on;
        
        figure('position',[0 0 1600 900]); hold on ; xticks([]); yticks([])
        
        %set(gca,'Color','k')
        tmp_p = pres_prob{lt}(:,:,sta_sm{lt}.staID(i_s));
        tmp_pt = pres_thr{lt}(:,:,sta_sm{lt}.staID(i_s));
        tmp_l = light_prob{lt}(:,:,sta_sm{lt}.staID(i_s));
        
        
        % Option 1
        mask = 0.3+0.7*double(tmp_pt);
        mask(mask_water{lt}) = 0;
        img_tmp = real2rgb(tmp_p,colormap);
        img_tmp = img_tmp.*mask;
        
        imagesc(lon{lt},lat{lt},img_tmp)
        borders('countries','w')
        % imagesc(lon{lt},lat{lt},tmp_p,'AlphaData',mask); 
        c_axis=caxis();
        contour(lon{lt},lat{lt},tmp_l,3,'color',[255, 252, 49]/255,'linewidth',2);
        

        
        tt = datestr(sta_sm{lt}.start(i_s),'dd-mmm HH:MM') + " to " + datestr(sta_sm{lt}.end(i_s),'dd-mmm HH:MM');
        tt = tt+ " (pres.: "+ num2str(pres_n{lt}(i_s)) + " hrs, light: "+ num2str(sta_sm{lt}.twlNb(i_s) + " twls)");
        
        if i_s>1
            tt = tt+ " | flight: "+num2str(round(hours(sta_sm{lt}.actEffort(max(1,i_s-1))))) + " hrs"  ;
            % tt = {tt{:} num2str(sta_sm{lt}.actNb(max(1,i_s-1))) + "act";
            if exist('gr','var')
                path_short = [lon{lt}(gr{lt}.sp.lon) lat{lt}(gr{lt}.sp.lat)];
                tmpd = reshape(lldistkm([path_short(i_s-1,2) path_short(i_s-1,1)],[gLAT(:) gLON(:)]),size(gLAT))./hours(sta_sm{lt}.actEffort(i_s-1));
                plot(path_short(i_s-1,1), path_short(i_s-1,2),'.w','linewidth',2,'MarkerSize',40)
                plot(path_short(i_s-1,1), path_short(i_s-1,2),'.g','linewidth',2,'MarkerSize',30)
                plot(path_short(i_s,1), path_short(i_s,2),'or','linewidth',2,'MarkerSize',12)
            else
                tmpd = reshape(lldistkm([gLAT(id) gLON(id)],[gLAT(:) gLON(:)]),size(gLAT))./hours(sta_sm{lt}.actEffort(i_s-1));
                plot(gLON(id), gLAT(id),'.w','linewidth',2,'MarkerSize',40)
                plot(gLON(id), gLAT(id),'.g','linewidth',2,'MarkerSize',30)
                [~,id]=max(tmp_p(:).*tmp_pt(:).*tmp_l(:).*mvt_pdf(tmpd(:)));
                plot(gLON(id), gLAT(id),'or','linewidth',2,'MarkerSize',12)
            end
            contour(lon{lt},lat{lt},tmpd,[7 17]/1000*60*60,'-g','linewidth',2); % ,'ShowText','on'
            
        else
            [~,id]=max(tmp_p(:).*tmp_pt(:).*tmp_l(:));
            plot(gLON(id), gLAT(id),'.w','MarkerSize',30)
        end

        % t=text(.01,.8,tt,'Units','normalized','FontSize',40,'FontWeight','bold','Color','w','VerticalAlignment','top');

        title(tt)

        caxis(c_axis); xticklabels('');yticklabels('');
        
        plot(raw{lt}.calib.lon,raw{lt}.calib.lat,'.w','MarkerSize',40)
        plot(raw{lt}.calib.lon,raw{lt}.calib.lat,'.r','MarkerSize',30)
        
        axis equal; axis([min(lon{lt}) max(lon{lt}) min(lat{lt}) max(lat{lt}) ]);
        
%         x0=[raw{lt}.calib.lon-4 raw{lt}.calib.lat-4];
%         dx = 1./[lldistkm(x0,x0+[1 0]) lldistkm(x0,x0+[0 1])]*100;
%         h=plot([x0(1) x0(1)], [x0(2) x0(2)+dx(2)], '-r', 'LineWidth', 2);%label(h,'100km','location','middle','slope')
%         h=plot([x0(1) x0(1)+dx(1)], [x0(2) x0(2)],'-r', 'LineWidth', 2); label(h,'100km')
%         
        exportgraphics(gcf,['combined_map_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '_' num2str(i_s) '.png'],'Resolution',300)
        %keyboard
        % close all
        
    end
     % keyboard
    % exportgraphics(gcf,['combined_map_48h_' tblLog.CommonName{lt} '_' raw{lt}.GDL_ID '.png'],'Resolution',300)
end
