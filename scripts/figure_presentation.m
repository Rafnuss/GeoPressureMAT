
% startup

load('../data/DEM/DEM_lt9.mat')

lt=10;

%% Get pressure

id_lon = find(glon>=tblLog.bndy_W(lt) & glon<=tblLog.bndy_E(lt));
id_lat = find(glat>=tblLog.bndy_S(lt) & glat<=tblLog.bndy_N(lt));
id_t = find(raw{lt}.pressure.date(1)<=spttime & spttime <= raw{lt}.pressure.date(end));
splt = permute(flip(ncread(file,'sp',[id_lon(1) numel(glat)-id_lat(end)+1 1 id_t(1)],[numel(id_lon) id_lat(end)-id_lat(1)+1 1 id_t(end)-id_t(1)+1]),2)/100,[2 1 4 3]);

%%
figure('position',[0 0 1600 900]); hold on; xticks(''); yticks('')
set(gca,'position',[0 0 1 1],'units','normalized')

m = imagesc(lon{lt},lat{lt},splt(:,:,1));
borders('countries','k')
axis equal; axis([min(lon{lt}) max(lon{lt}) min(lat{lt}) max(lat{lt}) ]); 

t=text(.01,1,datestr(spttime(id_t(1)),'dd-mmm-yyyy'),'Units','normalized','FontSize',40,'FontWeight','bold','Color','k','VerticalAlignment','top');

% colormap(cmocean('thermal')); caxis([750 1050])
colormap(cmocean('diff'));caxis([-20 20])

v = VideoWriter('pressure_surface_diff.mp4','MPEG-4');
v.Quality = 95;
v.FrameRate = 24*6;
open(v);
for i_t=1:min(v.FrameRate*20, numel(id_t))
    % m.CData = splt(:,:,i_t);
    m.CData = splt(:,:,i_t)-mean(splt,3);
    t.String = datestr(spttime(id_t(i_t)),'dd-mmm-yyyy');
    drawnow
    writeVideo(v,getframe(gcf));
end
close(v);


%% Pressure at the wrong place
i_s=1;
id_tgr = find(sta{lt}.start(i_s)<spttime(id_t) & spttime(id_t) < sta{lt}.end(i_s));
id_tge = sta{lt}.start(i_s)<raw{lt}.pressure.date & raw{lt}.pressure.date < sta{lt}.end(i_s);
pres_ge = movmean(movmedian(raw{lt}.pressure.obs(id_tge),3),dt*2+1);
id_t_gre = ismember(raw{lt}.pressure.date(id_tge),spttime(id_t(id_tgr)));
pres_gr = pres_ge(id_t_gre);


spltij = reshape(splt(60,100,id_tgr),1,[]);
spltij = reshape(splt(50,60,id_tgr),1,[]);
pres_gr_att = altitude2pressure(tblLog.DEMAttached(lt)-calibGeoDEM(lt,1), spltij);
pres_gr_min = altitude2pressure(calibGeoDEM(lt,2)-calibGeoDEM(lt,1)+dh_margin(1), spltij);
pres_gr_max = altitude2pressure(calibGeoDEM(lt,3)-calibGeoDEM(lt,1)+dh_margin(2), spltij);

figure; hold on;
plot(raw{lt}.pressure.date(id_tge),raw{lt}.pressure.obsWithOutliars(id_tge),'color',[.4 .4 .4]);  hold on
plot(spttime(id_t(id_tgr)),pres_gr,'color',col(1,:),'linewidth',2)
%plot(spttime(id_tgr),sp{lt}(id_tgr),'color',col(3,:))
plot(spttime(id_t(id_tgr)),pres_gr_att,'color',col(3,:),'linewidth',2)
plot(spttime(id_t(id_tgr)),pres_gr_min,'color',col(4,:),'linewidth',2)
plot(spttime(id_t(id_tgr)),pres_gr_max,'color',col(5,:),'linewidth',2)

figure; hold on;
borders('countries','k')
plot(lon{lt}(100),lat{lt}(60),'x')
plot(lon{lt}(60),lat{lt}(50),'x')
axis equal; axis([min(lon{lt}) max(lon{lt}) min(lat{lt}) max(lat{lt}) ]); 



%% 

[~,km] = classifyActivityTrainsetRaf(raw{lt},[]);
figure; hold on;
plot(raw{lt}.acceleration.date, raw{lt}.acceleration.act,'-k')
plot(raw{lt}.acceleration.date(km==2), raw{lt}.acceleration.act(km==2),'.')
plot(raw{lt}.acceleration.date(km==3), raw{lt}.acceleration.act(km==3),'.')