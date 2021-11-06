function out = getElevation(twlt, known_coord, varargin)
%GETELEVATION calculate the sun elevation angle for light measurements at a known location
%   getElevation 

if nargin>2
    plotit = varargin{1};
else
    plotit = true;
end

% This is the original method. I don't fully understand why going through the trouble of estimating both a zenith angle and a distribution of time
% error... 
% % Compute the zenith angle corresponding to each twlilight
z = refracted(zenith(twlt.Twilight, known_coord));
% 
% % Compute the time difference between measured and computed twilight at calibration site.
% % The reference zenith angle z0 might need to be adjusted in case a mesured twilight occurs before the computed twilight
% z0 = max(z)-0.01;
% cond = true;
% while (cond)
%     z0 = z0 + 0.01;
%     twl_t = twilight(twl.Twilight, known_coord(1), known_coord(2), twl.Rise, z0);
%     twl_dev = (twl.Rise*2-1) .* minutes(twl.Twilight-twl_t);
%     cond = ~all(twl_dev >= 0);
% end
% 
% % Fit a gamma distribution
% fit = fitdist(twl_dev,'gamma');

% Here is a more simple alternative assuming a zenith of 96 and allowing below twl_dev below zero
% twl_t = twilight(twlt.Twilight, known_coord(1), known_coord(2), twlt.Rise);
% twl_dev = (twlt.Rise*2-1) .* minutes(twlt.Twilight-twl_t);

if any(strcmp('isOutliar',twlt.Properties.VariableNames))
    twlt_isOutliar = twlt.isOutliar;
else 
    twlt_isOutliar=false(size(twl_dev));
end

if nargin>3
    bins = varargin{2};
    edges = [bins-diff(bins(1:2))/2 bins(end)+diff(bins(1:2))/2];
else
    [~,edges] = histcounts(z(~twlt_isOutliar));
    bins = edges(1:end-1)+diff(edges)/2;
end
NRise = histcounts(z(twlt.Rise&~twlt_isOutliar),edges);
NSet = histcounts(z(~twlt.Rise&~twlt_isOutliar),edges);


if plotit
    figure('position',[0 0 800 800]); 
    subplot(2,1,1); hold on;
    plot(twlt.Twilight(twlt.Rise&~twlt_isOutliar), z(twlt.Rise&~twlt_isOutliar),'.','markerSize',10)
    plot(twlt.Twilight(~twlt.Rise&~twlt_isOutliar), z(~twlt.Rise&~twlt_isOutliar),'.','markerSize',10)
    plot(twlt.Twilight(twlt_isOutliar), z(twlt_isOutliar),'.','color',[.2 .2 .2],'markerSize',10)
    if (sum(~twlt_isOutliar)>0) yline(median(z(~twlt_isOutliar)),'--k'); end
    ylabel('Zenith');legend('Sunrise','Sunset','Outliar')
    subplot(2,1,2); hold on;
    bar(bins, [NRise' NSet'],1,'stacked','FaceAlpha',0.6);
    ylabel('Number of Twilight')
    % yyaxis right; plot(0:bins(end),fit.pdf(0:bins(end)),'r','linewidth',2); ylabel('PDF')
    xlabel('Zenith'); 
    legend('Sunrise','Sunset')
end



out = table(median(z(~twlt_isOutliar)),bins,NRise,NSet,z','VariableNames',{'medZ','bins','NRise','NSet','z'});
end