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
out.z = refracted(zenith(twlt.Twilight, known_coord));
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
    out.bins = varargin{2};
    edges = [out.bins-diff(out.bins(1:2))/2 out.bins(end)+diff(out.bins(1:2))/2];
else
    [~,edges] = histcounts(out.z(~twlt_isOutliar));
    out.bins = edges(1:end-1)+diff(edges)/2;
end
out.NRise = histcounts(out.z(twlt.Rise&~twlt_isOutliar),edges);
out.NSet = histcounts(out.z(~twlt.Rise&~twlt_isOutliar),edges);


if plotit
    figure('position',[0 0 800 800]); 
    subplot(2,1,1); hold on;
    plot(twlt.Twilight(twlt.Rise&~twlt_isOutliar), out.z(twlt.Rise&~twlt_isOutliar),'.','markerSize',10)
    plot(twlt.Twilight(~twlt.Rise&~twlt_isOutliar), out.z(~twlt.Rise&~twlt_isOutliar),'.','markerSize',10)
    plot(twlt.Twilight(twlt_isOutliar), out.z(twlt_isOutliar),'.','color',[.2 .2 .2],'markerSize',10)
    if (sum(~twlt_isOutliar)>0); yline(median(out.z(~twlt_isOutliar)),'--k'); end
    ylabel('Zenith');legend('Suout.NRise','Suout.NSet','Outliar')
    subplot(2,1,2); hold on;
    bar(out.bins, [out.NRise' out.NSet'],1,'stacked','FaceAlpha',0.6);
    ylabel('Number of Twilight')
    % yyaxis right; plot(0:out.bins(end),fit.pdf(0:out.bins(end)),'r','linewidth',2); ylabel('PDF')
    xlabel('Zenith'); 
    legend('Suout.NRise','Suout.NSet')
end

out.medZ = median(out.z(~twlt_isOutliar));
end