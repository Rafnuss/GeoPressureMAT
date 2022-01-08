function out = findHEZenith(twl, varargin)
%FINDHEZENITH estimate location from consecutive twilights
%   This function estimates the location given the times at which the
%   

if nargin>1
    plotit = varargin{1};
else
    plotit = true;
end

z = 89:0.25:99;
crd = coord(twl, z, 0.08);

[~,id] = min(nanstd(crd.lat));
out = z(id);

if plotit
    figure;
    subplot(2,1,1);hold on;
    plot(twl.Twilight(1:2:end), crd.lat','color',[.8 .8 .8])
    plot(twl.Twilight(1:2:end), nanmean(crd.lat,2),'k');
    plot(twl.Twilight(1:2:end), nanmean(crd.lat,2)+nanstd(crd.lat,[],2),'--k')
    plot(twl.Twilight(1:2:end), nanmean(crd.lat,2)-nanstd(crd.lat,[],2),'--k')
    subplot(2,1,2);hold on;
    plot(z,nanstd(crd.lat))
    plot(out,min(nanstd(crd.lat)),'or')
end
end
