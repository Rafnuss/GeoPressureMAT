function [d_grid,t_grid,OBS]=sensorIMG(ax,date,obs,varargin)
%SENSORIMG Plot light or sensor image
%   Detailed explanation goes here

if numel(unique(diff(date)))~=1
    warning('Not regular grid!')
end

if nargin>3
    shiftK = varargin{1};
else
    shiftK = 0;
end

if isnat(date)
    return
end

d = dateshift(date,'start','day');
t = date - d;

t_grid = sort(unique(t));
d_grid = sort(unique(d));

[D,T] = meshgrid(d_grid,t_grid);
[~,Locb] = ismember([datenum(t) datenum(d)], [datenum(T(:)) datenum(D(:))],'rows');
OBS = nan(size(T));
OBS(Locb) = obs;
%OBS = OBS';

% OBS = reshape(obs,numel(t_grid),numel(d_grid));
K = round(size(OBS,2)*shiftK);
imagesc(ax,datenum(d_grid),datenum(t_grid),circshift(OBS,K));

if nargin>4
    twl = varargin{2};
    if ~any(strcmp('isOutliar',twl.Properties.VariableNames)), twl.isOutliar=true(size(twl.Twilight)); end
    x1=twl.Twilight(twl.Rise&~twl.isOutliar); 
    y1=dateshift(x1,'start','day');
    plot(ax,datenum(y1),mod(datenum(x1-y1)+K*days(diff(t_grid(1:2))),1),'.','MarkerSize',20)
    x2=twl.Twilight(~twl.Rise&~twl.isOutliar); 
    y2=dateshift(x2,'start','day');
    plot(ax,datenum(y2),mod(datenum(x2-y2)+K*days(diff(t_grid(1:2))),1),'.','MarkerSize',20)
    x3=twl.Twilight(twl.isOutliar); 
    y3=dateshift(x3,'start','day');
    plot(ax,datenum(y3),datenum(x3-y3),'.','MarkerSize',20,'Color',[.2 .2 .2])
end

if nargin>5
    calib = varargin{3};
    xline(ax,datenum(calib.first_period(end)),'--w','LineWidth',2)
    if ~any(isnat(calib.second_period))
        xline(ax,datenum(calib.second_period(1)),'--w','LineWidth',2)
    end
    plot(ax,datenum(y1),mod(datenum(twilight(x1, calib.lon, calib.lat, ones(size(x1)))-y1)+K*days(diff(t_grid(1:2))),1),'-','LineWidth',2);
    plot(ax,datenum(y2),mod(datenum(twilight(x2, calib.lon, calib.lat, zeros(size(x2)))-y2)+K*days(diff(t_grid(1:2))),1),'-','LineWidth',2);
end

datetick('x'); datetick('y')
view(2); axis tight;
colormap(bone)

end

