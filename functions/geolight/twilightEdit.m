function twl = twilightEdit(twl,varargin)

% p = inputParser;
% addOptional(p,'offset',17);
% addOptional(p,'window',4);
% addOptional(p,'outlierMins',45);
% addOptional(p,'stationaryMins',15);
% addOptional(p,'plot',true);
% parse(p,varargin{:});
% 
% names = fieldnames(p.Results);
% for i=1:length(names)
%     eval([names{i} '=' num2str(p.Results.(names{i})) ';']);
% end
% 
% 
% day = twilights.Twilight;
% hour = hourOffset(as.hour(twilights.Twilight), offset);
% sunr = twilights.Twilight(twilights.Rise);
% suns = twilights.Twilight(~twilights.Rise);
% 
% tmp = datenum(sunr-dateshift(sunr,'start','day'));
% v = [ones(round(window/2),1);0;ones(round(window/2),1)]/window;
% M = conv(tmp,v/4,'same');
% abs(tmp-M)*24*60>outlierMins;

if nargin>1
    plotit = varargin{1};
else
    plotit = true;
end


twl_day = dateshift(twl.Twilight,'start','day');
twl_hour = twl.Twilight-twl_day;

twl_isOutliar = false(size(twl.Twilight));

twl_isOutliar(twl.Rise) = isoutlierRaf(twl_hour(twl.Rise),twl_day(twl.Rise));
twl_isOutliar(~twl.Rise) = isoutlierRaf(twl_hour(~twl.Rise),twl_day(~twl.Rise));

if plotit
    figure('position',[0 0 1000 800]);
    subplot(2,1,1); hold on; ylabel('Sunrise'); box on;
    plot(twl_day(twl.Rise), twl_hour(twl.Rise),'-o')
    plot(twl_day(twl.Rise & ~twl_isOutliar), twl_hour(twl.Rise & ~twl_isOutliar),'-o')
    legend('Is Outliar', 'cleaned');
    subplot(2,1,2); hold on; ylabel('Sunset'); box on;
    plot(twl_day(~twl.Rise), twl_hour(~twl.Rise),'-o')
    plot(twl_day(~twl.Rise & ~twl_isOutliar), twl_hour(~twl.Rise & ~twl_isOutliar),'-o')
    legend('Is Outliar', 'cleaned');
end
% twl = twl(~twl_isOuliar,:);
twl.isOutliar = twl_isOutliar;
end


function isO = isoutlierRaf(y,x)

x=datenum(x)-datenum(x(1)); y=datenum(y);

% First removal of the obvious errors
isO = isoutlier(y,'movmedian',13,'ThresholdFactor',5);
%isO=false(size(y));

% figure('position',[0 0 1000 300]); hold on;
% plot(x,y,'-o'); plot(x(~isO),y(~isO),'-o');

% Remove trend
ydt=nan(size(y));
ydt(~isO) = detrend(y(~isO),15,'SamplePoints',x(~isO));

% Test: accounting for the displacement between staionary period
% coeff= 1/24*2;% minutes of displacement to change
% dx = minutes(sta.actDuration(twl.staID(twl.Rise)))*coeff;
% dx(~diff(twl.staID(twl.Rise)))=0;
% x2=cumsum([0; diff(x)]+dx);
% ydt = detrend(y,20,'Continuous',true,'SamplePoints',x2);



% % Method Grubbs
% isO = isoutlier(ydt);
% isO = isoutlier(ydt,'movmedian',51,'ThresholdFactor',3,'SamplePoints',x);
% isO(~isO) = isoutlier(ydt(~isO),'movmedian',7,'ThresholdFactor',3,'SamplePoints',x(~isO));

isO = isO|isoutlier(ydt,'quartiles');

% 
% isO = isoutlier(ydt,'gesd');

% figure('position',[0 0 1000 300]); hold on;
% % plot(x,y,'-o')
% plot(x,ydt,'-o')
% plot(x(~isO),ydt(~isO),'-o')



end
