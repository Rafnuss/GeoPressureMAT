function [activity,km] = classifyActivityRaf(dta,varargin)
%CLASSIFYACTIVITY uses activity data to classify migratory flapping flight.
%   Adapted from classifyFLAP

% Define threshold, otherwise automatic
if (nargin>1)
    threshold = varargin{1};
else
    threshold = 'auto';
end
% Define duration separing high activity and migration
if (nargin>2)
    d_max = varargin{2};
else
    d_max = duration(0,10,0);
end
% Define if detailed plot and calibration should be shown
if (nargin>3)
    plotit = varargin{3};
else
    plotit = 0;
end

% Get acceleration data
act=double(dta.act);
date=dta.date;

%% Part 1: Classification of Activity
if strcmp(threshold,'auto')
    km=zeros(size(act));
    % Classify as 1 or 2 according to low or high activity found by
    % clustering
    [km(act>0),C] = kmeans(act(act>0), 2,'Start',[0;100]);
    % Threadhols is between the two cluster
    threshold=mean(C);
else
    % If threshold provided, then compute directly the classes 1,2 and 3
    km = double(act>threshold);
    km(act>0) = km(act>0)+1;
end

% Show the classification
if plotit <= 1
    figure('position',[0 0 900 600]); 
    subplot(2,1,1);hold on; 
    histogram(act(km==0)); 
    histogram(act(km==1)); 
    histogram(act(km==2)); 
    set(gca,'YScale','log')
    xline(threshold)
    legend('No Activity','Low Activity','High Activity');
    xlabel('Activity'); ylabel('Histogram')
    
     subplot(2,1,2); hold on; 
    plot(dta.date, act,'-k')
    plot(dta.date(km==1), act(km==1),'.')
    plot(dta.date(km==2), act(km==2),'.')
end


%% Part 2 Computation of activity duration

% time resolution
res = diff(dta.date(1:2));

% Creat an activity ID of the classification (conitnous acitivity)
km_activity = [1; cumsum(diff(km)~=0)+1];

% Compute the duration per activity and activity unique ID (to find the class
% associated
activity=table();
[~,activity.id]=unique(km_activity);
activity.km = km(activity.id);
activity.duration = groupcounts(km_activity)*res;
%activity.act_max  = splitapply(@max,act,km_activity);
%activity.act_min  = splitapply(@min,act,km_activity);
%activity.act_mean = splitapply(@mean,act,km_activity);
%activity.date_max = splitapply(@max,date,km_activity);
%activity.date_min = splitapply(@min,date,km_activity);
% activity.date_mean= splitapply(@mean,date,km_activity);
if (isfield(dta, 'isDay'))
    % isDay -> 0: only night -> 1: only day
    activity.isDay = splitapply(@(x) mean(x) ,dta.isDay,km_activity);
end


% Show the high activity and their (cumulative) duration
if plotit <= 2 
    figure('position',[0 0 900 600]); hold on;
    histogram(activity.duration(activity.km==0),0:res:duration(12,0,0))
    histogram(activity.duration(activity.km==1),0:res:duration(12,0,0))
    histogram(activity.duration(activity.km==2),0:res:duration(12,0,0))
    legend('No Activity','Low Activity','High Activity');
    set(gca,'YScale','log');
    xlabel('Duration of Activity'); ylabel('Histogram')
    
%     figure('position',[0 0 900 300]); 
%     subplot(1,3,1);hold on
%     plot(activity.act_mean(activity.km==1),activity.duration(activity.km==1),'.')
%     plot(activity.act_mean(activity.km==2),activity.duration(activity.km==2),'.')
%     xlabel('Mean Activity Level');
%     ylabel('Activity duration')
%     xline(threshold); yline(d_max); box on;
%     subplot(1,3,2);hold on
%     plot(activity.act_min(activity.km==1),activity.duration(activity.km==1),'.')
%     plot(activity.act_min(activity.km==2),activity.duration(activity.km==2),'.')
%     xlabel('Min Activity Level');
%     ylabel('Activity duration')
%     xline(threshold); yline(d_max); box on;
%     subplot(1,3,3);hold on
%     plot(activity.act_max(activity.km==1),activity.duration(activity.km==1),'.')
%     plot(activity.act_max(activity.km==2),activity.duration(activity.km==2),'.')
%     xlabel('Max Activity Level');
%     ylabel('Activity duration')
%     xline(threshold); yline(d_max); box on;
    
    figure('position',[0 0 900 300]); 
    dd = duration(0,5,0):res:duration(0,15,0);
    b = nan(numel(dd)+1,3);
    for i=1:numel(dd)
        b(i,:) = histcounts(activity.isDay(activity.km==2 &activity.duration == dd(i)),[-eps eps 1-eps 1+eps]);
    end
    b(i+1,:) = histcounts(activity.isDay(activity.km==2 &activity.duration > dd(i) ),[-eps eps 1-eps 1+eps]);
    bar([minutes(dd) minutes(dd(end)+res)], b./sum(b,2))
    tmp = get(gca,'XTickLabel')+"min"; tmp(end)=">"+tmp{end};
    set(gca,'XTickLabel',tmp)
    legend('Only Night','Night and day','only day')
end

%% Part 3: Classify Migratory activity

% Migration is difined by:
% - Hight activity (km=2)
% - duration above d_max
% - happen solely during the night;
% - (we could add during migration period maybe?)
activity.km(activity.km==2&activity.duration>d_max)=3;

km(activity.km(km_activity)==3)=3;

if plotit <= 3
    figure('position',[0 0 900 300]);  hold on; box on;
    tmp = act; tmp(km<2)=nan;  
    tmp2 = date; tmp2(km<2)=NaT;
    tmp3 = minutes(activity.duration(km_activity));
    tmp3(tmp3>60)=60;
    plot(tmp2, tmp,'-k')
    plot(date, act,'.','color',[.95 .95 .95], 'MarkerSize',20)
    plot(date(km==2), act(km==2),'.','color',[.8 .8 .8], 'MarkerSize',20)
    plot(date(km==2 & tmp3>5), act(km==2 & tmp3>5),'.','color',[.6 .6 .6], 'MarkerSize',20)
    scatter(date(km==3), act(km==3),[], tmp3(km==3),'o','filled')
    yline(threshold,'--k')
    axis tight; ylim([0 100]);
    legend('Continous high activity','Low actvitiy (<19)','High activity of 5min','High activity 10 min',['Migration >' num2str(minutes(d_max)) ' min+ during the night'])
end

end

