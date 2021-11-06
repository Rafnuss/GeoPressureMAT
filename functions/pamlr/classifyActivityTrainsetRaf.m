function [activity,km2] = classifyActivityTrainsetRaf(raw,km)
%CLASSIFYACTIVITY uses activity data to classify migratory flapping flight.
%   Adapted from classifyFLAP

filename = ['../data/labels/activity_label/' raw.GDL_ID '_act_pres'];

if ~exist([filename '-labeled.csv'])>0
    % Create table of activity labeled (1: small activity, 2: high activity, 3: migration)
    label = strings(numel(raw.acceleration.act ),1);
    label(km~=0) = num2str(km(km~=0));
%     label(km==1) = "low activity";
%     label(km==2) = "high activity";
%     label(km==3) = "migration";
%     label(km==4) = "outliar";

    % Align end date for pressure, acceleration and light
    max_date = max([ raw.pressure.date(end) raw.acceleration.date(end)]);

    raw.pressure.date = raw.pressure.date(1):diff(raw.pressure.date(1:2)):max_date;
    raw.pressure.obs = [raw.pressure.obs ; -1*ones(numel(raw.pressure.date)-numel(raw.pressure.obs),1)];
    if ~isnat(raw.acceleration.date)
        raw.acceleration.date = raw.acceleration.date(1):diff(raw.acceleration.date(1:2)):max_date;
        label = [label;strings(numel(raw.acceleration.date)-numel(raw.acceleration.act),1)];
    else
        raw.acceleration.date=raw.pressure.date';
        label_acc=ones(size(raw.acceleration.date));
        label_acc(1:4)=[1 2 3 4 ];
        label = num2str(label_acc);
        raw.acceleration.act = zeros(size(raw.acceleration.date));
    end
    
    raw.acceleration.pit = [raw.acceleration.pit ; -1*ones(numel(raw.acceleration.date)-numel(raw.acceleration.pit),1)];
    raw.acceleration.act = [raw.acceleration.act ; -1*ones(numel(raw.acceleration.date)-numel(raw.acceleration.act),1)];

    
    T =  table(...
        [repmat("act",numel(raw.acceleration.date),1) ; repmat("pres",numel(raw.pressure.date),1)],...
        [datestr(raw.acceleration.date,'yyyy-mm-ddTHH:MM:SS.000Z') ; datestr(raw.pressure.date,'yyyy-mm-ddTHH:MM:SS.000Z')],...
        [raw.acceleration.act ; raw.pressure.obs],...
        [label; strings(numel(raw.pressure.date),1)],...
        'VariableNames',{'series','timestamp','value','label'});

    writetable(T,[filename '.csv'])

    %web(url)
    keyboard
end

% Read labled file
T2 = readtable([filename '-labeled.csv']);
km2 = T2.label(strcmp(T2.series,'act'));
km2(isnan(km2))=0;
tmp = T2.timestamp(strcmp(T2.series,'act'));
datekm2 = datetime(cellfun(@(x) [x(1:10) ' ' x(12:16)],tmp,'UniformOutput',false));
tmp_act = T2.value(strcmp(T2.series,'act'));
tmp_act(tmp_act==-1)=nan;

%% Similar as Part 2 of classifyActivityRaf

% Creat an activity ID of the classification (conitnous acitivity)
km_activity = [1; cumsum(diff(km2)~=0)+1];

% Compute the duration per activity and activity unique ID (to find the class
% associated
activity=table();
[~,activity.id]=unique(km_activity);
activity.km = km2(activity.id);
activity.duration = (groupcounts(km_activity)+1)*diff(datekm2(1:2));
% activity.act_max  = splitapply(@max,raw.acceleration.act,km_activity);
% activity.act_min  = splitapply(@min,raw.acceleration.act,km_activity);
activity.act_sum = splitapply(@sum,tmp_act,km_activity);
activity.date_max = splitapply(@max,datekm2,km_activity);
activity.date_min = splitapply(@min,datekm2,km_activity);
% if (isfield(raw.acceleration, 'isDay'))
%     % isDay -> 0: only night -> 1: only day
%     activity.isDay = splitapply(@(x) mean(x) ,raw.acceleration.isDay,km_activity);
% end


end