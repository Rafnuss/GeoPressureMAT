function twl = twilightEditTrainset(raw,twl,varargin)

if nargin>2
    plotit = varargin{1};
else
    plotit = true;
end

filename = ['Label/twilight_label/' raw.GDL_ID '_twl'];

twl_day = dateshift(twl.Twilight,'start','day');
twl_hour = minutes(twl.Twilight-twl_day);

% Write file
if ~exist([filename '-labeled.csv'])>0

    [G, act_day] = findgroups(dateshift(raw.acceleration.date,'start','day'));
    act_cum_acc = splitapply(@sum,raw.acceleration.act,G);

    T =  table(...
        ["Rise_" + twl.Rise ; repmat("act_cum",numel(act_day),1)],...
        [datestr(twl_day,'yyyy-mm-ddTHH:MM:SS.000Z') ; datestr(act_day,'yyyy-mm-ddTHH:MM:SS.000Z')],...
        [twl_hour ; act_cum_acc],...
        ["Outliar_" + twl.isOutliar; strings(numel(act_day),1)],...
        'VariableNames',{'series','timestamp','value','label'});

    writetable(T,[filename '.csv'])
    
    % https://github.com/Rafnuss/trainset
    % web('http://localhost:5000/')
    % https://trainset.herokuapp.com/
    keyboard
end

% Read labled file
T2 = readtable([filename '-labeled.csv']);
twl.isOutliar = contains(T2.label(contains(T2.series,'Rise_')),'Outliar_true');

if plotit
    figure('position',[0 0 1000 800]);
    subplot(2,1,1); hold on; ylabel('Sunrise'); box on;
    plot(twl_day(twl.Rise), twl_hour(twl.Rise),'-o')
    plot(twl_day(twl.Rise & ~twl.isOutliar), twl_hour(twl.Rise & ~twl.isOutliar),'-o')
    legend('Is Outliar', 'cleaned');
    subplot(2,1,2); hold on; ylabel('Sunset'); box on;
    plot(twl_day(~twl.Rise), twl_hour(~twl.Rise),'-o')
    plot(twl_day(~twl.Rise & ~twl.isOutliar), twl_hour(~twl.Rise & ~twl.isOutliar),'-o')
    legend('Is Outliar', 'cleaned');
end


end