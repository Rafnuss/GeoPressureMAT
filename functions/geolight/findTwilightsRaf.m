function out = findTwilightsRaf(tagdata, varargin)
%FINDTWILIGHTSRAF search for pairs of twilights spanning night.
%   Rewrite of FINDTWILIGHTS. Assume certain structure. see below
if nargin>1 && ~isempty(varargin{1})
    threshold = varargin{1};
else
    % by default the threshold is the minimum value recorded by the
    % geolocator
    threshold = max(mink(unique(tagdata.obs),2));
end



% get temporal resolution
res = diff(tagdata.date(1:2));
assert(all(diff(tagdata.date)==res),'Date need to be equaly spaced');

% pad with time to start at 00:00
pad_s = (dateshift(tagdata.date(1),'start','day'):res:(tagdata.date(1)-res))'; 
pad_e = (tagdata.date(end):res:(dateshift(tagdata.date(end),'end','day')-2*res))';
date = [pad_s ; tagdata.date ; pad_e];
light = [nan(size(pad_s)) ; tagdata.obs ; nan(size(pad_e))];

% Reshape light to be in day x hour format
light_r = reshape(light,1/days(res),[]);
date_r = reshape(date,1/days(res),[]);

if nargin>2
    shiftK = varargin{2};
    K = round(size(light_r,2)*shiftK);
    light_r = circshift(light_r, K);
    date_r = circshift(date_r, K);
end

% Compute exceedance of light
l = light_r >= threshold;

% Find the first light
[~,id_sr]=max(l);
sr = date_r(sub2ind(size(l),id_sr,1:size(l,2)));

% And the last light
[~,id_ss]=max(flipud(l));
ss = date_r(sub2ind(size(l),size(l,1)-id_ss+1,1:size(l,2)));

tmp = [sr; ss];
tmp = tmp(:);

% Export
% out = table(tmp(1:end-1), tmp(2:end), [repmat([1;2],numel(sr)-1,1);1], 'VariableNames',{'tFirst','tSecond','type'});
out = table(tmp, repmat([true;false],numel(sr),1), 'VariableNames',{'Twilight','Rise'});
end
