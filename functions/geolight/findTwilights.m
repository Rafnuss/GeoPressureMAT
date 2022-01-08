function out = findTwilights(tagdata, threshold)
%FINDTWILIGHTS search for pairs of twilights spanning night.
%   Search for sunset, sunrise pairs that correspond to a given light
%   threshold. Given a set of times (include) known to fall in the night,
%   findTwilights determines the twilights that span these times, and
%   computes the corresponding midnights. It then searches for periods of
%   darkness that lie approximately 24 hours from these midnights,
%   repeating the process until no new twilight pairs are found. If
%   interleave=TRUE, the sunrise and sunset times are interleaved
%   andreturned as a single sequence of twilights, otherwise sunset and
%   sunrise times are returned separately. The function
%   interleave.twilights takes a dataframe of separate sunset and sunrise
%   times and interleaves them to form a sequence of twilight times.

date = tagdata.date;
light = tagdata.obstrans;

% Compute exceedance of light
l = light >= threshold;

% Compute when exceeding the threshold: -1 is pseudo-sunset (i.e., getting from 
% below threashold to above (inversely for +1).
f = diff(l);
a = find(f == -1); 
b = find(f == 1);

% We force the start with a sunrise and delete the first sunset if before.
if b(1) < a(1)
    b = b(2:end);
end
a = a(1:numel(b));


% Keep only if night duration is greater than and less than 
keepMax = 1; % hours in day
keepMin = 5/24; % hours in days
% WHY USING +1 ??? keep = (date(b + 1) - date(a)) < keepMax & (date(b + 1) - date(a)) > keepMin; 
keep = (date(b) - date(a)) < keepMax & (date(b) - date(a)) > keepMin; % WHY USING +1 ??? This 
% Note that we are deleting both a and b when the duration threashold are
% not met. We could also delete a(1:end-1) and b(2:end) or a(1:end-1) and
% b(2:end). Because we force the start of a and b to be a sunset (l. 17),
% this order will force to keep the night as the SHORTEST period without
% any light. Therefore assuming that going below threshold is possibly by
% mistake but not going above. 
% light 
a = a(keep);
b = b(keep);


% Second filtering which keep only the nights which cover a certain timing
% (e.g., 00:00).
% mDate = date(a) + (date(b) - date(a))/2;
% keep = false(size(a));
% include = '';
% extend = 0; extend = 60 * extend;
% add = containsAny(date(a), date(b), include) & ~keep;
% while (any(add))
%     keep = keep | add;
%     mid = c(mDate(add) - 86400, mDate(add) + 86400);
%     add = containsAny(date(a) - extend, date(b) + extend, mid) & ~keep;
% end
% a = a(keep);
% b = b(keep);

% Figure
% figure; hold on
% plot(dateshift(date(a),'start','day'),date(a)-dateshift(date(a),'start','day'),'.')
% plot(dateshift(date(b),'start','day'),date(b)-dateshift(date(b),'start','day'),'.')

% Adjust the sunrise and sunset time with a linear interpolation of the
% light and the threashold. 
ss = date(a) + (threshold - light(a))/(light(a + 1) - light(a)) * (date(a + 1) - date(a));
sr = date(b) + (threshold - light(b))/(light(b + 1) - light(b)) * (date(b + 1) - date(b));

% Export
out = table(reshape([ss sr]',[],1), repmat([false; true], numel(ss),1), 'VariableNames',{'Twilight','Rise'});
end

% function f = containsAny(a, b, x)
%     f = logical(length(a));
%     for k = x
%         f = f | (a <= x(k) & b >= x(k));
%     end
% end
