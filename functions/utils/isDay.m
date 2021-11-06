function isday = isDay(date,twl)
% ISDAY 
% 
% 

% force horitzontal vector 
date = date(:)';
tmp = datenum(date) - datenum(twl.Twilight);
tmp(tmp<0)=nan; % remove when date are before twilight.
[~,id] = min( tmp );
isday = twl.Rise(id);

% Add padding 
% isday = logical(movmax(isday,3));

end