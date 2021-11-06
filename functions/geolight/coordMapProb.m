function prob_r = coordMapProb(twlt, gEf, lon, lat)
%COORDMAP estimate probability of each location from twilights

% twlt = twlt(~twlt.isOutliar,:);

% Compute full grid
[LON, LAT] = meshgrid(lon,lat);

% Compute the zenith angle at each point on the grid for all twilight
z = refracted(zenith(twlt.Twilight, [LON(:) LAT(:)]));

% Convert the difference of time into probability using the Gamm fitted function
% prob = gEf.pdfR(z);
% prob(~twlt.Rise,:) = gEf.pdfS(z(~twlt.Rise,:));
prob = gEf.pdf(z);

% Reshapre the map of probability and normalize it
prob_r = reshape(prob',numel(lat),numel(lon),[]);

% This is very wront to do! I keep it to remind me not to use it!
% prob_r = prob_r ./sum(sum(prob_r,1),2);

end
