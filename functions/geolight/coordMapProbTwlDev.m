function prob_r = coordMapProbTwlDev(twlt, gEf, lon, lat)
%COORDMAP estimate probability of each location from twilights

% twlt = twlt(~twlt.isOutliar,:);

% Compute full grid
[LON, LAT] = meshgrid(lon,lat);
    
% Compute the twilight at all coordinate of the grid for the date of each twl
twl_t = twilight(twlt.Twilight, LON(:), LAT(:), twlt.Rise);

% twlt.Twilight2 = twlt.Twilight + (~twlt.Rise*2-1) .* gEf.twl_dev/60/24;


% Compute the time difference between the measured and the computed on the gride
twl_dev = (~twlt.Rise*2-1) .* minutes(twl_t-twlt.Twilight);

% Convert the difference of time into probability using the Gamm fitted function
%gEf_pdf_tmp = @(x) double(abs(x)<1);
%prob = gEf_pdf_tmp(twl_dev);
% prob = gEf.pdf(twl_dev-1*(twlt.Rise*2-1)-2);
prob = gEf.pdfR(twl_dev);
prob(~twlt.Rise,:) = gEf.pdfS(twl_dev(~twlt.Rise,:));

% Reshapre the map of probability and normalize it
prob_r = reshape(prob',numel(lat),numel(lon),[]);
prob_r = prob_r ./sum(sum(prob_r,1),2);

end
