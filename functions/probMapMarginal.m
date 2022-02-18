function [map, path] = probMapMarginal(prob_map, mvt_pdf, lon, lat, thr_speed, flight_duration, shortestPath)

%% NOT REALLY WORKING, concept yes, but no real check 

sz = size(prob_map);

% % Normalize the probability map for each night
% prob_map = prob_map ./ sum(prob_map,[1 2]);
% 
% % Set threashold corresponding to all nodes corresponding to
% % thr_prob_percentile of the probability for each stationary period
% tmp = sort(reshape(prob_map,[],grt.snds(3)));
% thr_prob = tmp(sub2ind(size(tmp),sum(cumsum(tmp)<=(1-thr_prob_percentile)),1:grt.snds(3)));


mask = all(prob_map<0.0001,3);

nlm = sum(~mask(:));

prob_map_r = reshape(prob_map(repmat(~mask,1,1,sz(3))),[],sz(3));

[LON,LAT] = meshgrid(lon, lat);
Ddist = squareform(pdist([LAT(~mask) LON(~mask)],@lldistkm));

% number of nodes in the 3d grid
szr = [nlm sz(3)];
n = prod(szr);

% Shortest path
if shortestPath
    mapWeightDJ = -log(prob_map_r(:,1));
    mapSourceDJ = nan(1,n);
    mapSourceDJ(1:nlm) = 1:nlm;
end

%%  Forward map
mapF = nan(1,n);
mapF(1:nlm) = prob_map_r(:,1);

for i_s=1:(sz(3)-1)
    % Compute the transition probability of all position to all based on
    % the movement model
    speed = Ddist/flight_duration(i_s);
    cond = speed<thr_speed & prob_map_r(:,i_s+1)>0 & mapF((1:nlm)+(i_s-1)*nlm)>0;
    prob_mvt = zeros(size(Ddist));
    prob_mvt(cond) = mvt_pdf(speed(cond));
    
    % Compute the transitional probability
    prob = prob_mvt .* prob_map_r(:,i_s+1);
    
    % Compute map forward at the next step
    mapF((1:nlm)+(i_s-1+1)*nlm) = mapF((1:nlm)+(i_s-1)*nlm) * prob';
    
    % compute DJ
    if shortestPath
        [mapWeightDJ, pos] = min(-log(prob') + mapWeightDJ',[],1);
        mapSourceDJ((1:nlm)+(i_s-1+1)*nlm) = pos;
    end
end



% backward map
mapB = nan(1,n);
mapB(end-nlm+1:end) = prob_map_r(:,end);

for i_s=(sz(3)-1):-1:1
    speed = Ddist/flight_duration(i_s);
    cond = speed<thr_speed & prob_map_r(:,i_s)>0 & mapB((1:nlm)+(i_s-1+1)*nlm)>0;
    prob_mvt = zeros(size(Ddist));
    prob_mvt(cond) = mvt_pdf(speed(cond));
    
    prob = prob_mvt .* prob_map_r(:,i_s);
    
    mapB((1:nlm)+(i_s-1)*nlm) = mapB((1:nlm)+(i_s-1+1)*nlm)*prob';
end


%%
% Combine forward and backward
map_r = mapF.*mapB;

% reshape map
map = nan(sz);
map(repmat(~mask,1,1,sz(3)))=map_r;



%% Shortest path
if shortestPath
    mapSourceDJ = reshape(mapSourceDJ,[],sz(3));
    path_r=nan(1,sz(3));
    [~,path_r(end)] = min(mapWeightDJ);
    for i_s=(sz(3)-1):-1:1
        path_r(i_s)=mapSourceDJ(path_r(i_s+1),i_s+1);
    end
    tmp = nan(sz(1),sz(2));
    tmp(:)=1:numel(tmp);
    tmp=tmp(~mask);
    path=tmp(path_r);
else
    path = [];
end

end