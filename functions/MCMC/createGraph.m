function [gr,nds] = createGraph(prob_map,lat,lon,calib,actEffort,thr_prob_percentile)
% Create the graph (vectors of sink and source index in the matrix) from a
% proabibilty map and a threashold of ground_speed


% Start with building up the standard gr structure 
gr.lat = lat;
gr.lon = lon;
gr.snds = [numel(gr.lat), numel(gr.lon), size(prob_map,3)];
gr.actEffort = max(1,hours(actEffort));


% Create the index matrix of all coordinate (space-time) of the grid
% idx=reshape(1:numel(prob_map),size(prob_map));

% Normalize the probability map for each night
prob_map = prob_map ./ sum(prob_map,[1 2]);

% Set threashold corresponding to all nodes corresponding to
% thr_prob_percentile of the probability for each stationary period
tmp = sort(reshape(prob_map,[],gr.snds(3)));
thr_prob = tmp(sub2ind(size(tmp),sum(cumsum(tmp)<=(1-thr_prob_percentile)),1:gr.snds(3)));
nds = prob_map>=reshape(thr_prob,1,1,[]);

% Set first and last one to true only at the known location.
[~, tmp1] = min(abs(gr.lat(:)-calib.lat));
[~, tmp2] = min(abs(gr.lon(:)-calib.lon));
if any(isnat(calib.second_period))
    nds(:,:,1)=false;
    nds(tmp1,tmp2,1)=true;
else
    nds(:,:,[1 end])=false;
    nds(tmp1,tmp2,[1 end])=true;
end

assert(all(sum(nds,[1 2])>0),['No possible location at stationary period: ' num2str(find(sum(nds,[1 2])==0)')])

% Create the source S and target T
S=cell(gr.snds(3)-1,1);
T=cell(gr.snds(3)-1,1);
for i_s = 1:gr.snds(3)-1
    % Get index of the source and target according to the mask
    [S{i_s}, T{i_s}]= meshgrid(find(nds(:,:,i_s))+(i_s-1)*prod(gr.snds(1:2)), find(nds(:,:,i_s+1))+i_s*prod(gr.snds(1:2)));
end

% Convert the cells to matrix
gr.s = cell2mat(cellfun(@(x) x(:),S,'UniformOutput',false));
gr.t = cell2mat(cellfun(@(x) x(:),T,'UniformOutput',false));

[Slat,Slon,St]=ind2sub(gr.snds,gr.s);
[Tlat,Tlon,~]=ind2sub(gr.snds,gr.t);

% tmp1 = lldistkm([lon(Slon) lat(Slat)], [lon(Tlon) lat(Slat)],'pythagoran');
% tmp2 = lldistkm([lon(Slon) lat(Slat)], [lon(Slon) lat(Tlat)],'pythagoran');
% tmp =[tmp1 tmp2]./actEffortHr(St);
% 
% % Add the sign
% gr.gs = sum([sign(lon(Tlon)-lon(Slon)) sign(lat(Tlat)-lat(Slat))].*tmp .* [1 1i],2);
resolution=0.25*111;
gr.gs = resolution.*((Tlon-Slon).*cos(pi/180*lat(floor((Tlat+Slat)/2)))+1i.*(Tlat-Slat))./gr.actEffort(St);

end