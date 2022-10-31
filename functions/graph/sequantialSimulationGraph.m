function psim = sequantialSimulationGraph(nj,gr)

% number of nodes in the 3d grid
n=prod(gr.snds(1:3));
nll=prod(gr.snds(1:2));

% matrix of forward transition
% trans is O*T in paper.
trans = sparse(gr.s,gr.t, gr.po(gr.t) .* double(gr.pt),n,n);

mapF = cell(gr.snds(3),1);
mapF{1} = sparse(1,double(gr.firstNodes),1,1,n);
for i_sta=1:(gr.snds(3)-1)
    mapF{i_sta+1} = mapF{i_sta} * trans;
end

% Initialize path
path = nan(gr.snds(3),nj);
mapB = sparse(1,double(gr.lastNodes),1,1,n);
M = mapF{end}.*mapB;
tmp = cumsum(M,2);
path(end,:) = sum( (rand(nj,1).*tmp(:,end))> tmp,2)+1;

transT = sparse(gr.t, gr.s, double(gr.pt),n,n);

% loop through the states
for i_sta=(gr.snds(3)-1):-1:1

    % 
    M = mapF{i_sta} .* (sparse(1:nj,path(i_sta+1,:),1,nj,n) * transT);

    % Sampling
    tmp = cumsum(M(:,nll*(i_sta-1)+(1:nll)),2);
    ids = sum( (rand(nj,1).*tmp(:,end))> tmp,2)+1;

    %
    path(i_sta,:) = ids + nll*(i_sta-1);

end

psim.path = uint32(path);
% too much memory
% [psim.lon, psim.lat, ~] = path2lonlat(psim.path,gr);
end