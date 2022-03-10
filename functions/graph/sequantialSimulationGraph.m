function psim = sequantialSimulationGraph(nj,gr)

% number of nodes in the 3d grid
n=prod(gr.snds(1:3));
nll=prod(gr.snds(1:2));

% initialize path
path = nan(gr.snds(3),nj);
path(1,:) = gr.s(1);

% Generating the transition matrices from one sta to the next
[~,~,source_sta]=ind2sub(gr.snds,gr.s);

mapB = cell(max(source_sta),1);
mapB{gr.snds(3)} = sparse(1,double(gr.lastNodes'),1,1,n);
for i_sta=(gr.snds(3)-1):-1:1
    id = source_sta==i_sta;
    mapB{i_sta} = mapB{i_sta+1} * sparse(gr.t(id),gr.s(id),double(gr.p(id)),n,n);
end


% loop through the simulated path
for i_sta=2:gr.snds(3)

    % get index of all edges from this stationary period
    id = source_sta==i_sta-1;

    % create the local transF (only edges from previous sta to next sta
    tranF = sparse(gr.s(id),gr.t(id),double(gr.p(id)),n,n);

    %
    mapF = sparse(1:nj,path(i_sta-1,:),1,nj,n) * tranF;

    M = mapF(:,nll*(i_sta-1)+(1:nll)) .* mapB{i_sta}(:,nll*(i_sta-1)+(1:nll));
    %[~,~,c]=ind2sub(gr.snds,gr.s(id_prev))

    % Sampling
    %     ids=nan(nj,1);
    %     for i_j=1:nj
    %         ids(i_j) = randsample(size(M,2),1,true,M(i_j,:));
    %     end
    tmp = cumsum(M,2);
    ids = sum( (rand(nj,1).*tmp(:,end))> tmp,2)+1;

    %
    path(i_sta,:) = ids + nll*(i_sta-1);

end

psim.path = uint32(path);
% too much memory
% [psim.lon, psim.lat, ~] = path2lonlat(psim.path,gr);