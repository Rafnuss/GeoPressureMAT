function psim = sequantialSimulationGraph(nj,gr)

% number of nodes in the 3d grid
n=prod(gr.snds(1:3));
nll=prod(gr.snds(1:2));

% initialize path
path = nan(gr.snds(3),nj);
path(1,:) = gr.s(1);
% gr.lastNodes


% Generating the transition matrices from one sta to the next
[~,~,source_sta]=ind2sub(gr.snds,gr.s);
transF = cell(max(source_sta),1);
transB = cell(max(source_sta),1);
for i_sta = 1:max(source_sta)
    id = source_sta==i_sta;
    transF{i_sta} = sparse(gr.s(id),gr.t(id),gr.p(id),n,n);
    transB{i_sta} = sparse(gr.t(id),gr.s(id),gr.p(id),n,n);
end

% generate the simulation path 
sim_path = randperm(gr.snds(3));

% loop through the simulated path
for i_sim_path=1:numel(sim_path)
    if ~isnan(path(sim_path(i_sim_path),1))
        continue
    end

    i_prev = find(~isnan(path(1:sim_path(i_sim_path),1)),1,'last');
    i_next = find(~isnan(path(sim_path(i_sim_path):end,1)),1,'first')+sim_path(i_sim_path)-1;

    transF_tmp = sparse(1:nj,path(i_prev,:),1,nj,n) ;
    for i_c = (i_prev):(sim_path(i_sim_path)-1)
        transF_tmp = transF_tmp*transF{i_c};
    end
    if isempty(i_next)
        transB_tmp = sparse(repmat((1:nj)',1,numel(gr.lastNodes)),repmat(gr.lastNodes',nj,1),1,nj,n);
        i_next = gr.snds(3);
    else
        transB_tmp = sparse(1:nj,path(i_next,:),1,nj,n) ;
    end
    for i_c = (i_next-1):-1:sim_path(i_sim_path)
        transB_tmp = transB_tmp*transB{i_c};
    end

    M = transF_tmp .* transB_tmp;
    %[~,~,c]=ind2sub(gr.snds,find(M(1,:)>0))


    A=M(:,nll*(sim_path(i_sim_path)-1)+(1:nll));
    ids=nan(nj,1);
    for i_j=1:nj
        ids(i_j) = randsample(size(A,2),1,true,A(i_j,:));
    end
%     A=M(:,nll*(sim_path(i_sim_path)-1)+(1:nll));
%     tmp=cumsum(A,2);
%     sum( (rand(nj,1).*tmp(:,end))> tmp,2);

    % 
    path(sim_path(i_sim_path),:) = ids + nll*(sim_path(i_sim_path)-1);
    
end

psim.path = path;
[psim.lon, psim.lat, ~] = path2lonlat(psim.path,gr);