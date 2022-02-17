function psim = gibbsGraph(nj,path0,fixPath,gr)

% Build the graph
G = digraph(gr.s,gr.t,gr.p);

% Initiate the simulated nodes (path
path0=path0(:);
path = repmat(path0,1,nj);

ss = find(~fixPath);

tic
for i_p=2:nj
    for i_ss=1:numel(ss)
        
        % get statioanry id
        i_s=ss(i_ss);
        
        % Edges and Nodes going OUT from the PREVIOUS nodes sampled (current path)
        [eidOut,nidOut] = outedges(G,path(i_s-1,i_p));
        % idOut = find(gr.s==path(i_s-1,i_p));

        if i_s<size(path,1)
            % Edge and Nodes going IN the NEXT nodes from the previous path
            [eidIn,nidIn] = inedges(G,path(i_s+1,i_p-1));
            % idIn = find(gr.t==path(i_s+1,i_p-1));
            
            % Find nodes that are presents in both (going out from previous and coming in in the next)
            [iOut, iIn] = ismember(nidOut, nidIn);
            % [iOut, iIn] = ismember(gr.t(idOut), gr.s(idIn));

            % Possible nodes
            nid = nidOut(iOut);
            eid = eidIn(iIn(iIn~=0));
            % nid = gr.t(idOut(iOut));

            % Probability of the possible nodes (products of edge in and edge out
            prob = G.Edges.Weight(eidOut(iOut)) .* G.Edges.Weight(eid);
            % prob = gr.p(idOut(iOut)) .* gr.p(idIn(iIn(iIn~=0)));
        
        else
            nid = nidOut;
            prob = G.Edges.Weight(eidOut);
            % prob = gr.p(idOut);
        end
        
        % Random sample the nodes
        ids = randsample(numel(prob),1,true,prob);
        path(i_s,i_p) = nid(ids);

        assert(~isnan(nid(ids)))
    end
end

% Add path to the structure to return
psim.path = path;

% Lat and lon
[psim.lon, psim.lat, ~] = path2lonlat(psim.path,gr);

% Add windspeed to simulated path
path_edge = reshape(findedge(G,path(1:end-1,:),path(2:end,:)),gr.snds(3)-1,nj);
psim.ws = gr.ws(path_edge);
psim.gs = gr.gs(path_edge);
psim.as = gr.as(path_edge);
end