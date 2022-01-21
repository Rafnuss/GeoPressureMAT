function sp = shortestPathGraph(gr)

% Build the grpah
G = digraph(gr.s,gr.t,-log(gr.p));

% search shortest path
[sp.path,~,edgepath] = shortestpath(G,gr.s(1),gr.t(end));

% get lat and lon of shortest path
[sp.lon, sp.lat, ~] = path2lonlat(sp.path,gr);

% as well as edge information
sp.gs = gr.gs(edgepath);
sp.as = gr.as(edgepath);
sp.ws = gr.ws(edgepath);
sp.edge = edgepath;

if false
    figure; hold on; borders('countries','k');
    plot(sp.lon, sp.lat,'-o')
end

end