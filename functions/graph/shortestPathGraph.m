function sp = shortestPathGraph(gr)

G = digraph(string(gr.s),string(gr.t),-log(gr.p));
[p_s,~,edgepath] = shortestpath(G,string(gr.s(1)),string(gr.t(end)));
[sp.lat,sp.lon,~]=ind2sub(gr.snds,double(p_s));
sp.gs = gr.gs(edgepath);
sp.as = gr.as(edgepath);

%figure; hold on; borders('countries','k');
%plot(gr.lon(sp.lon), gr.lat(sp.lat),'-o')

end