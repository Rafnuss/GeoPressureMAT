function gr = filterGraph(gr,quantity,thr)
% Create the graph (vectors of sink and source index in the matrix) from a
% proabibilty map and a threashold of ground_speed


% Filter the graph based on groundspeed threashold
id = find(abs(gr.(quantity))<thr);

% Find all nodes in path between the calibration location
% filter from begining to end
NodesInPath = [nearest(digraph(gr.s(id),gr.t(id)), gr.s(1),Inf); gr.s(1); gr.t(end)];
sel2=ismember(gr.s(id),NodesInPath) & ismember(gr.t(id),NodesInPath);
% and filter from end to begining
% small trick in case there is not a single end node. We create an extra
% node and connect the nodes from the last stationary period to it. Then
% check if in path of the newly created node. 
[~,~,Tt]=ind2sub(gr.snds,gr.t(id));
unTEndS = unique(gr.t(id(Tt==max(Tt))));
NodesInPath = [nearest(digraph([gr.t(id); ones(size(unTEndS))],[gr.s(id) ; unTEndS]), 1, Inf); gr.s(1)];
sel1=ismember(gr.s(id),NodesInPath) & ismember(gr.t(id),NodesInPath);

id = id(sel2&sel1);

assert(numel(id)>0,'no nodes left on the graph')

gr.s=gr.s(id);
gr.t=gr.t(id);

for varName=["ws", "as", "gs"]
    if isfield(gr,varName)
        gr.(varName) = gr.(varName)(id);
    end
end

end