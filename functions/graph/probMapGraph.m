function [M, E] = probMapGraph(gr)

% number of nodes in the 3d grid
n=prod(gr.snds(1:3));

% matrix of forward transition
% trans is O*T in paper. 
trans = sparse(gr.s,gr.t, gr.po(gr.t) .* double(gr.pt),n,n);

% forward map
mapF = sparse(1,n);
% backward map
mapB = sparse(n,1);

for i_s=1:gr.snds(3)-1
    mapF(1,gr.firstNodes)=1;
    mapF = mapF*trans;
    
    mapB(gr.lastNodes,1)=1;
    mapB = trans*mapB;
end

mapF(1,gr.s(1))=1;
mapB(gr.lastNodes,1)=1;

map = mapF.*mapB';


M = reshape(full(map),gr.snds);
% M = M./sum(M,3);

if nargout>1
    E = diag(mapF) * trans * diag(mapB);
    E = full(E(sub2ind([n n],gr.s,gr.t)));
end
    

end