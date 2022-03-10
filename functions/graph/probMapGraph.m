function M = probMapGraph(gr)

% number of nodes in the 3d grid
n=prod(gr.snds(1:3));

% matrix of forward transition
transF = sparse(gr.s,gr.t,double(gr.p),n,n);

% matrix of backward transition
transB = sparse(gr.t,gr.s,double(gr.p),n,n);

% forward map
mapF = sparse(1,n);
% backward map
mapB = sparse(1,n);

for i_s=1:gr.snds(3)-1
    mapF(1,gr.s(1))=1;
    mapF = mapF*transF;
    
    mapB(1,gr.lastNodes)=1;
    mapB = mapB*transB;
end

mapF(1,gr.s(1))=1;
mapB(1,gr.lastNodes)=1;
map = mapF.*mapB;


M = reshape(full(map),gr.snds);
% M = M./sum(M,3);

end