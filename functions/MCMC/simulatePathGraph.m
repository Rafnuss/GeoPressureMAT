function M = simulatePathGraph(gr)

% number of nodes in the 3d grid
n=prod(gr.snds(1:3));

% matrix of forward transition
transF = sparse(gr.s,gr.t,gr.p,n,n);

% matrix of backward transition
transB = sparse(gr.t,gr.s,gr.p,n,n);


% forward map
mapF = sparse(1,n);
% backward map
mapB = sparse(1,n);

for i_s=1:gr.snds(3)-1
    mapF(1,gr.s(1))=1;
    mapF = mapF*transF;
    
    mapB(1,gr.t(end))=1;
    mapB = mapB*transB;
end

map = mapF.*mapB;
map(1,[gr.s(1) gr.t(end)])=1;

M = reshape(full(map),gr.snds);
M = M./nansum(M,[1 2]);

map2=map;
map2(1,gr.s(1))=0;
map2 = map * mat;

mat = map2 \ map;

% tmp = transF.* spfun(@(x) 1./x, sparse(M(:)));
tmp = transF .* sparse(M(:));

% forward map
mapF = sparse(1,n);

for i_s=1:gr.snds(3)-1
    mapF(1,gr.s(1))=1;
    mapF = mapF*tmp;
end
mapF = reshape(full(mapF),gr.snds);
mapF = mapF./nansum(mapF,[1 2]);


figure; subplot(1,2,1); imagesc(M(:,:,4))
 subplot(1,2,2); imagesc(mapF(:,:,4))

end