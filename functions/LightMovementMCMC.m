function [path,f,c,acc] = LightMovementMCMC(nj,path0,fixPath,prob,step,szG)

path0=path0(:);

% number of stationary period
ns = numel(path0);

% initialized empty path
path = repmat(path0,1,nj);


% Compute first prob for path0
f=nan(1,nj);
f(1)=prob(path(:,1));
assert(~isnan(f(1)),'Problem')

% current candidate
c = nan(1,nj);
c([1 2])=1;
acc=false(1,nj);


for j=2:nj
    % New Candidate
    [path_lat,path_lon] = ind2sub(szG,path(~fixPath,c(j)));

    path_lat_next = round(path_lat + randn(sum(~fixPath),1) .* step(~fixPath,1));
    path_lat_next(path_lat_next<1) = -path_lat_next(path_lat_next<1)+1;
    path_lat_next(path_lat_next>szG(1)) = szG(1)-(path_lat_next(path_lat_next>szG(1))-szG(1));

    path_lon_next = round(path_lon + randn(sum(~fixPath),1) .* step(~fixPath,2));
    path_lon_next(path_lon_next<1) = -path_lon_next(path_lon_next<1)+1;
    path_lon_next(path_lon_next>szG(2)) = szG(2)-(path_lon_next(path_lon_next>szG(2))-szG(2));

    path(~fixPath,j) = sub2ind(szG,path_lat_next,path_lon_next);

    % Compute probability
    % prob_mvt(path(:,j))
    f(j)=prob(path(:,j));

    % Accept/reject
    if j>1
        if f(j)/f(c(j))>rand
            c(j+1)=j;
            acc(j)=true;
        else
            c(j+1) = c(j); 
        end
    end
    
end

end