function [path,prob_path] = Gibbs(nj,path0,fixPath,prob,prob_mvt)

path0 = path0(:);
fixPath = fixPath(:);

% initialized empty path
path = repmat(path0,1,nj);

% log prob of each path
prob_path = zeros(1,nj);

% % Compute first prob for path0
% f=nan(1,nj);
% f(1)=prob(path(:,1));
% assert(~isnan(f(1)),'Problem')
% 
% % current candidate
% c = nan(1,nj);
% c([1 2])=1;
% acc=false(1,nj);

% find stationary period index to simulation (i.e., not fixed)
ss = find(~fixPath);

sprob = sort(prob);
tmp=sum(cumsum(sprob)<0.01);
idm = sprob(sub2ind(size(sprob),tmp+1,1:size(sprob,2)));

for j=2:nj
    for i_ss=1:numel(ss)
        % get statioanry id
        i_s=ss(i_ss);
        
        % prob static
        id=prob(:,i_s)>=idm(i_s);
        
        prob_next = ones(size(prob(:,i_s)));
        if i_s~=numel(path0)
            prob_next(id) = prob_mvt(path(i_s+1,j-1),i_s,id);
            prob_next(~id) = 0;
        end
        
        prob_prev = ones(size(prob(:,i_s)));
        if i_s~=1
            prob_prev(id) = prob_mvt(path(i_s-1,j),i_s-1,id);
            prob_prev(~id) = 0;
        end
        
        prob_tmp = prob(:,i_s) .* prob_next .* prob_prev;
        
        % Sample prob
        
        % Option 1: work, but slow
        % path(i_s,j) = randsample(numel(prob_tmp),1,true,prob_tmp);
        
        % Option 2: manual, 10-20% faster
        id_0 = find(prob_tmp>0&~isnan(prob_tmp));
        [prob_tmp_s,prob_tmp_id] = sort(prob_tmp(id_0));
        prob_tmp_s_cum = cumsum(prob_tmp_s)./sum(prob_tmp_s(:));
        [~,id_sampled] = min(abs(prob_tmp_s_cum-rand));
        path(i_s,j) = id_0(prob_tmp_id(id_sampled));
        prob_path(j) = prob_path(j)+log(prob(path(i_s,j),i_s))+log(prob_prev(path(i_s,j)));
        if i_ss==numel(ss)
            prob_path(j) = prob_path(j)+log(prob_next(path(i_s,j)));
        end
    end
end

end