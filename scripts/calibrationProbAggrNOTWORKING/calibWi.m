function q = calibWi(x,y,N,Ntest,block_time_i_b,id_calib,wi)

q=nan(Ntest,1);
for ii=1:Ntest
    ii_t = round(rand*N)+(1:block_time_i_b*24);
    if max(ii_t)<N
        e = x(:,:,ii_t)-y(:,:,ii_t); 
        e = e - mean(e,3);
        se = -sum(e.^2,3,'omitnan');
        % prob = exp(-sum(e.^2,3,'omitnan') .* (alpha + (1-alpha)*1/size(e,3)));
        prob = exp( se.* wi);
        % figure; hold on; imagesc(lon{lt},lat{lt},prob); plot(gLON(id_calib),gLAT(id_calib),'rx')
        q(ii) = sum(prob(prob<prob(id_calib)))./sum(prob(:));
    end
end
%h = chi2gof(reshape(q(lt,(i_s==1)+1,i_b,:),1,[]),'CDF',makedist('Uniform'))

end