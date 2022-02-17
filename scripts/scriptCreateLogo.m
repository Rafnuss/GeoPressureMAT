


file='../data/ECMWF/surface_pressure.nc';
spttime = datetime(double(ncread(file,'time'))/24 + datenum('1900-01-01 00:00:00'),'convertFrom','datenum');
glon=double(ncread(file,'longitude'));
glat=flip(double(ncread(file,'latitude')));

splt = permute(flip(ncread(file,'sp',[1 1 1 1000],[numel(glon) numel(glat) 1 2000]),2)/100,[2 1 4 3]);
splt = splt-mean(splt,3);

K = (1/(9^2))*ones(9);
figure;colormap(brewermap(9,'spectral'))
for i=1900:1:10000
 contourf(conv2(splt(250:end,:,i),K,'valid'))
 title(i)
 drawnow()
 pause(.5)
end

i=1914;
figure;colormap(brewermap(9,'spectral'))
contourf(conv2(splt(:,:,i),K,'valid'))
axis([68.8816  224.2733  276.7695  337])
caxis([ -2.1085 14])
%contourf(conv2(splt(280:end,1:180,i),K,'valid'))