addpath(genpath('../functions'))
addpath('~/Documents/GitHub/Flight-Matlab/functions/')
addpath('~/Documents/GitHub/Flight-Matlab/data/')

avonet = readtable("avonet_ebird.csv");

% GPS 15-20
% GPS with recovery 2.5
% pressure sensor 0.﻿37 
% light 0.32 - 1.3
% icarus 5

[.3 .5 1 2 5 20]/.05


figure; tiledlayout(2,1,'TileSpacing','tight','Padding','tight')
nexttile; hold on; box on; grid on;
ksdensity(log(avonet.Mass))

x_axis=[1 6 10 20 40 100 400];
xticks(log(x_axis))
yticks([])
xticklabels(x_axis)
xlim(log([1 400]))

nexttile; hold on; box on; grid on;
[f,xi]=ksdensity(log(avonet.Mass),'function','cdf');
plot(xi,1-f)
yti = fliplr(interp1(xi,1-f,log(x_axis)));
yticks(yti)
yticklabels(round(yti*100))
xticks(log(x_axis))
xticklabels(x_axis)

xticks(log(x_axis))
xticklabels(x_axis)
xlim(log([1 400]))
ylim([0 1])
xticklabels(x_axis*.05)

