function twl_gl = export2GeoLight(twl)

d = diff(twl.Twilight);
id = find(d < 1 & d~=0);
twl_gl = table(twl.Twilight(id), twl.Twilight(id+1), -twl.Rise(id)+2,'VariableNames',{'tFirst','tSecond','type'});
end