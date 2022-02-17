function out = geolight_convert(twl_gl) 

  tm = [twl_gl.tFirst ; twl_gl.tSecond];
  [tm,itmunique] = unique(tm);
  rise = [twl_gl.type == 1 ; twl_gl.type ~= 1];
  rise = rise(itmunique,:);
  [tm,itmsort] = sort(tm);
  rise = rise(itmsort,:);
  
  out = table(tm, rise, 'VariableNames',{'Twilight','Rise'});
end