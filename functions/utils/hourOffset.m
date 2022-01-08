function out = hourOffset(hr, varagin) 
   if nargin>1
       offset=varagin{1};
   else
       offset=0;
   end
  out = mod(hr - offset,24) + mod(offset,24);
end

