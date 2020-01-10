function [out] = w2q(in)
% computing FMT from w to q. use algortihm mtoq
% in = w vector
% out = q vector

lm = length(in);
natoms = round(log2(lm)); 
lin = log(in);
if 2^natoms == lm 
   if min(in)>0
      out = exp(sum(lin) - mtoq(lin));
      %		out = exp(-mtoq(log(in))); this is not OK, paper on matrix misses it.
   else
      'accident in wtoq: one of the weigths are non positive'
   end
else
   'ACCIDENT in wtoq: length of input vector not OK: should be a power of 2'
end