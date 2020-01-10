function [out] = b2pl(in)
% compute pl from b
% in = b
% out = pl

lm = length(in);
natoms = round(log2(lm)); 		
if 2^natoms == lm 
in = in(lm)-fliplr(in);
in(1) = 0;
out = in;
else
	'ACCIDENT in btopl: length of input vector not OK: should be a power of 2'
end
