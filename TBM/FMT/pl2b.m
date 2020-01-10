function [out] = pl2b(in)
% compute b from pl
% in = pl
% out = b

lm = length(in);
natoms = round(log2(lm)); 		
if 2^natoms == lm 
in = 1-fliplr(in);
out = in;
else
	'ACCIDENT in btopl: length of input vector not OK: should be a power of 2'
end
