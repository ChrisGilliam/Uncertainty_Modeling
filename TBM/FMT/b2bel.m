function [out] = b2bel(in)
% computing bel from b
% in = b
% out = bel

if in(1) == 1
	'ACCIDENT in btobel: you try to normalize with a 1 on m(¯)'
else
	
lm = length(in);
natoms = round(log2(lm)); 

if 2^natoms == lm 
in = (in-in(1));
in(1) = 0;
out = in;
else
	'ACCIDENT in btobel: length of input vector not OK: should be a power of 2'
end
end
