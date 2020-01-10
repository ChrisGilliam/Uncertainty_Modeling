function [out] = pl2q(in)
% computing FMT from pl to q. 
% in = pl vector
% out = q vector

out = abs(b2m(in));
out(1)=1;
